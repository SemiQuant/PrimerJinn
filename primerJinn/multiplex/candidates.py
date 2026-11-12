"""Candidate primers and primer pairs for one target.

Every primer position and length in the window around a target is enumerated,
filtered on the NEB-scale Tm window and simple sequence rules, then scored
with physical thermodynamics at the design annealing temperature:

    dG_eff = dG_duplex(Ta) + hairpin cost + cost of opening the template site

The template-site cost comes from one windowed ViennaRNA fold of the genomic
context per strand (early cycles). Short-listed pairs are then folded as
amplicon strands (later cycles) for site accessibility and 3' self-priming.

Coordinates: 1-based inclusive genome coordinates on the top strand.
A forward primer occupies [start, end] and its 3' end is `end`. A reverse
primer occupies [start, end] on the top strand, its sequence is the reverse
complement of that footprint, and its 3' end is `start`.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import primer3

from . import thermo
from .config import Config

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


@dataclass
class Primer:
    seq: str
    side: str          # 'F' or 'R'
    seqid: str
    start: int
    end: int
    tm: float                  # design scale (cfg.tm_scale)
    gc: float
    cost: float = 0.0
    tm_neb: float = 0.0        # NEB Q5 calculator scale, 500 nM
    dg_dup: float = 0.0        # physical duplex dG at Ta, kcal/mol
    hairpin_dg: float = 0.0    # primer3 thal hairpin dG at Ta
    hairpin_cost: float = 0.0
    self_d3: float = 0.0       # 3'-anchored homodimer dG at 37 C
    self_dany: float = 0.0
    site_open_genomic: float = 0.0
    dg_eff: float = 0.0
    tm_phys: float = 0.0
    dangle_dh: float = 0.0     # template dangling ends (kcal/mol, cal/K/mol)
    dangle_ds: float = 0.0
    ext5: int = 0              # bases added at the 5' end by binding balance (shown lowercase)
    n_variants: int = 0        # known strain variants under the primer (none in the 3' end)
    tm_core: Optional[float] = None       # design-scale Tm before extension
    tm_neb_core: Optional[float] = None
    notes: List[str] = field(default_factory=list)

    @property
    def three_prime(self) -> int:
        return self.end if self.side == "F" else self.start

    @property
    def length(self) -> int:
        return len(self.seq)

    @property
    def key(self) -> Tuple:
        return (self.seq, self.side, self.seqid, self.start)


@dataclass
class Pair:
    target_id: str
    fwd: Primer
    rev: Primer
    length: int
    gc: float
    cost: float = 0.0
    parts: Dict[str, float] = field(default_factory=dict)
    site_open_amp_f: float = 0.0
    site_open_amp_r: float = 0.0
    self_prime_top: float = 0.0
    self_prime_bottom: float = 0.0
    pair_d3: float = 0.0
    pair_dany: float = 0.0
    flags: List[str] = field(default_factory=list)

    @property
    def start(self) -> int:
        return self.fwd.start

    @property
    def end(self) -> int:
        return self.rev.end

    @property
    def insert(self) -> Tuple[int, int]:
        """Bases between the primers (what the reads cover after primer trimming)."""
        return self.fwd.end + 1, self.rev.start - 1


# ----------------------------------------------------------------------------
# small helpers
# ----------------------------------------------------------------------------
_DINUC = re.compile(r"((?:[ACGT]{2}))\1{4,}")


def max_run(s: str) -> int:
    best = run = 1
    for a, b in zip(s, s[1:]):
        run = run + 1 if a == b else 1
        best = max(best, run)
    return best


def gc_frac(s: str) -> float:
    return (s.count("G") + s.count("C")) / len(s) if s else 0.0


def dimer_cost(d3: float, dany: float, cfg: Config) -> float:
    """Pairwise dimer cost (DESIGN.md section 6)."""
    c = max(0.0, cfg.d3_soft - d3) ** 2 + 0.25 * max(0.0, cfg.dany_soft - dany) ** 2
    if d3 <= cfg.d3_hard or dany <= cfg.dany_hard:
        c += 1000.0
    return c


def hairpin_cost(seq: str, side: str, T: float, cfg: Config) -> Tuple[float, float]:
    """(cost kcal/mol, hairpin dG) of a primer's own hairpin at T. With adapter
    tails the tailed oligo is also scored, counting only what goes beyond the
    tail's own hairpin. One function for design, final tables and balance."""
    kw = cfg.thal_kw(T)

    def one(s):
        hp = primer3.calc_hairpin(s[-60:], **kw)
        if not hp.structure_found:
            return 0.0, 0.0
        dg = hp.dg / 1000.0
        return thermo.structure_cost(dg, T), dg

    cost, dg = one(seq)
    tail = cfg.tail_fwd if side == "F" else cfg.tail_rev
    if tail:
        base = one(tail)[0] if len(tail) >= 4 else 0.0
        c_t, dg_t = one(tail + seq)
        if c_t - base > cost:
            cost, dg = c_t - base, dg_t
    return cost, dg


def tailed(p: Primer, cfg: Config) -> str:
    tail = cfg.tail_fwd if p.side == "F" else cfg.tail_rev
    return tail + p.seq


def thal_pair(a: str, b: str, cfg: Config) -> Tuple[float, float]:
    """(d3, dany) in kcal/mol at 37 C. thal cannot align two sequences that are both
    longer than 60 nt, so tails are trimmed from the 5' end to fit."""
    if len(a) > 60 and len(b) > 60:
        a, b = a[-60:], b[-60:]
    kw = cfg.thal_kw(37.0)
    e1 = primer3.calc_end_stability(a, b, **kw).dg / 1000.0
    e2 = primer3.calc_end_stability(b, a, **kw).dg / 1000.0
    h = primer3.calc_heterodimer(a, b, **kw)
    dany = h.dg / 1000.0 if h.structure_found else 0.0
    return min(e1, e2, 0.0), min(dany, 0.0)


def template_dangles(p: "Primer", sub) -> Tuple[float, float]:
    """Template dangling-end dH/dS for primer p; sub(a, b) returns the top strand,
    1-based inclusive (returns '' outside the sequence)."""
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    if p.side == "F":   # template = bottom strand; complement of the top-strand bases
        b3, a3 = sub(p.end + 1, p.end + 1), sub(p.end, p.end)
        a5, b5 = sub(p.start, p.start), sub(p.start - 1, p.start - 1)
        b3, a3, a5, b5 = (comp.get(x, "") for x in (b3, a3, a5, b5))
    else:               # template = top strand
        b3, a3 = sub(p.start - 1, p.start - 1), sub(p.start, p.start)
        a5, b5 = sub(p.end, p.end), sub(p.end + 1, p.end + 1)
    return thermo.dangling_end_terms(b3, a3, a5, b5)


# ----------------------------------------------------------------------------
# accessibility of every primer site from one windowed fold per strand
# ----------------------------------------------------------------------------
def _plfold_up(seq: str, celsius: float, max_len: int, cfg: Config, window: int = 200, span: int = 150):
    """up[i][L] = P(window of length L ending at 1-based i is unpaired)."""
    import RNA

    from . import structure

    md = structure._md(celsius, cfg.struct_cond())
    md.window_size = min(window, len(seq))
    md.max_bp_span = min(span, md.window_size)
    fc = RNA.fold_compound(seq, md, RNA.OPTION_WINDOW)
    out: Dict[int, list] = {}

    def cb(v, v_size, i, maxsize, what, data):
        if what & RNA.PROBS_WINDOW_UP:
            data[i] = list(v)

    if not fc.probs_window(max_len, RNA.PROBS_WINDOW_UP, cb, out):
        raise RuntimeError("ViennaRNA probs_window failed")
    return out


def _dg_open(up, end1: int, length: int, celsius: float) -> float:
    row = up.get(end1)
    p = row[length] if row is not None and length < len(row) and row[length] is not None else 1e-300
    return -thermo.R_CAL / 1000.0 * (celsius + thermo.T0) * math.log(max(p, 1e-300))


def _amplicon_up(seq: str, celsius: float, max_len: int, cfg: Config):
    from . import structure

    return structure.unpaired_windows(seq, celsius, max_len, cfg.struct_cond())


# ----------------------------------------------------------------------------
# the per-target worker
# ----------------------------------------------------------------------------
def _side_reason(side: str, tm_hits: int, rej: Dict[str, int], cfg: Config) -> str:
    if tm_hits == 0:
        return f"no {side} primer within Tm {cfg.tm:.1f} +/- {cfg.tm_window:g} C"
    labels = (("gc", f"GC outside {cfg.primer_gc_min:.2f}-0.80"), ("poly", f"run > {cfg.max_poly}"),
              ("dinuc", "dinucleotide repeat"), ("3gc", "5 G/C in last 5"), ("ambiguous", "ambiguous base"),
              ("variant3", f"known variant in the last {cfg.variant_3prime_nt} nt"))
    parts = [f"{rej[k]} {lbl}" for k, lbl in labels if rej.get(k)]
    return f"{tm_hits} {side} primers in the Tm window but all failed sequence rules ({', '.join(parts)})"


def design_target(target, region: str, region_start: int, seq_len: int, cfg: Config,
                  keep_pairs: int = 120) -> dict:
    """Candidate pairs for one target, best first, WITHOUT the amplicon-strand
    folds (see amplicon_fold_job / apply_amplicon_fold, run across all cores).

    region is the genome sequence (upper case) of [region_start, region_start+len-1].
    Returns {'target': id, 'pairs': [Pair], 'stats': {...}, 'fail': reason or None}.
    """
    t = target
    if t.id in cfg.amp_by_target:            # a larger size class (--one_pool merged target)
        import copy as _copy
        cfg = _copy.copy(cfg)
        cfg.amp_min, cfg.amp_opt, cfg.amp_max = cfg.amp_by_target[t.id]
    T0, W, Ta = cfg.tm, cfg.tm_window, cfg.ta
    kmin, kmax = cfg.primer_min, cfg.primer_max
    rs = region_start
    n = len(region)
    stats: Dict[str, int] = {}
    rej = {"F": {}, "R": {}}

    def sub(a: int, b: int) -> str:  # genome 1-based inclusive
        return region[a - rs:b - rs + 1]

    gc_target = min(max(cfg.genome_gc, 0.30), 0.65)

    # geometry first, so the failure names the real cause
    room_f = t.start - 1 - cfg.margin          # bases available for a forward primer
    room_r = seq_len - (t.end + cfg.margin)    # bases available for a reverse primer
    if room_f < kmin or room_r < kmin:
        side = "start" if room_f < kmin else "end"
        return {"target": t.id, "pairs": [], "stats": stats,
                "fail": f"too close to the {side} of {t.seqid} (needs at least {kmin + cfg.margin} bp outside "
                        f"the target for a primer)"}


    def primer_rules(seq: str, side: str) -> Optional[Tuple[float, List[str]]]:
        """Hard filters and cheap cost terms. None = rejected."""
        r = rej[side]
        if any(c not in "ACGT" for c in seq):
            r["ambiguous"] = r.get("ambiguous", 0) + 1
            return None
        g = gc_frac(seq)
        if g < cfg.primer_gc_min or g > 0.80:
            r["gc"] = r.get("gc", 0) + 1
            return None
        if max_run(seq) > cfg.max_poly:
            r["poly"] = r.get("poly", 0) + 1
            return None
        if _DINUC.search(seq):
            r["dinuc"] = r.get("dinuc", 0) + 1
            return None
        end5 = seq[-5:]
        gc3 = end5.count("G") + end5.count("C")
        if gc3 == 5:
            r["3gc"] = r.get("3gc", 0) + 1
            return None
        cost = 0.5 * ((g - gc_target) / 0.10) ** 2
        notes = []
        if gc3 > 3:
            cost += 1.0
            notes.append("4 G/C in last 5")
        if gc3 == 0:
            cost += 0.5
            notes.append("no G/C in last 5")
        return cost, notes

    def collect(found, side, make):
        """Apply the rules to every length in the Tm window, keep the 2 cheapest."""
        scored = []
        for tm, pos, seq in found:
            r = primer_rules(seq, side)
            if r is None:
                continue
            p = make(pos, seq, tm)
            # Tm is a coarse filter (the window); binding at Ta (dG_eff) decides
            # (SantaLucia 2007, myth 3: match the amount bound at Ta, not Tm)
            cost, notes = r[0] + ((tm - T0) / 1.5) ** 2, list(r[1])
            if cfg.variants is not None:
                hits = cfg.variants.primer_hits(side, p.seqid, p.start, p.end)
                if any(h < cfg.variant_3prime_nt for h in hits):
                    rej[side]["variant3"] = rej[side].get("variant3", 0) + 1
                    continue
                if hits:
                    p.n_variants = len(hits)
                    cost += cfg.variant_penalty * len(hits)
                    notes.append(f"{len(hits)} known variant(s) under the primer, "
                                 f"{min(hits) + 1} nt or more from the 3' end")
            p.cost, p.notes = cost, notes
            scored.append((cost, len(scored), p))
        return [p for _, _, p in sorted(scored)[:2]]

    # -------- enumerate forward primers (3' end before the target)
    fwd: List[Primer] = []
    tm_hits = {"F": 0, "R": 0}
    f3_max = t.start - 1 - cfg.margin
    f3_min = t.end + cfg.margin + 2 * kmin - cfg.amp_max
    for e in range(max(f3_min, rs + kmin - 1), f3_max + 1):
        found = []
        for k in range(kmin, kmax + 1):
            s5 = e - k + 1
            if s5 < rs:
                break
            seq = sub(s5, e)
            if "N" in seq:
                break
            tm = cfg.primer_tm(seq)
            if tm > T0 + W:
                break
            if tm >= T0 - W:
                found.append((tm, s5, seq))
        tm_hits["F"] += len(found)
        fwd += collect(found, "F", lambda s5, seq, tm, e=e: Primer(seq, "F", t.seqid, s5, e, tm, gc_frac(seq),
                                                                     tm_neb=thermo.tm_neb(seq)))
    # -------- enumerate reverse primers (3' end after the target)
    rev: List[Primer] = []
    r3_min = t.end + 1 + cfg.margin
    r3_max = t.start - cfg.margin - 2 * kmin + cfg.amp_max
    for s3 in range(r3_min, min(r3_max, rs + n - kmin) + 1):
        found = []
        for k in range(kmin, kmax + 1):
            e5 = s3 + k - 1
            if e5 >= rs + n:
                break
            seq = revcomp(sub(s3, e5))
            if "N" in seq:
                break
            tm = cfg.primer_tm(seq)
            if tm > T0 + W:
                break
            if tm >= T0 - W:
                found.append((tm, e5, seq))
        tm_hits["R"] += len(found)
        rev += collect(found, "R", lambda e5, seq, tm, s3=s3: Primer(seq, "R", t.seqid, s3, e5, tm, gc_frac(seq),
                                                                       tm_neb=thermo.tm_neb(seq)))
    stats["fwd_enumerated"], stats["rev_enumerated"] = len(fwd), len(rev)
    stats["fwd_tm_window"], stats["rev_tm_window"] = tm_hits["F"], tm_hits["R"]
    for side, d in rej.items():
        for k, v in d.items():
            stats[f"{side}_rej_{k}"] = v
    if not fwd or not rev:
        why = [_side_reason(nm, tm_hits[sd], rej[sd], cfg) for nm, sd, lst in
               (("forward", "F", fwd), ("reverse", "R", rev)) if not lst]
        return {"target": t.id, "pairs": [], "stats": stats, "fail": "; ".join(why)}

    # -------- physical thermodynamics, hairpin (gene-specific and tailed), homodimer
    phys = cfg.phys()

    survivors = {"F": [], "R": []}
    lost = {"F": {}, "R": {}}
    for p in fwd + rev:
        r = thermo.analyze(p.seq, phys)
        p.tm_phys = r.tm_C
        p.dangle_dh, p.dangle_ds = template_dangles(p, sub)
        # duplex plus the template's dangling ends (SantaLucia and Hicks 2004 Table 3)
        p.dg_dup = r.dG(Ta) + p.dangle_dh - (Ta + thermo.T0) * p.dangle_ds / 1000.0
        p.hairpin_cost, p.hairpin_dg = hairpin_cost(p.seq, p.side, Ta, cfg)
        full = tailed(p, cfg)
        if p.hairpin_cost > cfg.hairpin_reject:
            lost[p.side]["hairpin"] = lost[p.side].get("hairpin", 0) + 1
            continue
        p.self_d3, p.self_dany = thal_pair(full, full, cfg)
        c = dimer_cost(p.self_d3, p.self_dany, cfg)
        if c >= 1000:
            lost[p.side]["self-dimer"] = lost[p.side].get("self-dimer", 0) + 1
            continue
        p.cost += c + (max(0.0, p.hairpin_cost - cfg.hairpin_penalty) / 0.5) ** 2
        survivors[p.side].append(p)
    fwd, rev = survivors["F"], survivors["R"]
    if not fwd or not rev:
        why = [f"every {nm} primer has a strong hairpin or self-dimer ({lost[sd]})"
               for nm, sd, lst in (("forward", "F", fwd), ("reverse", "R", rev)) if not lst]
        return {"target": t.id, "pairs": [], "stats": stats, "fail": "; ".join(why)}

    # -------- template-site accessibility in genomic context (one fold per strand)
    up_top = _plfold_up(region, Ta, kmax, cfg)
    up_bot = _plfold_up(revcomp(region), Ta, kmax, cfg)
    survivors = {"F": [], "R": []}
    for p in fwd + rev:
        if p.side == "F":
            # binds the bottom strand; top index a..b maps to bottom n-1-b..n-1-a
            p.site_open_genomic = _dg_open(up_bot, n - (p.start - rs), p.length, Ta)
        else:
            p.site_open_genomic = _dg_open(up_top, p.end - rs + 1, p.length, Ta)
        if p.site_open_genomic > cfg.site_open_reject:
            stats["rej_site_structure"] = stats.get("rej_site_structure", 0) + 1
            continue
        # dG_eff carries both costs into the pool-uniformity term; the extra penalty
        # applies only above the threshold so they are not charged twice
        p.cost += (max(0.0, p.site_open_genomic - cfg.site_open_penalty) / cfg.site_open_penalty) ** 2
        p.dg_eff = p.dg_dup + p.hairpin_cost + p.site_open_genomic
        survivors[p.side].append(p)
    fwd, rev = survivors["F"], survivors["R"]
    if not fwd or not rev:
        side = "forward" if not fwd else "reverse"
        return {"target": t.id, "pairs": [], "stats": stats,
                "fail": f"template structure blocks every {side} primer site (opening cost > "
                        f"{cfg.site_open_reject} kcal/mol at {Ta:.1f} C)"}

    # -------- pair up
    raw = []
    gc_amp_target = min(max(cfg.genome_gc, 0.30), 0.60)

    def amp_gc_cost(g):
        d = max(0.0, abs(g - gc_amp_target) - 0.08)
        return (d / 0.05) ** 2

    for f in fwd:
        for r in rev:
            L = r.end - f.start + 1
            if L < cfg.amp_min or L > cfg.amp_max:
                continue
            if f.end >= t.start or r.start <= t.end:
                continue
            amp_seq = sub(f.start, r.end)
            if "N" in amp_seq:   # an amplicon across an N gap has unknown length and content
                stats["rej_amplicon_N"] = stats.get("rej_amplicon_N", 0) + 1
                continue
            g = gc_frac(amp_seq)
            pre = (f.cost + r.cost + ((L - cfg.amp_opt) / 25.0) ** 2
                   + amp_gc_cost(g) + 0.5 * ((f.tm - r.tm) / 1.0) ** 2)
            raw.append((pre, f, r, L, g))
    stats["pairs_in_range"] = len(raw)
    if not raw:
        why = f"no primer pair gives an amplicon of {cfg.amp_min} to {cfg.amp_max} bp"
        if stats.get("rej_amplicon_N"):
            why += " without crossing unknown (N) sequence"
        return {"target": t.id, "pairs": [], "stats": stats, "fail": why}
    raw.sort(key=lambda x: x[0])
    use: Dict[Tuple, int] = {}
    shortlist = []
    for item in raw:
        _, f, r, _, _ = item
        if use.get(f.key, 0) >= 6 or use.get(r.key, 0) >= 6:
            continue
        use[f.key] = use.get(f.key, 0) + 1
        use[r.key] = use.get(r.key, 0) + 1
        shortlist.append(item)
        if len(shortlist) >= max(2 * keep_pairs, 200):
            break

    # -------- within-pair dimer
    pairs: List[Pair] = []
    for pre, f, r, L, g in shortlist:
        d3, dany = thal_pair(tailed(f, cfg), tailed(r, cfg), cfg)
        dc = dimer_cost(d3, dany, cfg)
        if dc >= 1000:
            stats["rej_pair_dimer"] = stats.get("rej_pair_dimer", 0) + 1
            continue
        pr = Pair(t.id, f, r, L, g)
        pr.pair_d3, pr.pair_dany = d3, dany
        pr.parts = {"fwd": f.cost, "rev": r.cost, "length": ((L - cfg.amp_opt) / 25.0) ** 2,
                    "amp_gc": amp_gc_cost(g), "tm_pair": 0.5 * ((f.tm - r.tm) / 1.0) ** 2, "pair_dimer": dc}
        pr.cost = pre + dc
        pairs.append(pr)
    if not pairs:
        return {"target": t.id, "pairs": [], "stats": stats, "fail": "every pair forms a strong primer dimer"}
    pairs.sort(key=lambda p: p.cost)
    pairs = pairs[:keep_pairs]
    return {"target": t.id, "pairs": pairs, "stats": stats, "fail": None, "region": (region, rs)}


# ----------------------------------------------------------------------------
# amplicon-strand structure (later cycles), one job per pair, run across all cores
# ----------------------------------------------------------------------------
def amplicon_fold_job(args):
    """Fold both strands of one amplicon (gene-specific part, no tails: tails are
    the same for every amplicon and their panhandle would mask template effects).
    Returns (site_open_f, site_open_r, self_prime_top, self_prime_bottom)."""
    amp, kf, kr, Ta, cfg = args
    bot = revcomp(amp)
    m = len(amp)
    up_t = _amplicon_up(amp, Ta, kr, cfg)
    up_b = _amplicon_up(bot, Ta, kf, cfg)
    return (_dg_open(up_b, m, kf, Ta), _dg_open(up_t, m, kr, Ta), 1.0 - up_t[m][1], 1.0 - up_b[m][1])


def apply_amplicon_fold(pr: Pair, res, cfg: Config) -> None:
    pr.site_open_amp_f, pr.site_open_amp_r, pr.self_prime_top, pr.self_prime_bottom = res
    extra = 0.0
    pen = cfg.site_open_penalty

    def over(x):
        return (max(0.0, x - pen) / pen) ** 2

    for so, gen in ((pr.site_open_amp_f, pr.fwd.site_open_genomic), (pr.site_open_amp_r, pr.rev.site_open_genomic)):
        if so > gen:  # the larger of the two contexts counts
            extra += over(so) - over(gen)
    for sp, name in ((pr.self_prime_top, "top"), (pr.self_prime_bottom, "bottom")):
        if sp > 0.1:
            extra += 2.0 * ((sp - 0.1) / 0.2) ** 2
        if sp > cfg.self_prime_flag:
            pr.flags.append(f"{name} strand 3' end paired with p={sp:.2f} (self-priming risk)")
    pr.parts["amplicon_structure"] = extra
    pr.cost += extra
