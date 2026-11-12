"""Binding balance by 5' extension.

Primers are chosen with their Tm inside the window, but a primer whose template
site is folded binds more weakly at the annealing temperature than the rest of
the pool (higher dG_eff). Such a primer keeps its core (and its 3' end, so
coverage and priming position do not change) and gains template-matching bases
at its 5' end until its dG_eff at Ta is close to the pool median. The added
bases are shown in lowercase, so the core still reads as a Tm-matched primer.

This is the Visual-OMP idea of matching the amount bound at Ta rather than Tm
(SantaLucia 2007, Methods Mol Biol 402:3, myth 3), applied only where it is
needed and kept visible.
"""
from __future__ import annotations

from typing import Dict, List

from . import candidates as cand
from . import thermo


def display(p: cand.Primer) -> str:
    """Primer sequence with the added 5' bases in lowercase."""
    return p.seq[:p.ext5].lower() + p.seq[p.ext5:]


def _extended(p: cand.Primer, k: int, store) -> tuple:
    """(sequence, start, end) of p with k template bases added at its 5' end, or None."""
    L = store.lengths[p.seqid]
    if p.side == "F":
        s = p.start - k
        if s < 1:
            return None
        seq = store.region(p.seqid, s, p.end)[0]
        return seq, s, p.end
    e = p.end + k
    if e > L:
        return None
    seq = cand.revcomp(store.region(p.seqid, p.start, e)[0])
    return seq, p.start, e


def _site_up(p: cand.Primer, max_len: int, T: float, cfg, store, pad: int = 300):
    """Accessibility table of p's template strand around its site, and a function
    giving the 1-based window end for a site extended by k bases."""
    region, rs = store.region(p.seqid, p.start - pad, p.end + pad)
    n = len(region)
    if p.side == "F":
        up = cand._plfold_up(cand.revcomp(region), T, max_len, cfg)
        return up, (lambda k: n - (p.start - k - rs))
    up = cand._plfold_up(region, T, max_len, cfg)
    return up, (lambda k: p.end + k - rs + 1)


def _dg_eff(seq: str, side: str, start: int, end: int, seqid: str, T: float, cfg, store, up, end1) -> tuple:
    """(dG_eff, hairpin cost, site opening) at T for a primer at [start, end]."""
    r = thermo.analyze(seq, cfg.phys())
    tmp = cand.Primer(seq, side, seqid, start, end, 0.0, 0.0)

    def sub(a, b):
        if a < 1 or b > store.lengths[seqid]:
            return ""
        return store.region(seqid, a, b)[0]

    dh, ds = cand.template_dangles(tmp, sub)
    dup = r.dG(T) + dh - (T + thermo.T0) * ds / 1000.0
    hc = cand.hairpin_cost(seq, side, T, cfg)[0]
    so = cand._dg_open(up, end1, len(seq), T)
    return dup + hc + so, hc, so


def extension_options(p: cand.Primer, T: float, cfg, store, max_k: int):
    """Every valid 5' extension of p (k = 1..max_k) with its (dG_eff, hairpin cost,
    site opening) at T, ignoring dimers. Stops at an N or a run longer than max_poly."""
    if max_k <= 0:
        return []
    up, end1_of = _site_up(p, p.length + max_k, T, cfg, store)
    out = []
    for k in range(1, max_k + 1):
        ext = _extended(p, k, store)
        if ext is None:
            break
        seq, s, e = ext
        if any(c not in "ACGT" for c in seq) or cand.max_run(seq) > cfg.max_poly:
            break
        if cfg.variants is not None:   # never grow over a known strain variant
            new = (s, p.start - 1) if p.side == "F" else (p.end + 1, e)
            if cfg.variants.within(p.seqid, *new):
                break
        dg, hc, so = _dg_eff(seq, p.side, s, e, p.seqid, T, cfg, store, up, end1_of(k))
        if hc > cfg.hairpin_reject:
            continue
        out.append((k, seq, s, e, dg, hc, so))
    return out


def best_extension(p: cand.Primer, target: float, T: float, cfg, store, max_k: int):
    """The extension that brings dG_eff closest to target, if it improves on no
    extension by at least 0.2 kcal/mol; None otherwise."""
    base = abs(p.dg_eff - target)
    best = None
    for opt in extension_options(p, T, cfg, store, max_k):
        d = abs(opt[4] - target)
        if d < (best[0] if best else base) - 1e-9:
            best = (d, opt)
        if opt[4] <= target:
            break
    if best is None or best[0] > base - 0.2:
        return None
    return best[1]


def extension_job(args):
    """Process-pool worker for selection-time balancing: (primer, target, T, cfg, max_k)."""
    p, target, T, cfg, max_k = args
    from .design import _G
    try:
        return best_extension(p, target, T, cfg, _G["store"], max_k)
    except Exception:
        return None


def apply_extension(p: cand.Primer, opt, cfg, store) -> None:
    """Turn p into its extended form (in place)."""
    k, seq, s, e, dg, hc, so = opt
    if p.tm_core is None:
        p.tm_core, p.tm_neb_core = p.tm, p.tm_neb
    p.seq, p.start, p.end = seq, s, e
    p.ext5 += k
    p.tm = cfg.primer_tm(seq)
    p.tm_neb = thermo.tm_neb(seq)
    p.tm_phys = thermo.analyze(seq, cfg.phys()).tm_C
    p.gc = cand.gc_frac(seq)
    p.hairpin_cost = hc
    p.site_open_genomic = so
    p.dangle_dh, p.dangle_ds = cand.template_dangles(
        p, lambda a, b: store.region(p.seqid, a, b)[0] if a >= 1 and b <= store.lengths[p.seqid] else "")
    p.dg_eff = dg
    p.self_d3, p.self_dany = cand.thal_pair(cand.tailed(p, cfg), cand.tailed(p, cfg), cfg)
    p.notes = list(p.notes) + [f"{p.ext5} nt added at the 5' end (lowercase) to balance binding"]


def balance(chosen: Dict[str, cand.Pair], pools: Dict[str, int], ta: dict, cfg, store,
            names: Dict[str, Dict[str, str]]) -> List[str]:
    """Extend weak primers in place. Returns one note per extended primer."""
    T = ta["rec_grid"]
    at = ta["at_recommended"]
    vals = sorted(v["dg_eff"] for v in at.values())
    if not vals:
        return []
    median = vals[len(vals) // 2] if len(vals) % 2 else 0.5 * (vals[len(vals) // 2 - 1] + vals[len(vals) // 2])
    target = median
    notes: List[str] = []

    pool_members: Dict[int, List[cand.Primer]] = {}
    for tid, pr in chosen.items():
        pool_members.setdefault(pools[tid], []).extend([pr.fwd, pr.rev])

    # weakest first, so later checks see the extended partners
    weak = []
    for tid, pr in chosen.items():
        for side, p in (("F", pr.fwd), ("R", pr.rev)):
            v = at.get(names[tid][side])
            if v and v["dg_eff"] > target + cfg.balance_tol:
                weak.append((v["dg_eff"], tid, side, p, pr))
    weak.sort(key=lambda x: -x[0])

    for dg0, tid, side, p, pr in weak:
        room = cfg.amp(tid)[2] - pr.length
        max_k = min(cfg.max_ext, room)
        if max_k <= 0:
            notes.append(f"{names[tid][side]}: not extended (dG_eff {dg0:.2f}): the amplicon is already at "
                         f"the length limit")
            continue
        up, end1_of = _site_up(p, p.length + max_k, T, cfg, store)
        best = (abs(dg0 - target), 0, None, None, None, dg0, None, None)   # no extension
        why = set()
        dimer_with: Dict[str, float] = {}
        baseline: Dict[int, tuple] = {}
        for k in range(1, max_k + 1):
            ext = _extended(p, k, store)
            if ext is None:
                break
            seq, s, e = ext
            if any(c not in "ACGT" for c in seq) or cand.max_run(seq) > cfg.max_poly:
                why.add("the next template bases are a long run or N")
                break   # extending further cannot remove a run or an N
            dg, hc, so = _dg_eff(seq, side, s, e, p.seqid, T, cfg, store, up, end1_of(k))
            if hc > cfg.hairpin_reject:
                why.add("the extended primer forms a hairpin")
                continue
            full = (cfg.tail_fwd if side == "F" else cfg.tail_rev) + seq
            ok = True
            for q in pool_members[pools[tid]]:
                other = q if q is not p else None
                qfull = cand.tailed(q, cfg) if other is not None else full
                d3, dany = cand.thal_pair(full, qfull, cfg)
                # refuse only a dimer the extension creates or worsens: a pool relaxed to the
                # -7 cap may already hold a -6.5 dimer that the extension does not touch
                b3, bany = baseline.setdefault(id(q), cand.thal_pair(
                    cand.tailed(p, cfg), cand.tailed(q, cfg) if other is not None else cand.tailed(p, cfg), cfg))
                if cand.dimer_cost(d3, dany, cfg) >= 1000 and (d3 < b3 - 0.3 or dany < bany - 0.3
                                                               or cand.dimer_cost(b3, bany, cfg) < 1000):
                    ok = False
                    who = "itself" if other is None else next(
                        (names[t][sd] for t, prx in chosen.items() for sd, qq in (("F", prx.fwd), ("R", prx.rev))
                         if qq is q), "another primer")
                    dimer_with[who] = min(d3, dimer_with.get(who, 0.0))
                    break
            if not ok:
                continue
            score = abs(dg - target)
            if score < best[0] - 1e-9:
                best = (score, k, seq, s, e, dg, hc, so)
            if dg <= target:
                break   # reached the pool median; longer only overshoots
        if best[1] == 0 or best[0] > abs(dg0 - target) - 0.2:
            # extension does not bring it meaningfully closer to the pool median
            why |= {f"a strong dimer with {w} (3' dG down to {d:.1f} kcal/mol)" for w, d in dimer_with.items()}
            reason = "; ".join(sorted(why)) or "one more base overshoots the pool median"
            notes.append(f"{names[tid][side]}: not extended (dG_eff {dg0:.2f}, median {target:.2f}): {reason}")
            continue
        _, k, seq, s, e, dg, hc, so = best
        apply_extension(p, (k, seq, s, e, dg, hc, so), cfg, store)
        pr.length = pr.rev.end - pr.fwd.start + 1
        pr.gc = cand.gc_frac(store.region(pr.fwd.seqid, pr.fwd.start, pr.rev.end)[0])
        notes.append(f"{names[tid][side]}: +{k} nt at the 5' end, dG_eff at {T:g} C {dg0:.2f} -> {dg:.2f} "
                     f"kcal/mol (pool median {target:.2f})")
    return notes
