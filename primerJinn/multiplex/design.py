"""Multiplex pool design: from genome positions to a primer pool.

    result = design_pool(sites, genome_meta, cfg)

Steps (docs/DESIGN.md):
  1. group positions into targets (exact dynamic programme)
  2. per target, in parallel: enumerate, filter and score primers and pairs,
     including template and amplicon secondary structure
  3. uniform binding: penalise pairs whose dG_eff at Ta is far from the pool median
  4. genome-wide specificity of every candidate primer (seed search + thal)
  5. pools by graph colouring; primer-dimer and product costs within pools
  6. choose one pair per target (simulated annealing)
  7. final checks: coverage, products, dimers, annealing temperature, amplicon structure
"""
from __future__ import annotations

import copy
import itertools
import math
import multiprocessing as mp
import statistics
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from . import candidates as cand
from .config import Config
from .optimize import Problem, conflict_graph, dsatur
from .targets import Site, Target, group_sites, span_limits


class DesignError(ValueError):
    """A design that cannot be produced from the inputs (reported without a traceback)."""


def _log(msg: str, verbose: bool = True):
    if verbose:
        print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _ctx():
    return mp.get_context("fork") if "fork" in mp.get_all_start_methods() else mp.get_context()


# ----------------------------------------------------------------------------
# sequence access
# ----------------------------------------------------------------------------
class SeqStore:
    """Upper-case sequences by id, loaded lazily from an indexed FASTA."""

    def __init__(self, fasta: str, fai: Optional[str] = None):
        import pysam

        self._fa = pysam.FastaFile(fasta, filepath_index=fai) if fai else pysam.FastaFile(fasta)
        self.lengths = dict(zip(self._fa.references, self._fa.lengths))
        self._cache: Dict[str, str] = {}

    def seq(self, seqid: str) -> str:
        if seqid not in self._cache:
            self._cache[seqid] = self._fa.fetch(seqid).upper()
        return self._cache[seqid]

    def region(self, seqid: str, start: int, end: int) -> Tuple[str, int]:
        start = max(1, start)
        end = min(self.lengths[seqid], end)
        return self.seq(seqid)[start - 1:end], start


# ----------------------------------------------------------------------------
# results
# ----------------------------------------------------------------------------
@dataclass
class Design:
    cfg: Config
    sites: List[Site]
    targets: List[Target]
    pools: Dict[str, int]                     # target id -> pool (1-based)
    chosen: Dict[str, cand.Pair]              # target id -> pair
    failed: Dict[str, str]                    # target id -> reason
    relaxed: Dict[str, str]                   # target id -> relaxation used
    coverage_problems: list
    products: List[dict]                      # unintended products within pools
    dimers: List[dict]                        # all within-pool primer pairs
    ta: dict                                  # annealing temperature advice
    amplicon_reports: Dict[str, object]
    offtarget_sites: Dict[str, List[dict]]    # primer name -> off-target sites
    energy: float = 0.0
    method: str = ""
    timings: Dict[str, float] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)
    names: Dict[str, Dict[str, str]] = field(default_factory=dict)   # target id -> {'F': name, 'R': name}
    balance_notes: List[str] = field(default_factory=list)
    dropped_sites: List[Site] = field(default_factory=list)   # requested but N in the reference
    passes: int = 1
    genome: dict = field(default_factory=dict)


# ----------------------------------------------------------------------------
# step 2: candidates per target (process pool)
# ----------------------------------------------------------------------------
_G: dict = {}


def _worker(args):
    target, cfg, keep = args
    store: SeqStore = _G["store"]
    pad = cfg.amp_max + 250
    region, rs = store.region(target.seqid, target.start - pad, target.end + pad)
    try:
        return cand.design_target(target, region, rs, store.lengths[target.seqid], cfg, keep_pairs=keep)
    except Exception as exc:  # report, do not kill the run
        return {"target": target.id, "pairs": [], "stats": {}, "fail": f"error: {exc!r}"}


def _relaxed(cfg: Config) -> Config:
    c = copy.copy(cfg)
    c.tm_window = cfg.tm_window * 2
    c.site_open_reject = 99.0
    c.hairpin_reject = 2.0
    c.max_poly = cfg.max_poly + 1
    c.primer_gc_min = max(0.08, cfg.primer_gc_min - 0.05)
    return c


def _fold_jobs(results, cfg):
    """Fold the amplicon strands of every short-listed pair of every target in one
    process pool (all cores busy even for one target), then re-rank the pairs."""
    jobs, where = [], []
    for r in results:
        if not r["pairs"]:
            continue
        region, rs = r.pop("region")
        for pr in r["pairs"]:
            amp = region[pr.fwd.start - rs:pr.rev.end - rs + 1]
            jobs.append((amp, pr.fwd.length, pr.rev.length, cfg.ta, cfg))
            where.append(pr)
    if not jobs:
        return
    if cfg.threads <= 1 or len(jobs) < 4:
        out = [cand.amplicon_fold_job(j) for j in jobs]
    else:
        with ProcessPoolExecutor(min(cfg.threads, len(jobs)), mp_context=_ctx()) as ex:
            out = list(ex.map(cand.amplicon_fold_job, jobs, chunksize=4))
    for pr, res in zip(where, out):
        cand.apply_amplicon_fold(pr, res, cfg)
    for r in results:
        r["pairs"].sort(key=lambda p: p.cost)


def _run_candidates(targets, cfg, store, keep, verbose):
    # load sequences in the parent: forked workers must not share a pysam file handle
    for sid in {t.seqid for t in targets}:
        store.seq(sid)
    _G["store"] = store
    jobs = [(t, cfg, keep) for t in targets]
    if cfg.threads <= 1 or len(jobs) == 1:
        results = [_worker(j) for j in jobs]
    else:
        with ProcessPoolExecutor(min(cfg.threads, len(jobs)), mp_context=_ctx()) as ex:
            results = list(ex.map(_worker, jobs))
    _fold_jobs(results, cfg)
    return results


def _candidates_with_fallback(targets, cfg, store, verbose):
    """Candidates for every target. A failing target with several sites (or a long
    range) is split in two; anything still failing is retried with relaxed limits.
    Returns (targets used, results by id, notes by id, failures by id)."""
    keep = max(cfg.candidates_per_target + 20, 60)

    def run(ts, c=cfg):
        return {r["target"]: r for r in _run_candidates(ts, c, store, keep, verbose)} if ts else {}

    results: Dict[str, dict] = {}
    notes: Dict[str, str] = {}
    failed: Dict[str, str] = {}
    final: List[Target] = []

    def accept(t, r, note=None):
        results[t.id] = r
        final.append(t)
        if note:
            notes[t.id] = note

    r1 = run(targets)
    stage2: List[Target] = []
    parent: Dict[str, Target] = {}
    for t in targets:
        if r1[t.id]["pairs"]:
            accept(t, r1[t.id])
        elif len(t.sites) > 1 or t.span > 60:
            hard = max(1, -(-t.span // 2))   # two halves, not three thirds
            subs = group_sites([Site(t.seqid, max(s.start, t.start), min(s.end, t.end), s.name) for s in t.sites],
                               span_comfort=hard, span_hard=hard, min_gap=1)
            for k, sub in enumerate(subs, 1):
                sub.id = f"{t.id}.{k}"
                sub.origins = t.origins
                parent[sub.id] = t
            _log(f"{t.id}: {r1[t.id]['fail']}; splitting into {len(subs)} targets", verbose)
            stage2.extend(subs)
        else:
            stage2.append(t)
    r2 = run([s for s in stage2 if s.id in parent])
    stage3: List[Target] = []
    for s in stage2:
        if s.id in r2 and r2[s.id]["pairs"]:
            accept(s, r2[s.id], f"split from {parent[s.id].id} ({r1[parent[s.id].id]['fail']})")
        else:
            stage3.append(s)
    rc = _relaxed(cfg)
    r3 = run(stage3, rc)
    for s in stage3:
        why = (r2.get(s.id) or r1.get(s.id) or r3[s.id])["fail"]
        if r3[s.id]["pairs"]:
            accept(s, r3[s.id], f"relaxed limits: Tm window +/-{rc.tm_window:g} C, runs <= {rc.max_poly}, "
                                f"primer GC >= {rc.primer_gc_min:.2f}, site structure and hairpin limits loosened "
                                f"({why})")
        else:
            failed[s.id] = (f"{r3[s.id]['fail']} (also with relaxed limits: Tm +/-{rc.tm_window:g} C, runs <= "
                            f"{rc.max_poly}, GC >= {rc.primer_gc_min:.2f}; strict stage: {why})")
    final.sort(key=lambda t: (t.seqid, t.start))
    return final, results, notes, failed


# ----------------------------------------------------------------------------
# step 4: specificity
# ----------------------------------------------------------------------------
def _specificity(pairs_by_target, cfg, genome_fasta, verbose):
    """Returns (seq index, per-seq intended sites, scored sites, products, spec params, genome)."""
    from . import specificity as spec

    seqs: List[str] = []
    index: Dict[str, int] = {}
    intended: Dict[int, set] = {}
    g = spec.load_genome(genome_fasta)
    cidx = {n: i for i, n in enumerate(g.names)}
    for pairs in pairs_by_target.values():
        for pr in pairs:
            for p in (pr.fwd, pr.rev):
                if p.seq not in index:
                    index[p.seq] = len(seqs)
                    seqs.append(p.seq)
                strand = 1 if p.side == "F" else -1
                intended.setdefault(index[p.seq], set()).add((cidx[p.seqid], strand, p.start - 1))
    params = spec.Params(ta_c=cfg.ta, mv_conc=cfg.mono_mM, dv_conc=cfg.mg_mM, dntp_conc=cfg.dntp_mM,
                         dna_conc=cfg.primer_nM, polymerase="proofreading" if cfg.proofreading else "taq",
                         max_size=cfg.offtarget_max_product)
    t0 = time.time()
    sites = spec.find_sites_seed(seqs, g, params, cfg.threads)
    scored = spec.score_sites(sites, seqs, g, params, cfg.threads)
    prods = spec.enumerate_products(scored, seqs, g, params, names=[str(i) for i in range(len(seqs))])
    _log(f"specificity: {len(seqs)} primers, {len(sites)} candidate sites, {len(scored)} scored, "
         f"{len(prods)} products ({time.time() - t0:.0f} s)", verbose)
    return seqs, index, intended, scored, prods, params, g


def classify_product(d, i, j, pr1, pr2, seqs, contig_name):
    """Kind of a product between primer index i (+ site) and j (- site), held by
    candidate pairs pr1 (contains seqs[i]) and pr2 (contains seqs[j]).

    'designed'  the amplicon of pr1 itself (pr1 is pr2)
    'repeat'    the designed pair amplifying a second copy of the same length (repeated locus)
    'outer'     F of one designed amplicon with R of another, both at their designed sites
    'other'     anything else (off-target sites, wrong pairings)
    """
    def at_f(pr):
        return (pr.fwd.seq == seqs[i] and pr.fwd.seqid == d["contig"] and abs(d["start"] - (pr.fwd.start - 1)) <= 5)

    def at_r(pr):
        return (pr.rev.seq == seqs[j] and pr.rev.seqid == d["contig"] and abs(d["end"] - pr.rev.end) <= 5)

    if pr1 is pr2:
        if at_f(pr1) and at_r(pr1):
            return "designed"
        if {seqs[i], seqs[j]} == {pr1.fwd.seq, pr1.rev.seq} and abs(d["size"] - pr1.length) <= 5:
            return "repeat"
        return "other"
    if at_f(pr1) and at_r(pr2):
        return "outer"
    return "other"


def _is_intended_locus(site, intended_set, tol=20):
    for ci, strand, start in intended_set:
        if site["contig"] == ci and site["strand"] == strand and abs(site["start"] - start) <= tol:
            return True
    return False


# ----------------------------------------------------------------------------
# step 5: pairwise dimer costs
# ----------------------------------------------------------------------------
def _dimer_chunk(args):
    tasks, cfg = args
    return [cand.thal_pair(a, b, cfg) for a, b in tasks]


def _pairwise_dimers(keys_pairs: List[Tuple[Tuple[str, str], Tuple[str, str]]], cfg: Config) -> Dict:
    """{(keyA, keyB): (d3, dany)} for (seq, side) keys, tails included."""
    def full(k):
        return (cfg.tail_fwd if k[1] == "F" else cfg.tail_rev) + k[0]

    tasks = [(full(a), full(b)) for a, b in keys_pairs]
    if not tasks:
        return {}
    if cfg.threads <= 1 or len(tasks) < 2000:
        res = _dimer_chunk((tasks, cfg))
    else:
        n = max(500, math.ceil(len(tasks) / (cfg.threads * 4)))
        chunks = [(tasks[i:i + n], cfg) for i in range(0, len(tasks), n)]
        res = []
        with ProcessPoolExecutor(cfg.threads, mp_context=_ctx()) as ex:
            for r in ex.map(_dimer_chunk, chunks):
                res.extend(r)
    return {kp: r for kp, r in zip(keys_pairs, res)}


# ----------------------------------------------------------------------------
# main entry
# ----------------------------------------------------------------------------
def design_pool(sites: Sequence[Site], genome_meta: dict, cfg: Config, verbose: bool = True) -> Design:
    """Design a pool. For Taq-type presets the recommended Ta can only be computed
    after design; when it differs from the design Ta by more than 1 C the whole
    design is re-scored once at the recommended Ta, because binding uniformity
    at one temperature does not carry over to another (SantaLucia 2007, Fig 5)."""
    ta_given = cfg.ta is not None
    work = copy.deepcopy(cfg)          # the caller's Config is never modified
    passes, earlier = 0, 0.0
    tried = set()
    while True:
        passes += 1
        d = _design_once(list(sites), genome_meta, work, verbose)
        rec = float(d.ta["recommended"])
        if ta_given or work.polymerase == "q5" or abs(rec - work.ta) <= 1.0 or passes >= 3 or rec in tried:
            break
        tried.add(work.ta)
        _log(f"recommended Ta {rec:g} C differs from the scoring Ta {work.ta:.1f} C; re-scoring at {rec:g} C",
             verbose)
        earlier += d.timings.get("total", 0.0)
        work = copy.deepcopy(work)
        work.ta = rec
    d.timings["earlier_passes"] = earlier
    d.passes = passes
    return d


def _design_once(sites: Sequence[Site], genome_meta: dict, cfg: Config, verbose: bool = True) -> Design:
    timings: Dict[str, float] = {}
    t_all = time.time()
    store = SeqStore(genome_meta["fasta"], genome_meta.get("fai"))
    # bare-integer advice is a note, not a design warning (older cached metadata kept it in warnings)
    warnings: List[str] = [w for w in genome_meta.get("warnings", []) if "bare integer positions" not in w]
    gc = (genome_meta.get("gc_percent") or 50.0) / 100.0
    cfg.resolve(gc)
    _log(f"genome GC {100 * gc:.1f}%: target Tm {cfg.tm:.1f} C ({cfg.tm_scale} scale) +/- {cfg.tm_window}, primers "
         f"{cfg.primer_min}-{cfg.primer_max} nt, amplicons {cfg.amp_min}-{cfg.amp_max} bp, design Ta {cfg.ta:.1f} C, "
         f"{cfg.polymerase}, {cfg.primer_nM:g} nM each, {cfg.mono_mM:g} mM mono, {cfg.mg_mM:g} mM Mg", verbose)

    # 0. requested bases must be real sequence
    keep_sites, dropped_sites = [], []
    for s_ in sites:
        ref = store.region(s_.seqid, s_.start, s_.end)[0]
        bad = sum(1 for c in ref if c not in "ACGT")
        if bad == len(ref):
            warnings.append(f"position {s_.label()} is N (unknown) in the reference and cannot be genotyped; "
                            f"left out")
            dropped_sites.append(s_)
            continue
        if bad:
            warnings.append(f"{bad} of {len(ref)} bases of {s_.label()} are N (unknown) in the reference")
        keep_sites.append(s_)
    if not keep_sites:
        raise DesignError("every requested position is an unknown (N) base in the reference")
    sites = keep_sites

    # 1. targets
    # comfortable span from a realistic primer length, hard limit from the shortest
    # allowed primer (failing merged targets are split later anyway)
    comfort = span_limits(cfg.amp_opt, cfg.amp_max, primer_len=min(cfg.primer_max, cfg.primer_min + 4),
                          margin=cfg.margin, flex=20)[0]
    _, hard, gap = span_limits(cfg.amp_opt, cfg.amp_max, primer_len=cfg.primer_min, margin=cfg.margin, flex=20)
    hard = max(hard, comfort)
    targets = group_sites(sites, comfort, hard, gap)
    _log(f"{len(sites)} positions -> {len(targets)} targets", verbose)
    for t in targets:
        L = store.lengths[t.seqid]
        if t.start <= cfg.amp_max or t.end > L - cfg.amp_max:
            warnings.append(f"target {t.id} lies within {cfg.amp_max} bp of an end of {t.seqid}; amplicons cannot "
                            f"cross the end, and if {t.seqid} is circular, products across the origin are not "
                            f"checked")

    # 2. candidates
    t0 = time.time()
    targets, results, relaxed, failed = _candidates_with_fallback(targets, cfg, store, verbose)
    timings["candidates"] = time.time() - t0
    _log(f"candidates for {len(targets)} targets ({timings['candidates']:.0f} s); failed: {len(failed)}", verbose)
    for tid, why in failed.items():
        warnings.append(f"target {tid} has no primer pair: {why}")
    if not targets:
        raise DesignError("no target produced any primer pair: " + "; ".join(f"{k}: {v}" for k, v in failed.items()))

    # 3. uniform binding across the pool
    pairs_by_target: Dict[str, List[cand.Pair]] = {t.id: results[t.id]["pairs"] for t in targets}
    dgs = [p.dg_eff for prs in pairs_by_target.values() for pr in prs for p in (pr.fwd, pr.rev)]
    dg_med = statistics.median(dgs)
    for prs in pairs_by_target.values():
        for pr in prs:
            u = ((pr.fwd.dg_eff - dg_med) / 0.5) ** 2 + ((pr.rev.dg_eff - dg_med) / 0.5) ** 2
            pr.parts["dg_uniformity"] = u
            pr.cost += u
        prs.sort(key=lambda p: p.cost)
        del prs[cfg.candidates_per_target:]

    # 4. specificity
    spec_data = None
    unary_extra: Dict[Tuple[str, int], float] = {}
    pair_products: Dict[Tuple[str, int, str, int], float] = {}
    offtarget_by_seq: Dict[int, List[dict]] = {}
    if cfg.specificity:
        t0 = time.time()
        seqs, index, intended, scored, prods, sparams, sgenome = _specificity(pairs_by_target, cfg,
                                                                              genome_meta["fasta"], verbose)
        timings["specificity"] = time.time() - t0
        spec_data = (seqs, index, intended, scored, sparams, sgenome)
        for s in scored:
            if s["primable"] and not _is_intended_locus(s, intended.get(s["primer"], ())):
                offtarget_by_seq.setdefault(s["primer"], []).append(s)
        tmap_t = {t.id: t for t in targets}
        # which candidate pairs contain each primer sequence
        holders: Dict[int, List[Tuple[str, int]]] = {}
        for tid, prs in pairs_by_target.items():
            for a, pr in enumerate(prs):
                for p in (pr.fwd, pr.rev):
                    holders.setdefault(index[p.seq], []).append((tid, a))
                off = offtarget_by_seq.get(index[pr.fwd.seq], []) + offtarget_by_seq.get(index[pr.rev.seq], [])
                if off:
                    c = min(2.0 * len(off) + 10.0 * sum(1 for s in off if s["score"] >= 0.5), 30.0)
                    unary_extra[(tid, a)] = unary_extra.get((tid, a), 0.0) + c
                    pr.parts["offtarget_sites"] = c
        for d in prods:
            i, j = int(d["fwd"]), int(d["rev"])
            likely = d["call"] == "likely"
            for (t1, a1) in holders.get(i, []):
                for (t2, a2) in holders.get(j, []):
                    if t1 == t2 and a1 != a2:
                        continue   # never chosen together
                    pr1, pr2 = pairs_by_target[t1][a1], pairs_by_target[t2][a2]
                    kind = classify_product(d, i, j, pr1, pr2, seqs, d["contig"])
                    if kind == "designed":
                        continue
                    if kind == "repeat":
                        # second copy of the locus: avoid when possible, but no primer choice
                        # inside a multi-copy locus (rRNA, IS elements) can remove it
                        pen = 50.0 if likely else 1.0
                    elif kind == "outer":
                        # tiles of one requested range: outer products are the normal ARTIC pattern
                        same_range = bool(tmap_t[t1].origins & tmap_t[t2].origins)
                        short = d["size"] < cfg.outer_min and not same_range
                        pen = (1000.0 if short else 2.0) if likely else 0.5
                    else:
                        pen = 1000.0 if likely else 3.0
                    if t1 == t2:
                        unary_extra[(t1, a1)] = unary_extra.get((t1, a1), 0.0) + pen
                        pr1.parts["products"] = pr1.parts.get("products", 0.0) + pen
                    else:
                        k = (t1, a1, t2, a2)
                        pair_products[k] = pair_products.get(k, 0.0) + pen

    # 5. pools and pairwise costs
    t0 = time.time()
    graph = conflict_graph(targets, cfg)
    colour = dsatur(graph)
    pools = {targets[i].id: colour[i] + 1 for i in range(len(targets))}
    n_pools = max(pools.values())
    _log(f"{n_pools} pool(s): " + ", ".join(
        f"pool {p}: {sum(1 for v in pools.values() if v == p)} amplicons" for p in range(1, n_pools + 1)), verbose)
    tindex = {t.id: i for i, t in enumerate(targets)}
    need = set()
    for t, u in itertools.combinations(targets, 2):
        if pools[t.id] != pools[u.id]:
            continue
        kt = {(p.seq, p.side) for pr in pairs_by_target[t.id] for p in (pr.fwd, pr.rev)}
        ku = {(p.seq, p.side) for pr in pairs_by_target[u.id] for p in (pr.fwd, pr.rev)}
        for a in kt:
            for b in ku:
                need.add((a, b) if a <= b else (b, a))
    dimer = _pairwise_dimers(sorted(need), cfg)
    timings["pairwise_dimers"] = time.time() - t0
    _log(f"{len(dimer)} cross-target primer pairs scored for dimers ({timings['pairwise_dimers']:.0f} s)", verbose)

    def g(a, b):
        k = (a, b) if a <= b else (b, a)
        d3, dany = dimer[k]
        return cand.dimer_cost(d3, dany, cfg)

    unary = []
    for t in targets:
        unary.append(np.array([pr.cost + unary_extra.get((t.id, a), 0.0)
                               for a, pr in enumerate(pairs_by_target[t.id])]))
    pairwise = {}
    for t, u in itertools.combinations(targets, 2):
        if pools[t.id] != pools[u.id]:
            continue
        P, Q = pairs_by_target[t.id], pairs_by_target[u.id]
        m = np.zeros((len(P), len(Q)))
        for a, pa in enumerate(P):
            ka = [(pa.fwd.seq, "F"), (pa.rev.seq, "R")]
            for b, pb in enumerate(Q):
                kb = [(pb.fwd.seq, "F"), (pb.rev.seq, "R")]
                m[a, b] = sum(g(x, y) for x in ka for y in kb) + pair_products.get((t.id, a, u.id, b), 0.0) \
                    + pair_products.get((u.id, b, t.id, a), 0.0)
        pairwise[(tindex[t.id], tindex[u.id])] = m

    # 6. optimise
    t0 = time.time()
    prob = Problem(unary, pairwise)
    x, energy, method = prob.solve(restarts=cfg.sa_restarts, seed=cfg.seed)
    def hard_terms(x, d3cap):
        dims, prods_hard = [], 0
        for (ti, ui) in pairwise:
            pa = pairs_by_target[targets[ti].id][x[ti]]
            pb = pairs_by_target[targets[ui].id][x[ui]]
            for a in ((pa.fwd.seq, "F"), (pa.rev.seq, "R")):
                for b in ((pb.fwd.seq, "F"), (pb.rev.seq, "R")):
                    d3, dany = dimer[(a, b) if a <= b else (b, a)]
                    if d3 <= d3cap or dany <= cfg.dany_hard:
                        dims.append((a[0], b[0], d3, dany))
            if pair_products.get((targets[ti].id, x[ti], targets[ui].id, x[ui]), 0) \
                    + pair_products.get((targets[ui].id, x[ui], targets[ti].id, x[ti]), 0) >= 1000:
                prods_hard += 1
        for i, t in enumerate(targets):
            if unary_extra.get((t.id, x[i]), 0) >= 1000:
                prods_hard += 1
        return dims, prods_hard

    dims, _ = hard_terms(x, cfg.d3_hard)
    if any(d3 <= cfg.d3_hard for _, _, d3, _ in dims) and cfg.d3_hard > -7.0:
        cfg2 = copy.copy(cfg)
        cfg2.d3_hard = -7.0

        def g2(a, b):
            k = (a, b) if a <= b else (b, a)
            return cand.dimer_cost(*dimer[k], cfg2)
        pw2 = {}
        for (ti, ui), m in pairwise.items():
            P, Q = pairs_by_target[targets[ti].id], pairs_by_target[targets[ui].id]
            m2 = np.zeros_like(m)
            for a, pa in enumerate(P):
                for b, pb in enumerate(Q):
                    m2[a, b] = sum(g2(xx, yy) for xx in [(pa.fwd.seq, "F"), (pa.rev.seq, "R")]
                                   for yy in [(pb.fwd.seq, "F"), (pb.rev.seq, "R")]) \
                        + pair_products.get((targets[ti].id, a, targets[ui].id, b), 0.0) \
                        + pair_products.get((targets[ui].id, b, targets[ti].id, a), 0.0)
            pw2[(ti, ui)] = m2
        x2, e2, method2 = Problem(unary, pw2).solve(restarts=cfg.sa_restarts, seed=cfg.seed)
        if e2 < energy:
            n_mid = sum(1 for _, _, d3, _ in hard_terms(x2, -6.0)[0] if -7.0 < d3 <= -6.0)
            x, energy, method = x2, e2, method2
            pairwise = pw2
            warnings.append(f"no choice met the -6 kcal/mol 3'-dimer cap; relaxed it to -7 "
                            f"({n_mid} primer pair(s) between -6 and -7 kcal/mol)")
    timings["optimise"] = time.time() - t0
    _log(f"optimised ({method}): energy {energy:.1f} ({timings['optimise']:.1f} s)", verbose)
    dims, nprod = hard_terms(x, cfg.d3_hard if energy < 1000 or cfg.d3_hard <= -7 else -7.0)
    if dims:
        warnings.append(f"the chosen pool still has {len(dims)} primer pair(s) past the dimer limits "
                        f"(see the dimers table)")
    if nprod:
        warnings.append(f"the chosen pool still has {nprod} likely unintended product(s) that no primer "
                        f"choice avoided (see the products table)")
    chosen = {t.id: pairs_by_target[t.id][x[i]] for i, t in enumerate(targets)}

    # drop amplicons that coverage does not need (split or relaxed targets can overlap)
    from .targets import coverage_gaps
    order = sorted(chosen, key=lambda tid: -chosen[tid].cost)
    for tid in order:
        rest = [(pr.fwd.seqid, *pr.insert) for k, pr in chosen.items() if k != tid]
        if len(chosen) > 1 and not coverage_gaps(sites, rest):
            _log(f"{tid} is not needed for coverage; dropped", verbose)
            del chosen[tid]
    targets = [t for t in targets if t.id in chosen]
    used = sorted({pools[t.id] for t in targets})
    renum = {p: k + 1 for k, p in enumerate(used)}
    pools = {t.id: renum[pools[t.id]] for t in targets}

    # 7. final checks
    from . import finalize

    t_fin = time.time()
    fin = finalize.finish(targets, chosen, pools, sites, cfg, store, spec_data, dimer, offtarget_by_seq, verbose)
    timings["final_checks"] = time.time() - t_fin
    balance_notes: List[str] = []
    if cfg.balance:
        from .balance import balance

        t_bal = time.time()
        balance_notes = balance(chosen, pools, fin["ta"], cfg, store, fin["names"])
        if any("not extended" not in n_ for n_ in balance_notes):
            _log(f"binding balance: {sum('not extended' not in n_ for n_ in balance_notes)} primer(s) extended "
                 f"at the 5' end", verbose)
            spec2, off2 = spec_data, offtarget_by_seq
            if cfg.specificity:
                single = {tid: [pr] for tid, pr in chosen.items()}
                seqs2, index2, intended2, scored2, _, sp2, g2 = _specificity(single, cfg, genome_meta["fasta"], verbose)
                spec2 = (seqs2, index2, intended2, scored2, sp2, g2)
                off2 = {}
                for s_ in scored2:
                    if s_["primable"] and not _is_intended_locus(s_, intended2.get(s_["primer"], ())):
                        off2.setdefault(s_["primer"], []).append(s_)
            fin = finalize.finish(targets, chosen, pools, sites, cfg, store, spec2, dimer, off2, verbose)
        timings["balance"] = time.time() - t_bal
    warnings.extend(fin["warnings"])
    if abs(fin["ta"]["recommended"] - cfg.ta) > 1.0:
        warnings.append(f"primers, structure and off-target sites were scored at {cfg.ta:.1f} C but the recommended "
                        f"Ta is {fin['ta']['recommended']:g} C; re-run with --ta {fin['ta']['recommended']:g} to "
                        f"score at that temperature")
    timings["total"] = time.time() - t_all
    _log(f"done in {timings['total']:.0f} s", verbose)
    return Design(cfg=cfg, sites=list(sites), targets=targets, pools=pools, chosen=chosen, failed=failed,
                  relaxed=relaxed, coverage_problems=fin["coverage"], products=fin["products"],
                  dimers=fin["dimers"], ta=fin["ta"], amplicon_reports=fin["amplicons"],
                  offtarget_sites=fin["offtarget"], energy=energy, method=method, timings=timings,
                  warnings=warnings, names=fin["names"], balance_notes=balance_notes, dropped_sites=dropped_sites,
                  genome={k: genome_meta.get(k) for k in ("query", "accession", "organism", "assembly_name",
                                                         "fasta", "fai", "gc_percent", "total_length",
                                                         "n_sequences")})
