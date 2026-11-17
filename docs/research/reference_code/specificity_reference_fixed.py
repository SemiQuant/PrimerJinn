#!/usr/bin/env python3
"""
specificity_reference.py: genome-wide specificity check for a multiplex primer pool.

Reference implementation (research prototype, tested on M. tuberculosis H37Rv, 4.4 Mb).

Pipeline
  1. find_sites_seed()   3'-anchored k-mer seed table (three 8-mers at 0, 5, 10 nt from the 3' end, 1 mismatch
                         each) + one numba pass over the genome; each hit gets an optimistic nearest-neighbour
                         dG at Ta (SantaLucia 2004 stacks + thal salt term) and is kept if dG_est <= est_max
     find_sites_blast()  alternative: blastn-short, hits end-filled to full ungapped footprints
     find_sites_brute()  exhaustive Hamming scan (ground truth for tests, bacteria only)
  2. score_sites()       primer3-py thal (calc_heterodimer = ANY, calc_end_stability = END1)
                         at the annealing temperature -> dG, fraction bound, extension factor.
                         Runs in a fork process pool: primer3-py 2.0.3 thal crashes (SIGBUS) under threads.
  3. enumerate_products() shifted duplicates merged per locus, then every '+' site x every '-' site
                         downstream within max_size,
                         from ANY primers in the pool (cross-pair products included)

Coordinates are 0-based. A site is the ungapped footprint of the primer on the genome.
  strand=+1: genome[start:start+L] resembles the primer; 3' end at start+L-1; extends rightwards
  strand=-1: genome[start:start+L] resembles revcomp(primer); 3' end at start; extends leftwards

Units: temperatures in degrees Celsius, dG in kcal/mol, concentrations mM (salts) and nM (oligo).
"""
from __future__ import annotations

import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field

import numba as nb
import numpy as np
import primer3

R_KCAL = 1.9872e-3  # kcal / (mol K)
_ENC = np.full(256, 4, dtype=np.uint8)
for _i, _c in enumerate("ACGT"):
    _ENC[ord(_c)] = _i
    _ENC[ord(_c.lower())] = _i
_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def encode(s: str) -> np.ndarray:
    return _ENC[np.frombuffer(s.encode("ascii"), dtype=np.uint8)]


# --------------------------------------------------------------------------- genome
@dataclass
class Genome:
    names: list
    seqs: list            # python str, upper case
    codes: list           # np.uint8 arrays, 0..3 = ACGT, 4 = other

    @property
    def size(self) -> int:
        return sum(len(s) for s in self.seqs)


def load_genome(fasta: str) -> Genome:
    names, seqs, cur = [], [], []
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur).upper())
                    cur = []
                names.append(line[1:].split()[0])
            else:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur).upper())
    return Genome(names, seqs, [encode(s) for s in seqs])


# --------------------------------------------------------------------------- parameters
@dataclass
class Params:
    # reaction conditions (primer3-py thal arguments)
    ta_c: float = 60.0            # annealing temperature used for dG / fraction bound
    mv_conc: float = 50.0         # mM monovalent
    dv_conc: float = 2.0          # mM divalent
    dntp_conc: float = 0.8        # mM dNTP
    dna_conc: float = 250.0       # nM, per primer
    conc_divisor: float = 1.0     # 1 = primer in excess (PCR); 4 = primer3 Tm convention
    polymerase: str = "taq"       # "taq" (no 3'->5' exo) or "proofreading" (Q5, Phusion, KAPA HiFi)
    # candidate search
    seeds: tuple = ((0, 8, 1, 0), (5, 8, 1, 0), (10, 8, 1, 0))  # (offset from 3' end, k, max mismatches, protected 3' bases)
    max_mm: int = 99              # optional hard cap on ungapped mismatches (off by default)
    max_mm_frac: float = 1.0      # optional cap on the mismatch fraction (off by default)
    est_max: float = -5.0         # keep a candidate only if its optimistic numba dG estimate at Ta <= this (kcal/mol)
    keep_min_score: float = 1e-6  # score_sites drops sites scoring below this (saves memory on big runs)
    flank: int = 2                # template bases added each side for dangling ends
    # calls
    site_min_score: float = 0.01      # 'primable' site
    product_likely: float = 0.01      # 'likely' product (score = score1 * score2)
    product_possible: float = 1e-4    # 'possible' product
    min_size: int = 40
    locus_merge_nt: int = 20       # sites of one primer/strand whose 3' ends lie within this distance are one locus
    max_size: int = 2000
    circular: bool = False         # True only for complete circular replicons (set per run)
    # extension model for non-proofreading polymerases (heuristic factors, see report section 5)
    taq_term_severe: float = 0.01     # 3'-terminal A:A, A:G, G:A, G:G, C:C (primer:template)
    taq_term_mild: float = 0.3        # any other 3'-terminal mismatch
    taq_pos_factor: tuple = (0.3, 0.4, 0.8, 0.9)   # single mismatch at positions 2, 3, 4, 5 from the 3' end
    taq_two_in_last5: float = 0.01    # cap when two or more mismatches lie within the last 5 nt



# --------------------------------------------------------------------------- fast dG prefilter (numba)
# SantaLucia & Hicks 2004 Watson-Crick nearest-neighbour stacks (dH kcal/mol, dS cal/K/mol), top strand 5'XY3'.
_NN = {"AA": (-7.6, -21.3), "TT": (-7.6, -21.3), "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
       "CA": (-8.5, -22.7), "TG": (-8.5, -22.7), "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
       "CT": (-7.8, -21.0), "AG": (-7.8, -21.0), "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
       "CG": (-10.6, -27.2), "GC": (-9.8, -24.4), "GG": (-8.0, -19.9), "CC": (-8.0, -19.9)}
NN_DH = np.zeros(16); NN_DS = np.zeros(16)
for _k, (_h, _s) in _NN.items():
    NN_DH["ACGT".index(_k[0]) * 4 + "ACGT".index(_k[1])] = _h
    NN_DS["ACGT".index(_k[0]) * 4 + "ACGT".index(_k[1])] = _s


def stack_dg_table(params) -> np.ndarray:
    """dG (kcal/mol) of each matched stack at Ta including the primer3 thal salt term per stack:
    0.368 * ln(Na_eq / 1000) cal/K/mol with Na_eq = mv + 120 sqrt(max(0, dv - dntp)) (thal.c saltCorrectS)."""
    T = params.ta_c + 273.15
    na = params.mv_conc + 120.0 * math.sqrt(max(0.0, params.dv_conc - params.dntp_conc))
    salt = 0.368 * math.log(na / 1000.0)
    return (NN_DH - T * (NN_DS + salt) / 1000.0).astype(np.float64)


def init_dg(params) -> float:
    return 0.2 - (params.ta_c + 273.15) * (-5.7) / 1000.0      # thal dplx_init_H = 200 cal, dplx_init_S = -5.7


@nb.njit(cache=True, nogil=True)
def _mm_and_est(g, s, pc, L, stack_dg, init):
    """Ungapped mismatch count and an optimistic duplex dG: initiation + every stack whose two pairs both
    match (mismatch loops cost nothing). Used only to discard hopeless candidates before thal."""
    mm = 0
    est = init
    prev_ok = False
    for i in range(L):
        ok = g[s + i] == pc[i]
        if not ok:
            mm += 1
        elif prev_ok:
            est += stack_dg[pc[i - 1] * 4 + pc[i]]
        prev_ok = ok
    return mm, est

# --------------------------------------------------------------------------- seed search (numba)
def _seed_variants(code_list, k, max_mm, protected_idx):
    """All sequences within max_mm substitutions of a k-mer (list of ints 0..3).
    protected_idx: set of indices (0-based within the k-mer) that may not be substituted."""
    out = {tuple(code_list)}
    frontier = {tuple(code_list)}
    for _ in range(max_mm):
        nxt = set()
        for v in frontier:
            for i in range(k):
                if i in protected_idx:
                    continue
                for b in range(4):
                    if b != v[i]:
                        w = v[:i] + (b,) + v[i + 1:]
                        if w not in out:
                            nxt.add(w)
        out |= nxt
        frontier = nxt
    return out


def build_seed_table(primers, seeds):
    """Return dict k -> (start_index[4^k+1], entry_primer, entry_strand, entry_shift).
    entry_shift = (site start) - (genome position of the first base of the seed)."""
    by_k = {}
    for pi, p in enumerate(primers):
        L = len(p)
        fw = [int(x) for x in encode(p)]
        rc = [int(x) for x in encode(revcomp(p))]
        for (off, k, mm, prot) in seeds:
            if off + k > L:
                continue
            a = L - off - k                       # seed start within the primer (5'->3')
            sub = fw[a:a + k]
            if 4 in sub:
                continue
            # protected: the last `prot` bases of the primer (only relevant when off < prot)
            prot_idx = {i for i in range(k) if (L - 1 - (a + i)) < prot}
            ent = by_k.setdefault(k, [])
            for v in _seed_variants(sub, k, mm, prot_idx):
                code = 0
                for b in v:
                    code = (code << 2) | b
                ent.append((code, pi, 1, -a))
                # reverse-complement seed for '-' sites: in rc(primer) the seed starts at `off`
                rcode = 0
                for b in reversed(v):
                    rcode = (rcode << 2) | (3 - b)
                ent.append((rcode, pi, -1, -off))
    tables = {}
    for k, ent in by_k.items():
        arr = np.array(sorted(set(ent)), dtype=np.int64).reshape(-1, 4)
        counts = np.bincount(arr[:, 0], minlength=4 ** k)
        start = np.zeros(4 ** k + 1, dtype=np.int64)
        np.cumsum(counts, out=start[1:])
        tables[k] = (start, arr[:, 1].astype(np.int32), arr[:, 2].astype(np.int8), arr[:, 3].astype(np.int32))
    return tables


@nb.njit(cache=True, nogil=True)
def _count_mm(g, s, pc, L, max_mm):
    """Ungapped mismatches of primer codes pc against g[s:s+L]; returns max_mm+1 on early exit."""
    mm = 0
    for i in range(L):
        if g[s + i] != pc[i]:
            mm += 1
            if mm > max_mm:
                return mm
    return mm


@nb.njit(cache=True, nogil=True)
def _scan_seeds(g, lo, hi, k, start, e_primer, e_strand, e_shift, fw, rc, plen, max_mm, max_mm_frac,
                stack_dg, init, est_max):
    """Scan genome codes g[lo:hi] (seed start positions) for seed hits; keep a hit if its ungapped
    mismatches are within the caps AND its optimistic dG estimate is <= est_max (kcal/mol)."""
    n = g.shape[0]
    cap = 1024
    out = np.empty((cap, 5), dtype=np.int64)
    cnt = 0
    mask = (1 << (2 * k)) - 1
    code = 0
    valid = 0
    end = min(hi + k - 1, n)
    for pos in range(lo, end):
        b = g[pos]
        if b > 3:
            valid = 0
            code = 0
            continue
        code = ((code << 2) | b) & mask
        valid += 1
        if valid < k:
            continue
        sp = pos - k + 1                       # seed start
        for e in range(start[code], start[code + 1]):
            pi = e_primer[e]
            L = plen[pi]
            s = sp + e_shift[e]
            if s < 0 or s + L > n:
                continue
            lim = min(max_mm, int(max_mm_frac * L))
            if e_strand[e] == 1:
                mm, est = _mm_and_est(g, s, fw[pi], L, stack_dg, init)
            else:
                mm, est = _mm_and_est(g, s, rc[pi], L, stack_dg, init)
            if mm <= lim and est <= est_max:
                if cnt == cap:
                    cap *= 2
                    new = np.empty((cap, 5), dtype=np.int64)
                    new[:cnt] = out[:cnt]
                    out = new
                out[cnt, 0] = pi
                out[cnt, 1] = e_strand[e]
                out[cnt, 2] = s
                out[cnt, 3] = mm
                out[cnt, 4] = int(round(est * 100))
                cnt += 1
    return out[:cnt]


@nb.njit(cache=True, nogil=True)
def _scan_brute(g, lo, hi, fw, rc, plen, max_mm):
    n = g.shape[0]
    cap = 4096
    out = np.empty((cap, 4), dtype=np.int64)
    cnt = 0
    P = plen.shape[0]
    for s in range(lo, hi):
        for pi in range(P):
            L = plen[pi]
            if s + L > n:
                continue
            for strand in (1, -1):
                if strand == 1:
                    mm = _count_mm(g, s, fw[pi], L, max_mm)
                else:
                    mm = _count_mm(g, s, rc[pi], L, max_mm)
                if mm <= max_mm:
                    if cnt == cap:
                        cap *= 2
                        new = np.empty((cap, 4), dtype=np.int64)
                        new[:cnt] = out[:cnt]
                        out = new
                    out[cnt, 0] = pi
                    out[cnt, 1] = strand
                    out[cnt, 2] = s
                    out[cnt, 3] = mm
                    cnt += 1
    return out[:cnt]


def _primer_arrays(primers):
    P = len(primers)
    Lmax = max(len(p) for p in primers)
    fw = np.full((P, Lmax), 9, dtype=np.uint8)
    rc = np.full((P, Lmax), 9, dtype=np.uint8)
    plen = np.zeros(P, dtype=np.int32)
    for i, p in enumerate(primers):
        L = len(p)
        plen[i] = L
        fw[i, :L] = encode(p)
        rc[i, :L] = encode(revcomp(p))
    return fw, rc, plen


def _chunks(n, nchunk):
    step = max(1, math.ceil(n / nchunk))
    return [(a, min(n, a + step)) for a in range(0, n, step)]


def _dedupe(arr):
    """rows: contig, primer, strand, start, mm -> unique sorted."""
    if len(arr) == 0:
        return np.empty((0, 5), dtype=np.int64)
    return np.unique(arr, axis=0)


def find_sites_seed(primers, genome: Genome, params: Params, threads: int = 0):
    """3'-anchored seed-table scan. Returns int64 array [contig, primer, strand, start, mm]."""
    threads = threads or os.cpu_count() or 1
    fw, rc, plen = _primer_arrays(primers)
    tables = build_seed_table(primers, params.seeds)
    sdg = stack_dg_table(params)
    ini = init_dg(params)
    res = []
    with ThreadPoolExecutor(threads) as ex:
        futs = []
        for ci, g in enumerate(genome.codes):
            for k, (start, ep, es, esh) in tables.items():
                for lo, hi in _chunks(len(g), threads * 2):
                    futs.append((ci, ex.submit(_scan_seeds, g, lo, hi, k, start, ep, es, esh, fw, rc, plen,
                                               params.max_mm, params.max_mm_frac, sdg, ini, params.est_max)))
        for ci, f in futs:
            a = f.result()[:, :4]
            if len(a):
                res.append(np.column_stack([np.full(len(a), ci, dtype=np.int64), a]))
    return _dedupe(np.vstack(res)) if res else np.empty((0, 5), dtype=np.int64)


def find_sites_brute(primers, genome: Genome, max_mm: int = 7, threads: int = 0):
    """Exhaustive ungapped scan (ground truth). O(genome x primers)."""
    threads = threads or os.cpu_count() or 1
    fw, rc, plen = _primer_arrays(primers)
    res = []
    with ThreadPoolExecutor(threads) as ex:
        futs = []
        for ci, g in enumerate(genome.codes):
            for lo, hi in _chunks(len(g), threads * 4):
                futs.append((ci, ex.submit(_scan_brute, g, lo, hi, fw, rc, plen, max_mm)))
        for ci, f in futs:
            a = f.result()
            if len(a):
                res.append(np.column_stack([np.full(len(a), ci, dtype=np.int64), a]))
    return _dedupe(np.vstack(res)) if res else np.empty((0, 5), dtype=np.int64)


# --------------------------------------------------------------------------- blast search
BLAST_DEFAULTS = dict(task="blastn-short", word_size=7, evalue=1000, reward=1, penalty=-1,
                      gapopen=2, gapextend=2, dust="no", soft_masking="false", strand="both",
                      max_target_seqs=100000, max_hsps=100000)


def find_sites_blast(primers, genome_fasta: str, genome: Genome, params: Params, threads: int = 0,
                     blast_opts: dict | None = None, db: str | None = None, workdir: str | None = None):
    """blastn-short search; every HSP is anchored on its 3' (query-end) side and end-filled to the
    full ungapped primer footprint (the ThermoAlign 'thermoalignment' idea), then verified by the
    same ungapped mismatch count as the seed method. Returns [contig, primer, strand, start, mm]."""
    if not shutil.which("blastn"):
        raise RuntimeError("blastn not found")
    threads = threads or os.cpu_count() or 1
    opts = dict(BLAST_DEFAULTS)
    opts.update(blast_opts or {})
    own = workdir is None
    workdir = workdir or tempfile.mkdtemp(prefix="spec_blast_")
    try:
        if db is None:
            db = os.path.join(workdir, "db")
            subprocess.run(["makeblastdb", "-in", genome_fasta, "-dbtype", "nucl", "-out", db],
                           check=True, capture_output=True)
        qf = os.path.join(workdir, "q.fa")
        with open(qf, "w") as fh:
            for i, p in enumerate(primers):
                fh.write(f">{i}\n{p}\n")
        cmd = ["blastn", "-query", qf, "-db", db, "-num_threads", str(threads),
               "-outfmt", "6 qseqid sseqid qstart qend sstart send"]
        for k, v in opts.items():
            cmd += [f"-{k}", str(v)]
        r = subprocess.run(cmd, check=True, capture_output=True, text=True)
    finally:
        if own:
            shutil.rmtree(workdir, ignore_errors=True)
    cidx = {n: i for i, n in enumerate(genome.names)}
    rows = []
    for line in r.stdout.splitlines():
        q, s, qs, qe, ss, se = line.split("\t")
        pi = int(q); qs = int(qs); qe = int(qe); ss = int(ss); se = int(se)
        L = len(primers[pi])
        if ss <= se:      # plus: primer 3' end lies at genome (1-based) se + (L - qe)
            start = (se + (L - qe)) - L           # 0-based footprint start
            strand = 1
        else:             # minus: subject coords decrease along the query; 3' end at se - (L - qe)
            start = (se - (L - qe)) - 1
            strand = -1
        rows.append((cidx[s], pi, strand, start))
    if not rows:
        return np.empty((0, 5), dtype=np.int64)
    rows = np.unique(np.array(rows, dtype=np.int64), axis=0)
    fw, rc, plen = _primer_arrays(primers)
    out = []
    for ci, pi, strand, start in rows:
        g = genome.codes[ci]
        L = int(plen[pi])
        if start < 0 or start + L > len(g):
            continue
        lim = min(params.max_mm, int(params.max_mm_frac * L))
        mm = _count_mm(g, start, fw[pi] if strand == 1 else rc[pi], L, lim)
        if mm <= lim:
            out.append((ci, pi, strand, start, mm))
    return _dedupe(np.array(out, dtype=np.int64)) if out else np.empty((0, 5), dtype=np.int64)


# --------------------------------------------------------------------------- thermodynamic scoring


def frac_bound(dg_kcal: float, temp_c: float, dna_conc_nm: float, divisor: float = 1.0) -> float:
    """Two-state fraction of template sites occupied by a primer at concentration C:
        theta = 1 / (1 + exp(dG / RT) / (C / divisor))
    divisor = 1: primer in large excess over template (the PCR case; SantaLucia 1998 notes that
                 C_T is replaced by C_A - C_B / 2 when strands are at unequal concentration).
    divisor = 4: the primer3 convention (equal strand concentrations); theta = 0.5 exactly at the
                 Tm that primer3 thal reports for the same duplex."""
    c = dna_conc_nm * 1e-9 / divisor
    x = dg_kcal / (R_KCAL * (temp_c + 273.15))
    if x > 500:
        return 0.0
    return 1.0 / (1.0 + math.exp(x) / c)


# (primer base, template base). Huang et al. 1992: A.A ~1e-6; G.A, A.G, G.G, C.C < 1e-6 relative efficiency
_SEVERE_TERMINAL = {("A", "A"), ("A", "G"), ("G", "A"), ("G", "G"), ("C", "C")}


def extension_factor(mm_pos, primer3_base, template_base_opposite, params: Params) -> float:
    """Heuristic polymerase factor (0..1) from the mismatch positions counted from the 3' end (1 = terminal).
    primer3_base: 3'-terminal primer base. template_base_opposite: template base facing it.
    Proofreading enzymes remove mismatched 3' bases (primer editing), so the factor is 1."""
    if params.polymerase == "proofreading":
        return 1.0
    last5 = [p for p in mm_pos if p <= 5]
    if not last5:
        return 1.0
    f = 1.0
    for p in last5:
        if p == 1:
            severe = (primer3_base, template_base_opposite) in _SEVERE_TERMINAL
            f *= params.taq_term_severe if severe else params.taq_term_mild
        else:
            f *= params.taq_pos_factor[p - 2]
    if len(last5) >= 2:
        f = min(f, params.taq_two_in_last5)
    return f


def _thal_chunk(tasks, kw):
    """Worker: tasks = [(primer, template, need_end)] -> [(dg_any, tm_any, dg_end)] (kcal/mol, C).
    Runs in a separate process: primer3-py 2.0.3 thal crashes (SIGBUS) when called from threads."""
    out = []
    for p, t, need_end in tasks:
        het = primer3.calc_heterodimer(p, t, **kw)
        dg_end = primer3.calc_end_stability(p, t, **kw).dg / 1000.0 if need_end else float("nan")
        out.append((het.dg / 1000.0, het.tm, dg_end))
    return out


def _run_thal(tasks, kw, threads):
    if threads <= 1 or len(tasks) < 4000:
        return _thal_chunk(tasks, kw)
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor
    n = max(1000, math.ceil(len(tasks) / (threads * 4)))
    chunks = [tasks[i:i + n] for i in range(0, len(tasks), n)]
    res = []
    # 'fork' avoids re-importing the caller's __main__ (spawn does, which breaks unguarded scripts);
    # workers only call primer3, so forking a process that has numba / thread pools is safe here.
    ctx = mp.get_context("fork") if "fork" in mp.get_all_start_methods() else mp.get_context()
    with ProcessPoolExecutor(threads, mp_context=ctx) as ex:
        for r in ex.map(_thal_chunk, chunks, [kw] * len(chunks)):
            res.extend(r)
    return res


def score_sites(sites: np.ndarray, primers, genome: Genome, params: Params, threads: int = 0):
    """Thermodynamic score per site. Returns a list of dicts (one per site)."""
    threads = threads or os.cpu_count() or 1
    kw = dict(mv_conc=params.mv_conc, dv_conc=params.dv_conc, dntp_conc=params.dntp_conc,
              dna_conc=params.dna_conc, temp_c=params.ta_c)
    fl = params.flank
    meta, tasks = [], []
    for ci, pi, strand, start, mm in sites:
        ci = int(ci); pi = int(pi); strand = int(strand); start = int(start)
        p = primers[pi]
        L = len(p)
        seq = genome.seqs[ci]
        a = max(0, start - fl)
        b = min(len(seq), start + L + fl)
        window = seq[a:b]
        foot = seq[start:start + L]
        if strand == 1:
            template = revcomp(window)            # strand the primer anneals to, 5'->3'
            same = foot                           # same sense as primer
        else:
            template = window
            same = revcomp(foot)
        if any(c not in "ACGT" for c in template):
            continue
        mm_pos = sorted(L - i for i in range(L) if p[i] != same[i])       # 1 = 3' terminal
        term_mm = bool(mm_pos) and mm_pos[0] == 1
        # proofreading: a frayed / mismatched 3' end is trimmed, binding of the rest decides.
        # Taq with a terminal mismatch: END1 would report a slipped alignment, so use ANY binding
        # and let extension_factor carry the terminal-mismatch penalty.
        need_end = not (params.polymerase == "proofreading" or term_mm)
        meta.append((ci, pi, strand, start, L, mm_pos, p[-1], revcomp(same[-1])))
        tasks.append((p, template, need_end))
    thal = _run_thal(tasks, kw, threads)
    out = []
    for (ci, pi, strand, start, L, mm_pos, p3, tmpl_base), (dg_any, tm_any, dg_end) in zip(meta, thal):
        f_any = frac_bound(dg_any, params.ta_c, params.dna_conc, params.conc_divisor)
        f_use = f_any if dg_end != dg_end else frac_bound(dg_end, params.ta_c, params.dna_conc, params.conc_divisor)
        ext = extension_factor(mm_pos, p3, tmpl_base, params)
        score = f_use * ext
        if score < params.keep_min_score:
            continue
        out.append(dict(contig=ci, primer=pi, strand=strand, start=start, length=L,
                        end3=(start + L - 1) if strand == 1 else start,
                        mm=len(mm_pos), mm_last5=sum(1 for x in mm_pos if x <= 5),
                        mm_pos=",".join(map(str, mm_pos)),
                        dg_any=round(dg_any, 2), dg_end=(round(dg_end, 2) if dg_end == dg_end else None),
                        tm_any=round(tm_any, 1), f_any=f_any, f_use=f_use, ext=ext, score=score,
                        primable=score >= params.site_min_score))
    return out


# --------------------------------------------------------------------------- products
def merge_loci(scored, merge_nt):
    """Collapse shifted duplicates: thal on a footprint a few nt away from a real binding site finds the
    real (shifted) duplex again, so one locus yields several high-scoring 'sites'. Keep the best per locus."""
    by = {}
    for s in scored:
        by.setdefault((s["contig"], s["primer"], s["strand"]), []).append(s)
    out = []
    for v in by.values():
        v.sort(key=lambda s: s["end3"])
        best = v[0]; first = v[0]["end3"]
        for s in v[1:]:
            # FIX (verification 2026-11-11): window anchored on the locus' first 3' end, not on the previous
            # site, so sites every <= merge_nt across a tandem repeat are not chained into one locus
            if s["end3"] - first <= merge_nt:
                if s["score"] > best["score"] or (s["score"] == best["score"] and s["mm"] < best["mm"]):
                    best = s
            else:
                out.append(best); best = s; first = s["end3"]
        out.append(best)
    return out


def enumerate_products(scored, primers, genome: Genome, params: Params, names=None, pairs=None,
                       only_primable=True, targets=None):
    """All products between a '+' site and a downstream '-' site (any primers). score = s1 * s2.
    pairs: optional set of (i, j) primer index tuples that are designed pairs (order free).
    targets: optional {(i, j): (contig_index, start, end)}, designed amplicon coordinates (0-based, end
             exclusive, +-5 nt). With targets only the product at those coordinates is 'intended'. Without
             targets every perfect product of a designed pair is 'intended', but if a pair has more than one
             such product all of them get extra_copy=True.
    FIX (verification 2026-11-11):
      1. Sites enter enumeration at score >= product_possible (was: >= site_min_score), else the 'possible'
         tier silently missed products of one strong site and one site scoring 1e-4 to 0.01.
      2. A designed pair amplifying a second genomic copy or a long product between repeated sites was
         labelled 'intended' and escaped the 'reject any unintended likely product' rule."""
    names = names or [str(i) for i in range(len(primers))]
    pairset = set()
    for a, b in (pairs or []):
        pairset.add((a, b)); pairset.add((b, a))
    tgt = {}
    for (a, b), loc in (targets or {}).items():
        tgt[(a, b)] = loc; tgt[(b, a)] = loc
    min_site = params.product_possible if only_primable else 0.0
    scored = [s for s in scored if s["score"] >= min_site]
    if params.locus_merge_nt > 0:
        scored = merge_loci(scored, params.locus_merge_nt)
    prods = []
    for ci in range(len(genome.seqs)):
        glen = len(genome.seqs[ci])
        S = [s for s in scored if s["contig"] == ci]
        plus = sorted((s for s in S if s["strand"] == 1), key=lambda s: s["start"])
        minus = sorted((s for s in S if s["strand"] == -1), key=lambda s: s["start"] + s["length"])
        m_end = np.array([s["start"] + s["length"] for s in minus], dtype=np.int64)
        for p in plus:
            shifts = (0, glen) if params.circular else (0,)
            for sh in shifts:
                lo = p["start"] + params.min_size - sh
                hi = p["start"] + params.max_size - sh
                i0 = np.searchsorted(m_end, lo, side="left")
                i1 = np.searchsorted(m_end, hi, side="right")
                for m in minus[i0:i1]:
                    end = m["start"] + m["length"] + sh
                    # 3' ends must not cross: '-' primer footprint must start at or after '+' footprint start
                    if m["start"] + sh < p["start"]:
                        continue
                    size = end - p["start"]
                    sc = p["score"] * m["score"]
                    intended = ((p["primer"], m["primer"]) in pairset and p["mm"] == 0 and m["mm"] == 0)
                    if intended and tgt:
                        loc = tgt.get((p["primer"], m["primer"]))
                        intended = (loc is not None and loc[0] == ci and abs(loc[1] - p["start"]) <= 5
                                    and abs(loc[2] - (end % glen if sh else end)) <= 5)
                    if sc >= params.product_likely:
                        call = "likely"
                    elif sc >= params.product_possible:
                        call = "possible"
                    else:
                        continue            # FIX: do not store the many 'unlikely' products (memory)
                    prods.append(dict(contig=genome.names[ci], start=p["start"], end=end % glen if sh else end,
                                      size=size, fwd=names[p["primer"]], rev=names[m["primer"]],
                                      fwd_mm=p["mm"], rev_mm=m["mm"], fwd_score=p["score"], rev_score=m["score"],
                                      score=sc, call=call, intended=intended,
                                      same_primer=p["primer"] == m["primer"], wraps_origin=bool(sh),
                                      pair_key=tuple(sorted((int(p["primer"]), int(m["primer"])))),
                                      extra_copy=False))
    if not tgt:
        from collections import Counter
        n_int = Counter(d["pair_key"] for d in prods if d["intended"])
        for d in prods:
            if d["intended"] and n_int[d["pair_key"]] > 1:
                d["extra_copy"] = True
    prods.sort(key=lambda d: -d["score"])
    return prods


def check_pool(primers, genome_fasta, params: Params | None = None, names=None, pairs=None,
               method="seed", threads=0, genome: Genome | None = None):
    """Convenience wrapper: returns (scored_sites, products, timings)."""
    params = params or Params()
    t = {}
    t0 = time.perf_counter()
    genome = genome or load_genome(genome_fasta)
    t["load_genome"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    if method == "seed":
        sites = find_sites_seed(primers, genome, params, threads)
    elif method == "blast":
        sites = find_sites_blast(primers, genome_fasta, genome, params, threads)
    elif method == "brute":
        sites = find_sites_brute(primers, genome, params.max_mm, threads)
    else:
        raise ValueError(method)
    t["find_sites"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    scored = score_sites(sites, primers, genome, params, threads)
    t["score_sites"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    prods = enumerate_products(scored, primers, genome, params, names, pairs)
    t["products"] = time.perf_counter() - t0
    t["n_candidate_sites"] = len(sites)
    return scored, prods, t


# --------------------------------------------------------------------------- CLI
def _read_primers(path):
    names, seqs = [], []
    with open(path) as fh:
        txt = fh.read()
    if txt.lstrip().startswith(">"):
        for block in txt.split(">")[1:]:
            lines = block.strip().splitlines()
            names.append(lines[0].split()[0]); seqs.append("".join(lines[1:]).upper())
    else:
        for i, line in enumerate(l for l in txt.splitlines() if l.strip()):
            f = line.split()
            if len(f) == 1:
                names.append(f"primer{i + 1}"); seqs.append(f[0].upper())
            elif set(f[0].upper()) <= set("ACGTN"):
                seqs.append(f[0].upper()); names.append(f[1])
            else:
                names.append(f[0]); seqs.append(f[1].upper())
    return names, seqs


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="Genome-wide specificity check for a primer pool")
    ap.add_argument("primers", help="FASTA, or text with one primer per line (seq [name] or name seq)")
    ap.add_argument("genome", help="reference FASTA")
    ap.add_argument("--ta", type=float, default=60.0)
    ap.add_argument("--polymerase", choices=["taq", "proofreading"], default="taq")
    ap.add_argument("--method", choices=["seed", "blast", "brute"], default="seed")
    ap.add_argument("--max-size", type=int, default=2000)
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--out", default="specificity")
    a = ap.parse_args(argv)
    names, seqs = _read_primers(a.primers)
    prm = Params(ta_c=a.ta, polymerase=a.polymerase, max_size=a.max_size)
    scored, prods, t = check_pool(seqs, a.genome, prm, names=names, method=a.method, threads=a.threads)
    with open(a.out + "_sites.tsv", "w") as fh:
        cols = ["primer", "contig", "strand", "start", "length", "mm", "mm_last5", "mm_pos", "dg_any", "dg_end",
                "tm_any", "f_any", "ext", "score", "primable"]
        fh.write("\t".join(cols) + "\n")
        for s in scored:
            row = dict(s); row["primer"] = names[s["primer"]]
            fh.write("\t".join(f"{row[c]:.4g}" if isinstance(row[c], float) else str(row[c]) for c in cols) + "\n")
    with open(a.out + "_products.tsv", "w") as fh:
        cols = ["contig", "start", "end", "size", "fwd", "rev", "fwd_mm", "rev_mm", "fwd_score", "rev_score",
                "score", "call", "same_primer", "wraps_origin"]
        fh.write("\t".join(cols) + "\n")
        for p in prods:
            fh.write("\t".join(f"{p[c]:.4g}" if isinstance(p[c], float) else str(p[c]) for c in cols) + "\n")
    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in t.items()}, file=sys.stderr)
    print(f"sites scored {len(scored)}, primable {sum(s['primable'] for s in scored)}, "
          f"products likely {sum(p['call'] == 'likely' for p in prods)}, "
          f"possible {sum(p['call'] == 'possible' for p in prods)}", file=sys.stderr)


if __name__ == "__main__":
    main()
