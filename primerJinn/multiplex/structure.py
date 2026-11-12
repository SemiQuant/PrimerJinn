"""
structure_reference.py

Reference implementation of single-stranded DNA secondary-structure metrics for
multiplex PCR primer design (PrimerJinn rebuild). Engine: ViennaRNA >= 2.7.0 with
DNA parameters. NUPACK 4 is an optional cross-check only (academic licence).

Tested with ViennaRNA 2.7.0, primer3-py 2.0.3, nupack 4.0.2, Python 3.13.

Conventions
-----------
* Energies in kcal/mol, temperatures in degrees Celsius, concentrations in mol/L
  unless a name ends in _mM / _nM.
* Sequences are 5'->3'. Positions passed to the public functions are 0-based,
  end-exclusive (Python slices) on the strand that is being folded.
* "Template strand" of a primer = the strand the primer hybridises to. For an
  amplicon written 5'->3' (top strand):
    - the forward primer equals the first len(F) nt of the top strand and binds
      the LAST len(F) nt of the bottom strand (reverse complement);
    - the reverse primer equals the first len(R) nt of the bottom strand and binds
      the LAST len(R) nt of the top strand.
  Template nucleotide j of the site (5'->3' on the template) pairs with primer
  nucleotide len-1-j, so the primer 3' end pairs with the 5'-most nt of the site.

IMPORTANT parameter note (measured, see secondary_structure.md)
---------------------------------------------------------------
The DNA set "Mathews2004" bundled with ViennaRNA 2.7.0 gives non-physical results
above about 55 C (small interior loops become stabilising, ensemble free energy of
a 300 nt strand gets MORE negative from 60 C to 90 C). The bundled "Mathews1999"
DNA set behaves physically from 37 to about 90 C (not above: non-monotonic near 99 C)
and is the default here, BUT its agreement with NUPACK dna04 depends on GC content
(verified 2026-11-11, verify_ss2/v2.out, PCR salt, 300 nt):
  * 60 to 70% GC: G_ens within about 2 kcal/mol of NUPACK (ss3/bench2.out, MAE 1.8 to 2.1);
  * 40 to 55% GC: Vienna 4 to 9 kcal/mol LESS stable than NUPACK;
  * about 80% GC: Vienna 3 to 14 kcal/mol MORE stable than NUPACK.
On a 20-bp AT-rich duplex at 65 C, M1999 gives -3.0 kcal/mol vs -4.6 (primer3 thal,
SantaLucia) and -7.1 (NUPACK): M1999 helices are weaker than SantaLucia for AT-rich
stacks. Treat absolute energies as GC-biased; compare primers within a pool only.

Other API facts verified on ViennaRNA 2.7.0
* RNA.params_load_DNA_Mathews1999() / ..._Mathews2004() return 1 and also set
  helical_rise 3.4 and backbone_length 6.76 for the salt model, but leave
  md.noGU = 0: set md.noGU = 1 yourself (otherwise G-T pairs appear in DNA folds).
* md.salt is monovalent salt in mol/L, default 1.021. There is no Mg2+ term, so
  Mg is folded in as Na_eq = [Na+] + 3.3*sqrt([Mg2+]_free) (NUPACK's DNA rule,
  verified numerically to 5 digits).
* RNA.pfl_fold_up() uses the GLOBAL default model (37 C), not a md object: use
  fold_compound(..., RNA.OPTION_WINDOW).probs_window(..., RNA.PROBS_WINDOW_UP, cb).
* fold_compound.eval_structure() segfaults on an unbalanced dot-bracket string
  (SIGSEGV, exit -11) and silently returns 100000 for a length mismatch: use
  safe_eval_structure() below.
* RNA.pfl_fold_up() can be made to honour T/noGU/salt only via the RNA.cvar
  globals (verified), which is not thread-safe; use probs_window.
* primer3-py 2.0.3 calc_tm() defaults to formamide_conc = 0.8 mol/L (argdefaults.py),
  which lowers every Tm by about 2 C versus primer3_core (default 0.0). Always pass
  formamide_conc=0.0: see primer_tm() below.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import RNA
import primer3

R_KCAL = 0.0019872  # kcal/(mol K)
_PARAM_SET = None


# --------------------------------------------------------------------------
# conditions
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Conditions:
    """PCR buffer. Defaults: 50 mM monovalent, 2.0 mM Mg2+, 0.8 mM total dNTP,
    200 nM of each primer."""
    mono_mM: float = 50.0
    mg_mM: float = 2.0
    dntp_mM: float = 0.8
    primer_nM: float = 200.0
    dmso_pct: float = 0.0          # % v/v
    dmso_c_per_pct: float = 0.6    # Tm depression per % DMSO (heuristic, 0.5 to 0.75)
    param_set: str = "Mathews1999"  # ViennaRNA DNA parameter set
    dangles: int = 2

    @property
    def mg_free_M(self) -> float:
        # dNTPs chelate Mg2+ about 1:1 (von Ahsen 2001; same rule as primer3)
        return max(self.mg_mM - self.dntp_mM, 0.0) / 1000.0

    @property
    def na_eq_M(self) -> float:
        """Sodium equivalent used for FOLDING. Na_eq = [mono] + 3.3*sqrt([Mg_free]),
        molar units (the conversion NUPACK 4 applies to DNA; verified numerically).
        ViennaRNA's own salt model takes monovalent salt only."""
        return self.mono_mM / 1000.0 + 3.3 * math.sqrt(self.mg_free_M)

    def t_eff(self, celsius: float) -> float:
        """DMSO handled as an effective temperature shift (crude heuristic)."""
        return celsius + self.dmso_pct * self.dmso_c_per_pct


def load_params(name: str = "Mathews1999", force: bool = False) -> None:
    """Load a ViennaRNA DNA parameter set (process-global state). The cache can go
    stale if other code calls RNA.params_load_*; pass force=True after that."""
    global _PARAM_SET
    if _PARAM_SET == name and not force:
        return
    ok = {"Mathews1999": RNA.params_load_DNA_Mathews1999,
          "Mathews2004": RNA.params_load_DNA_Mathews2004}[name]()
    if not ok:
        raise RuntimeError(f"could not load ViennaRNA DNA parameters {name}")
    _PARAM_SET = name


def _md(celsius: float, cond: Conditions, window: int = 0) -> "RNA.md":
    load_params(cond.param_set)
    md = RNA.md()
    md.temperature = cond.t_eff(celsius)
    md.noGU = 1                    # no G-T wobble pairs for DNA
    md.dangles = cond.dangles
    md.salt = cond.na_eq_M         # ViennaRNA reference salt is 1.021 M
    md.helical_rise = RNA.MODEL_HELICAL_RISE_DNA        # 3.4 A (set by params_load too)
    md.backbone_length = RNA.MODEL_BACKBONE_LENGTH_DNA  # 6.76 A
    if window:
        md.window_size = window
        md.max_bp_span = window
    return md


def safe_eval_structure(seq: str, structure: str, md) -> float:
    """eval_structure with input validation (ViennaRNA 2.7.0 segfaults on unbalanced
    dot-bracket strings and returns 100000 on a length mismatch)."""
    if len(seq) != len(structure):
        raise ValueError("sequence and structure lengths differ")
    depth = 0
    for ch in structure:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth < 0:
                raise ValueError("unbalanced dot-bracket")
        elif ch != ".":
            raise ValueError(f"unsupported character {ch!r}")
    if depth:
        raise ValueError("unbalanced dot-bracket")
    return RNA.fold_compound(seq.upper(), md).eval_structure(structure)


def primer_tm(primer: str, cond: "Conditions" = None) -> float:
    """Nearest-neighbour primer Tm via primer3-py with formamide explicitly 0
    (primer3-py 2.0.3 defaults to 0.8 mol/L formamide, about -2 C)."""
    cond = cond or Conditions()
    return primer3.calc_tm(primer, mv_conc=cond.mono_mM, dv_conc=cond.mg_mM,
                           dntp_conc=cond.dntp_mM, dna_conc=cond.primer_nM,
                           dmso_conc=cond.dmso_pct, dmso_fact=cond.dmso_c_per_pct,
                           formamide_conc=0.0)


def revcomp(seq: str) -> str:
    return seq.upper().translate(str.maketrans("ACGT", "TGCA"))[::-1]


# --------------------------------------------------------------------------
# folding primitives
# --------------------------------------------------------------------------
@dataclass
class FoldStats:
    celsius: float
    mfe: float
    structure: str
    g_ensemble: float
    frac_paired: float            # ensemble mean fraction of paired bases
    p_unpaired: List[float]       # per nucleotide, 0-based


def fold_stats(seq: str, celsius: float, cond: Conditions = Conditions()) -> FoldStats:
    seq = seq.upper()
    fc = RNA.fold_compound(seq, _md(celsius, cond))
    ss, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    _, g = fc.pf()
    bpp = fc.bpp()
    n = len(seq)
    paired = [0.0] * (n + 1)
    for i in range(1, n + 1):
        row = bpp[i]
        for j in range(i + 1, n + 1):
            p = row[j]
            if p > 1e-9:
                paired[i] += p
                paired[j] += p
    pu = [max(0.0, 1.0 - paired[i]) for i in range(1, n + 1)]
    return FoldStats(celsius, mfe, ss, g, 1.0 - sum(pu) / n, pu)


def unpaired_windows(seq: str, celsius: float, max_len: int,
                     cond: Conditions = Conditions()) -> Dict[int, List[Optional[float]]]:
    """Exact probability that a window is completely unpaired, for every window
    of length 1..max_len, from one partition-function pass (RNAplfold algorithm,
    window = whole sequence). Returns up[i][L] = P(nt i-L+1..i unpaired), i 1-based
    END position of the window."""
    seq = seq.upper()
    n = len(seq)
    md = _md(celsius, cond, window=n)
    fc = RNA.fold_compound(seq, md, RNA.OPTION_WINDOW)
    out: Dict[int, List[Optional[float]]] = {}

    def cb(v, v_size, i, maxsize, what, data):
        if what & RNA.PROBS_WINDOW_UP:
            data[i] = list(v)

    if not fc.probs_window(max_len, RNA.PROBS_WINDOW_UP, cb, out):
        raise RuntimeError("probs_window failed")
    return out


def dg_open(p_open: float, celsius: float) -> float:
    """Free energy cost of opening a window, dG_open = -RT ln P(window unpaired)."""
    return -R_KCAL * (celsius + 273.15) * math.log(max(p_open, 1e-300))


def dg_open_constrained(seq: str, start: int, end: int, celsius: float,
                        cond: Conditions = Conditions()) -> float:
    """Same quantity by the definition dG_open = G(window forced unpaired) - G(ensemble).
    Slower (two partition functions); used as a test of unpaired_windows()."""
    seq = seq.upper()
    md = _md(celsius, cond)
    fc = RNA.fold_compound(seq, md)
    _, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    _, g = fc.pf()
    fc2 = RNA.fold_compound(seq, md)
    for i in range(start + 1, end + 1):
        fc2.hc_add_up(i, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
    fc2.exp_params_rescale(mfe)
    _, g2 = fc2.pf()
    return g2 - g


# --------------------------------------------------------------------------
# primer : template duplex (nearest neighbour, primer3 thal)
# --------------------------------------------------------------------------
def duplex_dh_ds(primer: str, cond: Conditions = Conditions()) -> Tuple[float, float]:
    """dH (kcal/mol) and dS (kcal/mol/K) of primer + perfect complement, salt
    corrected by primer3 (SantaLucia 1998 NN parameters, monovalent + divalent +
    dNTP). Swap this for the project's main duplex model if a better one is used."""
    primer = primer[-60:]   # thal needs one strand <= 60 nt; longer duplexes are saturated anyway
    r = primer3.calc_heterodimer(primer, revcomp(primer), mv_conc=cond.mono_mM,
                                 dv_conc=cond.mg_mM, dntp_conc=cond.dntp_mM,
                                 dna_conc=cond.primer_nM, temp_c=60.0)
    return r.dh / 1000.0, r.ds / 1000.0


def primer_self_fold_dg(primer: str, celsius: float, cond: Conditions = Conditions()) -> float:
    """Ensemble free energy of the primer alone (<= 0). Opening cost is its negative."""
    fc = RNA.fold_compound(primer.upper(), _md(celsius, cond))
    _, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    return min(fc.pf()[1], 0.0)


@dataclass
class SiteMetrics:
    celsius: float
    p_open_site: float        # whole binding site unpaired on the template
    dg_open_site: float
    p_open_seed3: float       # template nt opposite the primer 3'-terminal 8 nt
    p_nucleation: float       # best 8-nt toehold anywhere in the site
    dg_duplex: float          # primer + open complement
    dg_primer_fold: float     # <= 0
    dg_eff: float             # dg_duplex + dg_open_site - dg_primer_fold
    f_bound: float            # equilibrium fraction of sites bound, with structure
    f_bound_nostruct: float   # same, ignoring template and primer structure
    p_primed_kinetic: float   # 1 - exp(-k_on [P] p_nucleation t_anneal), capped by f_bound


def site_metrics(template: str, site_start: int, site_end: int, primer: str,
                 celsius: float, cond: Conditions = Conditions(),
                 seed: int = 8, k_on: float = 1.0e6, t_anneal_s: float = 30.0,
                 _up: Optional[dict] = None) -> SiteMetrics:
    """Coupled-equilibrium metrics of one primer on one template strand.

    template[site_start:site_end] must be the reverse complement of the primer
    (or of its template-specific 3' part)."""
    L = site_end - site_start
    up = _up if _up is not None else unpaired_windows(template, celsius, L, cond)
    p_site = up[site_end][L]
    # primer 3' end pairs with the 5'-most nucleotides of the site
    p_seed3 = up[site_start + seed][seed]
    p_nuc = max(up[e][seed] for e in range(site_start + seed, site_end + 1))
    t_k = cond.t_eff(celsius) + 273.15
    dh, ds = duplex_dh_ds(primer, cond)
    dg_dup = dh - t_k * ds
    dg_pf = primer_self_fold_dg(primer, celsius, cond)
    dgo = -R_KCAL * t_k * math.log(max(p_site, 1e-300))
    dg_eff = dg_dup + dgo - dg_pf
    conc = cond.primer_nM * 1e-9

    def frac(dg):
        x = -dg / (R_KCAL * t_k) + math.log(conc)
        return 1.0 / (1.0 + math.exp(-x)) if x > -700 else 0.0

    f = frac(dg_eff)
    kin = 1.0 - math.exp(-k_on * conc * p_nuc * t_anneal_s)
    return SiteMetrics(celsius, p_site, dgo, p_seed3, p_nuc, dg_dup, dg_pf, dg_eff,
                       f, frac(dg_dup), min(f, kin))


# --------------------------------------------------------------------------
# amplicon level
# --------------------------------------------------------------------------
G4_RE = re.compile(r"(?=(G{3,}[ACGT]{1,7}?G{3,}[ACGT]{1,7}?G{3,}[ACGT]{1,7}?G{3,}))")


def g4_motifs(seq: str) -> List[Tuple[int, int, str]]:
    """Canonical G-quadruplex motif G3+ N1-7 G3+ N1-7 G3+ N1-7 G3+ on both strands.
    Returns (start, end, strand) in top-strand coordinates."""
    seq = seq.upper()
    hits = []
    for strand, s in (("+", seq), ("-", revcomp(seq))):
        last_end = -1
        for m in G4_RE.finditer(s):
            a, b = m.start(1), m.end(1)
            if a < last_end:
                continue
            last_end = b
            hits.append((a, b, strand) if strand == "+" else (len(s) - b, len(s) - a, strand))
    return hits


def duplex_tm_long(seq: str, cond: Conditions = Conditions()) -> float:
    """Empirical melting temperature of a long duplex (PCR product), exactly as
    primer3 long_seq_tm() (oligotm.c) computes PRIMER_PRODUCT_TM:
      Tm = 81.5 + 16.6 log10([Na_eq]/1000) + 41 * GC_fraction - 600/N - dmso*factor
      [Na_eq] (mM) = mono + 120 * sqrt(Mg - dNTP)   (primer3 divalent_to_monovalent)
    Accuracy is a few degrees; use it to rank amplicons and set denaturation advice,
    not as an absolute melting point."""
    seq = seq.upper()
    gc = (seq.count("G") + seq.count("C")) / len(seq)
    mg = cond.mg_mM if cond.mg_mM > 0 else 0.0
    dntp = cond.dntp_mM if mg > 0 else 0.0
    na_mM = cond.mono_mM + 120.0 * math.sqrt(max(mg, dntp) - dntp)
    return (81.5 + 16.6 * math.log10(na_mM / 1000.0) + 41.0 * gc - 600.0 / len(seq)
            - cond.dmso_pct * cond.dmso_c_per_pct)


def site_accessibility_context(context: str, site_start: int, site_end: int, celsius: float,
                               cond: Conditions = Conditions(), span: int = 150,
                               window: int = 200) -> float:
    """P(site unpaired) on a long template (genomic context, early cycles), RNAplfold
    style: local folding with base-pair span <= span averaged over windows of length
    `window`. On 10 M. tuberculosis sites (65 C, 164 mM Na_eq) W=200, L=150 on a
    900-nt genomic context reproduced the full-context value with r = 0.98 (log10 P).
    Positions are 0-based, end-exclusive on `context`."""
    seq = context.upper()
    md = _md(celsius, cond)
    md.window_size = min(window, len(seq))
    md.max_bp_span = min(span, md.window_size)
    fc = RNA.fold_compound(seq, md, RNA.OPTION_WINDOW)
    out: Dict[int, List[float]] = {}
    L = site_end - site_start

    def cb(v, v_size, i, maxsize, what, data):
        if what & RNA.PROBS_WINDOW_UP:
            data[i] = list(v)

    fc.probs_window(L, RNA.PROBS_WINDOW_UP, cb, out)
    return out[site_end][L]


def local_melting_tm(sub: str, cond: Conditions = Conditions(),
                     grid: Sequence[float] = tuple(range(40, 101)), engine: str = "vienna") -> float:
    """Ensemble melting temperature of a short local structure: the lowest grid
    temperature at which the open chain is the majority state,
    P(open) = exp(G_ensemble / RT) >= 0.5, i.e. G_ensemble >= -RT ln 2.
    Computed with the same salt as the PCR (no two-point extrapolation).
    Returns 101.0 if the structure is still folded at 100 C.
    Vienna values above about 90 C are unreliable (M1999 is non-monotonic near
    99 C): read any Vienna result > 90 C as "> 90 C", not as a number.
    engine='nupack' gives the same quantity with NUPACK dna04 (optional cross-check)."""
    sub = sub.upper()
    for t in grid:
        c0 = Conditions(**{**cond.__dict__, "dmso_pct": 0.0})
        if engine == "nupack":
            import nupack
            mdl = nupack.Model(material="dna", celsius=t, sodium=min(max(c0.mono_mM / 1000.0, 0.05), 1.1),
                               magnesium=min(c0.mg_free_M, 0.2))
            g = float(nupack.pfunc([sub], model=mdl)[1])
        else:
            fc = RNA.fold_compound(sub, _md(t, c0))
            _, m = fc.mfe()
            fc.exp_params_rescale(m)
            g = fc.pf()[1]
        if g >= -R_KCAL * (t + 273.15) * math.log(2.0):
            return t - cond.dmso_pct * cond.dmso_c_per_pct
    return 101.0


def self_priming_prob(strand: str, celsius: float, cond: Conditions = Conditions(),
                      k: int = 1) -> float:
    """k = 1: probability that the 3'-terminal nucleotide is paired inside an
    intramolecular helix. k > 1: probability that AT LEAST ONE of the 3'-terminal k
    nt is paired (marginals do not give 'all k paired'), i.e. the strand can fold back and
    self-prime (the mechanism Gimpel et al. 2025, doi:10.1038/s41467-025-64221-4,
    identified as the main cause of poor amplification in multi-template PCR).
    Computed from the base-pair probability matrix: for k=1 this is 1 - p_unpaired(n)."""
    fs = fold_stats(strand, celsius, cond)
    if k == 1:
        return 1.0 - fs.p_unpaired[-1]
    up = unpaired_windows(strand, celsius, k, cond)
    # lower bound on 'all k paired' is not available from marginals; report P(any of the k paired)
    return 1.0 - up[len(strand)][k]


def max_window_gc(seq: str, w: int = 50) -> float:
    seq = seq.upper()
    if len(seq) <= w:
        return (seq.count("G") + seq.count("C")) / len(seq)
    best = 0.0
    c = sum(ch in "GC" for ch in seq[:w])
    best = c
    for i in range(w, len(seq)):
        c += (seq[i] in "GC") - (seq[i - w] in "GC")
        best = max(best, c)
    return best / w


@dataclass
class LocalHairpin:
    start: int
    end: int
    structure: str
    dg: float           # at the folding temperature
    tm: float           # two-state melting temperature of this structure (dG = 0)
    dg_72: float = 0.0  # dG of this structure at 72 C (extension)


def structure_tm(sub: str, structure: str, cond: Conditions = Conditions(),
                 lo: float = 20.0, hi: float = 110.0) -> float:
    """Two-state melting temperature of one specific structure: the temperature at
    which its free energy relative to the open chain is zero (bisection on
    eval_structure). Values above about 90 C are extrapolations of the parameters."""
    def g(t):
        return safe_eval_structure(sub, structure, _md(t, cond))
    if g(lo) >= 0:
        return lo
    if g(hi) < 0:
        return hi
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if g(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def most_stable_local_structure(seq: str, celsius: float, cond: Conditions = Conditions(),
                                span: int = 60) -> Optional[LocalHairpin]:
    """Most stable locally closed structure (base-pair span <= span nt) at `celsius`,
    with its ensemble melting temperature from local_melting_tm() (Vienna; values
    differed from NUPACK by 0 to 20 C on Mtb hairpins, see secondary_structure.md)."""
    seq = seq.upper()
    md = _md(celsius, cond, window=min(span, len(seq)))
    fc = RNA.fold_compound(seq, md, RNA.OPTION_WINDOW)
    hits: list = []

    def cb(start, end, structure, energy, data):
        data.append((energy, start, end, structure))

    fc.mfe_window_cb(cb, hits)
    if not hits:
        return None
    energy, start, end, structure = min(hits)
    sub = seq[start - 1:start - 1 + len(structure)]
    tm = structure_tm(sub, structure, cond)
    dg72 = safe_eval_structure(sub, structure, _md(72.0, cond))
    return LocalHairpin(start - 1, start - 1 + len(structure), structure, energy, tm, dg72)


@dataclass
class StrandReport:
    strand: str
    celsius: float
    mfe: float
    g_ensemble: float
    g_per_nt: float
    frac_paired: float


@dataclass
class AmpliconReport:
    length: int
    gc: float
    max_gc_50: float
    duplex_tm: float
    g4: List[Tuple[int, int, str]]
    strands: List[StrandReport]
    hairpin_top: Optional[LocalHairpin]
    hairpin_bottom: Optional[LocalHairpin]
    fwd: Dict[float, SiteMetrics] = field(default_factory=dict)
    rev: Dict[float, SiteMetrics] = field(default_factory=dict)
    flags: List[str] = field(default_factory=list)


def amplicon_report(amplicon: str, fwd_len: int, rev_len: int,
                    temps: Sequence[float] = (60.0, 72.0), ta: Optional[float] = None,
                    cond: Conditions = Conditions(),
                    fwd_ctx: Optional[Tuple[str, int, int]] = None,
                    rev_ctx: Optional[Tuple[str, int, int]] = None) -> AmpliconReport:
    """Structure report for one amplicon (top strand 5'->3', primers included).

    temps: temperatures for strand-level statistics (annealing and extension).
    ta:    annealing temperature for the primer-site metrics (default temps[0]).
    fwd_ctx / rev_ctx: optional (template_strand, site_start, site_end) giving a
        different folding context for the primer site (for example the genomic
        flank used in the first cycles). Default: the amplicon strand itself.
    """
    amp = amplicon.upper()
    bottom = revcomp(amp)
    n = len(amp)
    ta = temps[0] if ta is None else ta
    fwd, rev = amp[:fwd_len], bottom[:rev_len]
    rep = AmpliconReport(
        length=n, gc=(amp.count("G") + amp.count("C")) / n, max_gc_50=max_window_gc(amp),
        duplex_tm=duplex_tm_long(amp, cond), g4=g4_motifs(amp), strands=[],
        hairpin_top=most_stable_local_structure(amp, ta, cond),
        hairpin_bottom=most_stable_local_structure(bottom, ta, cond))
    for t in temps:
        for name, s in (("top", amp), ("bottom", bottom)):
            fs = fold_stats(s, t, cond)
            rep.strands.append(StrandReport(name, t, fs.mfe, fs.g_ensemble,
                                            fs.g_ensemble / n, fs.frac_paired))
    f_t, f_a, f_b = fwd_ctx if fwd_ctx else (bottom, n - fwd_len, n)
    r_t, r_a, r_b = rev_ctx if rev_ctx else (amp, n - rev_len, n)
    rep.fwd[ta] = site_metrics(f_t, f_a, f_b, fwd, ta, cond)
    rep.rev[ta] = site_metrics(r_t, r_a, r_b, rev, ta, cond)
    # flags (thresholds: see secondary_structure.md, section C)
    for nm, m in (("fwd", rep.fwd[ta]), ("rev", rep.rev[ta])):
        if m.dg_open_site > 3.0:
            rep.flags.append(f"{nm}: site dG_open {m.dg_open_site:.1f} > 3.0 kcal/mol at {ta:g} C (reject)")
        elif m.dg_open_site > 1.5:
            rep.flags.append(f"{nm}: site dG_open {m.dg_open_site:.1f} > 1.5 kcal/mol at {ta:g} C (penalise)")
        if m.p_nucleation < 0.1:
            rep.flags.append(f"{nm}: no accessible 8-nt toehold (p {m.p_nucleation:.3f} < 0.1)")
    # Local hairpins and product Tm are reported as numbers, not flags: in GC-rich
    # genomes nearly every amplicon has a hairpin stable at 72 C and no validated
    # threshold links them to dropout (see docs/DESIGN.md section 5).
    if rep.g4:
        rep.flags.append(f"{len(rep.g4)} G-quadruplex motif(s)")
    if rep.max_gc_50 >= 0.80:
        rep.flags.append(f"50-bp window with {100*rep.max_gc_50:.0f}% GC")
    return rep


# --------------------------------------------------------------------------
# annealing temperature recommendation
# --------------------------------------------------------------------------
@dataclass
class TaRecommendation:
    ta: Optional[float]
    t_low: Optional[float]      # lowest T at which every site has an accessible toehold
    t_high: Optional[float]     # highest T at which every primer is still >= f_min bound
    limiting_primer: Optional[str]
    limiting_site: Optional[str]
    table: Dict[float, Dict[str, SiteMetrics]]
    advice: List[str]


def recommend_ta(amplicons: Dict[str, Tuple[str, int, int]],
                 cond: Conditions = Conditions(),
                 grid: Sequence[float] = tuple(range(55, 73)),
                 f_min: float = 0.90, p_nuc_min: float = 0.10,
                 margin: float = 1.0) -> TaRecommendation:
    """amplicons: name -> (amplicon top strand, fwd primer length, rev primer length).

    T_high = highest grid temperature where min over primers of f_bound >= f_min.
    T_low  = lowest grid temperature where min over sites of p_nucleation >= p_nuc_min.
    Ta = T_high - margin if T_low <= T_high - margin, otherwise no single-temperature
    solution: the limiting amplicon should be redesigned or an additive used."""
    table: Dict[float, Dict[str, SiteMetrics]] = {}
    for t in grid:
        row = {}
        for name, (amp, fl, rl) in amplicons.items():
            amp = amp.upper()
            bottom = revcomp(amp)
            n = len(amp)
            up_b = unpaired_windows(bottom, t, fl, cond)
            up_t = unpaired_windows(amp, t, rl, cond)
            row[name + ":F"] = site_metrics(bottom, n - fl, n, amp[:fl], t, cond, _up=up_b)
            row[name + ":R"] = site_metrics(amp, n - rl, n, bottom[:rl], t, cond, _up=up_t)
        table[float(t)] = row
    ok_f = [t for t in table if min(m.f_bound for m in table[t].values()) >= f_min]
    ok_p = [t for t in table if min(m.p_nucleation for m in table[t].values()) >= p_nuc_min]
    t_high = max(ok_f) if ok_f else None
    t_low = min(ok_p) if ok_p else None
    advice: List[str] = []
    lim_primer = lim_site = None
    if t_high is not None:
        nxt = t_high + 1 if (t_high + 1) in table else t_high
        lim_primer = min(table[nxt], key=lambda k: table[nxt][k].f_bound)
    if t_low is not None:
        prv = t_low - 1 if (t_low - 1) in table else t_low
        lim_site = min(table[prv], key=lambda k: table[prv][k].p_nucleation)
    ta = None
    if t_high is None:
        advice.append("No grid temperature keeps every primer >= f_min bound: primers too weak for this grid, redesign or lower the grid.")
    elif t_low is None or t_low > t_high - margin:
        ta = t_high - margin
        advice.append("Template structure still blocks at least one site at the highest temperature the primers tolerate: "
                      f"redesign {lim_site}, or add 5% DMSO / 1 M betaine and re-run with dmso_pct set.")
    else:
        ta = t_high - margin
    gcs = [(a.upper().count("G") + a.upper().count("C")) / len(a) for a, _, _ in amplicons.values()]
    tms = [duplex_tm_long(a, cond) for a, _, _ in amplicons.values()]
    if max(gcs) >= 0.65 or max(tms) > 92.0:
        advice.append(f"GC-rich pool (max amplicon GC {100*max(gcs):.0f}%, max duplex Tm {max(tms):.1f} C): "
                      "use 98 C denaturation with a polymerase that tolerates it, and consider 3 to 5% DMSO or 1 M betaine.")
    if ta is not None and ta >= 68.0:
        advice.append("Ta >= 68 C: two-step PCR (combined anneal/extend at Ta) is possible.")
    return TaRecommendation(ta, t_low, t_high, lim_primer, lim_site, table, advice)


# --------------------------------------------------------------------------
# NUPACK cross-check (optional)
# --------------------------------------------------------------------------
def nupack_crosscheck(seq: str, celsius: float, cond: Conditions = Conditions(),
                      site: Optional[Tuple[int, int]] = None, n_sample: int = 20000,
                      ensemble: str = "stacking") -> Optional[dict]:
    """Same quantities from NUPACK 4 (dna04 parameters). Returns None when NUPACK is
    not installed. NUPACK has no hard-constraint partition function in its Python
    API, so the window probability is estimated by Boltzmann sampling (standard
    error sqrt(p(1-p)/n_sample))."""
    try:
        import nupack
        import numpy as np
    except ImportError:
        return None
    seq = seq.upper()
    model = nupack.Model(material="dna", celsius=cond.t_eff(celsius),
                         sodium=min(max(cond.mono_mM / 1000.0, 0.05), 1.1),
                         magnesium=min(cond.mg_free_M, 0.2), ensemble=ensemble)
    _, g = nupack.pfunc(strands=[seq], model=model)
    m = nupack.mfe(strands=[seq], model=model)[0]
    pairs = nupack.pairs(strands=[seq], model=model).to_array()
    out = {"mfe": float(m.energy), "structure": str(m.structure), "g_ensemble": float(g),
           "frac_paired": float(1.0 - np.diag(pairs).mean())}
    if site is not None:
        a, b = site
        samples = nupack.sample(strands=[seq], num_sample=n_sample, model=model)
        p = sum(1 for s in samples if set(str(s)[a:b]) == {"."}) / n_sample
        out["p_open_site"] = p
        out["dg_open_site"] = dg_open(p, cond.t_eff(celsius)) if p > 0 else float("inf")
    return out


# --------------------------------------------------------------------------
# self-test / demo
# --------------------------------------------------------------------------
def _demo(ref_fasta: str) -> None:
    from Bio import SeqIO
    g = str(next(SeqIO.parse(ref_fasta, "fasta")).seq).upper()
    cond = Conditions()
    print(f"Conditions: {cond}\nNa_eq for folding = {cond.na_eq_M*1000:.1f} mM")
    examples = {"rpoB_RRDR": (761155, 22, 22), "katG_S315": (2155168, 22, 22)}
    amps = {}
    for name, (pos, fl, rl) in examples.items():
        amp = g[pos - 150:pos + 150]
        amps[name] = (amp, fl, rl)
        rep = amplicon_report(amp, fl, rl, temps=(60.0, 65.0, 72.0), ta=65.0, cond=cond)
        print(f"\n=== {name}  NC_000962.3:{pos-149}-{pos+150}  len {rep.length}  GC {100*rep.gc:.1f}%  "
              f"max GC(50 bp) {100*rep.max_gc_50:.0f}%  duplex Tm {rep.duplex_tm:.1f} C  G4 motifs {rep.g4}")
        for s in rep.strands:
            print(f"  {s.strand:6s} {s.celsius:4.0f} C  MFE {s.mfe:7.2f}  G_ens {s.g_ensemble:7.2f}  "
                  f"({s.g_per_nt:+.3f}/nt)  paired {100*s.frac_paired:4.1f}%")
        for nm, hp in (("top", rep.hairpin_top), ("bottom", rep.hairpin_bottom)):
            if hp:
                print(f"  most stable local structure, {nm}: {hp.start}-{hp.end} dG(65 C) {hp.dg:.2f}  Tm {hp.tm:.1f} C  {hp.structure}")
        for nm, d in (("fwd", rep.fwd), ("rev", rep.rev)):
            m = d[65.0]
            print(f"  {nm} @65 C: p_open(site) {m.p_open_site:.3f} dG_open {m.dg_open_site:.2f}  p_open(3' seed) {m.p_open_seed3:.3f} "
                  f"p_nuc {m.p_nucleation:.3f}  dG_duplex {m.dg_duplex:.2f}  dG_primer_fold {m.dg_primer_fold:.2f}  "
                  f"dG_eff {m.dg_eff:.2f}  f_bound {m.f_bound:.3f} (no structure {m.f_bound_nostruct:.3f})  primed/cycle {m.p_primed_kinetic:.3f}")
        print("  flags:", rep.flags or "none")
        # internal consistency test: probs_window vs constrained partition function
        n = len(amp)
        a = dg_open_constrained(amp, n - rl, n, 65.0, cond)
        assert abs(a - rep.rev[65.0].dg_open_site) < 0.02, (a, rep.rev[65.0].dg_open_site)
        print(f"  test: rev-site dG_open by constrained pf {a:.3f} == probs_window {rep.rev[65.0].dg_open_site:.3f}")
        # early cycles: the reverse-primer site inside genomic DNA (top strand, +-300 nt)
        a0 = pos - 150
        ctx = g[a0 - 300:a0 + 600]
        p_full = unpaired_windows(ctx, 65.0, rl, cond)[600][rl]
        p_pl = site_accessibility_context(ctx, 600 - rl, 600, 65.0, cond)
        sp = [self_priming_prob(s, 65.0, cond) for s in (amp, revcomp(amp))]
        print(f"  P(3'-terminal nt paired, self-priming risk) @65 C: top {sp[0]:.3f}, bottom {sp[1]:.3f}")
        print(f"  rev site in genomic context @65 C: P(open) full fold {p_full:.3f}, plfold W200/L150 {p_pl:.3f} "
              f"(amplicon strand {rep.rev[65.0].p_open_site:.3f})")
        x = nupack_crosscheck(amp, 65.0, cond, site=(n - rl, n))
        if x:
            top65 = [s for s in rep.strands if s.strand == "top" and s.celsius == 65.0][0]
            print(f"  NUPACK cross-check top strand 65 C: G_ens {x['g_ensemble']:.2f} (Vienna {top65.g_ensemble:.2f})  "
                  f"paired {100*x['frac_paired']:.1f}% (Vienna {100*top65.frac_paired:.1f}%)  "
                  f"rev-site dG_open {x['dg_open_site']:.2f} (Vienna {rep.rev[65.0].dg_open_site:.2f})")
    rec = recommend_ta(amps, cond)
    print(f"\nRecommended Ta for the 2-amplicon pool: {rec.ta} C (T_low {rec.t_low}, T_high {rec.t_high}; "
          f"limiting primer {rec.limiting_primer}, limiting site {rec.limiting_site})")
    for a in rec.advice:
        print("  advice:", a)
    print("  T    min f_bound   min p_nuc")
    for t, row in rec.table.items():
        print(f"  {t:4.0f}   {min(m.f_bound for m in row.values()):.3f}        {min(m.p_nucleation for m in row.values()):.3f}")


if __name__ == "__main__":
    import sys
    _demo(sys.argv[1] if len(sys.argv) > 1 else
          "/Users/semiquant/Bioinformatics/Projects/Archived/PrimerJinn/example/ref.fasta")
