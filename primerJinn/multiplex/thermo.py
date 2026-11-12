"""Reference implementation of DNA/DNA duplex thermodynamics for PCR primers.

Pure Python (no numpy needed). Perfectly matched duplex of an oligo with its
complement (primer on its target site).

Model stack (default):
  * Nearest-neighbour (NN) parameters: SantaLucia and Hicks 2004
    (Annu Rev Biophys Biomol Struct 33:415-440), at 1 M NaCl.
    SantaLucia 1998 unified (PNAS 95:1460-1465) is available as an option.
  * Salt: Owczarzy et al. 2008 (Biochemistry 47:5336-5353): Eq 16 with the
    Table 2 coefficients (p5343), Eqs 18-20 for the mixed regime (p5347) and
    the Fig 9 decision tree on R = sqrt([Mg_free]) / [Mon] (p5348); Owczarzy
    et al. 2004 Eq 22 (restated as Eq 4, p5339) for the monovalent regime.
    Checked digit by digit against the paper and its worked examples.
  * Free Mg2+: 1:1 Mg:dNTP binding with Ka = 3e4 M^-1 (Owczarzy 2008 Eq 17,
    p5346) or simple subtraction Mg - dNTP (the paper's default, p5346, p5351;
    the two differ by about 0.03 C at PCR conditions).
  * Monovalent cations: [Mon] = [Na+] + [K+] + [Tris+] (Eq 3); Tris+ is about
    0.5 x total Tris at pH 8.3 (p5339: 10 mM Tris counts as about 5 mM).
  * Accuracy (out of sample, Owczarzy 2008 SI Table S2, 32 duplexes <= 30 bp):
    salt correction from a measured 1 M Tm MAE 0.90 C; full prediction MAE
    2.4 C (1.6 C, bias +1.5 C, in PCR-like buffer); the error is nearest-
    neighbour error, largest for homopolymers and pyrimidine-rich strands.
  * Template dangling ends (SantaLucia and Hicks 2004 Table 3, p426):
    dangling_end_terms() for a primer bound inside a longer template.
  * Additives: DMSO linear Tm shift (default 0.6 C per % v/v, primer3
    default, Musielski 1981; von Ahsen 2001 propose 0.75), formamide
    Blake and Delcourt 1996: dTm = (0.453 fGC - 2.88) per mol/L.
  * Salt and additive effects are applied to Tm; an effective entropy
    dS_eff = dH / Tm_final - R ln(C_eff) is then derived (dH assumed salt
    and additive independent, as in SantaLucia 1998 and Owczarzy 2004/2008),
    so dG(T) and the fraction bound are consistent with the reported Tm.

Units: dH in kcal/mol, dS in cal/(K mol), dG in kcal/mol, T in C unless a name
ends in _K. Concentrations are given in mM (ions) and nM (oligos).

Run `python thermo_reference.py` to execute the self-tests.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace

R_CAL = 1.98720425864083  # gas constant, cal/(K mol) (CODATA 8.314462618 J/(K mol) / 4.184)
T0 = 273.15

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}

# ---------------------------------------------------------------------------
# NN parameters at 1 M NaCl. key = 5'-XY-3' on the top strand / its complement
# (dH kcal/mol, dS cal/(K mol)). The 10 unique WC stacks.
# ---------------------------------------------------------------------------
_NN_2004 = {
    "AA": (-7.6, -21.3), "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "GT": (-8.4, -22.4), "CT": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "CG": (-10.6, -27.2), "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9),
}
# SantaLucia 1998 unified differs only in AA/TT and in how initiation is written.
_NN_1998 = dict(_NN_2004)
_NN_1998["AA"] = (-7.9, -22.2)

PARAMS = {
    # init: applied once per duplex; term_AT: per terminal A-T pair;
    # term_GC: per terminal G-C pair; sym: self-complementary duplex.
    "santalucia2004": {"nn": _NN_2004, "init": (0.2, -5.7), "term_AT": (2.2, 6.9),
                        "term_GC": (0.0, 0.0), "sym": (0.0, -1.4)},
    # 1998 Table 2: "initiation with terminal G.C" 0.1/-2.8 and "with terminal
    # A.T" 2.3/4.1, one term for each duplex end.
    "santalucia1998": {"nn": _NN_1998, "init": (0.0, 0.0), "term_AT": (2.3, 4.1),
                        "term_GC": (0.1, -2.8), "sym": (0.0, -1.4)},
}


def _canon(seq: str) -> str:
    s = seq.upper().replace("U", "T").strip()
    bad = set(s) - set("ACGT")
    if bad:
        raise ValueError(f"only A, C, G, T allowed, got {sorted(bad)}")
    if len(s) < 2:
        raise ValueError("sequence must be at least 2 nt")
    return s


def revcomp(seq: str) -> str:
    return "".join(COMP[b] for b in reversed(seq))


def _stack(nn: dict, xy: str):
    if xy in nn:
        return nn[xy]
    return nn[revcomp(xy)]  # e.g. TT/AA == AA/TT


def nn_dh_ds(seq: str, params: str = "santalucia2004"):
    """dH (kcal/mol) and dS (cal/K/mol) at 1 M Na+ for seq + its perfect complement."""
    s = _canon(seq)
    p = PARAMS[params]
    dh, ds = p["init"]
    for end in (s[0], s[-1]):
        h, e = p["term_AT"] if end in "AT" else p["term_GC"]
        dh += h
        ds += e
    for i in range(len(s) - 1):
        h, e = _stack(p["nn"], s[i:i + 2])
        dh += h
        ds += e
    if s == revcomp(s):
        dh += p["sym"][0]
        ds += p["sym"][1]
    return dh, ds


def gc_fraction(seq: str) -> float:
    s = _canon(seq)
    return (s.count("G") + s.count("C")) / len(s)


# ---------------------------------------------------------------------------
# Cation bookkeeping
# ---------------------------------------------------------------------------

def free_mg_mM(mg_mM: float, dntp_mM: float, method: str = "ka", ka: float = 3.0e4) -> float:
    """Free Mg2+ (mM). method 'ka': 1:1 binding, Ka in M^-1. 'subtract': Mg - dNTP."""
    if mg_mM <= 0:
        return 0.0
    if dntp_mM <= 0:
        return mg_mM
    if method == "subtract":
        return max(mg_mM - dntp_mM, 0.0)
    mg, dn = mg_mM * 1e-3, dntp_mM * 1e-3
    b = ka * dn - ka * mg + 1.0
    return 1e3 * (-b + math.sqrt(b * b + 4.0 * ka * mg)) / (2.0 * ka)


def monovalent_mM(na_mM: float = 0.0, k_mM: float = 0.0, tris_mM: float = 0.0,
                  tris_frac_protonated: float = 0.5) -> float:
    return na_mM + k_mM + tris_frac_protonated * tris_mM


# ---------------------------------------------------------------------------
# Salt corrections. Each returns Tm (K) at the target condition from Tm at 1 M Na+.
# ---------------------------------------------------------------------------

def owczarzy2004_inv_corr(fgc: float, mon_M: float) -> float:
    """Owczarzy 2004 Eq 22: 1/Tm(Na) - 1/Tm(1 M)."""
    ln = math.log(mon_M)
    return (4.29 * fgc - 3.95) * 1e-5 * ln + 9.40e-6 * ln * ln


def owczarzy2008_inv_corr(fgc: float, n_bp: int, mg_M: float, mon_M: float) -> tuple[float, str]:
    """Owczarzy 2008 Eq 16 with decision tree. Returns (1/Tm - 1/Tm(1M), regime)."""
    if mg_M <= 0:
        if mon_M <= 0:
            raise ValueError("no cations")
        return owczarzy2004_inv_corr(fgc, mon_M), "mono(no Mg)"
    a, b, c, d = 3.92e-5, -9.11e-6, 6.26e-5, 1.42e-5
    e, f, g = -4.82e-4, 5.25e-4, 8.31e-5
    regime = "Mg(R>=6)"
    if mon_M > 0:
        r = math.sqrt(mg_M) / mon_M
        if r < 0.22:
            return owczarzy2004_inv_corr(fgc, mon_M), "mono(R<0.22)"
        if r < 6.0:
            lm = math.log(mon_M)
            a = 3.92e-5 * (0.843 - 0.352 * math.sqrt(mon_M) * lm)
            d = 1.42e-5 * (1.279 - 4.03e-3 * lm - 8.03e-3 * lm * lm)
            g = 8.31e-5 * (0.486 - 0.258 * lm + 5.25e-3 * lm ** 3)
            regime = "mixed(0.22<=R<6)"
    lmg = math.log(mg_M)
    corr = (a + b * lmg + fgc * (c + d * lmg)
            + (1.0 / (2.0 * (n_bp - 1))) * (e + f * lmg + g * lmg * lmg))
    return corr, regime


# ---------------------------------------------------------------------------
# Main API
# ---------------------------------------------------------------------------

@dataclass
class Conditions:
    na_mM: float = 0.0
    k_mM: float = 50.0
    tris_mM: float = 0.0          # total Tris; 0.5 counted as Tris+
    mg_mM: float = 2.0            # total Mg2+
    dntp_mM: float = 0.8          # total dNTP (sum of the four)
    primer_nM: float = 200.0      # each primer
    target_nM: float = 0.0        # template strand; 0 = primer in large excess
    dmso_pct: float = 0.0         # % v/v
    dmso_factor: float = 0.6      # C per % v/v
    formamide_M: float = 0.0
    salt: str = "owczarzy2008"    # 'owczarzy2008' | 'owczarzy2004' | 'santalucia1998' | 'none'
    free_mg_method: str = "ka"    # 'ka' | 'subtract'
    params: str = "santalucia2004"

    @property
    def mon_mM(self) -> float:
        return monovalent_mM(self.na_mM, self.k_mM, self.tris_mM)


def c_eff_M(seq: str, cond: Conditions) -> float:
    """Concentration term in Tm = dH/(dS + R ln C_eff).

    Non-self-complementary, primer A0 >= target B0: C_eff = A0 - B0/2
    (equal strands: CT/4 with CT = 2*A0). Self-complementary: C_eff = CT.
    """
    s = _canon(seq)
    a0 = cond.primer_nM * 1e-9
    b0 = cond.target_nM * 1e-9
    if s == revcomp(s):
        return a0  # single species, CT = A0
    hi, lo = max(a0, b0), min(a0, b0)
    return hi - lo / 2.0


@dataclass
class Result:
    seq: str
    dH: float           # kcal/mol (1 M Na, salt independent)
    dS_1M: float        # cal/(K mol) at 1 M Na
    dS: float           # effective cal/(K mol) at the stated conditions
    tm_1M_C: float
    tm_C: float
    regime: str

    def dG(self, t_C: float) -> float:
        return self.dH - (t_C + T0) * self.dS / 1000.0


def analyze(seq: str, cond: Conditions | None = None) -> Result:
    cond = cond or Conditions()
    s = _canon(seq)
    dh, ds = nn_dh_ds(s, cond.params)
    ce = c_eff_M(s, cond)
    rlnc = R_CAL * math.log(ce)
    tm1_K = dh * 1000.0 / (ds + rlnc)
    fgc = gc_fraction(s)
    n = len(s)
    mg_free = free_mg_mM(cond.mg_mM, cond.dntp_mM, cond.free_mg_method)
    mon = cond.mon_mM
    if cond.salt == "owczarzy2008":
        corr, regime = owczarzy2008_inv_corr(fgc, n, mg_free * 1e-3, mon * 1e-3)
        tm_K = 1.0 / (1.0 / tm1_K + corr)
    elif cond.salt == "owczarzy2004":
        corr = owczarzy2004_inv_corr(fgc, mon * 1e-3)
        tm_K, regime = 1.0 / (1.0 / tm1_K + corr), "mono only"
    elif cond.salt == "santalucia1998":
        # von Ahsen 2001 Na-equivalent with Mg - dNTP subtraction, then
        # SantaLucia 1998 entropy correction 0.368 (N-1) ln[Na_eq].
        na_eq = mon + 120.0 * math.sqrt(max(cond.mg_mM - cond.dntp_mM, 0.0))
        ds_s = ds + 0.368 * (n - 1) * math.log(na_eq * 1e-3)
        tm_K, regime = dh * 1000.0 / (ds_s + rlnc), "Na_eq"
    elif cond.salt == "none":
        tm_K, regime = tm1_K, "1 M Na"
    else:
        raise ValueError(cond.salt)
    tm_K -= cond.dmso_factor * cond.dmso_pct
    tm_K += (0.453 * fgc - 2.88) * cond.formamide_M
    ds_eff = dh * 1000.0 / tm_K - rlnc
    return Result(s, dh, ds, ds_eff, tm1_K - T0, tm_K - T0, regime)


# SantaLucia and Hicks 2004 Table 3 (p426): terminal dangling ends next to a
# Watson-Crick pair, 1 M NaCl. Key (unpaired X, paired Y) with the unpaired base
# 5' of Y (5'-dangle, "XY/comp(Y)") or 3' of Y (3'-dangle, "YX/comp(Y)").
# Values (dH kcal/mol, dG37 kcal/mol).
_DANGLE5 = {  # 5'-XY-3' / comp(Y), X unpaired
    ("A", "A"): (0.2, -0.51), ("C", "A"): (0.6, -0.42), ("G", "A"): (-1.1, -0.62), ("T", "A"): (-6.9, -0.71),
    ("A", "C"): (-6.3, -0.96), ("C", "C"): (-4.4, -0.52), ("G", "C"): (-5.1, -0.72), ("T", "C"): (-4.0, -0.58),
    ("A", "G"): (-3.7, -0.58), ("C", "G"): (-4.0, -0.34), ("G", "G"): (-3.9, -0.56), ("T", "G"): (-4.9, -0.61),
    ("A", "T"): (-2.9, -0.50), ("C", "T"): (-4.1, -0.02), ("G", "T"): (-4.2, 0.48), ("T", "T"): (-0.2, -0.10),
}
_DANGLE3 = {  # 5'-YX-3' / comp(Y), X unpaired
    ("A", "A"): (-0.5, -0.12), ("C", "A"): (4.7, 0.28), ("G", "A"): (-4.1, -0.01), ("T", "A"): (-3.8, 0.13),
    ("A", "C"): (-5.9, -0.82), ("C", "C"): (-2.6, -0.31), ("G", "C"): (-3.2, -0.01), ("T", "C"): (-5.2, -0.52),
    ("A", "G"): (-2.1, -0.92), ("C", "G"): (-0.2, -0.23), ("G", "G"): (-3.9, -0.44), ("T", "G"): (-4.4, -0.35),
    ("A", "T"): (-0.7, -0.48), ("C", "T"): (4.4, -0.19), ("G", "T"): (-1.6, -0.50), ("T", "T"): (2.9, -0.29),
}


def dangling_end_terms(template_beyond_3p: str, template_at_3p: str,
                       template_at_5p: str, template_beyond_5p: str) -> tuple:
    """(dH kcal/mol, dS cal/K/mol) of the two template dangling ends around a bound
    primer (template bases, written 5'->3' on the TEMPLATE strand).

    Beyond the primer's 3' end the template continues toward its own 5' end, so
    that base is a 5'-dangle on the template (key (beyond, paired)); beyond the
    primer's 5' end it is a 3'-dangle (key (beyond, paired)). Missing or
    non-ACGT neighbours contribute nothing. Only the first dangling base counts
    (SantaLucia and Hicks 2004: further bases add little above 25 C)."""
    dh = ds = 0.0
    for table, x, y in ((_DANGLE5, template_beyond_3p, template_at_3p),
                        (_DANGLE3, template_beyond_5p, template_at_5p)):
        v = table.get((x, y))
        if v:
            h, g37 = v
            dh += h
            ds += (h - g37) * 1000.0 / 310.15
    return dh, ds


def fraction_bound(seq: str, t_C: float, cond: Conditions | None = None) -> float:
    """Fraction of the TARGET strand hybridised to the primer at t_C (exact 2-state).

    Uses K = exp(-dG(T)/RT) and A + B <-> AB with A0 = primer, B0 = target.
    If target_nM == 0 the primer-excess limit K A0 / (1 + K A0) is returned.
    For a self-complementary oligo the fraction of strands in duplex is returned.
    """
    cond = cond or Conditions()
    r = analyze(seq, cond)
    t_K = t_C + T0
    lnk = -r.dG(t_C) * 1000.0 / (R_CAL * t_K)
    a0 = cond.primer_nM * 1e-9
    b0 = cond.target_nM * 1e-9
    if r.seq == revcomp(r.seq):
        # 2A <-> A2, K = [A2]/[A]^2; fraction of strands in duplex
        k = math.exp(min(lnk, 700.0))
        free = (-1.0 + math.sqrt(1.0 + 8.0 * k * a0)) / (4.0 * k)
        return 1.0 - free / a0
    if b0 == 0:
        x = lnk + math.log(a0)
        return 1.0 / (1.0 + math.exp(-x)) if x > -700 else 0.0
    k = math.exp(min(lnk, 700.0))
    # K (A0 - x)(B0 - x) = x  ->  x^2 - (A0 + B0 + 1/K) x + A0 B0 = 0 (small root)
    s = a0 + b0 + 1.0 / k
    x = (s - math.sqrt(max(s * s - 4.0 * a0 * b0, 0.0))) / 2.0
    return x / b0


# ---------------------------------------------------------------------------
# NEB Tm calculator scale and helpers used by the designer
# ---------------------------------------------------------------------------

#: The NEB Tm calculator (Q5 setting) equals SantaLucia 1998 + Owczarzy 2004 at
#: 150 mM Na+, with C = primer concentration and Mg2+ ignored. Checked against
#: 240 NEB API values at 500 and 25 nM: max error 0.000 C (2026-11-11,
#: docs/research/reference_code/neb_calib.json). A vendor scale, not a physical Tm.
NEB_CONDITIONS = Conditions(na_mM=150.0, k_mM=0.0, tris_mM=0.0, mg_mM=0.0, dntp_mM=0.0,
                            primer_nM=500.0, target_nM=0.0, salt="owczarzy2004",
                            params="santalucia1998")


def tm_neb(seq: str, primer_nM: float = 500.0) -> float:
    """Tm on the NEB calculator scale (Q5), default 500 nM as NEB uses."""
    cond = NEB_CONDITIONS if primer_nM == 500.0 else replace(NEB_CONDITIONS, primer_nM=primer_nM)
    return analyze(seq, cond).tm_C


def theta_from_dg(dg_kcal: float, t_C: float, primer_nM: float) -> float:
    """Fraction of template sites bound by a primer in excess: 1 / (1 + exp(dG/RT) / C)."""
    x = dg_kcal * 1000.0 / (R_CAL * (t_C + T0)) - math.log(primer_nM * 1e-9)
    if x > 700:
        return 0.0
    if x < -700:
        return 1.0
    return 1.0 / (1.0 + math.exp(x))


def structure_cost(dg_kcal: float, t_C: float) -> float:
    """Free-energy cost (kcal/mol, >= 0) of a competing two-state structure with
    folding dG `dg_kcal`: RT ln(1 + exp(-dG/RT))."""
    rt = R_CAL * (t_C + T0) / 1000.0
    x = -dg_kcal / rt
    return rt * (x + math.log1p(math.exp(-x))) if x > 30 else rt * math.log1p(math.exp(x))


# ---------------------------------------------------------------------------
# Self-tests
# ---------------------------------------------------------------------------

def _close(a, b, tol, msg):
    assert abs(a - b) <= tol, f"{msg}: got {a}, expected {b} +/- {tol}"


def _self_test(verbose: bool = True) -> None:
    # 1. Owczarzy 2008 / US 8,055,451 Example 6: 25-mer, measured Tm 83.4 C in
    #    1 M Na+ (2 uM total strands), 0.5 mM MgCl2 + 10 mM Tris-HCl (5 mM Tris+).
    #    Published prediction 346.3 K = 73.2 C (measured 73.3 C).
    seq = "CAGTGGGCTCCTGGGCGTGCTGGTC"
    corr, regime = owczarzy2008_inv_corr(gc_fraction(seq), len(seq), 0.5e-3, 5e-3)
    tm = 1.0 / (1.0 / (83.4 + T0) + corr) - T0
    _close(gc_fraction(seq), 0.72, 1e-9, "fGC example 6")
    _close(tm, 73.2, 0.06, "Owczarzy 2008 worked example")
    assert regime.startswith("mixed"), regime

    # 1b. Owczarzy 2008 worked example (p5347-5348): AAGGCGAGTCAGGCTCAGTG,
    #     Tm(1 M Na+) 76.3 C, 1.5 mM Mg2+, 5 mM Tris+ (R = 7.75, Mg regime):
    #     correction 7.0537e-5 K^-1, predicted 67.9 C (measured 68.5 C).
    seq2 = "AAGGCGAGTCAGGCTCAGTG"
    corr2, regime2 = owczarzy2008_inv_corr(gc_fraction(seq2), len(seq2), 1.5e-3, 5e-3)
    _close(corr2, 7.0537e-5, 5e-9, "Owczarzy 2008 p5348 correction term")
    assert regime2 == "Mg(R>=6)", regime2
    _close(1.0 / (1.0 / (76.3 + T0) + corr2) - T0, 67.9, 0.06, "Owczarzy 2008 p5348 Tm")

    # 1c. SantaLucia and Hicks 2004 Table 3 extremes: AC/G -0.96, GT/A +0.48 (5'-dangles).
    _close(_DANGLE5[("A", "C")][1], -0.96, 1e-12, "dangle AC/G")
    _close(_DANGLE5[("G", "T")][1], 0.48, 1e-12, "dangle GT/A")

    # 2. Intermediate coefficients of Example 6 (a=3.82e-5, d=1.53e-5, g=8.91e-5 at Mon=5 mM).
    lm = math.log(0.005)
    _close(3.92e-5 * (0.843 - 0.352 * math.sqrt(0.005) * lm), 3.82e-5, 0.006e-5, "a(Mon)")
    _close(1.42e-5 * (1.279 - 4.03e-3 * lm - 8.03e-3 * lm * lm), 1.53e-5, 0.006e-5, "d(Mon)")
    _close(8.31e-5 * (0.486 - 0.258 * lm + 5.25e-3 * lm ** 3), 8.91e-5, 0.006e-5, "g(Mon)")

    # 3. 1998 vs 2004 parameter sets give the same dG37 within 0.1 kcal/mol for
    #    a sequence without AA/TT stacks, and SantaLucia 1998 Table 2 values.
    dh, ds = nn_dh_ds("CGTGACGTGACGGACT", "santalucia1998")
    dh4, ds4 = nn_dh_ds("CGTGACGTGACGGACT", "santalucia2004")
    _close(dh - 310.15 * ds / 1000, dh4 - 310.15 * ds4 / 1000, 0.1, "dG37 1998 vs 2004")

    # 4. Self-complementary symmetry handling.
    assert revcomp("CGCGAATTCGCG") == "CGCGAATTCGCG"
    dh_s, ds_s = nn_dh_ds("CGCGAATTCGCG")
    #    hand sum: CG GC CG GA AA AT TT TC CG GC CG + init
    hand_h = (-10.6 - 9.8 - 10.6 - 8.2 - 7.6 - 7.2 - 7.6 - 8.2 - 10.6 - 9.8 - 10.6) + 0.2
    hand_s = (-27.2 - 24.4 - 27.2 - 22.2 - 21.3 - 20.4 - 21.3 - 22.2 - 27.2 - 24.4 - 27.2) - 5.7 - 1.4
    _close(dh_s, hand_h, 1e-9, "dH Dickerson")
    _close(ds_s, hand_s, 1e-9, "dS Dickerson")

    # 5. Free Mg: Ka model, 2 mM Mg + 0.8 mM dNTP.
    fm = free_mg_mM(2.0, 0.8, "ka")
    assert 1.2 < fm < 1.25, fm
    _close(free_mg_mM(3.0, 0.8, "subtract"), 2.2, 1e-12, "subtract (patent Ex. 7)")

    # 6. Fraction bound is 0.5 at Tm by construction (primer excess and equal strands).
    c = Conditions()
    r = analyze("TACTTCCAGTGCTCAGCGTA", c)
    _close(fraction_bound("TACTTCCAGTGCTCAGCGTA", r.tm_C, c), 0.5, 1e-6, "theta(Tm) excess")
    c2 = Conditions(primer_nM=200, target_nM=200)
    r2 = analyze("TACTTCCAGTGCTCAGCGTA", c2)
    _close(fraction_bound("TACTTCCAGTGCTCAGCGTA", r2.tm_C, c2), 0.5, 1e-6, "theta(Tm) equal")
    # C_eff 200 nM (excess) vs 100 nM (equal strands, CT/4 with CT = 400 nM): ~1 C
    _close(r.tm_C - r2.tm_C, 1.0, 0.3, "excess vs CT/4 shift")
    assert fraction_bound("TACTTCCAGTGCTCAGCGTA", r.tm_C - 5, c) > 0.9

    # 7. Regimes: no Mg -> Owczarzy 2004; huge Mg with little K -> Mg regime.
    assert analyze("ACGTACGTACGTACGTAC", Conditions(mg_mM=0, dntp_mM=0)).regime.startswith("mono")
    assert analyze("ACGTACGTACGTACGTAC", Conditions(k_mM=1, mg_mM=50, dntp_mM=0)).regime == "Mg(R>=6)"

    # 8. The published model is piecewise and NOT continuous at R = 0.22 and
    #    R = 6. For PCR-like monovalent levels (10-100 mM) the jump is < 1 C
    #    for a 20-mer; at Mon <= 5 mM it can reach 3-16 C (outside PCR use).
    def tm_at(mon_mM, mg_mM):
        return analyze("TACTTCCAGTGCTCAGCGTA", Conditions(k_mM=mon_mM, mg_mM=mg_mM, dntp_mM=0)).tm_C
    for mon in (10.0, 20.0, 50.0, 100.0):
        for rr in (0.22, 6.0):
            mg_x = (rr * mon * 1e-3) ** 2 * 1e3
            _close(tm_at(mon, mg_x * 0.9999), tm_at(mon, mg_x * 1.0001), 1.0, f"jump R={rr} Mon={mon}")

    if verbose:
        print("all self-tests passed")
        print(r)
        print("dG(60 C) =", round(r.dG(60.0), 3), "kcal/mol; theta(60 C) =",
              round(fraction_bound(r.seq, 60.0, c), 4))


if __name__ == "__main__":
    _self_test()
