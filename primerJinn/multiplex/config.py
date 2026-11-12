"""Design settings, polymerase presets, adapter tails and GC-dependent defaults.

Every number here is explained in docs/DESIGN.md.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Tuple

from . import thermo

TAILS = {
    "none": ("", ""),
    # Illumina TruSeq-style partial adapters, as in the original PrimerJinn
    "illumina": ("ACACTCTTTCCCTACACGACGCTCTTCCGATCT", "GACTGGAGTTCAGACGTGTGCTCTTCCGATCT"),
    # Illumina Nextera / 16S amplicon protocol overhangs
    "nextera": ("TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG", "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"),
}

#: mg_mM / dntp_mM: published final concentrations of the 1x reaction. Monovalent
#: content of these buffers is not published; 50 mM is assumed (override --mono_mM).
#: tm_scale: 'physical' = SantaLucia 2004 + Owczarzy 2008 in the stated buffer at the
#: per-primer concentration (right for Taq-type enzymes); 'neb' = the NEB Tm
#: calculator's Q5 scale (vendor-calibrated, used only for Q5).
POLYMERASES = {
    # KAPA2G Robust HotStart ReadyMix (Roche/KAPA TDS KR0382 v3.16): engineered Taq-type
    # enzyme, no 3'->5' proofreading, A-tailed products, 2 mM MgCl2 and 0.2 mM each dNTP
    # at 1x, anneal 55/60/65 C for primer Tm <55 / 55-65 / >65 C for 15 s, extend 72 C
    # 15 s/kb, denature 95 C (5 min initial for >70% GC), 5% DMSO for amplicons >65% GC.
    "kapa": dict(mg_mM=2.0, dntp_mM=0.8, proofreading=False, ta_rule="model", tm_scale="physical",
                 label="KAPA2G Robust HotStart ReadyMix"),
    "taq": dict(mg_mM=1.5, dntp_mM=0.8, proofreading=False, ta_rule="model", tm_scale="physical",
                label="standard Taq (1.5 mM MgCl2)"),
    "q5": dict(mg_mM=2.0, dntp_mM=0.8, proofreading=True, ta_rule="neb", tm_scale="neb",
               label="NEB Q5"),
    "hifi": dict(mg_mM=2.0, dntp_mM=0.8, proofreading=True, ta_rule="model", tm_scale="physical",
                 label="proofreading high-fidelity enzyme (KAPA HiFi, SuperFi II, Phusion-like)"),
}


def target_tm_from_gc(gc: float) -> float:
    """The user's rule on the NEB scale: about 65 C at 65% genome GC, 55 to 60 C otherwise."""
    if gc >= 0.60:
        return 65.0
    if gc >= 0.55:
        return 60.0 + (gc - 0.55) / 0.05 * 5.0
    if gc >= 0.30:
        return 60.0
    if gc >= 0.20:
        return 55.0 + (gc - 0.20) / 0.10 * 5.0
    return 55.0


def primer_length_range(gc: float) -> Tuple[int, int]:
    if gc < 0.30:
        return 20, 40
    if gc < 0.40:
        return 18, 34
    if gc < 0.60:
        return 18, 32
    return 17, 32


@dataclass
class Config:
    # targets and amplicons
    amp_min: int = 225
    amp_opt: int = 300
    amp_max: int = 375
    margin: int = 3                 # clear bases between a primer 3' end and the first required base
    # primers
    tm: Optional[float] = None      # target primer Tm on tm_scale; None = from genome GC
    tm_scale: Optional[str] = None  # 'physical' or 'neb'; None = polymerase preset
    tm_window: float = 1.5
    primer_min: Optional[int] = None
    primer_max: Optional[int] = None
    max_poly: Optional[int] = None  # longest allowed mononucleotide run; None = 4 (5 below 25% GC)
    primer_gc_min: Optional[float] = None  # None = min(0.20, max(0.10, genome GC - 0.08))
    tails: str = "none"
    tail_fwd: str = ""
    tail_rev: str = ""
    # reaction
    polymerase: str = "kapa"
    mono_mM: float = 50.0
    mg_mM: Optional[float] = None   # None = polymerase preset
    dntp_mM: float = 0.8            # total of the four dNTPs
    primer_nM: float = 100.0        # each primer in the pool
    dmso_pct: float = 0.0
    ta: Optional[float] = None      # design annealing temperature; None = derived
    balance: bool = True            # extend weak primers at the 5' end (lowercase) to even out binding
    balance_tol: float = 0.3        # kcal/mol weaker than the pool median before a primer is extended
    max_ext: int = 8                # most bases added to one primer
    # scoring thresholds (kcal/mol); uncalibrated starting points, see DESIGN.md
    hairpin_penalty: float = 0.5    # structure cost of the primer's own hairpin at Ta
    hairpin_reject: float = 1.0
    site_open_penalty: float = 1.5
    site_open_reject: float = 3.0
    d3_soft: float = -4.0
    d3_hard: float = -6.0
    dany_soft: float = -8.0
    dany_hard: float = -12.0
    self_prime_flag: float = 0.2
    # specificity
    specificity: bool = True
    offtarget_max_product: int = 2000
    # pools and search
    outer_min: int = 1000          # designed amplicons closer than this (outer product) go to separate pools
    candidates_per_target: int = 40
    sa_restarts: int = 12
    seed: int = 1
    threads: int = 0
    # filled in by resolve()
    genome_gc: float = field(default=0.5, repr=False)

    # ------------------------------------------------------------------
    def resolve(self, genome_gc: float) -> "Config":
        """Fill GC-dependent defaults. Returns self."""
        self.genome_gc = genome_gc
        if self.tm is None:
            self.tm = target_tm_from_gc(genome_gc)
        if self.max_poly is None:
            self.max_poly = 5 if genome_gc < 0.25 else 4
        if self.primer_gc_min is None:
            self.primer_gc_min = min(0.20, max(0.10, genome_gc - 0.08))
        lo, hi = primer_length_range(genome_gc)
        self.primer_min = self.primer_min or lo
        self.primer_max = self.primer_max or hi
        if self.polymerase not in POLYMERASES:
            raise ValueError(f"unknown polymerase '{self.polymerase}', choose from {sorted(POLYMERASES)}")
        if self.mg_mM is None:
            self.mg_mM = POLYMERASES[self.polymerase]["mg_mM"]
        if self.tm_scale is None:
            self.tm_scale = POLYMERASES[self.polymerase]["tm_scale"]
        if self.tm_scale not in ("physical", "neb"):
            raise ValueError("tm_scale must be 'physical' or 'neb'")
        if self.tails not in TAILS and self.tails != "custom":
            raise ValueError(f"unknown tails '{self.tails}'")
        self.tail_fwd = (self.tail_fwd or "").upper()
        self.tail_rev = (self.tail_rev or "").upper()
        if self.tails == "custom":
            if not (self.tail_fwd or self.tail_rev):
                raise ValueError("--tails custom needs --tail_fwd and/or --tail_rev")
            for nm, tl in (("--tail_fwd", self.tail_fwd), ("--tail_rev", self.tail_rev)):
                if set(tl) - set("ACGT"):
                    raise ValueError(f"{nm} may contain only A, C, G and T")
        else:
            if self.tail_fwd or self.tail_rev:
                raise ValueError("--tail_fwd/--tail_rev need --tails custom")
            self.tail_fwd, self.tail_rev = TAILS[self.tails]
        if not (self.amp_min <= self.amp_opt <= self.amp_max):
            raise ValueError("need amp_min <= amp_opt <= amp_max")
        if self.ta is None:
            self.ta = self.design_ta()
        import os
        self.threads = self.threads or max(1, (os.cpu_count() or 2) - 1)
        return self

    @property
    def proofreading(self) -> bool:
        return POLYMERASES[self.polymerase]["proofreading"]

    def design_ta(self) -> float:
        """Temperature at which binding, structure and specificity are scored.
        Q5: NEB rule, lowest Tm + 1 C, so about the target Tm. Taq-type enzymes on the
        physical scale: primers are about 90% bound near Tm - 3 C and template
        structure costs about 1 C more, so target - 4 C (refined after design by
        the structure-aware equilibrium model)."""
        if POLYMERASES[self.polymerase]["ta_rule"] == "neb":
            # NEB rule at the real primer concentration: offset of the 500 nM design
            # scale from the reaction concentration, for a 20-mer at genome GC
            n_gc = int(round(20 * min(max(self.genome_gc, 0.2), 0.8)))
            ref = ("GC" * 10)[:n_gc] + ("AT" * 10)[:20 - n_gc]
            off = thermo.tm_neb(ref, 500.0) - thermo.tm_neb(ref, self.primer_nM)
            return min(self.tm + 1.0 - self.tm_window - off - 0.6 * self.dmso_pct, 72.0)
        return self.tm - 4.0

    def primer_tm(self, seq: str) -> float:
        """Primer Tm on the design scale."""
        if self.tm_scale == "neb":
            return thermo.tm_neb(seq)
        return thermo.analyze(seq, self.phys()).tm_C

    def phys(self, primer_nM: Optional[float] = None) -> thermo.Conditions:
        return thermo.Conditions(k_mM=self.mono_mM, mg_mM=self.mg_mM, dntp_mM=self.dntp_mM,
                                 primer_nM=primer_nM or self.primer_nM, dmso_pct=self.dmso_pct)

    def thal_kw(self, temp_c: float) -> dict:
        return dict(mv_conc=self.mono_mM, dv_conc=self.mg_mM, dntp_conc=self.dntp_mM,
                    dna_conc=self.primer_nM, temp_c=temp_c)

    def struct_cond(self):
        from .structure import Conditions
        return Conditions(mono_mM=self.mono_mM, mg_mM=self.mg_mM, dntp_mM=self.dntp_mM,
                          primer_nM=self.primer_nM, dmso_pct=self.dmso_pct)
