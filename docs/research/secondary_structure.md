# Template and primer secondary structure for multiplex primer choice and annealing temperature

Scope: PrimerJinn multiplex rebuild. Topic: use single-stranded DNA structure of the template (genomic DNA in early cycles, amplicon strands later) and of the primers to choose primers and recommend the annealing temperature (Ta) and additives. Test organism: *M. tuberculosis* H37Rv, NC_000962.3, from `example/ref.fasta` (measured here: 4,411,532 bp, 65.61% GC).

Everything numeric below was either run here (scripts and outputs listed in section 6) or read from the cited primary source. Where I could not verify something, it says so.

Units: kcal/mol, degrees C, mol/L (M) or mM. "Na_eq" means sodium-equivalent monovalent salt.

---

## Key decisions (short version)

1. Engine: ViennaRNA 2.7.0 with the **Mathews1999 DNA** parameter set, `noGU=1`, `dangles=2`, `salt = Na_eq` (M). Do not use the bundled Mathews2004 DNA set: it becomes more stable as temperature rises (measured below). NUPACK 4 is an optional cross-check only (academic licence). Do not use seqfold (35 s per 300 nt fold, MFE only, 9 to 15 kcal/mol weaker than both other tools).
2. Site metric: `dG_open = -RT ln P(site unpaired)` from `fold_compound(..., OPTION_WINDOW).probs_window(L, PROBS_WINDOW_UP, cb)`. Effective binding `dG_eff = dG_duplex + dG_open(site) + dG_open(primer self-structure)`, fraction bound `f = K[P]/(1+K[P])`.
3. Fold the amplicon strand (later cycles dominate yield) and the genomic context with RNAplfold-style W=200, L=150 (early cycles). Use the worse of the two.
4. Thresholds (soft, to be calibrated): penalise `dG_open(site) > 1.5`, reject `> 3.0` kcal/mol at Ta; require an accessible 8-nt toehold `P >= 0.1`; Ta = highest temperature at which every primer keeps `f_bound >= 0.9`, minus 1 C.
5. GC-rich advice is rule-based: if any amplicon has duplex Tm > 92 C (primer3 product-Tm formula) or a 50-bp window >= 80% GC, recommend 98 C denaturation, longer denaturation, and 1 to 2 M betaine or 3 to 5% DMSO, and recompute with the DMSO Tm shift.
6. Be honest in the UI: model agreement on *site* accessibility is only moderate (Vienna vs NUPACK Spearman 0.74 to 0.83, per-site factor up to about 5). Use structure metrics to rank and filter, not as quantitative efficiency predictions.

---

## Part A. Evidence

### A1. Template structure at the primer site slows or blocks hybridisation

| Source | What was measured | Quantitative result |
|---|---|---|
| Gao, Wolf, Georgiadis 2006, NAR 34:3370. doi:10.1093/nar/gkl422 (PMC1488884, full text read) | 25-mer probe/target pairs with designed intramolecular structure, UV and SPR kinetics, 0.5 M NaCl | Solution k_on 12.0e5, 7.2e5, 3.0e5, 2.0e5 M^-1 s^-1 for 0, 3+3, 4+3, 4+4 intramolecular bp. "addition of three intramolecular base pairs ... slow hybridization by a factor of two" (about -6 kcal/mol combined stability). With >= 4 bp per strand, kinetics are no longer two-state. |
| Mann et al. 2009 (Pythia), NAR 37:e95. doi:10.1093/nar/gkp443 (PMC2715258; local copy `dl/pythia.txt`) | Primer design by chemical-equilibrium analysis of primer binding, primer folding, dimers, and template folding in the priming region | Coverage in RepeatMasked human sequence: median 89% vs 51% for Primer3; recall 97% vs 48% at 81% sensitivity. Shows that an equilibrium model with folding terms is a workable design criterion. It does not isolate the effect of template structure alone. |
| Gimpel et al. 2025, Nat Commun 16:9187. doi:10.1038/s41467-025-64221-4 (PMC12533003, full text read) | Synthetic multi-template amplicon pools, 1D-CNN on per-sequence efficiency | GC content and base composition did not explain poor amplifiers. The main mechanism was a motif near the 3' adapter complementary to the 5' adapter, giving **hairpin formation and self-priming** that "inhibits primer annealing". AUROC 0.88. |
| Benita et al. 2003, NAR 31:e99. doi:10.1093/nar/gng101 | 1,438 human exons, logistic regression on 27 features | Template GC was the main predictor; regionalised GC (21-bp sliding window, 61% GC threshold) gave specificity 84.3%, sensitivity 94.8%. Primer features mostly did not matter. |

What I did **not** find: a primary study giving a validated numeric rule like "template structure at the site with dG below X kcal/mol at Ta causes dropout". The commonly quoted thresholds (for example dG of template hairpins more negative than about -3 kcal/mol at Ta) appear in vendor notes and reviews, not in a validated dataset I could check. The thresholds in Part C are therefore derived from physics (fraction bound) and must be calibrated against PrimerJinn's own amplicon depth data.

### A2. GC-rich templates, incomplete denaturation, and additives

| Source | Finding (verified from text or abstract) |
|---|---|
| Aird et al. 2011, Genome Biol 12:R18. doi:10.1186/gb-2011-12-2-r18 (PMC3188800, full text read) | Fast ramping left "insufficient time above a critical threshold temperature, causing incomplete denaturation". Standard Phusion PCR depleted >70% GC fragments to <1/100 of mid-GC. Extending denaturation (30 s to 3 min initial, 10 s to 80 s per cycle) plus **2 M betaine** rescued high GC: 50-bp windows up to 94% GC got > half the mean coverage. Betaine caused slight loss at 10 to 40% GC. AccuPrime Taq HiFi with long denaturation and anneal/extend at 65 C gave the most even profile. |
| Rees et al. 1993, Biochemistry 32:137. doi:10.1021/bi00052a019 | Betaine is "isostabilizing" (AT and GC pairs equally stable) at about **5.2 M**; it does not change B-form or polyelectrolyte behaviour. At PCR concentrations (1 to 2 M) it only partly reduces the GC dependence. |
| Henke et al. 1997, NAR 25:3957. doi:10.1093/nar/25.19.3957 | Betaine improved amplification of GC-rich cDNA (PSMA variants, c-jun), attributed to reduced secondary structure. |
| Musso et al. 2006, J Mol Diagn 8:544. doi:10.2353/jmoldx.2006.060058 | Betaine + DMSO + 7-deaza-dGTP together were "essential" for 67 to 79% GC targets (RET, PHOX2B, HOXD). Exact concentrations: not verified (full text not retrievable as text). |
| McConlogue, Brow, Innis 1988, NAR 16:9869. doi:10.1093/nar/16.20.9869 | Original "structure-independent" PCR with 7-deaza-dGTP (abstract not available; title only verified). |
| Frey et al. 2008, Nat Protoc 3:1312. doi:10.1038/nprot.2008.112 | "Slowdown PCR" for >83% GC: 7-deaza-dGTP, ramp 2.5 C/s, cooling 1.5 C/s to Ta, 48 cycles. |
| Chester and Marshak 1993, Anal Biochem 209:284. doi:10.1006/abio.1993.1121 | DMSO lowers primer Tm and a DMSO gradient behaves as a temperature gradient. Per-percent value not in abstract. |
| primer3 manual, PRIMER_DMSO_FACTOR (local `dl/p3manual.htm`) | Tm(DMSO) = Tm - factor x %DMSO; default 0.6 C/% (Musielski 1981), alternatives 0.75 (von Ahsen 2001), 0.5 (Cullen 1976), 0.675 (Escara 1980). I read this in the manual, not the primary papers. |

Amplicon duplex melting (computed here, primer3 `long_seq_tm` formula from `src_p3/oligotm.c`, PCR buffer 50 mM K+, 2 mM Mg2+, 0.8 mM dNTP, so Na_eq = 50 + 120 x sqrt(1.2) = 181 mM):

```
Tm = 81.5 + 16.6 log10(Na_eq/1000) + 41 x GC_fraction - 600/N - 0.6 x %DMSO
rpoB RRDR 300 bp (66.7% GC): 94.5 C     katG S315 300 bp (64.3% GC): 93.6 C
```

At a 95 C denaturation step with a short hold, these amplicons are within 1 C of their midpoint, so a large part of each product may not fully melt; that matches Aird's mechanism. With 5% DMSO the estimate drops by 3 C. Caveat: this empirical formula is accurate only to a few C, and primer3-py 2.0.3 `calc_tm(..., max_nn_length=60)` returns about 2 C lower than the formula in the primer3 source for the same input (measured, cause not traced).

### A3. G-quadruplexes

| Source | Finding |
|---|---|
| Huppert and Balasubramanian 2005, NAR 33:2908. doi:10.1093/nar/gki609 (full text read) | Folding rule: d(G3+ N1-7 G3+ N1-7 G3+ N1-7 G3+) forms a quadruplex in 100 mM KCl. |
| Woodford, Howell, Usdin 1994, JBC 269:27029 (PMID 7929444); Weitzmann, Woodford, Usdin 1996, JBC 271:20958. doi:10.1074/jbc.271.34.20958 | K+-dependent DNA synthesis arrest on the G-rich template strand; arrest only when the G-rich strand is the template, eliminated by blocking guanine N7. |
| Perrone et al. 2017, Sci Rep 7:5743. doi:10.1038/s41598-017-05867-z (full text read) | Putative G4s in the Mtb genome; 45 genes with a putative G4 near the TSS (>= 3 G per island); four confirmed by CD in KCl. |

Relevance: PCR buffers contain about 50 mM K+, which stabilises G4. 7-deaza-dGTP removes N7 and weakens Hoogsteen pairing in the copies (not the original genomic template). Neither ViennaRNA's `md.gquad` (RNA-derived G4 energies) nor NUPACK models DNA G4 reliably, so treat G4 as a motif flag. The module scans both strands with the Huppert rule. Neither demo amplicon had a hit.

---

## Part B. Tools, verified by running them here

### B1. ViennaRNA 2.7.0 (Python `import RNA`)

Verified API facts (`ss3/api.py`, output reproduced):

```
RNA 2.7.0  nupack 4.0.2.0  seqfold 0.9.0
default md: T 37.0 salt 1.021 dangles 2 noGU 0 helical_rise 2.8 backbone_length 6.0
params_load_DNA_Mathews2004 -> 1
after M2004 load: noGU 0 helical_rise 3.4 backbone 6.76 salt 1.021
params_load_DNA_Mathews1999 -> 1
parameter_set strings exist: True True       (RNA.parameter_set_dna_mathews1999 / _2004)
params_load_from_string -> 1
noGU 0 MFE 65C -25.8 GT pairs 10
noGU 1 MFE 65C -22.29 GT pairs 0
Mg-like attributes on md: []
constraint: dG_open(last 22)= 1.322 P= 0.1399
probs_window: P(last 22 unpaired)= 0.1399  dG= 1.322
pfl_fold_up ... P22 at end 0.0013 (uses global md defaults: noGU,T)
Vienna mfe+pf+bpp 300nt: 0.086 s
Vienna probs_window up to 25 nt, 300nt: 0.093 s
```

Working recipe:

```python
import RNA
RNA.params_load_DNA_Mathews1999()          # returns 1; also sets helical_rise 3.4, backbone 6.76
md = RNA.md()
md.temperature = 65.0                      # C
md.noGU = 1                                # DNA load does NOT set this
md.dangles = 2
md.salt = 0.05 + 3.3 * (0.0012 ** 0.5)     # mol/L, monovalent only; reference 1.021 M
fc = RNA.fold_compound(seq, md)
ss, mfe = fc.mfe(); fc.exp_params_rescale(mfe); _, G_ens = fc.pf()
bpp = fc.bpp()                             # 1-based upper triangle
# p_unpaired[i] = 1 - sum_j bpp[min(i,j)][max(i,j)]
md.window_size = len(seq); md.max_bp_span = len(seq)
fw = RNA.fold_compound(seq, md, RNA.OPTION_WINDOW)
def cb(v, vsize, i, maxsize, what, data):
    if what & RNA.PROBS_WINDOW_UP: data[i] = list(v)   # data[i][L] = P(i-L+1..i unpaired)
up = {}; fw.probs_window(22, RNA.PROBS_WINDOW_UP, cb, up)
# equivalent definition: dG_open = G(fc with hc_add_up(k, CONSTRAINT_CONTEXT_ALL_LOOPS) for k in site) - G_ens
```

Pitfalls found:
- `RNA.pfl_fold_up()` ignores your `md` and uses the global default model (37 C), giving P = 0.0013 instead of 0.14 above. Use `probs_window`.
- `fold_compound.eval_structure()` **segfaults** on an unbalanced dot-bracket string (earlier script `ss2/t3.py`, exit 139). Validate structures first.
- Salt: `md.salt` is monovalent only, in M, default 1.021. The C source (`src/ViennaRNA/params/salt.c`, read on GitHub) corrects stacking, loops, multiloops and duplex initiation relative to the default salt and cites Chen and Znosko 2013 for the duplex-initiation fit. No Mg2+ model. Salt matters a lot: rpoB window at 65 C, G_ens -33.9 kcal/mol at 1.021 M, -24.1 at 0.2 M, -13.7 at 0.05 M.

**Mathews2004 DNA set is broken above about 55 C** (`ss2/t1.py`, rpoB window NC_000962.3:761006-761305, dangles 2, 1.021 M):

| T (C) | 37 | 50 | 55 | 60 | 65 | 70 | 80 | 90 | 99 |
|---|---|---|---|---|---|---|---|---|---|
| MFE, Vienna M2004 | -74.0 | -52.2 | -45.5 | -43.9 | -44.9 | -48.0 | -64.8 | -96.3 | -129.0 |
| MFE, Vienna M1999 | -68.0 | -43.1 | -35.0 | -27.7 | -22.3 | -17.3 | -10.9 | -8.1 | -10.2 |
| MFE, NUPACK dna04 stacking 1 M | -74.2 | -46.1 | -36.2 | -27.3 | -22.4 | -17.0 | -8.9 | -3.8 | -1.2 |
| G_ens, Vienna M1999 | -76.0 | -53.8 | -46.4 | -39.7 | -33.9 | -29.0 | -22.6 | -20.7 | -23.4 |
| G_ens, NUPACK stacking | -79.8 | -54.4 | -46.1 | -39.5 | -34.0 | -29.4 | -22.0 | -17.2 | -14.8 |

Cause, isolated with a 1x2 interior loop (`ss3/checks.py`): the loop cost (energy with the loop minus the same stem without it) at 37/55/72/90 C is 3.3/1.6/-0.02/-1.7 kcal/mol for M2004 (a loop becomes stabilising), 6.7/5.2/3.8/2.3 for M1999, and 3.7/3.8/3.9/4.1 for NUPACK dna04 (physically expected: entropic, grows with T). So M1999 also has a questionable loop temperature dependence, but it stays positive and its ensemble energies track NUPACK within about 2 kcal/mol from 55 to 72 C. M1999 is above 90 C also not monotonic (99 C row), so do not use Vienna above about 90 C.

### B2. NUPACK 4.0.2

- `nupack.Model(material='dna', celsius=T, sodium=Na, magnesium=Mg, ensemble='stacking')`; repr shows `dna04.2.json`. Ensembles accepted: `stacking` (default), `nostacking`, `some-nupack3`, `all-nupack3` (docs also list `dangle-stacking`, `coaxial-stacking`).
- Ranges enforced by assertion (tested): sodium 0.05 to 1.1 M, magnesium 0 to 0.2 M. Docs (docs.nupack.org/model): salt corrections are DNA-only; sodium is the sum of Na+, K+, NH4+; references SantaLucia 1998, SantaLucia and Hicks 2004, Peyret 2000, Koehler and Peyret 2005.
- Mg is converted exactly as Na_eq = Na + 3.3 sqrt(Mg) (M). Verified by root-finding: for Mg 0.5, 1.2, 2, 5 mM with Na 50 mM, the Na-only value that reproduces the ensemble energy was 0.12379, 0.16432, 0.19758, 0.28335 M, identical to the formula to 5 digits.
- Functions used: `pfunc([s], model)` returns (Q, G); `mfe([s], model)[0].energy/.structure`; `pairs([s], model).to_array()` (diagonal = P(unpaired)); `structure_energy([s], structure=..., model=...)`; `sample([s], num_sample, model)` (1000 samples of 300 nt in 0.05 s). No hard-constraint pfunc was found in the Python API, so window accessibility is estimated by sampling (SE = sqrt(p(1-p)/n)).
- Licence: NUPACK is free for academic non-commercial use only (stated on nupack.org licence terms; I did not re-read the full licence text in this session). Keep it an optional import.

### B3. seqfold 0.9.0

`seqfold.dg(seq, temp=65.0)` returns MFE only (no partition function, no salt argument). 300 nt fold took 27.8 to 35 s (vs 0.09 s Vienna). On 20 folds its MFE was 9 to 15 kcal/mol less negative than Vienna M1999 and NUPACK (for example rpoB + strand 65 C: seqfold -9.4, Vienna -25.3, NUPACK -26.1). Not suitable.

### B4. Benchmark: 10 Mtb loci x 2 strands x 5 temperatures (`ss3/bench2.py`, `ss3/bench2.out`)

Windows: 300 nt centred on NC_000962.3 positions 6735, 7576, 761125, 779240, 781750, 801268, 1461180, 1473290, 1673428, 2155168 (GC 56.7 to 68.7%). Site = 3'-terminal 22 nt of each strand (where the opposite primer binds on an amplicon strand); "P8" = the 8 site nt that pair with the primer's 3' end. Vienna M1999 noGU dangles 2; NUPACK stacking; NUPACK site probabilities from 2000 samples.

Reference salt (Vienna 1.021 M, NUPACK 1.0 M), means over 20 strands:

| T | MFE V99 / V04 / N | G_ens V99 / V04 / N | MAE G_ens V99-N / V04-N | paired % V99 / N | r(p_unpaired) | P22 open V99 / N | P8 open V99 / N |
|---|---|---|---|---|---|---|---|
| 55 | -32.3 / -38.9 / -31.9 | -41.3 / -50.8 / -41.9 | 1.96 / 8.99 | 42.9 / 42.0 | 0.84 | 0.074 / 0.056 | 0.20 / 0.17 |
| 60 | -26.2 / -36.7 / -25.1 | -35.6 / -50.1 / -36.0 | 1.85 / 14.1 | 40.8 / 38.2 | 0.83 | 0.100 / 0.068 | 0.23 / 0.19 |
| 65 | -20.9 / -37.5 / -19.1 | -30.6 / -52.1 / -30.9 | 1.84 / 21.1 | 38.6 / 34.8 | 0.81 | 0.128 / 0.089 | 0.26 / 0.23 |
| 68 | -18.2 / -39.9 / -15.7 | -27.9 / -54.7 / -28.2 | 1.90 / 26.5 | 37.3 / 32.7 | 0.80 | 0.145 / 0.101 | 0.28 / 0.24 |
| 72 | -14.7 / -44.6 / -11.8 | -24.8 / -59.8 / -25.0 | 2.04 / 34.8 | 35.5 / 29.9 | 0.78 | 0.168 / 0.114 | 0.30 / 0.26 |

PCR buffer (50 mM monovalent, 1.2 mM free Mg; Vienna salt = 0.1643 M; NUPACK sodium 0.05, magnesium 0.0012):

| T | G_ens V99 / N | MAE | paired % V99 / N | P22 V99 / N | P8 V99 / N | min P22 (V99) | strands with P22 < 0.1 (V99) |
|---|---|---|---|---|---|---|---|
| 55 | -29.2 / -29.5 | 1.94 | 37.2 / 34.2 | 0.19 / 0.11 | 0.33 / 0.25 | 0.0000 | 55% |
| 60 | -23.9 / -24.6 | 1.84 | 34.2 / 30.5 | 0.24 / 0.13 | 0.38 / 0.29 | 0.0001 | 55% |
| 65 | -19.4 / -20.6 | 1.95 | 31.0 / 26.6 | 0.28 / 0.17 | 0.41 / 0.32 | 0.0002 | 45% |
| 68 | -17.0 / -18.5 | 1.98 | 29.0 / 24.3 | 0.31 / 0.18 | 0.44 / 0.34 | 0.0005 | 35% |
| 72 | -14.3 / -16.2 | 2.14 | 26.3 / 21.3 | 0.35 / 0.21 | 0.46 / 0.37 | 0.0014 | 35% |

Spearman rank correlation of P22 between Vienna and NUPACK (PCR buffer): 0.83 (55 C), 0.80, 0.77, 0.77, 0.74 (72 C). Per-site disagreements can be large, for example 7576 minus strand P22 0.010 (Vienna) vs 0.233 (NUPACK), 1461180 plus strand 0.72 vs 0.14. Timing per 300 nt strand: Vienna mfe+pf+bpp+window 0.21 s, NUPACK mfe+pfunc+pairs 0.24 s (12 parallel workers on a shared machine).

How much structure do 300 nt GC-rich strands keep? In PCR salt, ensemble free energy is still -14 to -16 kcal/mol at 72 C and about a quarter of bases are paired; at 65 C about 30% are paired. Primer sites at the 3' end of amplicon strands were fully open with probability 0.28 (Vienna) or 0.17 (NUPACK) on average at 65 C, and 35 to 55% of the 20 strands had P22 < 0.1. Structure does not "melt away" at typical Ta for this genome.

### B5. Folding context: amplicon strand vs genomic DNA (`ss3/checks.py`)

P(22-nt reverse-primer site unpaired), Vienna M1999, 65 C, 164 mM:

| centre | amplicon strand (site at 3' end) | genomic +-300 nt (full fold) | plfold W200 L150 | plfold W240 L160 | 120 nt local only |
|---|---|---|---|---|---|
| 6735 | 0.665 | 0.422 | 0.613 | 0.609 | 0.668 |
| 7576 | 0.482 | 0.242 | 0.357 | 0.320 | 0.486 |
| 761125 | 0.121 | 0.036 | 0.100 | 0.098 | 0.126 |
| 779240 | 0.006 | 0.000 | 0.000 | 0.000 | 0.006 |
| 781750 | 0.086 | 0.026 | 0.054 | 0.049 | 0.071 |
| 801268 | 0.012 | 0.004 | 0.010 | 0.010 | 0.012 |
| 1461180 | 0.721 | 0.450 | 0.663 | 0.661 | 0.722 |
| 1473290 | 0.103 | 0.072 | 0.152 | 0.140 | 0.102 |
| 1673428 | 0.000 | 0.000 | 0.002 | 0.002 | 0.000 |
| 2155168 | 0.485 | 0.080 | 0.170 | 0.168 | 0.422 |

Correlation of log10 P: amplicon vs genomic r = 0.96; genomic full vs plfold W200/L150 r = 0.98; genomic vs 120-nt local r = 0.96. An internal site (genomic DNA) has partners on both sides and is less accessible than the same site at an amplicon 3' end (factor 1.4 to 6 here).

---

## Part C. Design

### C1. Site accessibility and effective primer binding

For primer P (length L) binding template strand T at site s (0-based, end exclusive; the primer 3' end pairs with the 5'-most site nt):

```
P_open(s)  = up[s_end][L]                         (probs_window)
dG_open(s) = -R T ln P_open(s)                    R = 0.0019872 kcal/mol/K
dG_pf      = min(G_ens(primer alone), 0)          primer self-structure, Vienna
dG_duplex  = dH - T dS                            primer3 thal heterodimer with perfect complement,
                                                  mv, dv, dNTP, primer conc of the PCR
dG_eff     = dG_duplex + dG_open(s) - dG_pf
f_bound    = 1 / (1 + exp(dG_eff/RT - ln[P]))     primer in excess, [P] in M
P_seed3    = P(template nt opposite the primer 3'-terminal 8 nt unpaired)
P_nuc      = max over 8-nt windows in the site of P(window unpaired)   (toehold)
```

This is the coupled-equilibrium model used by Pythia (Mann 2009) and by FISH probe models: structure costs are additive in free energy. It assumes the site must be fully open, which is conservative (a primer can invade a partly structured site via a toehold), so P_nuc is reported alongside.

Thresholds at Ta (soft; calibrate with PrimerJinn amplicon depth):
- penalise dG_open(site) > 1.5 kcal/mol (P_open < 0.11 at 65 C), reject > 3.0 kcal/mol (P_open < 0.012);
- require P_nuc >= 0.1;
- require f_bound >= 0.9 for every primer at the recommended Ta;
- in a multiplex score, use dG_eff (not dG_duplex) when equalising primer strength across the pool.

Context: compute dG_open on (a) the amplicon strand (site at the strand 3' end) and (b) genomic context with W=200, L=150; take the larger dG_open. Use (a) for ranking within the pool since amplicon strands dominate after about 5 to 10 cycles.

### C2. Amplicon-level structure score

Per strand, at Ta and at the extension temperature (72 C, or Ta for two-step):
- `G_ens` and `G_ens/N` (kcal/mol/nt) and fraction paired = 1 - mean(p_unpaired);
- most stable local structure (Vienna `mfe_window_cb`, span <= 60 nt) and its ensemble melting temperature `Tm_local` = lowest T where G_ens(sub) >= -RT ln 2 (open chain is the majority state), computed in PCR salt on a 1 C grid;
- self-priming risk = P(3'-terminal nt paired) on each amplicon strand (Gimpel 2025 mechanism);
- G4 motifs (Huppert rule) on both strands;
- max 50-bp GC and product Tm (primer3 formula).

Flags: Tm_local > extension temperature with dG < -3 kcal/mol (polymerase pause risk; no quantitative dropout model verified); self-priming > 0.2; any G4; product Tm > 92 C; 50-bp window >= 80% GC.

Measured on the two demo loci (65 C, PCR salt): rpoB inverted repeat `CAAGCGCCGACTGTCGGCGCTGGG` Tm_local 101 C (Vienna) / 101 C (NUPACK); katG hairpins Vienna 96 and 100 C vs NUPACK 91 and 80 C. So Vienna's local hairpin Tm can be 5 to 20 C higher than NUPACK; treat Tm_local as a ranking value. Part of the gap is salt: an 8-bp GC hairpin at 95 C is -1.73 kcal/mol in Vienna at 0.164 M but +0.39 in NUPACK at 0.164 M (both near -1.6 to -2.5 at 1 M), so Vienna's salt correction is weaker for small hairpins.

### C3. Annealing temperature rule and additive advice

On a 1 C grid (55 to 72 C):
1. T_high = highest T where min over all primers of f_bound >= 0.9 (with structure terms).
2. T_low = lowest T where min over all sites of P_nuc >= 0.1.
3. Ta = T_high - 1 C. If T_low > Ta, report the limiting site (lowest P_nuc at Ta) and suggest redesigning that primer or adding an additive.
4. GC advice: if max product Tm > 92 C or any 50-bp window >= 80% GC or max amplicon GC >= 65%: 98 C denaturation (polymerase permitting), longer denaturation holds (Aird 2011: 3 min initial, about 80 s per cycle rescued high GC on a fast cycler), and 1 to 2 M betaine or 3 to 5% DMSO; 7-deaza-dGTP for G4 hits or > 75 to 80% GC regions (note it can reduce staining and some downstream steps; not tested here). If DMSO is used, recompute everything with Tm shifted by 0.6 C per % (primer3 default, unverified primary value range 0.5 to 0.75).
5. If Ta >= 68 C, two-step cycling is possible.

Demo result (rpoB RRDR + katG S315 loci, primers = 22 terminal nt of each 300 nt window, not designed primers): Ta 63 C (T_low 56, T_high 64); limiting primer katG reverse (f_bound 0.87 at 65 C); limiting site katG forward (dG_open 3.5 kcal/mol at 65 C, flagged "reject").

### C4. Whole strand or local window, and kinetics

- Early cycles: template is genomic DNA. The site is internal, and after denaturation the full complementary genomic strand is present at the same (low) concentration. At typical 1e3 to 1e5 genome copies, reannealing of the two genomic strands is slow compared with primer binding at 200 nM (k_on about 1e6 M^-1 s^-1 for a 25-mer per Gao 2006, so the pseudo-first-order rate is about 0.2 to 1 s^-1). Intramolecular folding is faster still (unimolecular). So for genomic templates, local intramolecular structure is the relevant competitor; fold the context with RNAplfold-style windows (W=200, L=150 reproduced the 900-nt full fold with r = 0.98).
- Later cycles: template is the amplicon strand (about 300 nt) at rising concentration. The site is at the strand 3' end. Competitors: intramolecular folding (fold the whole strand, which is cheap at 300 nt) and reannealing with the complementary amplicon strand, which becomes dominant late (plateau). Product reannealing is not modelled here; it is the standard explanation for plateau and is not sequence-specific enough to inform primer choice.
- Kinetic consequence: Gao 2006 shows that even modest structure (3 bp per strand) halves k_on and >= 4 bp makes binding non-two-state. With annealing holds of 15 to 30 s, a site with P_nuc of 0.01 gives a pseudo-first-order binding of about 0.2 x 0.01 x 30 = 0.06 per cycle (heuristic: the module's `p_primed_kinetic` uses k_on 1e6 M^-1 s^-1, 30 s, 200 nM). This is a rough indicator, not a validated predictor.
- What is predictive: GC content and extreme local GC predict failure (Benita 2003, Aird 2011); self-priming hairpins predict poor amplification in synthetic pools (Gimpel 2025); equilibrium models with folding terms design primers that work in hard regions (Mann 2009). What is not established: a calibrated mapping from dG_open to per-amplicon depth in multiplex amplicon sequencing. The thresholds above are physically motivated starting points and should be fitted to PrimerJinn's own read-depth data.

### C5. Reference module

`structure_reference.py` (same folder). Main functions:

| function | returns |
|---|---|
| `Conditions(mono_mM=50, mg_mM=2, dntp_mM=0.8, primer_nM=200, dmso_pct=0)` | buffer; `na_eq_M` = mono + 3.3 sqrt(Mg - dNTP) |
| `fold_stats(seq, T)` | MFE, structure, G_ens, fraction paired, per-nt p_unpaired |
| `unpaired_windows(seq, T, max_len)` | exact P(window unpaired) for all windows (probs_window) |
| `dg_open_constrained(seq, a, b, T)` | same by hard constraints (test oracle) |
| `site_accessibility_context(ctx, a, b, T, span=150, window=200)` | genomic-context accessibility |
| `site_metrics(template, a, b, primer, T)` | dG_open, P_seed3, P_nuc, dG_duplex, dG_eff, f_bound |
| `most_stable_local_structure`, `local_melting_tm(sub, cond, engine='vienna'|'nupack')` | local hairpin and its Tm |
| `self_priming_prob(strand, T)` | P(3'-terminal nt paired) |
| `g4_motifs`, `duplex_tm_long`, `max_window_gc` | flags |
| `amplicon_report(amp, fwd_len, rev_len, temps, ta)` | everything above plus flags |
| `recommend_ta(amplicons)` | Ta, T_low, T_high, limiting primer/site, advice |
| `nupack_crosscheck(seq, T, site=(a,b))` | NUPACK G_ens, fraction paired, sampled site dG_open; None if NUPACK missing |

Tests (`test_structure_reference.py`, `8 passed in 2.1 s`): probs_window equals the constraint definition within 0.02 kcal/mol; G_ens is monotonic in T from 55 to 90 C; no G-T pairs; product Tm equals the primer3 source formula; Na_eq matches NUPACK's Mg conversion; f_bound bounds; G4 finder on both strands; local Tm drops with salt.

Demo output (`python3 structure_reference.py`, saved in `ss3/demo.out`), abridged:

```
Na_eq for folding = 164.3 mM
=== rpoB_RRDR  NC_000962.3:761006-761305  len 300  GC 66.7%  max GC(50 bp) 78%  duplex Tm 94.5 C  G4 motifs []
  top      65 C  MFE  -14.11  G_ens  -22.71  (-0.076/nt)  paired 33.9%
  bottom   65 C  MFE  -12.12  G_ens  -20.86  (-0.070/nt)  paired 34.2%
  top      72 C  MFE   -8.63  G_ens  -16.54  (-0.055/nt)  paired 27.7%
  most stable local structure, top: 135-159 dG(65 C) -5.28  Tm 101.0 C  (.(((((((.....))))))).).
  fwd @65 C: p_open(site) 0.352 dG_open 0.70  p_open(3' seed) 0.460 p_nuc 0.754  dG_duplex -13.23  dG_primer_fold -0.18  dG_eff -12.34  f_bound 0.950 (no structure 0.986)
  rev @65 C: p_open(site) 0.401 dG_open 0.61  ... dG_eff -15.33  f_bound 0.999
  test: rev-site dG_open by constrained pf 0.615 == probs_window 0.615
  P(3'-terminal nt paired, self-priming risk) @65 C: top 0.003, bottom 0.117
  rev site in genomic context @65 C: P(open) full fold 0.080, plfold W200/L150 0.167 (amplicon strand 0.401)
  NUPACK cross-check top strand 65 C: G_ens -23.13 (Vienna -22.71)  paired 28.9% (Vienna 33.9%)  rev-site dG_open 1.11 (Vienna 0.61)
=== katG_S315  NC_000962.3:2155019-2155318  len 300  GC 64.3%  max GC(50 bp) 80%  duplex Tm 93.6 C  G4 motifs []
  top      65 C  MFE   -8.73  G_ens  -17.12  (-0.057/nt)  paired 29.0%
  fwd @65 C: p_open(site) 0.005 dG_open 3.50  p_open(3' seed) 0.026 p_nuc 0.143  dG_duplex -17.83  dG_primer_fold -2.12  dG_eff -12.22  f_bound 0.940 (no structure 1.000)
  rev @65 C: p_open(site) 0.485 dG_open 0.49  ... dG_eff -11.64  f_bound 0.870 (no structure 0.948)
  flags: fwd site dG_open 3.5 > 3.0 (reject); local structure Tm 96 C; duplex Tm 93.6 C; 50-bp window 80% GC
  P(3'-terminal nt paired, self-priming risk) @65 C: top 0.020, bottom 0.238
  NUPACK cross-check top strand 65 C: G_ens -19.43 (Vienna -17.12)  paired 24.0% (Vienna 29.0%)  rev-site dG_open 1.08 (Vienna 0.49)
Recommended Ta for the 2-amplicon pool: 63.0 C (T_low 56.0, T_high 64.0; limiting primer katG_S315:R, limiting site katG_S315:F)
  advice: GC-rich pool (max amplicon GC 67%, max duplex Tm 94.5 C): use 98 C denaturation ..., consider 3 to 5% DMSO or 1 M betaine.
```

Runtime: about 20 s for the whole demo including an 18-temperature Ta scan of 4 sites; a single 300 nt amplicon report is under 1 s. For a pool of 100 amplicons the Ta scan is about 18 x 200 x 0.1 s = 6 min single-threaded; parallelise per amplicon or use a coarser grid (2 C) then refine.

---

## 5. Open questions

1. No calibration data: thresholds (1.5 / 3.0 kcal/mol, P_nuc 0.1, f_bound 0.9) need fitting to PrimerJinn amplicon depth from real Mtb runs.
2. Which model is closer to reality for site accessibility? Vienna and NUPACK disagree per site by up to about 5x; no experimental ground truth was found for DNA in PCR buffer at 60 to 72 C.
3. Vienna's DNA loop temperature dependence (both sets) and weaker small-hairpin salt correction vs NUPACK are not resolved; a patched Mathews2004 with purely entropic loops (`patchpar.py`, earlier attempt) had MAE 3.6 to 5.2 kcal/mol vs NUPACK, worse than M1999.
4. Effect of betaine on folding energies is not modelled by any tool here; DMSO is modelled only as a Tm shift.
5. primer3-py 2.0.3 `calc_tm` for long sequences is about 2 C lower than the primer3 C-source formula; cause not traced.
6. Exact additive concentrations in Musso 2006 and the per-percent DMSO value in the primary papers were not verified (only via the primer3 manual).

## 6. Files

- Report: `secondary_structure.md` (this file)
- Module and tests: `structure_reference.py`, `test_structure_reference.py`
- Scripts and outputs: `ss3/api.py`, `ss3/bench2.py`, `ss3/bench2.out`, `ss3/bench2.json`, `ss3/checks.py`, `ss3/demo.out`, `ss2/t1.py` + `ss2/t1.out` (parameter set vs temperature), `bench1.py` + `ss3/bench1.out` (patched M2004 variants)
- Papers (full text as text): `ss3/papers/PMC3188800.txt` (Aird), `PMC1488884.txt` (Gao), `PMC12533003.txt` (Gimpel), `PMC1140081.txt` (Huppert), `PMC5515968.txt` (Perrone); Pythia in `dl/pythia.txt`; primer3 source `src_p3/oligotm.c`, manual `dl/p3manual.htm`
