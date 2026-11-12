# Fact-check of ampseq_practice.md

Checked 2026-11-11. Every item was re-checked against the primary text (PMC full texts already on disk, re-searched by me), freshly downloaded source (GitHub), PubMed metadata, the live NEB Tm API (80 new oligos, no personal data sent), or code I ran with primer3-py 2.0.3 and Biopython 1.85 as independent implementations. My scripts and outputs: `research/verify_ampseq/` and the `_fixed.py` files in `research/ampseq_exp/`.

## Headline problems (fix before implementing)

1. **Hidden formamide in primer3-py.** primer3-py 2.0.3 `calc_tm()` defaults to `formamide_conc=0.8` (mol/L) (`primer3/argdefaults.py`). That lowers every Tm by about 2.1 C. All Tm numbers in the report computed with `calc_tm` carry this shift, and the "mv 290 mM" fit absorbed it. Proof: at 1 M Na and no salt correction, primer3-py returns 67.50 C for GAGGAACACCACTAGTACCG at 500 nM. My own SantaLucia 1998 implementation and Biopython both give 69.60 C, and primer3-py also gives 69.60 C once `formamide_conc=0` is set. Primer3Plus's own defaults use `PRIMER_FORMAMIDE_CONC 0.0`.
   - The NEB fit still works numerically, but only while that library default holds. A **better, explicit model** is `primer3.calc_tm(seq, mv_conc=195, dv_conc=2, dntp_conc=0.8, dna_conc=<nM>, salt_corrections_method='owczarzy', formamide_conc=0)`. Its RMSE against NEB is 0.35 C on the report's 240 oligos and 0.46 C on 80 fresh oligos (16 to 34 nt, 15 to 75% GC). The report's model gives 0.50 and 0.64 C on the same sets. Biopython `Tm_NN(DNA_NN3, Na=195, dnac1=C/4, dnac2=0, saltcorr=6)` reproduces it exactly (max difference <0.01 C).
   - `dv_conc` and `dntp_conc` do nothing in either fit. At these salt levels sqrt([Mg free])/[Na+] < 0.22, so primer3's Owczarzy 2008 scheme uses the monovalent-only branch (I checked: dv 0 and dv 2 give identical Tm). Raising dv to about 3 mM would switch branches and change Tm without warning.
   - "Owczarzy Mg correction" is the right citation (primer3 option 2 = Owczarzy 2008, per the primer3 manual). In this regime, though, Mg is not what does the work.
   - Section 0 row "primer3 at PrimalScheme conditions, median 59.1" is an artefact. With `formamide_conc=0`, primer3 reproduces the ARTIC tsv values to +0.33 C (SD 0.14) over all 196 primers. With the default 0.8 M it is 1.82 C low. So the correct figure is about 61.2 C, matching the tsv, not 59.1 C. Likewise "primer3, 50 mM Na+, 50 nM: 55.8 / 56.1" should be 57.9 / 58.3 C.
   - The same trap affects primalscheme3. Its `thermo.calc_tm` calls primer3-py `calc_tm` without a formamide argument (GitHub artic-network/primalscheme3, `primalscheme3/core/thermo.py`). With primer3-py 2.x it therefore reports Tm about 2 C lower than the values in the ARTIC V3 tsv. When you quote a "PrimalScheme 59.5 to 62.5 C" window, say which primer3-py version it assumes.
2. **Annealing recommendation contradicts the evidence and the report's own code.** Recommendation 6 says "lowest Tm at the real pool concentration + 1 C". For ARTIC V3 (lowest NEB-equivalent Tm at 15 nM about 59.8 C) that gives about 60.8 C, but ARTIC runs at 65 C and drops to 63 C only if amplicon 64 fails. `gc_mapping.py` in fact uses Tm at 500 nM + 1. The NEB +1 rule was set for 500 nM singleplex reactions. Use the 500 nM convention for Ta, and keep Tm at the pool concentration for dG and occupancy matching. Fixed in `ampseq_exp/gc_mapping_fixed.py`.
3. **Tailed-dimer statistic was biased.** `tails2.py` compared one untailed partner against the minimum over two tailed partners, and left primer p untailed. The TruSeq Read 2 tail was also 32 nt instead of 34 (missing the 5' GT). In a fair one-partner comparison with both primers tailed (`tails2_fixed.py`, 3,000 pairs), the fraction of 3'-anchored dG < -5 kcal/mol at 37 C is: 0.004 to 0.021 (35% GC), **0.016 to 0.061** (50% GC, the report said 0.015 to 0.084), and **0.096 to 0.204** (65% GC, the report said 0.095 to 0.216). The direction holds; the 50% GC effect is about 3.8-fold, not 5.6-fold.

## Claim by claim

| # | Claim | Verdict | Evidence / correction |
|---|---|---|---|
| 1 | primer3 mv 290 + Owczarzy reproduces NEB Q5 (RMSE 0.50, bias -0.04 at 500 nM, -0.09 at 25 nM); current PrimerJinn approximation bias -1.40, SD 1.30 | Corrected | Numbers reproduce from `neb_calib.json`. On 80 fresh oligos the RMSE is 0.64 C (bias -0.10), with a GC-dependent residual (+0.9 C at 10% GC, -1.0 C at 90% GC). The fit only works because of the hidden formamide 0.8 M default (see above). Use mv 195 with formamide 0: RMSE 0.35 / 0.46 C. Current repo approximation (`PCRinSilico.py` `calcTm(dv 2, mv 70, dna 3300)`) against fresh NEB data: bias -1.50, SD 1.73 at 500 nM; +2.72, SD 1.97 at 15 nM. |
| 2 | ARTIC V3 NEB Q5 median 66.5 C (64.1 to 69.4); tsv 60.0 to 62.5; 22 to 30 nt; GC 30 to 55% | Confirmed | Fresh GitHub download is identical to the saved tsv: n=196, 22 to 30 nt, GC 30.0 to 54.55%, tsv Tm 60.01 to 62.45. My corrected model gives a median of 66.4 C (64.1 to 69.0) over all 196. The tsv values were made without formamide (see headline 1). |
| 3 | PrimalScheme v1 Tm 59.5 to 62.5 (opt 61), mv 100, dv 2, dNTP 0.8, 15 nM; size 19 to 34 (17 to 30 high GC); GC 30 to 55 (40 to 65); homopolymer 5 | Confirmed | Freshly downloaded aresti/primalscheme `config.py` is byte-identical to the saved copy. primalscheme3 config (artic-network/primalscheme3 main) also matches the saved `ps3_config.py`: 19 to 36, high GC 17 to 30, GC 30 to 55 and 30 to 67, hairpin 55, dimer -26.0. |
| 4 | ARTIC Q5 at 15 nM, 98 C 15 s / 65 C 5 min, 25 to 35 cycles; Quick 2017 says 65 to 68 C | Confirmed | LoCost v3 text: 15 nM per primer; 98 C 30 s activation; 98 C 15 s / 65 C 5 min for 25 to 35 cycles. Quick 2017 quote found verbatim; its table gives 0.015 µM per primer and the text 0.02 µM. **Missed:** LoCost v3 says "If you see amplicon 64 dropout then decrease the annealing/extension temperature to 63 °C". The report lists 63 C as UNVERIFIED although that file is on disk. |
| 5 | primerJinn MTB 8-plex: 20-mers, mean NEB Tm 65 (64 to 67), 454 to 791 bp, Q5 35 cycles 98/10 s, 65/20 s, 72/30 s | Confirmed | PMC10714587 text and all 16 primers in the XML are 20 nt. NEB API values for pair 1 (64.25 and 64.48) are in the saved JSON. |
| 6 | NEB Ta = Tm_lower + 1 (Q5), - 5 (Taq); 500 nM (Q5), 200 nM (Taq); protocol text +3 C and two-step at Ta 72 C or higher | Confirmed, with an addition | Fresh API calls: Q5 Ta - Tm_lower = +0.95 to +1.04 across 40 pairs, **capped at 72 C** (pairs with high Tm gave Ta 72). Taq gives -4.96 to -5.03. The concentration notes appear only when the entered concentration differs from the recommendation. Protocol text quotes found verbatim in q5.txt. |
| 7 | Concentration 500 to 25 nM lowers NEB Tm by 3.58 C; primer3 20-mer drops 4.9 C from 500 to 15 nM | Confirmed | 3.58 C recomputed from the JSON; fresh data give 4.22 C for 500 to 15 nM. primer3 drop 59.02 to 54.09 = 4.93 C (formamide shifts both ends equally). |
| 8 | Gimpel 2025: 5% gives 2-fold after 12 cycles; about 2% of sequences at up to 80%; GC/base composition near random; CNN AUROC 0.88, AUPRC 0.44; 6-nt motif 94.0 to 78.8% | Confirmed, wording note | All numbers found in PMC12533003. AUROC 0.88/0.87 and AUPRC 0.42/0.44 are for GCall/GCfix. The report's formula ((1+E1)/(1+E2))^n does not give 2-fold if "5%" means 5 percentage points of E (1.85/1.90 to the 12th power = 0.73). Gimpel's efficiency is a relative per-cycle yield: 0.95^12 = 0.54, and 0.8^3 = 0.51 matches their "halving every 3 cycles". **Missed:** code is public at github.com/BorgwardtLab/PCR-bias (BSD-3); the report says this was not checked. |
| 9 | Illumina: Tm on gene-specific part (60 to 65 C), hairpin and dimer on full primer | Confirmed | Quote verbatim in 16s.txt, Part # 15044223 Rev. B; overhang sequences match. First PCR uses KAPA HiFi at 55 C, 25 cycles. |
| 10 | Tail raises Tm by about 20 C; tailed-partner dimer fraction 0.015 to 0.084 (50% GC), 0.095 to 0.216 (65% GC) | Tm: confirmed. Fractions: corrected | Tm rise reproduced: +18.7 to +19.0 C in the corrected Q5 model and +20.6 to +20.9 C in Biopython (mv 50, Mg 2); +23 C for an AT-rich 26-mer. Fair-comparison fractions are 0.016 to 0.061 and 0.096 to 0.204 (see headline 3). These are 37 C dG values, not values at Ta. |
| 11 | Aird 2011: plateau 11 to 56% GC on fast cycler; >65% GC to about 1/100; slow cycler 13 to 84%; long denaturation + 2 M betaine helped | Confirmed with nuance | All numbers verbatim in PMC3188800. Nuance: long denaturation alone or betaine alone did **not** recover extreme-GC loci. Together they did, but at the cost of low-GC loci (plateau shifted to 23 to 90%). The flattest profile came from AccuPrime Taq HiFi with 65 or 60 C annealing/extension. |
| 12 | 60 C extension needed for ≥90% AT; López-Barragán found better AT-rich coverage at 60 C vs 65 and 70 C | Confirmed (secondary), Su 1996 unverifiable | Quote verbatim in PMC3026866. Their gain began only at about 80% AT per 100 bp, with a maximum ratio of about 2.8. Anneal was 65 C, 19 cycles, library PCR (not amplicon multiplex). Su 1996 (PMID 8628694) is in PMC145803, but only as a scanned PDF; the text was not retrievable, so the claim is unverifiable at the primary source. |
| 13 | SADDLE dG about -11.5 (window -10.5 to -12.5); 90.7% to 4.9% dimers at median 45 nM; Olivar defaults 60 C, 0.18, -11.8 | Confirmed | Quotes verbatim in PMC9001684 (also 17 cycles, Vent exo-, 3 min at 60 C annealing) and PMC11282221. Missing: SADDLE's dG reference temperature and salt are not given in the main text, so -11.5 cannot be applied to a primer3 or NUPACK dG without checking SADDLE Supplementary S9. |
| 14 | Deeplex: 95 C 2 min; 35 cycles 94 C 1 min, 55 C 30 s, 72 C 1 min 30 s; enzyme and primers not disclosed | Confirmed | manual.txt (V5, March 2023), 24-plex; proprietary master mix; ramp 3.0 C/s block on a GeneAmp 9700. |
| 15 | 3' 5-mer end stability 0.86 (TATAT) to 6.86 (GCGCG); Primer3Plus default 9.0 never rejects | Confirmed | The primer3 manual states these values. I also ran primer3 `check_primers` and got END_STABILITY 6.86 / 0.86 / 4.02, independent of salt. Primer3Plus GitHub `default_settings.json` has PRIMER_MAX_END_STABILITY 9.0. The primer3 core default is 100.0. Note: table 1.6 in the report still lists "9 kcal/mol", which contradicts its own advice to use about 6.5. |

## Other items checked

- DOIs for Quick 2017, Gimpel 2025, Aird 2011, Olivar 2024, primerJinn 2023, López-Barragán 2011, Itokawa 2020, Colman 2015, Mariner-Llicer 2021 and NOMADS 2024 all match the PMC headers.
- Itokawa: 18_LEFT and 76_RIGHT were complementary over 10 nt at their 3' ends; Ta range 63.1 to 68.6 C. Confirmed.
- Colman 2015: Q5 2X + 5 µl of 5 M betaine in 25 µl (1 M); 100/400/100/200/400/400 nM; 25 cycles at 98/60/72 C; index PCR 6 cycles. Confirmed.
- Mariner-Llicer: KAPA HiFi HS, 30 cycles of 98 C 20 s, 65 C 15 s, 72 C 2 min. Confirmed.
- NOMADS: 0.6 µM total, doubled or halved at fixed total; dimer score < -6; 79% AT; "reduced annealing and extension temperatures". Confirmed; exact temperatures not stated.
- Andreson 2008 (PubMed abstract and text): binding-site count is the most important factor; GM1 = four factors; failure 17% to 6%. Confirmed. The data were human-genome PCRs, so for MTB the relevant genome for counting sites is host DNA.
- Benita 2003 (abstract): 61% GC threshold, 21 bp window, 84.3% specificity and 94.8% sensitivity, human exons. Confirmed. The report's "max(80%, genome GC + 20 pp)" rule is its own invention and has not been validated.
- Kebschull 2015: the quote matches the abstract. Confirmed.
- **Onodera and Melcher 2004: misrepresented.** AGG 3.27%, TTA 0.42% and the least-used list are correct. But the authors write that it is "unlikely" the frequencies are "due primarily to primer design software or literature recommendations". The report says the opposite ("partly driven by design software defaults, as the authors note"). The survivorship caveat is fair; the attribution is wrong.
- Phusion Ta in the saved NEB responses is Tm_lower + 3.5 to 3.8 C (not a constant +3.5). The report already marks it as unverified.
- Fraction-bound example (section 5.2) reproduced at mv 50. Under the Q5-equivalent salt (mv 195), occupancy at 60 C and 15 nM is 0.62, not 0.17. At 65 C it is 0.05, so the qualitative conclusion holds at 65 C only.
- WHO tNGS rapid communication: not checked.

## Recommendations: status

1. Tm call: replace with `mv_conc=195, formamide_conc=0` (explicit), and keep `salt_corrections_method='owczarzy'`. Label the result "NEB-Q5-equivalent Tm (500 nM)". Also report Tm at the real concentration.
2. Replace the current approximation: agreed (fresh data RMSE 2.3 to 3.4 C).
3. Tm_opt ramp (58 to 65 C): the author's own interpolation, not tested on any panel. Keep it as a default and say so in the output.
4 and 5. Length and GC windows: unvalidated proposals. They agree with the length-for-Tm simulation, but no published panel anchors them.
6. Ta: **change** to Tm_lower(500 nM convention) + 1, capped at 72 C, with a 63 C fallback note (ARTIC). Do not apply +1 to Tm at pool concentration.
7. Two-step default: ARTIC uses it at 65 C; the NEB protocol says two-step needs Ta of 72 C or higher. Present it as an option backed by practice, not by NEB.
8. 60 C extension: evidence supports it for ≥80 to 90% AT regions. A genome GC < 30% cut-off is a reasonable proxy but has not been tested.
9. Tails handling: confirmed by Illumina guidance.
10. Adapter-motif check: mechanism confirmed by Gimpel; the 4 to 6 nt / 20 nt window is the author's choice. Gimpel's code could be used directly.
11. Reject contiguous 3' overlap ≥6 nt: no source. Itokawa's failure case was 10 nt. Treat 6 nt as a tunable threshold.
12. Soft penalties for 3' rules: reasonable. Fix table 1.6 (9 kcal/mol to about 6.5).
13. Off-target count: supported by Andreson.
14. Amplicon GC window: unvalidated synthesis.
15. Rebalancing worksheet: supported by NOMADS and Colman practice.

## Points the report missed

- primer3-py 2.0.3 `calc_tm` defaults to formamide 0.8 M, and primalscheme3 inherits it. Every Tm in the tool must pass `formamide_conc=0` explicitly, or Tm conventions drift by about 2 C between library versions.
- `dv_conc`/`dntp_conc` are inert in the fitted regime, and there is a discontinuity at the Owczarzy 2008 branch point (R = 0.22).
- NEB caps Q5 Ta at 72 C.
- The ARTIC 63 C fallback is documented in LoCost v3.
- The GC-dependent residual of the NEB emulation (about ±1 C at GC extremes) matters for the P. falciparum and MTB ends of the range. Consider a GC term or a per-genome check against NEB.
- Gimpel code availability (BorgwardtLab/PCR-bias).
- The SADDLE dG target needs its reference temperature and salt before reuse.
- Dimer dG thresholds in `tails2` are at 37 C; the recommendation says Ta. The thresholds need re-deriving at Ta.

## Files

- `research/ampseq_exp/neb_calib_fixed.py` (+ `.out`): explicit-formamide NEB Q5 emulation, refit, fresh validation, Biopython cross-check.
- `research/ampseq_exp/gc_mapping_fixed.py` (+ `.out`): uses the fixed Tm model, computes Tm at pool concentration per primer, Ta from the 500 nM convention.
- `research/ampseq_exp/tails2_fixed.py` (+ `.out`): fair tailed-dimer comparison.
- `research/verify_ampseq/`: `neb_fresh.py`/`neb_fresh.json` (80 fresh oligos, NEB API v0.19.0), `refit.py`, `diag.py`, `nn_own.py`, `bio_cross*.py`, `gh/` (fresh GitHub copies of the PrimalScheme configs, ARTIC V3 tsv, primalscheme3 `thermo.py`, Primer3Plus defaults).
