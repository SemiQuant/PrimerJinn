# Evidence-based design targets for multiplex amplicon sequencing panels, by genome GC

Scope: primer Tm targets vs genome/local GC, coverage uniformity, amplicon GC window, adapter tails, cycling and primer concentration, 3' end rules. Every number below is either quoted from a primary source I read in this session (file paths given) or produced by code I ran (scripts and outputs under `research/ampseq_exp/` and `research/neb/`). Items I could not verify are marked **UNVERIFIED**.

Working folder: `/private/tmp/claude-501/-Users-semiquant--claude-spares-spare-n9Plyg/5e124d61-7974-4f46-b4a7-ff4fbdb28e2f/scratchpad/research/`

Reused from earlier attempts: PMC full texts in `ampseq/dl/`, `papers/`, `mx/`; PrimalScheme configs `ampseq/dl/primalscheme_config.py` (PrimalScheme v1, Quick and Smith 2020) and `ampseq/dl/ps3_config.py` (primalscheme3); Deeplex manual `ampseq/web/deeplex/manual.txt`; NEB Q5 protocols.io PDF `ampseq/web/q5pdf/q5.txt`; Illumina 16S guide `ampseq/web/16s.txt`; ARTIC LoCost `ampseq/web/locost.txt`; primer3 manual `ampseq/web/p3manual.htm`; NEB Tm API OpenAPI spec `tm_dl/reference.json`. The NEB product HTML pages saved earlier were Cloudflare challenge pages (no content) and were not used.

New in this session: `ampseq_exp/tails_conc.py`, `end5_tails.py`, `tails2.py`, `neb_calib.py`, `artic_tm.py`, `gc_mapping.py` with `.out` files; `neb/neb_api_results.jsonl`, `neb/artic_neb.txt`; `ampseq_exp/artic_dl/v3.tsv` (ARTIC V3 primer table from github.com/artic-network/artic-ncov2019).

---

## 0. The single most important finding: Tm numbers are convention dependent

The user's rule ("about 65 C for ~65% GC genomes, 55 to 60 C otherwise") mixes Tm conventions. The same primer gets very different "Tm" values depending on the salt model and the oligo concentration entered:

| Primer set | Tool / convention | Tm reported |
|---|---|---|
| primerJinn MTB pair 1 (GAGGAACACCACTAGTACCG / CTCGATGACTTTACGGCCAT) | NEB Q5 calculator, 500 nM | 64.2 / 64.5 C (API, run) |
| same | NEB Q5, 25 nM | 59.9 / 60.2 C |
| same | primer3, 50 mM Na+, 2 mM Mg, 0.8 mM dNTP, 50 nM | 55.8 / 56.1 C |
| same | current PrimerJinn "Q5 approximation" (mv 70, dv 2, dna 3300 nM) | 62.2 / 62.9 C |
| ARTIC V3, 196 primers | Value printed in the ARTIC tsv ("tm (use 65)") | 60.0 to 62.5, median 60.8 |
| ARTIC V3, random 30 | NEB Q5, 500 nM | median 66.5 (64.1 to 69.4) |
| ARTIC V3, random 30 | NEB-Q5-equivalent model at the real 15 nM | median 62.0 |
| ARTIC V3, random 30 | primer3 at PrimalScheme conditions (mv 100, dv 2, dNTP 0.8, 15 nM) | median 59.1 |

Sources: `neb/neb_api_results.jsonl`, `ampseq_exp/tails_conc.out`, `ampseq_exp/end5_tails.out`, `ampseq_exp/artic_tm.out`.

So in the NEB Q5 convention used by the primerJinn paper, SARS-CoV-2 ARTIC primers (38% GC genome) sit at about 66 C, the same place as the M. tuberculosis primers (65 C). The "60 to 63 C" quoted for ARTIC is primer3 at 100 mM monovalent salt. **The tool must fix one Tm convention, state it in output, and express every target in it.** Recommended convention below.

### 0.1 Reproducing the NEB Q5 calculator with primer3 (verified, run)

`ampseq_exp/neb_calib.py` sent 240 random oligos (18 to 32 nt, GC 20 to 70%) to the NEB Tm API (`https://tmapi.neb.com/tm`, prodcode `q5-0`, no email sent) at 500 nM and 25 nM, then fitted primer3-py 2.0.3 `calc_tm` (SantaLucia 1998 NN, Owczarzy 2008 Mg correction) over the monovalent concentration:

```
best fit at 500 nM (dv 2 mM, dNTP 0.8 mM, Owczarzy): rmse 0.50 mv 290 bias -0.04 sd 0.49
same mv at 25 nM vs NEB 25 nM: bias -0.09 sd 0.49
NEB mean shift 500->25 nM: 3.58 C
current PrimerJinn approx (mv70 dv2 dna3300, default santalucia salt) vs NEB 500 nM: bias -1.40 sd 1.30
```

Interpretation: `primer3.calc_tm(seq, mv_conc=290, dv_conc=2, dntp_conc=0.8, dna_conc=C_nM, salt_corrections_method='owczarzy')` reproduces NEB Q5 Tm to 0.5 C RMSE at both 500 and 25 nM. The 290 mM is an effective value that absorbs the undisclosed Q5 buffer, not a real salt concentration. It is 2.7 times better than the current PrimerJinn approximation (SD 1.30 C, bias -1.40 C). The concentration term behaves as NN theory predicts, so the tool can use the real per-primer pool concentration.

NEB API facts observed in responses (`neb/neb_api_results.jsonl`):
- Q5: "the rule used to determine Q5 annealing temperatures (Ta = Tm_lower+1°C)"; observed Ta 65.2 for Tm_lower 64.25. Notes also say "The recommended primer concentration for Q5 reactions is 500 nM."
- Taq: Ta = Tm_lower - 5 (observed 52.4 for Tm_lower 57.41). "The recommended primer concentration is 200 nM."
- Phusion: the API divides the entered concentration by 4 (0.5 entered, 0.125 used) and gave Ta = Tm_lower + 3.5 for 20-mers. Exact Phusion rule **UNVERIFIED**.
- The NEB Q5 protocol (protocols.io, doi:10.17504/protocols.io.cjtunm, 2015) says instead "a 10–30 seconds annealing step at 3°C above the Tm of the lower Tm primer" and "When primers with annealing temperatures ≥ 72°C are used, a 2-step thermocycling protocol ... is possible". The calculator (API v0.19.0, queried 2026-11-11) now uses +1 C. Use the calculator rule; note the protocol text is older.

---

## 1. Tm targets by GC: what published panels actually use

### 1.1 Mycobacterium tuberculosis (65.6% GC)

| Panel | Polymerase | Primers / Tm | Cycling | Primer conc | Source |
|---|---|---|---|---|---|
| primerJinn 8-plex (Limberis and Metcalfe 2023, BMC Bioinformatics 24:468, doi:10.1186/s12859-023-05609-1) | NEB Q5 HotStart MasterMix | 20-mers; default opt Tm 65 C (62 to 68), mean 65 C (64 to 67) NEB Tm; amplicons 454 to 791 bp | 35 cycles 98 C 10 s, 65 C 20 s, 72 C 30 s | not stated in text | `ampseq/dl/PMC10714587.txt` |
| Next Gen-RDST / SMOR 6-plex (Colman et al. 2015 PLoS One, doi:10.1371/journal.pone.0126626) | Q5 Hot Start 2X + 1 M betaine (5 µl of 5 M in 25 µl) | universal tails | 25 cycles 98 C 10 s, 60 C 15 s, 72 C 20 s; index PCR 6 cycles at 60 C | 100 to 400 nM per target, unequal (gyrA 100, eis 400, rpoB 100, katG 200, inhA 400, rrs 400 nM) | `ampseq/dl/PMC4430321.txt` |
| Deeplex Myc-TB 24-plex (GenoScreen; manual RUO V5, March 2023) | proprietary master mix | not disclosed | 95 C 2 min; 35 cycles 94 C 1 min, 55 C 30 s, 72 C 1 min 30 s | not disclosed | `ampseq/web/deeplex/manual.txt` |
| Nanopore MTB amplicon (Mariner-Llicer et al. 2021 Microb Genom, doi:10.1099/mgen.0.000740) | KAPA HiFi HS | 2.5 to 3.6 kb amplicons, concentrations optimised for equimolar output | 30 cycles 98 C 20 s, 65 C 15 s, 72 C 2 min | per-primer optimised (Table S1) | `ampseq/dl/PMC8767343.txt` |

WHO: the July 2023 rapid communication named three first-in-class tNGS products, Deeplex Myc-TB (GenoScreen), AmPORE-TB (Oxford Nanopore Diagnostics) and TBseq (ShengTing) (https://www.who.int/publications/i/item/9789240076372; https://tbksp.who.int/en/node/2704). WHO documents do not specify primer Tm or cycling. Product drug lists from secondary summaries were inconsistent and are **not used**.

Take-away: high-fidelity (Q5, KAPA) MTB panels anneal at 60 to 65 C with primers at about 65 C NEB-Q5 Tm. Deeplex anneals at 55 C with an undisclosed enzyme, so its primer Tm cannot be compared.

### 1.2 SARS-CoV-2, ARTIC / PrimalScheme (38% GC)

- PrimalScheme v1 config (`ampseq/dl/primalscheme_config.py`): `PRIMER_OPT_TM 61.0`, `MIN 59.5`, `MAX 62.5`, size `DEFAULT (19, 34, opt 22)`, `HIGH_GC (17, 30, opt 20)`, GC `DEFAULT 30 to 55 (opt 50)`, `HIGH_GC 40 to 65 (opt 60)`, `HIGH_GC_WARN_THRESHOLD = 55`, `PRIMER_MAX_HOMO 5`, `PRIMER_MAX_HAIRPIN_TH 47`, Tm conditions `MV_CONC 100, DV_CONC 2, DNA_CONC 15, DNTP_CONC 0.8`.
- primalscheme3 config (`ampseq/dl/ps3_config.py`): Tm 59.5 to 62.5; size 19 to 36 default, 17 to 30 high GC; GC 30 to 55 default, 30 to 67 high GC; hairpin Th max 55; `dimer_score -26.0`; same salt values.
- ARTIC V3 primers (`ampseq_exp/artic_tm.out`, computed from the official tsv): 196 primers, length 22 to 30 (median 24), primer GC 30 to 55% (median 45.5%, well above the 38% genome), NEB Q5 Tm at 500 nM median 66.5 C.
- Quick et al. 2017 Nat Protoc 12:1261 (doi:10.1038/nprot.2017.066, `ampseq/dl/PMC5902022.txt`): "Primers are designed with a high annealing temperature within a narrow range (65-68°C) which allows PCR to performed as a 2-step protocol (95°C denaturation, 65°C combined annealing and extension)". Reaction: Q5, 0.015 µM per primer in the table (text elsewhere says 0.02 µM), cycling 98 C 30 s then 98 C 15 s, 65 C 5 min, cycles 2 to 40. "lower primer concentrations and longer annealing/extension time were optimal."
- ARTIC LoCost v3 (protocols.io doi:10.17504/protocols.io.bp2l6n26rgqe/v3, `ampseq/web/locost.txt`): 15 nM per primer, 98 C 15 s, 65 C 5 min, 25 to 35 cycles.
- Itokawa et al. 2020 PLoS One (doi:10.1371/journal.pone.0239403, `ampseq/dl/PMC7500614.txt`): Ta gradient 63.1 to 68.6 C; optimum around 65 C; dimer-prone amplicons dropped as Ta fell.
- The "63 C" annealing used in some later ARTIC versions is **UNVERIFIED** here. One optimisation paper used 25 cycles at 65 C then 15 cycles at 62.5 C (Frontiers Med 2022, doi:10.3389/fmed.2022.836728).

### 1.3 Plasmodium falciparum (19% GC)

- Su et al. 1996 NAR 24:1574 (doi:10.1093/nar/24.8.1574): abstract unavailable in PubMed and PMC; the finding is quoted by López-Barragán et al. 2011 (doi:10.1016/j.molbiopara.2010.11.013, `ampseq/dl/PMC3026866.txt`): "Use of a 60°C extension temperature has been shown to be necessary in order to amplify regions with AT content 90% or higher because the DNA segments are already denatured at a 72°C extension temperature". López-Barragán compared 60, 65 and 70 C extension (Phusion master mix, anneal 65 C 30 s) and found better coverage of AT-rich regions at 60 C.
- NOMADS8/16 (Girgis et al. 2024 Nat Commun, doi:10.1038/s41467-024-45688-z, `ampseq/dl/PMC10869361.txt`): amplicons 79% AT; "reduced annealing and extension temperatures" citing Su 1996; total primer 0.6 µM split equimolar, then doubled or halved per primer from first-run coverage. Dimer screen score < -6 flagged. Exact temperatures not found in the text (**UNVERIFIED**). The `multiply` config (`ampseq/gh/default.json`) uses primer3 Tm 57 to 63 (opt 60), GC 20 to 80%, size 18 to 27.
- Oyola et al. 2012 BMC Genomics (doi:10.1186/1471-2164-13-1, `ampseq/dl/PMC3312816.txt`): for P. falciparum libraries KAPA HiFi with 60 mM TMAC gave the best AT-rich coverage; most other enzymes were inhibited by TMAC.
- MAD4HatTeR (Aranda-Díaz et al. 2025 Sci Rep, doi:10.1038/s41598-025-94716-5) used Paragon CleanPlex chemistry; primer rules are proprietary.

### 1.4 How many nucleotides does a target Tm cost at a given GC? (run)

From `ampseq_exp/neb_calib.out`: length needed for random sequence to reach the NEB-Q5 Tm (500 nM convention), shown as p10-median-p90 (search floor 14 nt, ceiling 45 nt):

| Q5 Tm | GC 20% | 30% | 40% | 50% | 60% | 65% | 70% |
|---|---|---|---|---|---|---|---|
| 58 | 22-29-38 | 19-24-30 | 16-21-25 | 14-18-22 | 14-15-19 | 14-14-17 | 14-14-17 |
| 60 | 25-34-43 | 21-26-35 | 16-21-27 | 14-18-23 | 14-16-20 | 14-15-18 | 14-14-17 |
| 62 | 28-38-45 | 21-28-37 | 18-23-30 | 15-19-25 | 14-17-21 | 14-16-19 | 14-14-18 |
| 65 | 34-45-45 | 23-33-44 | 19-25-34 | 16-22-28 | 14-18-23 | 14-17-21 | 14-16-19 |

This is the physical reason AT-rich panels lower Tm: holding 65 C at 20% GC needs 45+ nt primers. Above 40% GC, 65 C is reachable at 19 to 25 nt.

### 1.5 Proposed continuous mapping (my proposal, built from 1.1 to 1.4)

Convention: NEB-Q5 Tm at 500 nM, computed as in 0.1. The tool should also report Tm at the actual per-primer pool concentration.

```
Tm_opt(G)  = 58 + 7 * clip((G - 0.20)/0.20, 0, 1)      # 58 C at <=20% GC, 65 C at >=40% GC
window     = Tm_opt +/- 1.5 C  (pool-wide spread target <= 3 C)
primer GC  = [clip(g_loc - 0.15, 0.15, 0.45), clip(g_loc + 0.15, 0.50, 0.70)]   # g_loc = 100-300 bp local GC around the site
length     = <30% GC: 24-30-40 ; 30-40%: 20-25-34 ; 40-60%: 18-22-30 ; >60%: 16-19-25 (min-opt-max)
anneal     = Tm_lower(at entered conc) + 1 C for Q5 (NEB API), capped at 72; two-step if anneal >= ~65
extension  = 60 C if G < 0.30 (Su 1996; López-Barragán 2011), else 72 C (or combined with anneal in two-step)
```

Output of `ampseq_exp/gc_mapping.py`:

```
P. falciparum    GC 0.19: tm 58.0 (56.5-59.5), tm@25nM 54.4, primer_gc 0.15-0.50, length 24/30/40, anneal 59.0, ext 60
SARS-CoV-2       GC 0.38: tm 64.3 (62.8-65.8), tm@25nM 60.7, primer_gc 0.23-0.53, length 20/25/34, anneal 65.3, ext 72
human            GC 0.41: tm 65.0 (63.5-66.5), tm@25nM 61.4, primer_gc 0.26-0.56, length 18/22/30, anneal 66.0, ext 72
E. coli          GC 0.51: tm 65.0 (63.5-66.5), ...         primer_gc 0.36-0.66, length 18/22/30
M. tuberculosis  GC 0.66: tm 65.0 (63.5-66.5), tm@25nM 61.4, primer_gc 0.45-0.70, length 16/19/25, anneal 66.0, ext 72
```

Checks against evidence: ARTIC (38%) actual NEB Tm median 66.5 vs 64.3 proposed (within ARTIC's spread 64 to 69); primerJinn MTB 65 vs 65; P. falciparum lower Tm with 60 C extension matches Su/López-Barragán and NOMADS practice. The 58 C floor and the 20 to 40% ramp are my interpolation, **not directly measured** in a published panel. The user's "55 to 60 C otherwise" is reproduced only in a low-salt primer3 convention (section 0) and should not be applied in the Q5 convention, where it would put primers 6 to 8 C below the ARTIC practice.

### 1.6 Concrete defaults for three genome classes (Q5 convention, 500 nM)

| Class | Tm opt (window) | Length min/opt/max | Primer GC | Max 3' end dG (5 nt) | Anneal | Extension |
|---|---|---|---|---|---|---|
| High GC (>60%), e.g. MTB | 65 (63.5 to 66.5) | 16/19/25 | 45 to 70% | 9 kcal/mol | 66 C two-step, or 65 to 66 C + 72 C | 72 C |
| Mid (40 to 60%) | 65 (63.5 to 66.5) | 18/22/30 | g_loc +/- 15 pp, 30 to 65% typical | 9 | 65 to 66 C two-step | 72 C or combined |
| Low (<40%) | 58 to 64 by formula | 20 to 24 / 25 to 30 / 34 to 40 | 15 to 55% | 9 | Tm_lower + 1 | 60 C when <30% GC, else 72 C |

---

## 2. Determinants of coverage uniformity

1. **Small efficiency differences compound.** Gimpel et al. 2025 Nat Commun (doi:10.1038/s41467-025-64221-4, `ampseq/dl/PMC12533003.clean.txt`): "a template with an amplification efficiency just 5% below the average will be underrepresented by a factor of around two after only 12 PCR cycles". Relative abundance after n cycles is ((1+E1)/(1+E2))^n.
2. **Sequence-specific self-priming against the adapter, not GC, explains the worst amplifiers.** Gimpel: 12,000 random sequences with truncated TruSeq adapters, 90 cycles; about 2% of sequences amplified very poorly (efficiencies "as low as 80% relative to the population mean"). A pool fixed at 50% GC (GCfix) skewed as much as the unconstrained pool. A Lasso logistic regression on GC and base frequencies performed close to random. A 1D-CNN reached AUROC 0.88, AUPRC 0.44. CluMo found motifs (CGTG-derived, for example TCGTGT) complementary to the 5' end of the adapter that cause hairpin self-priming; a 6-nt motif complementary to the adapter dropped efficiency from 94.0 +/- 1.9% to 78.8 +/- 2.1% (5' side) or 81 +/- 4% (3' side). Model-based filtering cut required sequencing depth for 99% recovery fourfold. A conventional filter (GC 40 to 60%, homopolymer <= 4, mfold dG > -15 kcal/mol) removed 3,413 of 10,000 sequences but did not prevent the bias. Caveat: random synthetic inserts with fixed adapters; transfer to genomic amplicons with gene-specific primers is plausible but not shown.
3. **Primer dimers drive amplicon dropout in multiplex.** Itokawa 2020: ARTIC amplicons 18 and 76 dropped out because 18_LEFT and 76_RIGHT shared 10 nt 3' complementarity; replacing 12 primers (N1 set) evened coverage. SADDLE (Xie et al. 2022 Nat Commun, doi:10.1038/s41467-022-29500-4): dimer fraction in a 96-plex fell from 90.7% (naive) to 4.9% (optimised).
4. **Primer binding strength matching.** SADDLE: "primers that hybridize to their cognate templates with ΔG° ≈ −11.5 kcal/mol have the best tradeoff between amplification efficiency/uniformity and nonspecific hybridization", candidates truncated to dG between -10.5 and -12.5; Olivar's SADDLE defaults: 60 C, salinity 0.18, max dG -11.8 (`papers/PMC11282221.txt`). SADDLE used a median 45 nM per primer, 17 cycles, Vent (exo-).
5. **PCR stochasticity at low input.** Kebschull and Zador 2015 NAR (doi:10.1093/nar/gkv717): "PCR stochasticity is the major force skewing sequence representation" for low-copy pools; template switches and polymerase errors had little effect on distribution. Relevant for low-burden clinical samples (paucibacillary TB).
6. **GC and thermocycler.** Aird et al. 2011 Genome Biol (doi:10.1186/gb-2011-12-2-r18): 10 cycles of Phusion with standard Illumina cycling depleted >65% GC loci to about 1/100 and <12% GC to about 1/10; plateau 11 to 56% GC on a fast-ramping cycler, 13 to 84% on a slow cycler (2.2 C/s). Extending denaturation (30 s to 3 min initial; 10 s to 80 s per cycle) and 2 M betaine helped high GC; AccuPrime Taq HiFi with 65 C anneal/extension gave the flattest profile.
7. **Other predictors.** Andreson et al. 2008 NAR (doi:10.1093/nar/gkn290): over 80,000 PCRs, the best failure model (GM1) used the number of genomic binding sites of the primer plus primer GC%, cutting failure from 17% to 6%. Benita et al. 2003 NAR (doi:10.1093/nar/gng101): over 1,438 exons, template GC was the dominant predictor; regionalised GC (21 bp windows over a 61% GC threshold) improved prediction to 84.3% specificity, 94.8% sensitivity. pcrEfficiency (Mallona 2011, doi:10.1186/1471-2105-12-404): primer self-complementarity, dimers, GC imbalance between primers and primer GC were significant univariately.
8. **Amplicon length.** Within 150 to 450 bp, I found no quantitative per-bp efficiency coefficient in the papers above. Gimpel lists amplicon length among reported causes without a coefficient. **UNVERIFIED**: magnitude of length effect at 300 +/- 75 bp with 30 to 60 s extension. Keeping lengths within the user's 150 bp band and using long combined anneal/extension (ARTIC 5 min) is the practical mitigation.
9. **Empirical rebalancing works.** NOMADS doubled or halved individual primer concentrations after a first run at fixed total. Colman 2015 used 100 to 400 nM per pair. The tool should output a rebalancing worksheet (coverage ratio to new concentration).

---

## 3. Amplicon GC window and extreme GC

- Evidence: Aird plateau 11 to 56% GC (fast ramp, standard Phusion protocol) and up to 13 to 84% with slower ramp or longer denaturation; Benita 61% GC 21 bp window threshold.
- In high GC genomes every amplicon is above Aird's 56% plateau edge, so the relevant risk is local GC peaks, not mean GC. In MTB (65.6%), use relative limits.
- Proposed rules (my synthesis):
  - Amplicon mean GC within genome GC +/- 10 percentage points, hard limits 20% and 75%.
  - No 21 bp window above max(80%, genome GC + 20 pp); no 50 bp window below 10% GC unless the genome is <25% GC.
  - Penalise amplicon ensemble MFE at the annealing temperature (ViennaRNA DNA params) near primer sites; the user already plans this.
  - Recommend Q5 High GC Enhancer or 1 M betaine for high-GC panels (Colman used 1 M betaine with Q5 on MTB; NEB: enhancer "can improve reaction performance" for GC-rich templates); longer denaturation (Aird) for >70% GC amplicons. For <25% GC panels: 60 C extension; TMAC 60 mM only with KAPA HiFi (Oyola).

---

## 4. 5' adapter tails

Vendor rule, Illumina 16S Metagenomic Library Prep guide (Part 15044223 Rev. B, `ampseq/web/16s.txt`): "The locus‐specific portion of primer (not including overhang sequence) must have a melting temperature (Tm) of 60°–65°C ... For the Tm calculation only, the gene‐specific portion must be used in calculation. For hairpin and dimer calculations, the fully‐assembled primer sequence (including the overhang) should be used." Overhangs: forward TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG, reverse GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG; first PCR 55 C anneal with KAPA HiFi.

Run (`ampseq_exp/tails_conc.out`, primer3 mv 50, dv 2, dNTP 0.8, 50 nM):

```
TruSeq F + GAGGAACACCACTAGTACCG: gs 55.5 C, full 76.3 C, delta +20.8
TruSeq F + CTCGATGACTTTACGGCCAT: gs 56.2 C, full 76.4 C, delta +20.3
Nextera F + GAGGAACACCACTAGTACCG: gs 55.5 C, full 76.3 C, delta +20.8
tail-tail dimers at 60 C: TruF-TruR dG -3.55 (Tm 12.0 C); NxtF-NxtR dG -0.78
```

Mechanism: from cycle 3 onward, products carry the tail complement, so the full 53 to 54 nt primer binds with Tm about 20 C higher than the gene-specific part. All amplicons then anneal efficiently whatever their gene-specific Tm, which partly compresses early efficiency differences but also means early-cycle differences are locked in. This is why some protocols run a few cycles at a Ta set by the gene-specific Tm, then raise Ta (for example the Frontiers 2022 ARTIC variant used 65 then 62.5 C, the opposite direction; tail-specific step-up protocols are **UNVERIFIED** here).

Tails and dimers (run, `ampseq_exp/tails2.out`, 3,000 random 20-mer pairs per GC level, 3'-anchored dG from primer3 `calc_end_stability`, 37 C):

```
GC 0.35: p1 untailed -3.91 vs tailed partner -6.07; frac<-5: 0.003 vs 0.021
GC 0.5 : p1 untailed -5.60 vs tailed partner -7.68; frac<-5: 0.015 vs 0.084
GC 0.65: p1 untailed -8.20 vs tailed partner -9.19; frac<-5: 0.095 vs 0.216
```

Adding the 33 nt tail to the partner raises the fraction of pairs with a strong 3'-anchored dimer 2.3 to 7 fold. Gimpel's self-priming mechanism is the amplicon-level version of the same problem (insert motifs complementary to the adapter). **Decision: compute Tm on the gene-specific part only; compute all dimers, hairpins and amplicon self-priming checks on the full tailed oligo and tailed amplicon.** This matches Ultiplex, which checks hairpins and dimers with tags attached (`papers/PMC8600765.txt`).

---

## 5. Cycling and primer concentration

### 5.1 Annealing temperature rules
- Q5: Ta = Tm_lower + 1 C (NEB calculator API notes, 2026), or +3 C (NEB protocols.io text, 2015); two-step when Ta >= 72 C per the protocol, but ARTIC runs two-step at 65 C with 5 min holds and Q5.
- Taq: Ta = Tm_lower - 5 C (NEB API, observed).
- Touchdown PCR: Don et al. 1991 NAR 19:4008 (doi:10.1093/nar/19.14.4008), designed to suppress spurious priming. I found no multiplex-amplicon study comparing touchdown with two-step for coverage evenness (**gap**). Rationale against touchdown in high-plex pools: different primers become active in different early cycles, which converts Tm spread into abundance spread. Rationale for: suppresses mispriming in high-GC genomes. Recommend two-step as default, touchdown as an option.

### 5.2 Primer concentration
- Typical per-primer concentrations in high-plex pools: ARTIC 15 nM (LoCost), 15 to 20 nM (Quick 2017); SADDLE median 45 nM; NOMADS 0.6 µM total over the pool; Colman 100 to 400 nM in a 6-plex; NEB single-plex recommendations: Q5 500 nM, Taq 200 nM.
- Effect on Tm (NN theory, Tm = dH / (dS + R ln(C/4)) for non-self-complementary duplexes in primer3; verified numerically): primerJinn pair 1 primer3 Tm 59.0 C at 500 nM, 55.8 at 50 nM, 54.1 at 15 nM (drop 4.9 C). NEB Q5: 64.2 at 500 nM, 60.9 at 50 nM, 59.9 at 25 nM; mean drop 500 to 25 nM over 240 oligos 3.58 C.
- Fraction of template sites occupied at equilibrium, theta = K[P]/(1 + K[P]), K = exp(-dG(T)/RT), primer in excess (`tails_conc.out`, primer3 mv 50 conditions, so absolute numbers are low relative to Q5 buffer; relative effect holds):

```
GAGGAACACCACTAGTACCG T=60C ['15nM:0.168', '50nM:0.403', '500nM:0.871']
GAGGAACACCACTAGTACCG T=65C ['15nM:0.006', '50nM:0.020', '500nM:0.173']
```

At pool concentrations the duplex is far from saturated at Ta, so per-cycle priming is kinetically limited; this is why ARTIC uses 5 min combined anneal/extension. Design consequence: match primers on dG at the real Ta and real concentration (equivalently Tm at real concentration), not on Tm at 500 nM, when concentrations differ between primers. When concentrations are rebalanced (NOMADS, Colman), recompute.

---

## 6. 3' end rules: evidence vs folklore

| Rule | Evidence | Verdict |
|---|---|---|
| Avoid 3' complementarity between pool primers | Itokawa 2020 (10 nt 3' overlap caused dropout); SADDLE; PrimerROC | Strong. Hard constraint. |
| Limit 3' end stability (dG of last 5 nt) | primer3 PRIMER_MAX_END_STABILITY; Primer3Plus default 9.0 kcal/mol; range 0.86 (TATAT) to 6.86 (GCGCG) under SantaLucia 1998 (my implementation reproduced both values) | Mechanistically sound (mispriming tolerance); threshold is a convention. Note 9.0 does not bind at SantaLucia values (max 6.86); use about 6.5 if you want it to bite, or penalise continuously. |
| Max GC in last 5 nt | primer3 PRIMER_MAX_END_GC default 5 (off); commonly 3 | Convention. Use <= 3 for high-GC genomes, as a soft penalty. |
| GC clamp (1 to 2 S at 3' end) | Onodera and Melcher 2004 (doi:10.1016/j.mcp.2004.05.007): in 2,137 published viral PCR primers all top-16 3' triplets ended in S; AGG 3.27% vs TTA 0.42%; least used included TTA, TAA, CGA, ATT, CGT, GGG | Weak (survivorship of published primers, partly driven by design software defaults, as the authors note). Use as a mild preference, not a filter. |
| Avoid 3' terminal T | same review of recommendations; evidence only from usage frequency | Folklore-level. |
| Poly-X <= 4 to 5 | PrimalScheme MAX_HOMO 5; primer3 MAX_POLY_X 5; Gimpel conventional filter used <= 4 | Convention; also helps sequencing (homopolymer errors, especially nanopore). Use <= 4 for low-GC genomes where long A/T runs are common, relaxed to 5 if needed. |
| Dinucleotide repeats | No primary study found in this session | **UNVERIFIED**; treat as a low-complexity soft penalty (Olivar uses a sequence complexity risk term). |
| Avoid SNPs under 3' end | PrimalScheme mismatch penalties [3,2,1] by distance from 3'; multiply removed 11 SNP sites; Quick 2017 | Strong for pathogen panels; include a variant mask. |
| Genomic uniqueness of the 3' end | Andreson 2008: binding-site count was the most important failure predictor | Strong. Include 3'-anchored off-target count against host (human) and target genome. |

---

## 7. Open questions

1. Q5 annealing offset: calculator +1 C vs protocol +3 C; and actual Q5 buffer composition (undisclosed). A Ta gradient on the first pool should settle it.
2. Whether Deeplex's 55 C anneal implies low-Tm primers or just Taq-type chemistry; manual does not say.
3. No quantitative amplicon-length efficiency coefficient for 150 to 450 bp found.
4. Gimpel's CNN was trained on synthetic random inserts with fixed TruSeq stubs; untested on genomic amplicons. Code availability for direct use not checked.
5. The 20 to 40% GC ramp in the Tm mapping and the P. falciparum floor (58 C Q5 convention) need a wet-lab or literature anchor (NOMADS temperatures not found in text).
6. Two-step vs touchdown for evenness in high-plex: no head-to-head data found.
