# Multiplex designer: design decisions

This records what the new designer (`primerJinn/multiplex/`, command `primerjinn-multiplex`) does and why. Evidence is in `docs/research/` (one report per topic plus an independent fact-check, `*.verified.md`). Where a fact-check corrected a report, the corrected value is used here.

## 1. Inputs

- **Genome**: local FASTA, assembly accession (GCF_/GCA_), nucleotide accession, or organism name. Resolved through NCBI Datasets v2 and E-utilities with no API key (`genome.py`, from `docs/research/genome_fetch.md`). Ambiguous names stop with a list of choices. Cached under `~/.cache/primerjinn/genomes` (`PRIMERJINN_CACHE` overrides).
- **Positions**: 1-based inclusive. Bare integers, `seqid:pos`, `start-end`, tables with a header, XLSX, BED (0-based, converted). Bare integers on a genome with several chromosomes are refused.

## 2. Polymerase and Tm scale

Default enzyme: **KAPA2G Robust HotStart ReadyMix** (user's choice, 2026-11-11). From its data sheet (Roche/KAPA TDS KR0382 v3.16): engineered Taq-type enzyme without 3'-5' proofreading, A-tailed products, 2 mM MgCl2 and 0.2 mM of each dNTP at 1x, anneal 55/60/65 C for primer Tm <55/55-65/>65 C for 15 s, extend at 72 C for 15 s/kb, denature at 95 C (5 min initial for >70% GC), 5% DMSO for amplicons above 65% GC. Monovalent salt is not published; 50 mM is assumed.

Presets (`--polymerase`): `kapa` (default), `taq` (1.5 mM Mg), `q5`, `hifi`. For Taq-type enzymes the physical model describes primer binding directly, so the target Tm is on the **physical scale** (below). For Q5 the target is on the **NEB scale**, because Q5 carries a DNA-binding domain and NEB tells users to anneal above the primer Tm.

## 2a. Two Tm scales, kept apart

The fact-checks showed that "Tm" numbers from different tools differ by 2 to 6 C for the same primer, mostly because of hidden conventions. The tool therefore states which scale every number is on.

| Scale | Model | Used for |
|---|---|---|
| **NEB scale** (`tm_neb`) | SantaLucia 1998 NN + Owczarzy 2004 monovalent correction at 150 mM Na+, C = 500 nM, Mg ignored | The user-facing target Tm, the Tm window, the Q5 annealing rule. Reproduces the NEB Tm calculator (Q5) exactly: 240 primers at 500 and 25 nM, max error 0.000 C (checked 2026-11-11 against `docs/research/reference_code/neb_calib.json`). It is a vendor-calibrated scale, not a physical Tm. |
| **Physical** (`dg_dup`, `tm_phys`) | SantaLucia and Hicks 2004 NN + Owczarzy 2008 Mg2+/monovalent correction with the R decision tree, free Mg from 1:1 dNTP binding, C = per-primer concentration (primer in excess) | Free energy of binding at the annealing temperature, fraction bound, uniformity. Measured accuracy: MAE 1.25 C on 342 measured Tms (Owczarzy data, in-sample for salt). Matches Biopython `Tm_NN(DNA_NN4, saltcorr=7)` within 0.005 C. |

Things the research proved wrong or risky, now avoided:
- primer3-py 2.0.3 `calc_tm` defaults to 0.8 M formamide (about -2 C on every Tm). We never use primer3 `calc_tm`.
- primer3's Owczarzy branch drops the length term (C integer division). We compute Tm ourselves.
- The old PrimerJinn "Q5 approximation" (70 mM mono, 3300 nM) was off by -1.4 C (SD 1.3) against NEB.

## 3. Target Tm from genome GC

The user's rule, on the preset's scale (physical for KAPA/Taq, NEB for Q5; overridable with `--tm` and `--tm_scale`):

| Genome GC | Target Tm |
|---|---|
| >= 60% | 65 C |
| 55 to 60% | linear 60 to 65 C |
| 30 to 55% | 60 C |
| 20 to 30% | linear 55 to 60 C |
| < 20% | 55 C |

Candidates must lie within target +/- 1.5 C (`--tm_window`), so the pool Tm spread is at most 3 C. Primer length range by genome GC: < 30%: 20 to 40 nt; 30 to 40%: 18 to 34; 40 to 60%: 18 to 32; > 60%: 17 to 32. Primer GC must be at least min(0.20, genome GC - 0.08) (0.11 for P. falciparum) and runs at most 4 nt (5 below 25% genome GC); both are CLI options and are loosened once more for targets that still fail.

## 4. Binding strength at the annealing temperature (the modern part)

Tm matching alone leaves a 2 C spread in real binding on M. tuberculosis rpoB primers; matching the free energy at the annealing temperature (Ta) halves it (equilibrium report). So, for every primer at the design Ta:

    dG_eff = dG_dup(Ta)                     physical NN model, real buffer
           + RT ln(1 + exp(-dG_hairpin/RT))  primer hairpin cost (primer3 thal at Ta)
           + dG_open(site)                   cost of opening template structure at the binding site
    theta  = 1 / (1 + exp(dG_eff/RT) / C)    fraction of template sites bound, primer in excess

Pool uniformity is scored on dG_eff (spread target <= 1 kcal/mol), on top of the Tm window. theta is reported but not used as an absolute pass/fail because the monovalent content of commercial buffers (Q5, Phusion, KAPA) is not published, so absolute theta is uncertain by roughly 1 to 3 C equivalent.

Hairpin rule (fact-check corrected the report): penalise a cost above 0.5 kcal/mol, reject above 1.0 (a -1 kcal/mol hairpin costs 1.13 kcal/mol at 60 C). With adapter tails the hairpin is also scored on the tailed oligo, counting only what goes beyond the tail's own hairpin. Hairpin and site-opening costs enter dG_eff (and so the uniformity term) in full, and an extra penalty only above their thresholds, so they are not charged twice.

## 5. Secondary structure of templates and amplicons

Engine: ViennaRNA 2.7 with the **Mathews 1999 DNA** parameters (the bundled Mathews 2004 DNA set is non-physical above about 55 C), `noGU=1`, `dangles=2`, salt as Na_eq = mono + 3.3 sqrt(free Mg) (NUPACK's DNA rule). NUPACK is optional (academic licence) and only a cross-check.

1. **Primer site accessibility, genomic context** (early cycles): one RNAplfold-style windowed fold (window 200, span 150) per target region and strand at the design Ta gives P(site unpaired) for every candidate primer at once. dG_open = -RT ln P. Penalise above 1.5 kcal/mol, reject above 3.0.
2. **Amplicon strands** (later cycles), for the short-listed pairs: fold both strands at Ta (gene-specific amplicon without tails: the tails are the same on every amplicon and their panhandle would mask the template's own structure). These folds run across all cores for all targets at once. Recorded: site accessibility at the 3' end of each strand, probability that the 3'-terminal base is paired (self-priming, the main cause of poor amplifiers in Gimpel et al. 2025, doi 10.1038/s41467-025-64221-4), the strongest local hairpin, ensemble free energy, G-quadruplex motifs, highest 50 bp window GC, product duplex Tm. These feed the pair cost and the per-amplicon flags.
3. **Annealing temperature.** For Taq-type enzymes the recommended Ta is the model Ta: the highest temperature at which every primer is at least 90% bound on unfolded template (the classical Tm - 4 to 5 C), capped at the temperature where template and primer structure leave the weakest primer at 50% binding. Requiring 90% binding *with* structure was too strict: GC-rich templates fold more as the temperature falls, and the M. tuberculosis example dropped to 54 C against KAPA's 60 C; the capped rule gives 60 C. KAPA's own 55/60/65 C rule is shown alongside. For Q5 the NEB rule (lowest Tm at the reaction primer concentration + 1 C, minus 0.6 C per % DMSO) is used. If a primer site is still structured at the recommended Ta, the advice only suggests a higher Ta when the model shows that net binding improves (usually it does not; duplex binding weakens faster than structure opens), otherwise DMSO, betaine or redesign. GC-rich pools get DMSO advice (5% for KAPA) and longer initial denaturation.

All structure thresholds are physical starting points, not calibrated against sequencing depth. Vienna and NUPACK disagree per site by up to 23-fold and by 3 to 14 kcal/mol at extreme GC, so structure is used to rank and to flag, not as a precise prediction.

## 6. Primer dimers

For every pair of primers in a pool (including self and within-pair), with primer3 thal at 37 C and the real salt (published thresholds are all on the 37 C scale):

    d3   = min(end_stability(a, b), end_stability(b, a))   3'-anchored
    dany = heterodimer(a, b)
    cost = max(0, -4 - d3)^2 + 0.25 max(0, -8 - dany)^2 + 1000 if d3 <= -6 or dany <= -12

The -6 cap is relaxed to -7 only when the best choice actually contains a 3'-dimer past -6 (not when the remaining hard cost is a product). With adapter tails the full tailed oligos are scored (tails make d3 about 1 kcal/mol more negative). thal cannot score two sequences both over 60 nt; the tail is trimmed from the 5' end to fit.

## 7. Specificity

3'-anchored seed search (three 8-mers at 0, 5 and 10 nt from the 3' end, one mismatch each) over both strands of the whole genome with numba, then primer3 thal scoring at Ta, then a polymerase factor (proofreading enzymes trim 3' mismatches, so they get no mismatch protection; Taq gets published-range factors). Recall 96 to 98% of primable sites, 5 to 13 times faster than BLAST. Calls: a site is primable at score >= 0.01; a product (any + site and downstream - site from any two primers in the pool, 40 to 2000 bp) is likely at >= 0.01. Products are classified (contig-aware):

- **unintended** (off-target sites or wrong pairings): cost 1000 if likely, so the optimiser avoids them;
- **outer product** (forward primer of one designed amplicon with the reverse primer of another, both at their designed sites): hard only if shorter than `--outer_min` (1000 bp), otherwise a small cost, because a short extension suppresses them (ARTIC relies on this);
- **repeated locus** (the designed pair also amplifies a second copy of the same length, such as rRNA operons or IS elements): cost 50, so it is avoided when possible, and reported with an explanation when it is not (reads from the copies cannot be told apart).

Candidate primers with off-target primable sites are penalised.

## 8. Targets, amplicons, pools

- Positions are grouped into targets by an exact dynamic programme (`targets.py`): fewest amplicons, tight spans avoided, cuts placed where neighbours need not overlap. Longer ranges are tiled.
- Amplicon length: 225 to 375 bp (hard), optimum 300, cost ((L - 300)/25)^2. Every amplicon is then within 150 bp of every other.
- Pools: a conflict graph has an edge when two targets are on the same sequence and close enough that their amplicons must overlap, or would form an outer product under 1000 bp (`--outer_min`). Tiles of one requested range conflict only with their neighbours, so a tiled range alternates between 2 pools as in ARTIC. DSATUR colouring gives the pools; one pool when there are no edges.
- A target that finds no pair is split in two (not three) and retried; anything still failing is retried with relaxed limits. After optimisation, amplicons that coverage does not need are dropped. Requested positions that are N in the reference are left out with a warning.

## 9. Optimisation

Per target, up to 60 candidate pairs (own enumeration of every primer position and length, filtered and ranked, then paired). Objective = sum of per-pair costs + sum of within-pool pairwise costs:

- per pair: ((Tm_F - T0)/0.75)^2 + ((Tm_R - T0)/0.75)^2 + ((dG_eff - median)/0.5)^2 per primer + ((L - 300)/25)^2 + ((GC_amp - GC_genome)/0.08)^2 + structure costs + within-pair and self dimer cost + off-target cost
- pairwise (same pool only): dimer cost of all four primer combinations + 1000 per unintended likely product

Optimiser: exhaustive when the product of candidate counts is at most 1e5; otherwise 30 randomised greedy starts, then simulated annealing restarts (single-target moves, geometric cooling) and a final 1-opt descent. On Mtb benchmarks this matched a proven MILP optimum at 16 targets and beat greedy by 3 times at 100 targets.

The weights are uncalibrated. They should be tuned once read depth from a first sequenced pool is available.

## 10. Outputs

`<out>.xlsx` and TSVs: primers (order sheet), amplicons, pools, pairwise dimers, off-target products, PCR conditions; `<out>_amplicons.fasta`, `<out>.bed`, `<out>_report.md`.

## 11. Known limits

- Thresholds for dimers, structure and off-target products come from models and literature, not from our own wet-lab data.
- Equilibrium binding ignores kinetics; annealing steps may not reach equilibrium for structured GC-rich templates.
- Buffer monovalent content for commercial enzymes is unknown; physical numbers carry a 1 to 3 C uncertainty. A Ta gradient on the first pool is advised.
- Strain variation under primers is not handled (no SNP mask yet).
- Circular replicons are treated as linear: no amplicon across the origin, and products across it are not checked (targets near a sequence end get a warning).
