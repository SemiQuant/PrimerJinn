# Multiplex designer: design decisions

This records what the new designer (`primerJinn/multiplex/`, command `primerjinn-multiplex`) does and why. Evidence is in `docs/research/` (one report per topic plus an independent fact-check, `*.verified.md`). Where a fact-check corrected a report, the corrected value is used here.

## 1. Inputs

- **Genome**: local FASTA, assembly accession (GCF_/GCA_), nucleotide accession, or organism name. Resolved through NCBI Datasets v2 and E-utilities with no API key (`genome.py`, from `docs/research/genome_fetch.md`). Ambiguous names stop with a list of choices. Cached under `~/.cache/primerjinn/genomes` (`PRIMERJINN_CACHE` overrides).
- **Positions**: 1-based inclusive. Bare integers, `seqid:pos`, `start-end`, tables with a header, XLSX, BED (0-based, converted). Bare integers on a genome with several chromosomes are refused.

## 2. Polymerase and Tm scale

Default enzyme: a **generic Taq master mix** (preset `kapa`; the user's choice on 2026-11-11 was KAPA2G Robust HotStart ReadyMix, and the preset's numbers come from its data sheet, Roche/KAPA TDS KR0382 v3.16; the user asked on 2026-11-12 that user-facing text call it a generic Taq master mix): engineered Taq-type enzyme without 3'-5' proofreading, A-tailed products, 2 mM MgCl2 and 0.2 mM of each dNTP at 1x, anneal 55/60/65 C for primer Tm <55/55-65/>65 C for 15 s, extend at 72 C for 15 s/kb, denature at 95 C (5 min initial for >70% GC), 5% DMSO for amplicons above 65% GC. Monovalent salt is not published; 50 mM is assumed.

Presets (`--polymerase`): `kapa` (default), `taq` (1.5 mM Mg), `q5`, `hifi`. For Taq-type enzymes the physical model describes primer binding directly, so the target Tm is on the **physical scale** (below). For Q5 the target is on the **NEB scale**, because Q5 carries a DNA-binding domain and NEB tells users to anneal above the primer Tm.

## 2a. Two Tm scales, kept apart

The fact-checks showed that "Tm" numbers from different tools differ by 2 to 6 C for the same primer, mostly because of hidden conventions. The tool therefore states which scale every number is on.

| Scale | Model | Used for |
|---|---|---|
| **NEB scale** (`tm_neb`) | SantaLucia 1998 NN + Owczarzy 2004 monovalent correction at 150 mM Na+, C = 500 nM, Mg ignored | The user-facing target Tm, the Tm window, the Q5 annealing rule. Reproduces the NEB Tm calculator (Q5) exactly: 240 primers at 500 and 25 nM, max error 0.000 C (checked 2026-11-11 against `docs/research/reference_code/neb_calib.json`). It is a vendor-calibrated scale, not a physical Tm. |
| **Physical** (`dg_dup`, `tm_phys`) | SantaLucia and Hicks 2004 NN + Owczarzy 2008 Mg2+/monovalent correction with the R decision tree, free Mg from 1:1 dNTP binding (Ka 3e4 M^-1, paper Eq 17), C = per-primer concentration (primer in excess) | Free energy of binding at the annealing temperature, fraction bound, uniformity. Coefficients and both worked examples in Owczarzy 2008 reproduced exactly. **Out-of-sample accuracy** (Owczarzy 2008 SI Table S2, 32 duplexes up to 30 bp, 50 Tm values not used for fitting): salt correction from a measured 1 M Tm MAE 0.90 C (paper 0.89 C); full prediction MAE 2.4 C, and in PCR-like buffers (20 to 150 mM monovalent, 0.5 to 5 mM Mg2+) MAE 1.6 C with a bias of +1.5 C (over-predicted). The error is nearest-neighbour error, largest for homopolymers and pyrimidine-rich strands; primer3-py on the same set: MAE 2.7 C (santalucia) and 4.2 C (owczarzy). Matches Biopython `Tm_NN(DNA_NN4, saltcorr=7)` within 0.005 C. |

DMSO (`--dmso`) shifts the duplex Tm by 0.6 C per % and is applied as an effective temperature shift to the ViennaRNA folds and the primer3 hairpin and dimer terms. The target Tm is not changed, so a DMSO design uses longer primers at a similar Ta; to keep shorter primers and anneal lower instead, lower `--tm` by about 0.6 C per %.

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

    dG_eff = dG_dup(Ta)                     physical NN model, real buffer, plus the template's
                                            dangling ends (SantaLucia and Hicks 2004 Table 3)
           + RT ln(1 + exp(-dG_hairpin/RT))  primer hairpin cost (primer3 thal at Ta)
           + dG_open(site)                   cost of opening template structure at the binding site
    theta  = 1 / (1 + exp(dG_eff/RT) / C)    fraction of template sites bound, primer in excess

Pool uniformity is scored on dG_eff (spread target <= 1 kcal/mol), on top of the Tm window. The Tm window (+/- 1.5 C) is a coarse filter and the Tm-centring weight is small; dG_eff decides, following SantaLucia 2007 (Methods Mol Biol 402:3, myth 3: match the amount bound at Ta, not Tm; Fig 5). Template dangling ends (one base beyond each primer end; +0.48 to -0.96 kcal/mol each, mean about -0.75 kcal/mol for the pair) are included because SantaLucia 2007 (p9) and SantaLucia and Hicks 2004 (p426) say they cannot be neglected for a primer binding a longer target; checked against primer3 thal (mean difference 0.02 kcal/mol). For Taq-type presets the recommended Ta is only known after design, so when it differs from the scoring Ta by more than 1 C the design is re-scored at the recommended Ta (at most 3 passes, each on a fresh copy of the settings); the report states the scoring Ta and warns if they still differ by more than 1 C.

**Binding balance during selection (2026-11-12).** Before optimisation, every candidate primer more than 0.3 kcal/mol weaker than the median of all candidates gets its best 5' extension (the length that brings its dG_eff at the scoring Ta closest to the median, same rules as below, ignoring dimers), and each pair using it gets an extended variant with recomputed length, within-pair dimer, self-dimer, length and uniformity costs. The unextended pairs stay candidates, so the optimiser weighs every extension against dimers and products with the rest of the pool. On M. tuberculosis this let the weak pncA forward primer be extended (a different rpsL pair was chosen to avoid the dimer that had blocked it), and the dG_eff spread at Ta fell from 2.28 to 1.22 kcal/mol (middle half 0.59 to 0.39); the recommended Ta then matched the scoring Ta, so one pass sufficed (94 s instead of 220 s).

**Binding balance by 5' extension after selection (user's proposal, 2026-11-11), now a final fine-tune.** A few primers bind weakly at Ta because their template site is folded (site-opening cost 2 to 2.5 kcal/mol, for example the M. tuberculosis pncA tiles). After the design and the Ta are known, each primer more than 0.3 kcal/mol weaker than the pool median keeps its Tm-matched core and gains template-matching bases at its 5' end, shown in **lowercase**, up to 8 nt (`--max_ext`, off with `--no_balance`). The length chosen is the one that brings its dG_eff closest to the median, recomputed each step with the longer site's opening cost, dangling ends and hairpin. An extension is refused if it adds a hard dimer with any primer in its pool, a run longer than `max_poly`, an N, or pushes the amplicon past `amp_max`; the 3' end and the covered positions never change. Specificity and every final check are then re-run on the extended primers. The report lists every extended primer and every weak primer that could not be extended, with the reason. Both the core Tm (inside the window) and the full-length Tm are reported. On the M. tuberculosis example one added base is worth about 1.2 to 1.7 kcal/mol, so only primers well below the median gain from it: 2 were extended, and the middle half of the pool sits within about 0.6 kcal/mol; the full range stays wider because one primer's extensions all form a dimer past the cap.

Reference points from the literature, translated into this model: SADDLE's candidate window (-10.5 to -12.5 kcal/mol at 60 C, SantaLucia and Hicks 2004 stacks, "salinity" 0.18, initiation 2.09, no terminal AT penalty; Xie et al. 2022 SI section S9) equals dG_dup(60 C) of about -10.4 to -12.4 kcal/mol here in the KAPA buffer, a 2 kcal/mol window at 30 to 45 nM primers with a 3 min anneal. Pythia (Mann et al. 2009), the closest equilibrium method, was conservative: 17 of 20 primer pairs it rejected still amplified, which supports not using absolute theta as pass/fail. theta is reported but not used as an absolute pass/fail because the monovalent content of commercial buffers (Q5, Phusion, KAPA) is not published, so absolute theta is uncertain by roughly 1 to 3 C equivalent.

Hairpin rule (fact-check corrected the report): penalise a cost above 0.5 kcal/mol, reject above 1.0 (a -1 kcal/mol hairpin costs 1.13 kcal/mol at 60 C). With adapter tails the hairpin is also scored on the tailed oligo, counting only what goes beyond the tail's own hairpin. Hairpin and site-opening costs enter dG_eff (and so the uniformity term) in full, and an extra penalty only above their thresholds, so they are not charged twice.

## 5. Secondary structure of templates and amplicons

Engine: ViennaRNA 2.7 with the **Mathews 1999 DNA** parameters (the bundled Mathews 2004 DNA set is non-physical above about 55 C), `noGU=1`, `dangles=2`, salt as Na_eq = mono + 3.3 sqrt(free Mg) (NUPACK's DNA rule). NUPACK is optional (academic licence) and only a cross-check.

1. **Primer site accessibility, genomic context** (early cycles): one RNAplfold-style windowed fold (window 200, span 150) per target region and strand at the design Ta gives P(site unpaired) for every candidate primer at once. dG_open = -RT ln P. Penalise above 1.5 kcal/mol, reject above 3.0.
2. **Amplicon strands** (later cycles), for the short-listed pairs: fold both strands at Ta (gene-specific amplicon without tails: the tails are the same on every amplicon and their panhandle would mask the template's own structure). These folds run across all cores for all targets at once. Recorded: site accessibility at the 3' end of each strand, probability that the 3'-terminal base is paired (self-priming, the main cause of poor amplifiers in Gimpel et al. 2025, doi 10.1038/s41467-025-64221-4), the strongest local hairpin, ensemble free energy, G-quadruplex motifs, highest 50 bp window GC, product duplex Tm. These feed the pair cost and the per-amplicon flags.
3. **Annealing temperature.** For Taq-type enzymes the recommended Ta is the model Ta: the highest temperature at which every primer is at least 90% bound on unfolded template (about the lowest physical Tm - 3 to 4 C at 100 nM), capped at the temperature where template and primer structure leave the weakest primer at 50% binding. Requiring 90% binding *with* structure was too strict: GC-rich templates fold more as the temperature falls, and the M. tuberculosis example dropped to 54 C against KAPA's 60 C; the capped rule gives 60 C. KAPA's own 55/60/65 C rule is shown alongside. For Q5 the NEB rule (lowest Tm at the reaction primer concentration + 1 C, minus 0.6 C per % DMSO) is used. If a primer site is still structured at the recommended Ta, the advice only suggests a higher Ta when the model shows that net binding improves (usually it does not; duplex binding weakens faster than structure opens), otherwise DMSO, betaine or redesign. GC-rich pools get DMSO advice (5% for KAPA) and longer initial denaturation; if GC-rich amplicons still fail, the advice escalates to the Musso et al. 2006 mix (1.3 M betaine, 5% DMSO, 50 uM 7-deaza-dGTP), which was the only combination that gave specific 68 to 79% GC products with Taq. AT-rich advice (60 C extension; Su et al. 1996: ~90% AT stretches amplified at 60 C but not 65 or 72 C, while a 70% AT region extended at 70 C) fires only when a 100 bp window is at least 85% AT.

All structure thresholds are physical starting points, not calibrated against sequencing depth. The genomic site-opening cost uses a windowed fold (window 200, base-pair span 150). Against a full-context fold of 600 nt it is on average 0.7 kcal/mol lower (up to 1.2), and a 400/300 window adds about 0.5. Its ranking of sites agrees well with the full fold (r = 0.98 on log P), so read the 1.5 and 3.0 kcal/mol thresholds as relative to this setting. Vienna and NUPACK disagree per site by up to 23-fold and by 3 to 14 kcal/mol at extreme GC, so structure is used to rank and to flag, not as a precise prediction.

## 6. Primer dimers

For every pair of primers in a pool (including self and within-pair), with primer3 thal at 37 C and the real salt (published thresholds are all on the 37 C scale):

    d3   = min(end_stability(a, b), end_stability(b, a))   3'-anchored
    dany = heterodimer(a, b)
    cost = max(0, -4 - d3)^2 + 0.25 max(0, -8 - dany)^2 + 1000 if d3 <= -6 or dany <= -12

Wet-lab anchor (SADDLE, Xie et al. 2022, p3-4): up to 3 nt of 3' complementarity gave no significant dimers, 5 nt at the 3' end did, 10 nt at the 5' end did not. With the thresholds here, about 42% of random 5 nt 3' overlaps are penalised and none of the pure A/T ones. SantaLucia 2007 (myth 4) notes that many "primer dimers" are in fact short products from nearby mismatched genomic sites, which the specificity step covers.

The -6 cap is relaxed to -7 only when the best choice actually contains a 3'-dimer past -6 (not when the remaining hard cost is a product). With adapter tails the full tailed oligos are scored (tails make d3 about 1 kcal/mol more negative). thal cannot score two sequences both over 60 nt; the tail is trimmed from the 5' end to fit.

## 6a. Strain variants

`--variants` reads known variant positions (VCF, gzipped VCF with deletions covering their REF bases and an optional INFO AF filter, BED, or any position format `--positions` accepts). A primer with a known variant within `variant_3prime_nt` = 5 nt of its 3' end is rejected: 3'-terminal mismatches block or bias extension (Huang 1992; Lefever 2013; to be checked against their full texts). A variant further from the 3' end costs `variant_penalty` = 3 per variant (about one Tm-window unit), so a clean primer is preferred but a target is not lost. Balance extensions stop before a variant. On M. tuberculosis, masking the 3'-terminal base of all 30 primers of a design replaced every one of them with all 24 positions still covered.

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

Per target, 60 candidate pairs are folded and scored; after the uniformity term the best 40 (`--candidates`) go to the optimiser. Every primer position and length is enumerated; for each 3' end the two cheapest lengths inside the Tm window that pass the sequence rules are kept. Objective = sum of per-pair costs + sum of within-pool pairwise costs, exactly as coded (candidates.py, design.py):

- per primer: ((Tm - T0)/1.5)^2 + 0.5 ((GC - GC_t)/0.10)^2 with GC_t = clamp(genome GC, 0.30, 0.65); +1.0 for 4 G/C in the last 5 nt, +0.5 for none; hairpin ((max(0, cost - 0.5))/0.5)^2; site opening ((max(0, cost - 1.5))/1.5)^2; self-dimer cost. Hard rules: only ACGT, GC 0.20 (lower for AT-rich genomes) to 0.80, runs <= max_poly, no 5x dinucleotide repeat, not 5 G/C in the last 5 nt, hairpin cost <= 1.0, site opening <= 3.0.
- per pair: the two primer costs + ((L - 300)/25)^2 + ((max(0, |GC_amp - GC_a| - 0.08))/0.05)^2 with GC_a = clamp(genome GC, 0.30, 0.60) + 0.5 ((Tm_F - Tm_R)/1)^2 + within-pair dimer cost + amplicon-strand structure (site opening above the genomic value, 3' self-priming above 0.1) + ((dG_eff - pool median)/0.5)^2 for each primer + off-target-site cost (2 per primable site, +10 per site scoring >= 0.5, at most 30) + product costs.
- pairwise (same pool only): dimer cost of all four primer combinations + product costs (1000 per likely unintended product, 2 per tolerated outer product).

Optimiser (PRISM, Chen et al. 2026, uses the same one-pair-per-target constraint and a single-swap local search, which our final 1-opt descent includes; its guarantee is for a maximisation objective and does not carry over): exhaustive when the product of candidate counts is at most 1e5; otherwise 30 randomised greedy starts, then simulated annealing restarts (single-target moves, geometric cooling) and a final 1-opt descent. On Mtb benchmarks this matched a proven MILP optimum at 16 targets and beat greedy by 3 times at 100 targets.

The weights are uncalibrated. They should be tuned once read depth from a first sequenced pool is available.

## 10. Outputs

`<out>.xlsx` and TSVs: primers (with every metric), amplicons, pairwise dimers, products, off-target sites, binding curves, PCR conditions; `<out>_order_sheet.tsv` (forward name, forward sequence, reverse name, reverse sequence); `<out>_amplicons.fasta`, `<out>.bed`, `<out>_report.md`. With `--sequencing`: `<out>_rebalance.tsv`, `<out>_rebalance_targets.tsv`, `<out>_conc.tsv`, `<out>_feedback.tsv/.md`, and after a redesign `<out>_order_sheet_new.tsv`.

## 11. Known limits

- Thresholds for dimers, structure and off-target products come from models and literature, not from our own wet-lab data. The Taq 3'-mismatch factors in the specificity model have not yet been checked against the full texts of Huang 1992, Lefever 2013, Stadhouders 2010 and Rejali 2018 (supplied, not yet read: the reading agent hit a usage limit).
- Equilibrium binding ignores kinetics; annealing steps may not reach equilibrium for structured GC-rich templates.
- Buffer monovalent content for commercial enzymes is unknown; physical numbers carry a 1 to 3 C uncertainty. A Ta gradient on the first pool is advised.
- Strain variation under primers is not handled (no SNP mask yet).
- Circular replicons are treated as linear: no amplicon across the origin, and products across it are not checked (targets near a sequence end get a warning).

## 12. Feedback from sequencing (2026-11-12, simplified 2026-11-13)

`primerJinn/multiplex/feedback.py` reads one table of read counts per amplicon (`amplicon`, `reads`, optional `sample`; column aliases) produced by any pipeline, with rows named `unaligned` (also `off_target`, `dimer`) holding the reads that match no designed amplicon. The earlier dependence on the DrDx.TB run analysis tables (read fates, per-primer usage, suggestions) was dropped because that tool is private; the counts table is the whole interface. Per sample `fold = (reads + 0.5) / (median over amplicons with reads + 0.5)`, pooled over samples as the geometric mean; status good 0.5-2x, warn 0.2-5x, bad, drop (no reads). Amplicon rows are matched to targets by target id, primer name without `_F`/`_R`, or target name (case-insensitive; `A/B` names two merged amplicons).

Decision per target: **replace** when it has no reads (or no row); below `redesign_below` (0.2) it is weak, and the concentration used (a `conc_nM` column, `--conc` per primer, or the design default) decides: if `current / fold` (the undamped concentration that would reach the median) is within `conc_max` (500 nM) the primers are **raised**, otherwise the pair is replaced; a pair already at `conc_max` and under 0.5x is replaced too. Otherwise **rebalance**, suggestions never above `conc_max`. `--replace` / `--keep` override. Concentration rule for rebalanced primers: `new = current x clip((1 / fold) ^ alpha, 1 / cap, cap)` with alpha 0.5 and cap 4 (a damped step: a 4x over-amplifier is halved, not quartered, because read share is not linear in primer concentration). Warnings: more than 30 % of reads outside the designed amplicons, fewer than 200 reads per amplicon on average, total stock volume above half the reaction.

Redesign (`redesign.py`, `candidates.fixed_target`, `design._apply_previous`): kept targets enter the designer as a single fixed candidate pair built from the previous coordinates and scored by the same code paths as new candidates (duplex dG with template dangling ends, hairpin cost, self-dimer, template-site opening with the window widened to the longest kept primer, within-pair dimer), so uniformity, cross-target dimers, outer products and genome-wide specificity are all computed between kept and new primers. Kept primers are never extended by the balance passes and never dropped as redundant; a kept pair that no longer matches the genome fails loudly instead of being split or relaxed. Replaced targets exclude the 3' end positions of the failed primers and pay `avoid_penalty` (2.0) for overlapping their sites. Pool numbers of kept targets seed the DSATUR colouring. Target ids are deterministic for the same positions file, genome and amplicon options, so names and ids carry over; a kept id of the form `Txx.k` (a target the fallback had split) is rebuilt from the kept insert.

## 13. Long deletions: junction amplicons (2026-11-12)

A VCF deletion whose deleted length (END minus POS, POS being the anchor base) is above `--long_variant` (100 bp) becomes a junction site: the anchor base and the base after the deletion, with `gap` = the deleted interval (`genome.vcf_kind`, `parse_positions`). It is never grouped or merged (`targets.junction_targets`, ids `J01` ...; `targets.build_targets` is the one place that builds targets for the designer and the precheck). `candidates.design_target` enumerates the forward primer before the anchor base and the reverse primer after the base following the deletion (the usual margins apply to the junction bases), and computes the product length and sequence with the deleted bases removed (`Pair.gap`, `candidates.product_length`, `candidates.amplicon_seq`); every later consumer of the amplicon sequence (amplicon-strand folds, finalize's structure report, GC, balance, FASTA) goes through `amplicon_seq`. Template-site accessibility is still folded in the reference context; the mutant template differs only within a window of the junction, which the primers do not cover. `Pair.insert` stays genomic, so coverage of the junction site and the BED file are consistent. The wild-type product (primer to primer on the reference) is reported (`genomic_length`, `wildtype_product_bp`) and the advice explains that with a short extension it is not made, so absence of the product is wild type or dropout. The fallback split is skipped for junction targets. In one-pool mode junction targets are added after the layout; since they overlap amplicons inside the deletion they may force a second pool, which the design reports as a warning naming them rather than failing. `--vcf_snps_only` keeps only single-base records (`vcf_kind == "snp"`).
