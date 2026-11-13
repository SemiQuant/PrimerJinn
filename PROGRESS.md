# PrimerJinn multiplex designer: progress and handover

Read this first when picking the work up in a new session. Keep it current: update the status table and the log whenever a phase finishes or a decision changes.

Last updated: 2026-11-13 (session 3)

## Goal

Design a multiplex PCR primer pool for amplicon sequencing from two inputs: a reference genome name (or accession, or FASTA path) and a list of genome positions. The original PrimerJinn (`getMultiPrimerSet`, `PCRinSilico`) was the starting point; the new designer lives in `primerJinn/multiplex/` and runs as `primerjinn-multiplex` (or `python3 -m primerJinn.multiplex.cli`).

## Requirements (from the user)

| # | Requirement | How it is met |
|---|---|---|
| 1 | Input is a genome name and a list of positions | `--genome` takes a name, accession or FASTA (NCBI lookup, cached). `--positions` takes a file or inline list |
| 2 | All positions covered | Hard check: every position must lie between the primers of an amplicon. Failing targets are split in two, then retried with relaxed limits. Redundant amplicons are dropped |
| 3 | Amplicons ideally 300 bp | Optimum 300, cost ((L-300)/25)^2 |
| 4 | Amplicons within 150 bp of each other | Amplicon LENGTHS (confirmed by the user 2026-11-11): hard range 225 to 375 bp |
| 5 | Tm and GC similar | Tm window +/- 1.5 C, GC cost, and pool uniformity of binding free energy at Ta (dG_eff) |
| 6 | Target Tm from genome GC | 65 C at >= 60% GC, 60 C for 30 to 55%, ramps between, 55 C at <= 20%. `--tm` overrides |
| 7 | Modern, accurate thermodynamics | SantaLucia 2004 + Owczarzy 2008 (Mg, dNTP); NEB Q5 scale reproduced exactly for Q5; primer3 calc_tm avoided (hidden 0.8 M formamide default, integer-division bug) |
| 8 | Amplicon secondary structure informs design and PCR temperature | ViennaRNA (Mathews 1999 DNA parameters): template-site opening cost in genomic context, amplicon-strand site accessibility and 3' self-priming for every short-listed pair, structure-aware binding curves set the recommended Ta |
| 9 | Default polymerase: standard Taq, probably KAPA Robust ReadyMix (user, 2026-11-11); shown to users as a "generic Taq master mix" (user, 2026-11-12) | Preset `kapa` (key unchanged; label and advice say generic Taq master mix): 2 mM Mg, 0.8 mM dNTP, no proofreading, Tm on the physical scale, Ta from the model with the data sheet's 55/60/65 C rule shown alongside, cycling and 5% DMSO advice |
| 11 | Close the loop after sequencing: rebalance concentrations and redesign amplicons that fail, from a simple user-made read-count table (user, 2026-11-12; simplified to one counts file with an unaligned row on 2026-11-13 because the DrDx.TB run analysis is private) | `--sequencing counts.tsv` + `--previous` (feedback.py, redesign.py); README documents the format, including the unaligned rowget reads |
| 12 | GUI looks like the DrDx.Me dist apps, full DrDx.Me header (user, 2026-11-12) | gui.py restyled with the sq-dashboard tokens, bundled fonts and logo |
| 10 | Keep this document | This file |

## Status

| Phase | State | Notes |
|---|---|---|
| 1. Scout existing code and environment | Done | |
| 2. Research, fact-checked | Done | `docs/research/` (7 reports, each with an independent `*.verified.md`; `fact_check_digest.txt`) |
| 3. Design | Done | `docs/DESIGN.md` |
| 4. Implementation | Done | `primerJinn/multiplex/`, 83 tests pass (`python3 -m pytest tests`, about 2 min) |
| 5. Tests on real genomes | Done | M. tuberculosis, E. coli, S. aureus (with Illumina tails), P. falciparum (incl. AT-rich intergenic locus), SARS-CoV-2 spike tiling. All exit 0, all positions covered |
| 6. Independent review and fixes | Done, one round | 4 genome testers + 4 code reviewers, each finding checked by a separate skeptic; all confirmed findings fixed (see log). The fixes themselves were re-tested on all genomes but not re-reviewed by independent agents |
| 7. Feedback from sequencing | Implemented, not calibrated | `--sequencing counts.tsv` rebalances and redesigns from read counts per amplicon; thresholds (0.2x weak, 500 nM conc_max, alpha 0.5, cap 4) are judgement calls until a real pool has been through the loop |

## How to run

```bash
python3 -m primerJinn.multiplex.cli --genome example/ref.fasta --positions example/mtb_positions.tsv --out results/mtb
```

Outputs: `<out>.xlsx`, TSV tables, `<out>_amplicons.fasta`, `<out>.bed`, `<out>_report.md` (conditions, advice, warnings).

Results on 2026-11-11 (KAPA default, 100 nM each, 16 cores):

| Genome | Input | Result | Ta | Time |
|---|---|---|---|---|
| M. tuberculosis H37Rv (65.6% GC) | 24 resistance sites | 15 amplicons, 2 pools | 60 C (KAPA rule 60) | 116 s |
| E. coli MG1655 (50.8%) | 12 codons + 801 bp range | 11 amplicons, 2 pools | 56 C (KAPA 60) | 51 s |
| S. aureus NCTC 8325 (32.9%), Illumina tails | 17 codons | 6 amplicons, 1 pool | 57 C | 33 s |
| P. falciparum 3D7 (19.3%) | 15 codons on 5 chromosomes | 7 amplicons, 2 pools | 51 C | 43 s |
| SARS-CoV-2 (38%) | 3.8 kb spike range + 5 sites | 19 amplicons, 2 alternating pools (ARTIC-like) | 56 C | 86 s |

## Code map (`primerJinn/multiplex/`)

| Module | Role |
|---|---|
| `cli.py` | Command line, position parsing glue |
| `config.py` | Settings, polymerase presets (kapa, taq, q5, hifi), tails, GC-dependent defaults |
| `genome.py` | Genome name/accession to cached FASTA (NCBI, no key); position parser validated against the genome |
| `targets.py` | Group positions into targets (exact DP), coverage check |
| `thermo.py` | NN thermodynamics, salt corrections, NEB scale, fraction bound |
| `structure.py` | ViennaRNA DNA folding: accessibility, self-priming, hairpins, amplicon reports |
| `candidates.py` | Per-target enumeration and scoring of primers and pairs; amplicon fold jobs |
| `specificity.py` | Genome-wide seed search and thal scoring of off-target sites and products |
| `optimize.py` | Pool colouring (DSATUR, tiles alternate) and simulated annealing over pair choices |
| `design.py` | Orchestration, product classification, redundancy pass |
| `finalize.py` | Coverage, products, dimers, Ta advice, amplicon reports |
| `report.py` | TSV, XLSX, FASTA, BED, order sheets, markdown |
| `feedback.py` | Sequencing feedback: read the counts table (amplicon, reads, optional conc_nM and sample; unaligned row; column aliases), match rows to targets, decide keep/raise/lower/replace, damped concentrations, summary and rebalance files |
| `redesign.py` | Previous `_primers.tsv` -> fixed pairs for kept targets, 3' ends and sites to avoid for replaced ones (`Config.fixed/avoid_3prime/avoid_sites/prev_pools`) |
| `gui.py`, `static/` | Streamlit GUI in the DrDx.Me dashboard style (fonts, logo, style.css bundled; `.streamlit/config.toml` theme) |

## Decisions worth knowing

- Tm scale depends on the polymerase: physical (SantaLucia 2004 + Owczarzy 2008 in the stated buffer at the per-primer concentration) for Taq-type enzymes; NEB calculator scale for Q5 (Q5 has a DNA-binding domain, and NEB's Ta rule sits above Tm). Both values are always reported.
- Recommended Ta (Taq-type): every primer >= 90% bound on unfolded template, capped where structure leaves the weakest primer at 50%. Requiring 90% binding with structure was too strict (54 C vs KAPA's 60 C on M. tuberculosis).
- Product kinds: unintended (hard), outer product between designed amplicons (soft unless under 1000 bp and not tiles of one range), repeated locus (cost 50, explained).
- Hairpins in amplicons are reported as numbers, not flags (no validated threshold; almost every GC-rich amplicon has one stable at 72 C).
- Thresholds for dimers (-4 soft, -6 hard 3'-anchored dG at 37 C), template structure (1.5 penalise, 3.0 reject kcal/mol), off-target score (0.01) and the objective weights are model or literature starting points.

## Open questions for the user

- Answered 2026-11-11: "within 150 bp" means amplicon lengths; 100 nM per primer is right.
- Wet-lab pilot: "maybe sometime" (user). Until then thresholds stay model-based.
- Papers from `docs/papers_needed.md` were supplied as PDFs in `/Users/semiquant/Downloads/papers` (not copied into the repo).

## Next steps

0. Read the four Taq 3'-mismatch papers in `/Users/semiquant/Downloads/papers` and check the extension factors in `specificity.py` (Params.taq_*).
1. Wet-lab pilot: run the M. tuberculosis pool with a Ta gradient (57 to 63 C) and 15 s vs 60 s annealing; use read depth to calibrate the structure, dimer and uniformity weights.
2. Circular replicons (products across the origin).
2b. Run a real sequenced pool through `--sequencing counts.tsv` (with the `conc_nM` column) and tune the feedback thresholds (`--redesign_below` 0.2x, `--conc_max` 500 nM, alpha 0.5, cap 4) and the optimiser weights against the observed depth.
3. Optional second independent review of the fixes made after validation.
4. Fix the `PrimerJinn/` vs `primerJinn/` case mismatch in git before publishing.

## Repo notes

- git tracks the package as `PrimerJinn/` but on disk (case-insensitive macOS) it is `primerJinn/`; `setup.py` uses `primerJinn`. A fresh clone on Linux would not import. Fix with `git mv` when convenient.
- Nothing has been committed yet in this work.
- Validation scripts and outputs from session 1 were in the session scratchpad (temporary); the evidence that matters is summarised here and in `docs/research/`.

## Log

- 2026-11-12 (session 3, late): The WHO-catalogue VCF (3357 records: 795 SNPs, 2533 short indels, 29 long deletions, 1864 distinct positions) under `--one_pool` hit the layout error (katG and ethA tiles must overlap; cannot merge under 1000 bp) and the GUI's design then stopped (precheck errors block the design). Changes: `one_pool_layout` no longer raises; blocking targets stay separate and go to a second pool with a "NOT ONE POOL" note (precheck WARN, design warning naming them); `--one_pool_drop` drops them instead (listed as not designed); `--exclude` (CLI) and, in the GUI, a "Leave out targets" multiselect fed by the last Check positions plus a free-text names/loci box (`precheck.run(exclude=...)`, `_apply_exclude`: names, loci or target ids of a provisional grouping). The GUI keeps the last check result in session state so the table stays while the user picks. test_one_pool updated (fallback instead of raise), test_precheck exclude test; 80 tests pass. `OnePoolError` is kept but no longer raised by the layout.
  `--vcf_max_indel N` (GUI Advanced): VCF records whose alleles differ in length by more than N bp are left out (`genome.vcf_indel_length`; symbolic alleles count their span); WHO file with 50 bp, one pool, cap 2500: no conflicts, one tube.
  Pool forecast in the check (`precheck._pool_forecast`: DSATUR on the conflict graph, WARN naming the targets that need extra pools and why). Junction targets whose deletion lies inside a standard target are not made (`targets.build_targets`, note). User's WHO VCF with --one_pool: cap 1000 -> 7 pools forecast (8 junctions + 2 dense tiles); cap 2500 -> 29 targets, 3 pools, only the two katG deletions (2337 and 931 bp, katG split in two amplicons) still conflict.
  Also fixed: `optimize.conflict_graph` stopped scanning once the outer product from a target's START was long, so a 2.3 kb junction target did not conflict with tiles starting inside it (a 200 bp outer product in one tube on the WHO file); the scan now stops only when a later target can neither overlap nor give a short outer product. WHO file without one pool: 41 amplicons, 6 pools (nested katG junction amplicons conflict with every tile inside them), 440 s.

- 2026-11-12 (session 3, later): Long deletions and VCF SNP filter. Prompted by a WHO-catalogue VCF with 0.4 to 2.3 kb katG/ethA/fbiA deletions that `--one_pool` could not lay out (tiles of a deleted region cannot share a tube, and tiling only shows a deletion as dropout). `--long_variant N` (default 100): a VCF deletion longer than N becomes a junction amplicon (`Site.gap`, `Target.gap`, ids J01..; `genome.vcf_kind`, `targets.junction_targets/build_targets`): primers on both flanks, product length and sequence computed on the deleted allele (`Pair.gap`, `candidates.product_length/amplicon_seq`, used by folds, finalize, balance, FASTA), wild-type product size reported (`wildtype_product_bp`), advice on reading it, BED name `_junction`, redesign round-trips the `deletion` column. `--vcf_snps_only` keeps single-base records. Precheck reports VCF record kinds and what the options did. GUI: checkbox and threshold under Positions. E2E on example/ref.fasta with a synthetic VCF (2 SNPs, 420 bp and 900 bp <DEL>, a 12 bp indel): 5 amplicons, 1 pool, junction products 301 bp (721/1201 bp on wild type). Tests in test_targets.py and test_regressions.py. Not handled: long insertions/duplications/inversions (still ranges), nested deletions in one tube (each gets its own junction pair, and they conflict, so one-pool mode warns and uses a second pool). The user's VCF also had spans that do not match the HGVS names (c.794_907del over 2.3 kb); told the user to check the file.

- 2026-11-13: Missed positions are now listed with a reason: `Design.dropped_reasons`, `report.uncovered()` -> `<out>_uncovered.tsv` (also the coverage_problems), GUI expander "Positions not covered" and download. Before, one-pool-dropped sites were mislabelled as N in the reference.
- 2026-11-13: User saw "0 warnings" then 5 pools with one pool on: the GUI shows the last check result even after settings change. Now `ss.pc_args` records the args of the check and a warning says the check is stale when they differ. Junction targets are renumbered after dropping those inside a standard amplicon (J01.. contiguous) and the note names the covering target.
- 2026-11-13: Concentration-aware decision (user): the counts table takes an optional `conc_nM` column (or `--conc` per primer); a weak amplicon is raised when `current / fold` stays within `--conc_max` (500 nM, CLI and GUI Advanced) and redesigned otherwise; suggestions are capped at `--conc_max`. `Thresholds.conc_max`, `targets` table gains `current_nM`.
- 2026-11-13: Feedback input simplified (user: the DrDx.TB run analysis is private). `--sequencing` now takes one counts table (`amplicon`, `reads`, optional `sample`) with a row named `unaligned` (or `off_target`, `dimer`) for the reads outside the designed amplicons; per-primer off-target and dimer thresholds (`--offtarget_max`, `--dimer_max`), the zip/folder reader, sequence matching and the links to the analysis tool were removed from the CLI, GUI, README and DESIGN.md. `example/sequencing_run/` is now one `counts.tsv`. Tests: `tests/test_feedback.py` rewritten (9).
- 2026-11-12 (session 3): Sequencing feedback and redesign (input format since replaced, see 2026-11-13). `primerJinn/multiplex/feedback.py` reads the DrDx.TB run analysis output (`results.zip`, its `tables/` folder, or the TSVs; `amplicon_summary`, `read_fates`, `primer_usage`, `primer_suggestions`, `primer_dimers`; column aliases for other pipelines) and the sequenced design's `_primers.tsv` (+ optional `--conc`), matches units to targets by name (the job names amplicons after the forward primer minus `_F`, so ids round-trip; merged units `A/B` and `members`) or by sequence, pools folds over samples as the job does, and decides per target: replace (drop, design failure, < 0.2x median, >= 20% off-target or >= 5% dimer reads of a primer, `fix_primer`) or rebalance with the job's damped rule `(1/fold)^0.5` capped at 4x; writes `_rebalance.tsv`, `_rebalance_targets.tsv`, `_conc.tsv` (for the job's concentration slot), `_feedback.tsv/.md`. Redesign (`redesign.py`): kept targets become single fixed candidates (`candidates.fixed_target`, scored by the same code as new ones, incl. site opening with a widened window), never extended (`balance`, `_add_extended_variants`), never dropped as redundant, pools seeded into DSATUR (`optimize.dsatur(preset)`), split ids `Txx.k` rebuilt (`design._apply_previous`); replaced targets exclude the failed 3' ends (`Config.avoid_3prime`) and pay 2.0 for overlapping the failed sites. Every design now writes `_order_sheet.tsv` (4-column, the format the job and vendors take) and after a redesign `_order_sheet_new.tsv`. CLI: `--sequencing --previous --conc --replace --keep --no_redesign --redesign_below --offtarget_max --dimer_max --alpha --cap --reaction_ul --stock_uM`; `--genome/--positions` are optional for a rebalance-only run. Verified: 5-site MTB design -> synthetic run (one over, one weak, one dropout, one off-target primer) -> rebalance only (4 lowered, 3 to replace) -> redesign: kept primers byte-identical, pools kept, new 3' ends; order sheet and conc table parsed by the job's own `read_table`/`load_conc`. Real run tables (`example/sequencing_run/`, run_20260925) reproduce the job's suggested concentrations within 2%. Wording: "KAPA2G Robust HotStart ReadyMix" -> "generic Taq master mix" in labels, help, advice, report, README and DESIGN.md (preset key `kapa` and code unchanged). README: "After sequencing" section with the analysis link and the input format for other pipelines. GUI (gui.py rewritten by a subagent, 711 lines): DrDx.Me dashboard look from `sq-dashboard.css` tokens (dark slate, pale-green accent, Montserrat/Roboto/Roboto Mono bundled as woff2 in `primerJinn/multiplex/static/` with OFL.txt, DrDx.Me white wordmark, everything embedded as data URIs so the page stays offline), full header (logo, name, subtitle, All tools, OFFLINE pill), sidebar 'Inputs', stat tiles, sq-tables, pills, dark Altair theme; `.streamlit/config.toml` theme and `--theme.*` flags in `main()`; new 'From sequencing' sidebar section (results zip/tables, previous `_primers.tsv` defaulting to the session's last design, concentration table, replace/keep, thresholds in Advanced) with 'Rebalance' and 'Rebalance and redesign' buttons and a 'Sequencing feedback' results tab; `PRIMERJINN_GUI_RESULT=<prefix>` preloads a result for screenshots. Verified in headless Chrome at 1400 and 390 px against dist/pcropt.html; README screenshots retaken plus docs/img/gui_feedback.png. Not exercised: file uploads through the browser. Tests: test_feedback.py (8), test_redesign.py (4); 75 pass. Example regenerated (15 amplicons, 2 pools, 60 C) with `mtb_order_sheet.tsv`.

- 2026-11-12: `--one_pool` / `--one_pool_max` (CLI, GUI, input check). Targets that may not share a tube (overlapping tiles, outer product < --outer_min) are merged; only merged targets that need more room get larger size classes (same 150 bp spread, each ending at the longest need, 25 bp steps, cap 1000), the rest keep 225-375 (user preferred small amplicons plus a few big ones over shifting everything up). Per-target windows: Config.amp_by_target / cfg.amp(tid), used by candidates, balance, extension variants, optimiser; split fallback off in this mode; `size_class_bp` column and a 'pool layout' conditions row. MTB 24 sites: 9 x 225-375, Rv0678+pncA 475-625, gyrB+gyrA 775-925; 1 pool, 12 amplicons, all covered, Ta 60 C, no outer products, dG_eff spread 2.00 (two-pool design 1.22). 58 tests pass. Idea not done: suggested primer concentrations for long classes (user plans empirical rebalancing).

- 2026-11-12: VCF input and an input precheck. `--positions` (multiplex) and `getMultiPrimerSet --region_file` now take tables, BED or VCF (genome.read_vcf/is_vcf; one target per record covering REF; gVCF reference blocks skipped; END or |SVLEN| for symbolic DEL/DUP/INV/CNV). New primerJinn/multiplex/precheck.py runs before every design and alone with `--check` (both tools; GUI 'Check positions'): software, genome, every unreadable row, unknown names grouped with the likely cause (version mismatch, one-sequence FASTA), VCF REF vs genome bases and ##contig lengths, BED-style 2-bp tables, N gaps, sequence ends, grouping, variants, settings, output folder. getMultiPrimerSet: FASTA with uneven lines or gzip now read via a rewrapped/decompressed cache copy; regions keep the old overlap rule; a seqid column on a one-sequence FASTA maps to it (WARN). Two adversarial review rounds (4 lenses + 2 skeptics per finding, then a recheck) found 30 issues, all fixed with tests; 52 tests pass in base and the conda env. E2E: MTB design from a VCF, getMultiPrimerSet design from VCF regions on a gzipped uneven FASTA.

- 2026-11-12: environment.yml (conda-forge + bioconda: numba with matching llvmlite, ViennaRNA, pysam, primer3-py 2.0.3, BLAST+, streamlit, altair), then `pip install --no-deps -e .`. Prompted by the user's base env failing with a conda numba 0.64 over a pip llvmlite 0.44. Clean env tested: 41 pass (1 NUPACK reference test skipped), 3-site MTB design OK. The CLI now imports specificity at start-up and exits 2 with the fix in under a second if numba is broken (test in test_regressions.py). setup.py python_requires >= 3.10. Note: conda viennarna has no pip metadata, so plain `pip install .` inside the conda env would pull ViennaRNA from PyPI; hence --no-deps.

- 2026-11-12: README GUI screenshots in docs/img (24 MTB sites from example/mtb_positions.tsv: 15 amplicons, 2 pools, Ta 60 C), captured by driving headless Chrome over the DevTools protocol (Streamlit on a spare port). GUI: primer table scrolls inside its box; `use_container_width` replaced by `width="stretch"` (deprecated), so the gui extra needs streamlit >= 1.50.

- 2026-11-12: Streamlit GUI (primerJinn/multiplex/gui.py, entry point `primerjinn-gui`, extra `.[gui]`): check positions, design (runs the CLI, live log), results with metrics, advice, primer table with lowercase highlight, amplicon map and zoom, binding curves, dimers/products, downloads. Tested in the browser on 3 M. tuberculosis sites (3 amplicons, Ta 62 C). Helper tests in tests/test_gui.py; 34 tests pass. `.claude/launch.json` starts it on port 8517.

- 2026-11-12: Strain awareness: `--variants` (VCF/BED/TSV of positions) via primerJinn/multiplex/variants.py; primers with a known variant in the last 5 nt are rejected, others penalised, balance extensions stop at variants; report column and conditions row. Tests in tests/test_variants.py (31 tests pass). MTB check: masking all 30 primer 3' ends replaced every primer, coverage kept.

- 2026-11-12: Binding balance during selection: extended pair variants for weak primers are offered to the optimiser (design._add_extended_variants, balance.extension_options/best_extension/apply_extension). MTB: spread 2.28 -> 1.22 kcal/mol (middle half 0.39), pncA F now extended, one pass, 94 s. 28 tests pass.

- 2026-11-12: Full code review (4 Sonnet reviewers: newest code, docs vs code, numerical consistency, robustness; findings checked and fixed by hand). Fixed: crash when preset tails met the Ta re-score (Config.resolve not idempotent; design_pool now works on copies and iterates up to 3 passes until scoring and recommended Ta agree within 1 C); amplicons crossing N gaps; false 'mostly N' failures; N positions now count against coverage and exit code 1; Ta grid capped by --ta; one hairpin-cost function (tails included) for design, tables and balance; DMSO applied to hairpin/dimer terms; scale-aware design Ta for --tm_scale overrides; input validation; balance tail cap, stale GC and relaxed-dimer handling; timings; report wording (scoring Ta, specificity status, weakest primer, balance notes, amplicon flags, hifi advice); XLSX formula-safe cells; DESIGN.md formulas now match the code. Confirmed consistent: dG_eff identical across design/finalize/balance (0.003 kcal/mol), dangling-end orientation, output files vs genome in ~30 runs. 27 tests pass.

- 2026-11-11 (late night): Added binding balance by 5' extension (user's idea): weak primers keep their Tm-matched core and gain lowercase 5' bases until their dG_eff at Ta is near the pool median; refused when it creates a hard dimer, a long run or an over-long amplicon; specificity and final checks re-run. MTB: 2 primers extended, middle-50% dG_eff spread 0.59 kcal/mol, one primer (pncA F) blocked by a dimer. Report shows core and full Tm, lowercase sequences, and the reason for every primer not balanced. 23 tests pass. Workflow subagents unavailable until the weekly limit resets (Nov 17), so this was built and checked directly.

- 2026-11-11 (night): Read the supplied papers. Owczarzy 2008 salt model and both worked examples reproduced exactly; Ka 3e4 confirmed in the paper. Out-of-sample test on its SI Table S2 (32 duplexes <= 30 bp): salt part MAE 0.90 C, full prediction 2.4 C (PCR-like buffer 1.6 C, bias +1.5 C, nearest-neighbour error). Added template dangling ends (SantaLucia and Hicks 2004 Table 3, checked against primer3 thal), weakened Tm centring and added one re-score at the recommended Ta (SantaLucia 2007, myth 3), corrected AT-rich advice (Su 1996: 60 C extension, only for ~85-90% AT stretches) and added the Musso 2006 GC-rich escalation. SADDLE and Pythia anchors recorded in DESIGN.md. Not done: the Taq 3'-mismatch papers (Huang, Lefever, Stadhouders, Rejali) because the background agents hit the weekly usage limit; most skeptic checks of the paper findings also did not run, so the applied changes were verified by hand (code runs, paper text, primer3 cross-check). 22 tests pass; M. tuberculosis (15 amplicons, 2 pools, Ta 59 C) and E. coli (11 amplicons, Ta 55 C) re-run cleanly.

- 2026-11-11 (late): User set the default to standard Taq, probably KAPA Robust ReadyMix. Added the `kapa` preset from the data sheet (KR0382 v3.16), physical Tm scale for Taq-type enzymes, new structure-capped Ta rule. Independent validation (5 genomes, 4 code reviews, every finding verified by a skeptic) found: tails > 60 nt crash, tail panhandle masking sites, Q5 Ta at 500 nM instead of the reaction concentration, inconsistent pooling and product limits, no ARTIC-style tiling, AT-rich filter failures, 3-way split with redundant amplicon, false dimer-relaxation warning, double-counted costs, contig-blind designed-product check, slow per-target folding, read-only FASTA crash, thousands-separator split, N positions reported as covered. All fixed, regression tests added (21 tests), all genomes re-run.
- 2026-11-11: Implemented the designer and ran it on M. tuberculosis. Hairpin flags replaced by numbers, structure-aware Ta advice.
- 2026-11-11: Research run hit usage limits twice; restarted on Opus 5.5 and finished. NEB Q5 Tm reproduced exactly by SantaLucia 1998 + Owczarzy 2004 at 150 mM Na, C = primer concentration (240 primers, max error 0.000 C).
- 2026-11-10: Scouted repo and environment. Position parsing and target grouping written and tested.
