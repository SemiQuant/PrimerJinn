# PrimerJinn multiplex designer: progress and handover

Read this first when picking the work up in a new session. Keep it current: update the status table and the log whenever a phase finishes or a decision changes.

Last updated: 2026-11-11 (session 1)

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
| 9 | Default polymerase: standard Taq, probably KAPA Robust ReadyMix (user, 2026-11-11) | Preset `kapa` (default): 2 mM Mg, 0.8 mM dNTP, no proofreading, Tm on the physical scale, Ta from the model with KAPA's 55/60/65 C rule shown alongside, KAPA cycling and 5% DMSO advice |
| 10 | Keep this document | This file |

## Status

| Phase | State | Notes |
|---|---|---|
| 1. Scout existing code and environment | Done | |
| 2. Research, fact-checked | Done | `docs/research/` (7 reports, each with an independent `*.verified.md`; `fact_check_digest.txt`) |
| 3. Design | Done | `docs/DESIGN.md` |
| 4. Implementation | Done | `primerJinn/multiplex/`, 21 tests pass (`python3 -m pytest tests`) |
| 5. Tests on real genomes | Done | M. tuberculosis, E. coli, S. aureus (with Illumina tails), P. falciparum (incl. AT-rich intergenic locus), SARS-CoV-2 spike tiling. All exit 0, all positions covered |
| 6. Independent review and fixes | Done, one round | 4 genome testers + 4 code reviewers, each finding checked by a separate skeptic; all confirmed findings fixed (see log). The fixes themselves were re-tested on all genomes but not re-reviewed by independent agents |
| 7. Wet-lab calibration | Not started | Needs read depth from a first sequenced pool |

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
| `report.py` | TSV, XLSX, FASTA, BED, markdown |

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
2. Circular replicons (products across the origin), SNP masking under primers, primer rebalancing worksheet from observed coverage.
3. Optional second independent review of the fixes made after validation.
4. Fix the `PrimerJinn/` vs `primerJinn/` case mismatch in git before publishing.

## Repo notes

- git tracks the package as `PrimerJinn/` but on disk (case-insensitive macOS) it is `primerJinn/`; `setup.py` uses `primerJinn`. A fresh clone on Linux would not import. Fix with `git mv` when convenient.
- Nothing has been committed yet in this work.
- Validation scripts and outputs from session 1 were in the session scratchpad (temporary); the evidence that matters is summarised here and in `docs/research/`.

## Log

- 2026-11-11 (late night): Added binding balance by 5' extension (user's idea): weak primers keep their Tm-matched core and gain lowercase 5' bases until their dG_eff at Ta is near the pool median; refused when it creates a hard dimer, a long run or an over-long amplicon; specificity and final checks re-run. MTB: 2 primers extended, middle-50% dG_eff spread 0.59 kcal/mol, one primer (pncA F) blocked by a dimer. Report shows core and full Tm, lowercase sequences, and the reason for every primer not balanced. 23 tests pass. Workflow subagents unavailable until the weekly limit resets (Nov 17), so this was built and checked directly.

- 2026-11-11 (night): Read the supplied papers. Owczarzy 2008 salt model and both worked examples reproduced exactly; Ka 3e4 confirmed in the paper. Out-of-sample test on its SI Table S2 (32 duplexes <= 30 bp): salt part MAE 0.90 C, full prediction 2.4 C (PCR-like buffer 1.6 C, bias +1.5 C, nearest-neighbour error). Added template dangling ends (SantaLucia and Hicks 2004 Table 3, checked against primer3 thal), weakened Tm centring and added one re-score at the recommended Ta (SantaLucia 2007, myth 3), corrected AT-rich advice (Su 1996: 60 C extension, only for ~85-90% AT stretches) and added the Musso 2006 GC-rich escalation. SADDLE and Pythia anchors recorded in DESIGN.md. Not done: the Taq 3'-mismatch papers (Huang, Lefever, Stadhouders, Rejali) because the background agents hit the weekly usage limit; most skeptic checks of the paper findings also did not run, so the applied changes were verified by hand (code runs, paper text, primer3 cross-check). 22 tests pass; M. tuberculosis (15 amplicons, 2 pools, Ta 59 C) and E. coli (11 amplicons, Ta 55 C) re-run cleanly.

- 2026-11-11 (late): User set the default to standard Taq, probably KAPA Robust ReadyMix. Added the `kapa` preset from the data sheet (KR0382 v3.16), physical Tm scale for Taq-type enzymes, new structure-capped Ta rule. Independent validation (5 genomes, 4 code reviews, every finding verified by a skeptic) found: tails > 60 nt crash, tail panhandle masking sites, Q5 Ta at 500 nM instead of the reaction concentration, inconsistent pooling and product limits, no ARTIC-style tiling, AT-rich filter failures, 3-way split with redundant amplicon, false dimer-relaxation warning, double-counted costs, contig-blind designed-product check, slow per-target folding, read-only FASTA crash, thousands-separator split, N positions reported as covered. All fixed, regression tests added (21 tests), all genomes re-run.
- 2026-11-11: Implemented the designer and ran it on M. tuberculosis. Hairpin flags replaced by numbers, structure-aware Ta advice.
- 2026-11-11: Research run hit usage limits twice; restarted on Opus 5.5 and finished. NEB Q5 Tm reproduced exactly by SantaLucia 1998 + Owczarzy 2004 at 150 mM Na, C = primer concentration (240 primers, max error 0.000 C).
- 2026-11-10: Scouted repo and environment. Position parsing and target grouping written and tested.
