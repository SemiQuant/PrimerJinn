# State of the art in multiplex PCR primer pool design (2019 to 2026)

Status: complete (third attempt; sections 1 to 4.4 from earlier attempts were re-checked against cloned source by grep on 2026-11-11: Olivar weights and SA schedule, PrimalScheme3 dimer_score -26.0, varVAMP thresholds, NGS-PrimerPlex -6/-10 defaults; PRISM metadata confirmed via PubMed PMID 42378440). Sources are local copies under `research/papers`, `research/mx`, `research/src/*/repo` (cloned tool source).

## 1. Tool-by-tool extraction

Notation: "paper" = read from the PMC full text saved locally; "code" = read in the cloned repository at the commit given.

### 1.1 SADDLE (Xie et al. 2022, Nat Commun 13:1881, doi:10.1038/s41467-022-29500-4, PMC9001684)

- Candidates (paper): from "pivot" nucleotides that must be in the insert, proto-primers are placed with the 3' end just outside the pivot, then trimmed from the 3' end until primer-to-template dG° is between -10.5 and -12.5 kcal/mol (target about -11.5 kcal/mol, "best tradeoff between amplification efficiency/uniformity and nonspecific hybridization"). GC 0.25 to 0.75. All fP x rP combinations within the amplicon length window are enumerated as pair candidates.
- Dimer score (paper, Eq. 2): `Badness(pa,pb) = sum over reverse-complementary subsequences of length >= 4 of 2^len * 2^numGC / ((d1+1)(d2+1))`, d1, d2 = distance of the subsequence to each 3' end; len capped at 8. Rationale stated: up to 3 nt 3' complementarity gave no significant NTC dimers in qPCR; 10 nt at the 5' end gave none, 5 nt at the 3' end did. Authors deliberately avoided nearest-neighbour parameters, citing uncertain effective salt and poor extrapolation of dH/dS to high temperature. Hash-table evaluation is O(N x P).
- Loss: `L(S) = sum_{b>=a} Badness(pa,pb)` (includes self terms, precomputed).
- Optimiser: simulated annealing over one pair per target. Change 1 or more pairs; accept a worse set with `p = exp((L(Sg) - L(T)) / C(g))` for g < gt, then 0 (pure descent). Multiple restarts recommended. They stopped at 40,000 generations rather than 60,000 to avoid "over-training" an imperfect loss.
- Pools: single tube.
- Validation (paper): 96-plex (192 primers) dimer fraction of reads 90.7% (random set PS1) -> 39.6% (PS2) -> 4.9% (PS3); loss fell about 24-fold. 384-plex (768 primers): about 1% dimer reads. 96-plex design about 490 s for 60,000 iterations on a laptop (MATLAB). Vent (exo-) polymerase, 17 cycles, median 45 nM primer.

### 1.2 Olivar (Wang et al. 2024, Nat Commun 15:6306, doi:10.1038/s41467-024-49957-9, PMC11282221; code treangenlab/Olivar b632a2f)

- Candidates: per-nucleotide risk array = weighted sum of 4 components in [0,1] (SNP frequency, extreme GC, low complexity, non-specificity by BLAST), default weights 1.0 each (code `tiling`). Primer design regions (PDR, 40 nt) chosen to minimise risk; SADDLE-style generator (code `design.primer_generator`, defaults: 60 °C, Na 0.18 M, dG_max -11.8 kcal/mol, GC 0.2 to 0.75, max length 36).
- Dimer score (code, `PrimerSetBadnessFast`), which differs from the paper formula: k-mers of the probing primer's reverse complement are looked up in hash tables of every primer's last 4, 5, 6 nt (3'-anchored) and of all 7-mers and 8-mers (weighted by 1/(distance from 3' end + 1)). Weights END4=1, END5=4, END6=20, MIDDLE7=100, MIDDLE8=500, times 2^numGC, times 1/(j+1) where j is the position in the probing reverse complement, times primer concentration.
- Optimiser (code `tiling_helper.optimize`): per pool, random start; each step changes one target's pair; accept worse with probability `2^((L_old - L_new)/T)`. T0 = 1000 + 10*max(N/2 - 100, 0); linear cooling over NUMSTEPS = 10 + int(pool_size/10) steps of 1000 trials each, then the same number of zero-temperature steps.
- Pools: tiled amplicons alternate pool 1 and 2 (`'tube': i%2 + 1`), so adjacent overlapping amplicons are never in one tube.
- Validation (paper): SARS-CoV-2, 146 amplicons, 292 primers, designed in <30 min; SADDLE loss pool 1 49,692 -> 6,640, pool 2 95,747 -> 7,192 vs PrimalScheme 69,474 and 52,369; mapping rate about 90% (similar to ARTIC v4.1); fewer SNP-overlapping primers (4 vs 18) and predicted non-specific products (5 vs 27) than PrimalScheme; up to 3-fold higher mapping in wastewater.

### 1.3 PrimalScheme (Quick et al. 2017 Nat Protoc 12:1261, doi:10.1038/nprot.2017.066, PMC5902022) and PrimalScheme3 (code ChrisgKent/primalscheme3 f71ec1c, dimer kernel primalschemers 414a2b3)

- PrimalScheme v1: greedy left-to-right tiling, 2 pools alternating, Primer3 candidates (code primalscheme1 04837c1).
- PrimalScheme3 defaults (code `core/config.py`): 2 pools, amplicon 400 bp (+/-10%), primer 19 to 36 nt (high-GC mode 17 to 30), GC 30 to 55% (high-GC mode 30 to 67%), Tm 59.5 to 62.5 °C, hairpin Tm max 55 °C, homopolymer max 5, Na 100 mM, Mg 2 mM, dNTP 0.8 mM, oligo 15 nM, `dimer_score = -26.0`. Optional fraction-annealed filter at 65 °C (`calc_annealing` via primer3 heterodimer with the template, K = exp(-dG/RT)), target 30% +/- 5.
- Dimer score ("primaldimer", Rust): for every offset where at least one of the two 3'-terminal bases of seq1 pairs with seq2, score = SantaLucia NN dG sum over the aligned window + dangling-end terms + a 3'-extension term (-(3 for G/C, 2 for A/T)/(i+1) over the last 4 bases, -2 more if all 4 match) + fitted bonuses (longest match coefficient 2.33, match proportion coefficient 3.25, bubble penalty). Coefficients were fitted (BONUS_ARRAY). A pair "interacts" if any offset scores <= -26 in either direction. It is a binary filter, not a continuous loss.
- Panel mode (code `panel/panel_main.py`, `panel_classes.py`): greedy. Primer pairs are scored by the number of still-uncovered region bases inside the primer-trimmed amplicon x (1 - GC deviation), x10 if the pair bridges a coverage break; best pair is tried in each pool in round-robin and added to the first pool where it does not overlap another amplicon (`does_overlap`), does not interact (primaldimer <= -26) with any primer already there, and creates no new off-target product (k-mer match DB). No backtracking by default.

### 1.4 varVAMP (Fuchs et al. 2025, Nat Commun 16:5067, doi:10.1038/s41467-025-60175-9, PMC12126543; code jonas-fuchs/varVAMP 85d870b)

- Candidates: k-mers from a consensus with ambiguity codes; Primer3 Tm 56 to 63 °C (opt 60), GC 35 to 65%, length 18 to 24, hairpin Tm <= 47 °C, 1 to 3 G/C in the last 5 bases, penalty includes 3' mismatch weights (32, 16, 8, 4, 2) and a degeneracy penalty (code `default_config.py`).
- Dimer test (code `primers.is_dimer`): `primer3.calc_heterodimer(mv=100 mM, dv=2 mM, dNTP=0.8 mM, dna=15 nM)`; dimer if Tm > 35 °C or dG < -9 kcal/mol, or if the ASCII structure shows a perfect end overlap of more than 5 nt (extensible end-to-end overlap even when Tm and dG pass).
- Tiling: amplicons are nodes, overlaps are edges, Dijkstra shortest path over penalty costs; pools alternate. Dimers inside a pool are resolved afterwards by swapping in overlapping alternative primers (`check_and_solve_heterodimers`), not by a global optimiser.

### 1.5 MFEprimer-3.0 (Wang et al. 2019, NAR 47:W610, doi:10.1093/nar/gkz351, PMC6602485)

- QC tool, not a pool optimiser. k-mer binding-site search allowing mismatches except at the 3'-terminal base; dimer detection "modified from AutoDimer": slide two sequences one base at a time, keep alignments with score >= cut-off, report the best, and report dG. Exact default cut-offs for `mfeprimer dimer` could not be verified (the web docs list only global defaults: Na 50 mM, Mg 1.5 mM, dNTP 0.25 mM, oligo 50 nM, k = 9). The binary is closed source in the cloned repo (README only).

### 1.6 NGS-PrimerPlex (Kechin et al. 2020, PLoS Comput Biol 16:e1008468, doi:10.1371/journal.pcbi.1008468, PMC7802936; code aakechin/NGS-PrimerPlex e59627f)

- Candidates: Primer3 per target; SNP overlap check against dbSNP (frequency >= 0.05); non-target amplicons by alignment.
- Dimer test (code): any dimer `calcHeterodimer/calcHomodimer` dG must be >= -10 kcal/mol (`-minmultdimerdg2`), and a 3'-anchored dimer dG >= -6 kcal/mol (`-minmultdimerdg1`). The 3'-anchored dG (`calcThreeStrikeEndDimer`) finds the longest 3'-terminal run (>= 4 nt) of one primer whose reverse complement occurs in the other, then runs `calcHeterodimer` on that segment plus one flanking base. Defaults: Na 50 mM, Mg 3 mM, dNTP 0.8 mM, primer 250 nM. Neighbouring amplicons that would produce a short product between inner primers are tested with adapter-tailed sequences when adapters are given.
- Pooling (code `makeFinalMultiplexes`): graph where nodes are pairs and edges mean "compatible" (no dimer, no product overlap, no non-target product); pools are found as cliques (`networkx.find_cliques`, first maximal clique, trimmed to ceil(n/k) by degree), recursively for k pools. Contiguous regions are joined by shortest path.
- Validation: LRRK2 panel on 354 samples, median 97.4% of target bases at >= 30 reads (paper).

### 1.7 PMPrimer (Yang et al. 2023, Sci Rep 13:16313, doi:10.1038/s41598-023-43825-0, PMC10542359; code AGIScuipeng/PMPrimer a7ffccb)

- Designs conserved (degenerate) primers across diverse templates: Shannon entropy conserved regions, haplotype-based gap tolerance, Primer3. Multiplex evaluation is coverage and taxon specificity, maximum Tm difference 12 °C by default. The code has only a homodimer check (`primer3.calc_homodimer`); no cross-dimer pool optimisation was found. Limited relevance.

### 1.8 PrimerPooler (Brown et al. 2017, Biol Methods Protoc 2:bpx006, doi:10.1093/biomethods/bpx006, PMC6994079; code ssb22/PrimerPooler ee2b2fb)

- Not a designer: splits existing pairs into k pools. Interaction = worst of all alignments by match score or NN dG (Allawi and SantaLucia parameters from MPprimer), default 60 °C, Na 50 mM in the paper example.
- Objective (code `pool-split.c`): lexicographic minimax. Each pool's badness packs (worst score, count at worst, count at worst-1, count at worst-2) into one 64-bit integer, with score = -dG x 2. Moves swap a pair between pools; hill climbing from many random starts in parallel; overlapping amplicons are forced into different pools. Initial pool count suggested by single-pass packing at dG -7 kcal/mol.
- Validation: 1153 pairs into 3 pools (388, 389, 376), then 144 sub-pools of 6 to 9-plex; mean 3269 reads per target, 95% of targets >= 50 reads; worst residual interaction -6 kcal/mol at 60 °C, 50 mM Na.

### 1.9 Ultiplex (Yuan et al. 2021, BMC Genomics 22:839, doi:10.1186/s12864-021-08149-1, PMC8600765)

- Primer3 candidates; primers are checked together with their 5' tags: hairpin Tm > 45 °C removed, heterodimer (tag + primer, `primer3.calcHeterodimer`) Tm > 40 °C removed; last 7 nt of the 3' end must not hit SNPs or repeats; BLAST specificity; Tm difference between compatible pairs < 5 °C. Pairs are grouped into maximal mutually compatible groups.
- Validation: 295 targets submitted, 271 clustered into one compatible group; the experimental 108-pair panel gave mean 5000x, 99.21% alignment, 7 targets < 30x, and doubling primer concentration rescued weak targets (paper; the abstract says "108 primers", the results say 108 primer pairs).

### 1.10 openPrimeR (Kreer et al. 2020, J Immunol Methods 480:112752, doi:10.1016/j.jim.2020.112752, PMID 31991148; code matdoering/openPrimeR bbd319d)

Citation corrected in this pass: an earlier draft cited "BMC Bioinformatics 21:69, doi:10.1186/s12859-020-3311-y"; that DOI does not resolve in Crossref, and PubMed esearch `Kreer[au] AND openPrimeR` returns the J Immunol Methods paper.

- Aim is minimal primer sets covering many templates (immunoglobulin repertoires). Optimiser: greedy or ILP (lpSolveAPI) set cover; ILP maximises covered templates for a fixed set size k. Default design constraints (`A_Taq_PCR_design.xml`): cross-dimer dG >= -5 kcal/mol (relaxable to -7), self-dimer >= -5, Tm difference <= 5 °C, GC 0.4 to 0.6, length 18 to 22, GC clamp 1 to 3. Relevant mainly as the ILP precedent. (DOI from memory of the citation, not re-checked online.)

### 1.11 primerJinn (Limberis and Metcalfe 2023, BMC Bioinformatics 24:468, doi:10.1186/s12859-023-05609-1, PMC10714587)

- Candidates: Primer3 per region (product 400 to 800 bp default in the current code), Q5-like Tm.
- Selection: features per pair = (mean Tm, product size, Tm difference, within-pair heterodimer dG); agglomerative clustering (Ward, Euclidean); choose the cluster covering most regions and fill missing regions from the nearest cluster. Cross-pair dimers between different targets are not part of the selection features (`get_features` only scores forward vs reverse of the same pair). Current code calls `calc_heterodimer(temp_c=target_tm, mv=70 mM, dv=2 mM, dna=3300 nM)`.
- Validation: 8 M. tuberculosis drug-resistance regions, mean Tm 65 °C (64 to 67), mean amplicon 665 bp (454 to 791), Q5 Hot Start, 35 cycles, anneal 65 °C, 82% reads aligned to targets. Paper limitation: cannot fix overlapping input regions.

### 1.12 Other notable tools 2019 to 2026

- multiply (Hendry et al. 2024, Nat Commun 15:1394, doi:10.1038/s41467-024-45688-z, PMC10869361; code JasonAHendry/multiply 8505b40): target-panel design for P. falciparum and also ships a `tb-amr.ini`. Dimer score: ungapped, 5'-overhang (extensible) alignments with NN dG, double mismatch +0.2, plus an end bonus of -0.5 per complementary base in the last 4 nt of either 3' end ("like Johnston et al. 2019"). Cost = sum of individual costs (Primer3 penalty, BLAST mispriming count, SNP counts, weight 1 each) + sum of pairwise dimer scores (weight -1). Optimiser: randomised-order greedy repeated N = 10,000 times (also brute force, random).
- PRISM (Chen, Han, ... Boucher, Bioinformatics 2026 42(7), doi:10.1093/bioinformatics/btag478, PMID 42378440, PMC13378452, code github.com/yhhan19/PRISM-new): constrained submodular maximisation balancing coverage and SADDLE-style Badness with a local search having a constant-factor approximation guarantee; reports lower Badness than PrimalScheme, Olivar and primerJinn on viral genomes (abstract only available locally, details not verified).
- MultiPrime (2024, iMeta 3:e143, doi:10.1002/imt2.143, PMC10989836): degenerate/mismatch-tolerant tNGS primer design with a dimer examination step; reduced dimer percentage. Algorithmic dimer details not extracted.
- PrimerROC (Johnston et al. 2019, Sci Rep 9:209, doi:10.1038/s41598-018-36612-9, PMC6338771): the key empirical dimer study (below).
- PRIMEval (Conzemius et al. 2019, Sci Rep 9:19213, doi:10.1038/s41598-019-55883-4, PMC6917790): in silico evaluation of multiplex primer/probe sets with NN thermodynamics; QC only.
- ThermoPlex (2025, Biol Methods Protoc, doi:10.1093/biomethods/bpaf074, PMC12598146): ensemble NN thermodynamics and multi-reaction equilibrium for specificity; explicitly does not yet penalise primer-primer dimers (tested 4 to 5-plex).
- Gimpel et al. 2025 (Nat Commun, doi:10.1038/s41467-025-64221-4, PMC12533003): 1D-CNN predicts sequence-specific amplification efficiency in multi-template PCR (AUROC 0.88); identifies adapter-mediated self-priming as the main cause of poor amplification. Relevant to adapter-tailed primers.

PrimalScheme 2017 (paper, verified): default amplicon 400 nt (suggested 200 to 2000), default overlap 75 nt, Primer3 candidates, "alternate target genome regions" go to one of two pools "so that neighbouring amplicons do not overlap within the same pool (which would result in a short overlap product being generated preferentially)". Lower primer concentration and longer annealing/extension were found optimal.

## 2. Comparison table

| Tool | Candidate generation | Dimer model | Pool optimiser | Overlap handling | Validation |
|---|---|---|---|---|---|
| SADDLE 2022 | 3'-trimmed to template dG -10.5 to -12.5 kcal/mol | Badness = sum 2^len 2^GC / ((d1+1)(d2+1)), len 4 to 8 | Simulated annealing, Metropolis, then descent | single tube | 96-plex dimer reads 90.7% -> 4.9%; 384-plex about 1% |
| Olivar 2024 | risk array (SNP, GC, complexity, BLAST), SADDLE generator dG <= -11.8 at 60 °C | Olivar Badness: 3'-anchored 4/5/6-mers weight 1/4/20, 7/8-mers 100/500, x 2^GC / (j+1) | SA, accept 2^(dL/T), linear cooling, then greedy | 2 pools alternating | SARS-CoV-2 146 amplicons, about 90% mapping |
| PrimalScheme3 | k-mer digestion of MSA, Tm 59.5 to 62.5, GC 30 to 55 (67 high-GC) | primaldimer: NN dG + 3' extension + fitted bonuses, interacts if <= -26 | greedy, first feasible pool | overlap and dimer and mispriming checks per pool | ARTIC lineage (not re-verified here) |
| varVAMP 2025 | degenerate k-mers, Primer3 penalty | primer3 heterodimer Tm > 35 °C or dG < -9 kcal/mol (37 °C), or > 5 nt perfect end overlap | Dijkstra tiling, post hoc dimer swaps | 2 pools alternating | viral panels (paper) |
| MFEprimer-3 | QC only | AutoDimer-like score + NN dG | none | none | n/a |
| NGS-PrimerPlex 2020 | Primer3, dbSNP, alignment | primer3 dG any >= -10, 3'-anchored dG >= -6 kcal/mol (37 °C, 50 mM Na, 3 mM Mg) | clique search on compatibility graph | graph edges exclude overlaps | LRRK2, 354 samples, 97.4% bases >= 30x |
| PrimerPooler 2017 | none | worst alignment, NN dG | hill climbing, lexicographic minimax, random restarts | overlaps forced apart | 1153 pairs, 95% targets >= 50 reads |
| Ultiplex 2021 | Primer3 + BLAST | primer3 heterodimer of tag+primer Tm > 40 °C rejected; hairpin Tm > 45 °C | maximal compatible groups | via compatibility | 108 pairs, mean 5000x, 7 targets < 30x |
| openPrimeR 2020 | tree/naive | dG >= -5 kcal/mol (relax -7) | ILP set cover or greedy | n/a | repertoire primers |
| primerJinn 2023 | Primer3 | primer3 heterodimer dG within pair | Ward clustering | none (stated limitation) | 8-plex Mtb, 82% on target |
| multiply 2024 | Primer3 (several settings) | PrimerROC-like ungapped 5'-overhang NN score with 3'-end bonus | randomised greedy x 10,000 | n/a | P. falciparum 8 and 16-plex |
| PRISM 2026 | n/a | Badness | submodular local search with guarantee | tiling | in silico only (abstract) |

## 3. What the empirical literature says about dimer thresholds

- SADDLE qPCR pilot (paper): <= 3 nt of 3' complementarity gave no significant dimer even in NTC; 5 nt at the 3' end did; 10 nt at the 5' end did not.
- PrimerROC (Johnston 2019, paper): over 300 primer pairs, plain dG (any structure) was a weak classifier; best plain-dG AUC across temperature/salt was 0.85 (37 °C, 0.5 M Na), versus 0.94 for their 3'-anchored, extensible-only score. Stable structure at a single 3' end was sufficient to give abundant dimer artefacts. Their length/GC-adjusted score had a dimer-free cut-off of -4.43 (their score units, close to kcal/mol, at their conditions). Non-extensible (5'-5') structures were ignored by design.
- Thresholds used in code (verified): NGS-PrimerPlex 3'-anchored dG >= -6 and any dG >= -10 kcal/mol; varVAMP any dG >= -9 kcal/mol or dimer Tm <= 35 °C; PrimerPooler packing -7, best residual -6 at 60 °C; openPrimeR -5 (relax to -7); Primer3 default PRIMER_PAIR_MAX_COMPL_END_TH = 47 °C (structure Tm, not dG; p3manual). Note: primer3-py `calc_heterodimer` and `calc_end_stability` default to temp_c = 37 °C, so every dG threshold above is a 37 °C value; varVAMP and NGS-PrimerPlex call them without temp_c.
- Adapter tails: NGS-PrimerPlex and Ultiplex score dimers with the tails attached. Gimpel et al. 2025 identified adapter-mediated self-priming as the main cause of poor amplification in multi-template PCR.

## 4. Experiments run here (all code and outputs in `research/multiplex/`)

Environment: primer3-py 2.0.3, scipy 1.15.2 (HiGHS MILP), Python 3.13.2, macOS arm64, 16 cores. Reference H37Rv `example/ref.fasta` from the repo (read only), positions from `example/mtb_positions.tsv` (24 targets).

### 4.1 What `calc_end_stability` measures (exp1_semantics.py, exp1_out.txt)

`calc_end_stability(x, y)` returns the best structure in which the 3'-terminal base of x is paired to y (both given 5'->3'). If y carries rc(x's last 6 nt) at its 5' end (x can extend), end(x,y) = -2.51 and end(y,x) = +5.97 kcal/mol at 65 °C; a 5'-end-only structure gives +1.52 while `calc_heterodimer` still reports -0.82. So the 3'-anchored score is `d3(a,b) = min(end(a,b), end(b,a))`, and it ignores non-extensible structures, which matches the PrimerROC finding. Cost: about 71 to 78 us per call single-threaded (40,000 calls in 2.85 to 3.11 s).

### 4.2 dG of a clean 3'-3' overlap of k nt (exp1b_scan.py, exp1b_out.txt; 250 nM, 0.8 mM dNTP)

| k (GC-rich 3' end ...GGCGC) | 37 °C, 50 mM Na, 2 mM Mg | 65 °C same salt |
|---|---|---|
| 4 | -4.90 | -1.93 |
| 5 | -6.33 | -2.75 |
| 6 | -8.75 | -4.33 |
| 7 | -10.54 | -5.51 |

For an AT/GC mixed 3' end (...ACTAGTACCG) k = 5, 7, 9 gave -3.72, -6.77, -8.62 at 37 °C. At 65 °C even k = 12 gives only -4.83, and the dimer "Tm" for k = 5 GC-rich is 12.8 °C, yet SADDLE observed dimers with 5 nt. Dimers form during set-up and early cycles and need only transient binding for extension, so equilibrium occupancy at Ta is not the right metric for dimers. Use 37 °C dG (the scale every published threshold uses) for dimer scoring, and reserve Ta-based fraction bound for primer-template binding. A -6 kcal/mol 3'-anchored cut-off (NGS-PrimerPlex) at 2 mM Mg corresponds to about 5 nt of GC-rich or 7 nt of mixed 3' complementarity, consistent with SADDLE's 5-nt observation.

### 4.3 Candidate generation and interaction statistics for the 24 Mtb targets (exp2_candidates.py, exp2_out.txt)

Positions merged when the span is <= 250 bp; ranges longer than one amplicon can cover (315 bp, = 375 - 2 x (25 + 5)) were tiled. 24 targets -> 16 windows (rpoB 5 codons -> 1 window of 183 bp; Rv0678 and pncA -> 2 tiles each). Merging embB M306 and G406 (span 302 bp) gave zero Primer3 pairs (all rejected on product size or Tm difference), so the merge limit must leave room: 250 bp worked. Primer3 with Tm 62 to 68 (opt 65), 18 to 25 nt, GC 40 to 75%, product 225 to 375 (opt 300), 100 pairs per window took 2.9 s and produced only 737 unique primers (pairs share primers heavily). All 271,953 primer-primer combinations x 3 primer3 calls took 3.4 s on 16 processes.

Fraction of primer pairs at or below a 3'-anchored dG37: -4: 12.3%, -5: 5.2%, -6: 2.4%, -7: 1.0%, -8: 0.44%, -9: 0.19%. Any-structure heterodimer dG37: -6: 16.5%, -8: 3.1%, -10: 0.52%, -12: 0.10%.

### 4.4 Optimisers on the same objective (exp3_optimise.py)

Objective (section 5) with hard thresholds d3 <= -7 or dh <= -12 (+1000 each). Pools from DSATUR colouring of a proximity graph (edge if windows are < 1,500 bp apart): 2 pools (gyrB vs gyrA, the Rv0678, pncA and embB tiles split).

K = 30 pairs per window (exp3_out_K30.txt):

| Method | Energy | Hard violations | Worst d3 (kcal/mol) | Tm range °C | Length range | Time |
|---|---|---|---|---|---|---|
| 1 random draw | 1134.7 | 1 | -8.63 | 63.4 to 66.0 | 267 to 368 | 0 s |
| greedy, 200 random orders | 46.7 | 0 | -4.80 | 64.5 to 65.7 | 295 to 369 | 0.01 s |
| SA 8 x 20,000 + ICM descent | 46.5 | 0 | -4.48 | 64.5 to 65.7 | 294 to 369 | 0.44 s |
| exact MILP (HiGHS, 480 x + 25,539 y vars) | 46.5 (optimal) | 0 | -4.48 | same | same | 157.7 s |

K = 100 (exp3_out_K100.txt): greedy 46.5, SA 46.3 (0.49 s); MILP not run. SA reached the proven optimum at K = 30 in all 8 restarts. The longest amplicons (369 bp) come from the 296-bp tiles, so tile spans should be shorter to keep lengths near 300.

Note on the MILP: it is the standard linearisation of a quadratic assignment (y_ij >= x_i + x_j - 1 for every positive pairwise cost). It proved optimality but took 158 s for 16 windows x 30 candidates, against 0.44 s for SA. It is useful as a test oracle, not as the production optimiser.

Hard cap -6 instead of -7 (exp6_hard6.py, exp6_out.txt, K = 100): a random draw has 10 hard violations (E 10,179.7), greedy and SA both reach 0 violations with worst d3 -4.48 to -4.45 kcal/mol and E 46.3 (SA, 0.74 s). So the stricter NGS-PrimerPlex value of -6 costs nothing on this 24-target Mtb panel.

### 4.5 Scale test, one pool, random single positions on H37Rv (exp4_scale.py)

Same objective, hard 3'-anchored threshold -7 kcal/mol, K = 100 Primer3 pairs per target.

| N targets | Primers | Interactions scored (3 primer3 calls each, 16 processes) | Greedy x50 | SA 4 restarts + ICM |
|---|---|---|---|---|
| 50 | 2,828 | 4,000,206 in 48.1 s | E 298.7, worst d3 -6.77, 1.0 s | E 214.8, worst d3 -6.16, 2.3 s |
| 100 | 5,555 | 15,431,790 in 182.2 s | E 2292.1, 1 hard violation, worst d3 -7.75, 4.3 s | E 816.6, 0 violations, worst d3 -6.98, 6.0 s |

SA beat greedy by 28% at N = 50 and by 64% at N = 100 (where greedy left a hard violation), and the restarts disagreed (N = 50: 226, 221, 226, 215; N = 100: 868, 817, 995, 850), so at this size the landscape is rugged and more restarts or an iterated local search pay off. At 100-plex in one tube the best design only just met the -7 cap (worst -6.98), so the -6 cap is probably not reachable at 100-plex with 100 candidates per target: expect to relax to -7, raise K, or split into two pools at that scale. Pair scoring time grows with P^2 (48 s at 2,828 primers, 182 s at 5,555). Scoring all primer pairs is the dominant cost (O(P^2), P = unique primers); it can be cut by scoring only primers that appear in the top K pairs and by a cheap k-mer prefilter (only call primer3 when the two primers share a reverse-complementary 4-mer at a 3' end), see 5.2.

### 4.6 Cheap scores versus primer3 3'-anchored dG, and adapter tails (exp5_scores.py, exp5_out.txt)

On 2,996 random primer pairs from the Mtb candidates:

| Score | Spearman with -d3 (37 °C) | Recall of pairs with d3 <= -6 when flagging the same number of top-scoring pairs (n = 83) |
|---|---|---|
| SADDLE paper Badness (Eq. 2) | 0.55 | 0.39 |
| Olivar code Badness | 0.26 | 0.16 |
| primer3 any-structure heterodimer dG | 0.44 | 0.39 |

So the k-mer Badness scores are not good surrogates for a 3'-anchored nearest-neighbour dG; they are fast proxies for ranking huge sets. With 50 to 100 targets the primer3 calls are affordable (section 4.5), so there is no need for the proxy except as a prefilter.

Illumina overhang tails (TF = TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG on forward, TR = GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG on reverse): attaching them changed the 3'-anchored dG37 by a median of -1.07 kcal/mol (5th percentile -4.92), and 51% of pairs became more stable by more than 1 kcal/mol. Any-structure dG changed by median -1.29. Tail vs tail: TF/TR heterodimer -4.12, TF self -3.06, TR self -5.87 kcal/mol; tail hairpin Tm 34.8 and 42.1 °C. Conclusion: score dimers on the full tailed oligos when the user will order tailed primers; the tail-tail terms are a constant offset common to every pair and should be excluded from the comparison (they cannot be designed away).

## 5. Recommended algorithm for the new primerJinn designer

### 5.1 What to adopt, and from whom

| Decision | Adopt | Source of the idea |
|---|---|---|
| Candidate generation | Primer3 `design_primers` per window, 50 to 200 pairs, product 225 to 375 bp (opt 300), Tm window centred on the genome-derived target | primerJinn, NGS-PrimerPlex, PrimalScheme |
| Primer-template affinity | Uniform fraction bound at Ta (or dG to template) in the unary cost, not only Tm | SADDLE (dG -10.5 to -12.5), Olivar, PrimalScheme3 (`calc_annealing`, 30% at 65 °C) |
| Dimer model | primer3 thermodynamic alignment: 3'-anchored dG (`calc_end_stability` both ways) as the main term, any-structure `calc_heterodimer` as secondary | PrimerROC (3' anchoring, extensible only), NGS-PrimerPlex (-6 and -10 kcal/mol), varVAMP (end-overlap rule) |
| Objective | Sum of unary costs + sum of pairwise dimer and off-target costs, soft quadratic with hard caps | SADDLE (sum over pairs), PrimerPooler (worst case matters) |
| Optimiser | Repeated randomised greedy start, then simulated annealing with single-target moves and 1-opt descent, several restarts; optional MILP oracle for tests | SADDLE, Olivar, multiply |
| Pools | Conflict graph over windows, DSATUR colouring; single pool whenever the graph has no edges | PrimalScheme, PrimerPooler, NGS-PrimerPlex |

### 5.2 Dimer scoring (exact calls)

Conditions: the user's real reaction, defaults for a Q5 or KAPA-type mix: `mv_conc` 50 mM (monovalent, K+ and Na+ together), `dv_conc` 2 mM Mg2+, `dntp_conc` 0.8 mM (0.2 mM each), `dna_conc` = per-primer concentration in nM (for example 50 to 250). Score at `temp_c=37` because every published dG threshold (NGS-PrimerPlex, varVAMP, openPrimeR, PrimerROC) is a primer3 or NN value at 37 °C, and section 4.2 shows that at 60 to 65 °C even 12-nt 3' overlaps look harmless (dG > -5) although 5-nt 3' overlaps make dimers in practice (SADDLE). Dimers are made during set-up, ramping and the first cycles, so a low-temperature scale is the right one.

```python
C = dict(mv_conc=50, dv_conc=2.0, dntp_conc=0.8, dna_conc=250, temp_c=37)
def d3(a, b):     # 3'-anchored, extensible structures only (kcal/mol)
    return min(primer3.calc_end_stability(a, b, **C).dg,
               primer3.calc_end_stability(b, a, **C).dg) / 1000
def dany(a, b):   # best structure of any kind
    return primer3.calc_heterodimer(a, b, **C).dg / 1000
# homodimers: d3(a, a), dany via calc_homodimer(a)
# tailed primers: call on tail+primer, i.e. a = TF + F, b = TR + R
```

Thresholds (kcal/mol, 37 °C, 2 mM Mg2+):

| Quantity | Soft start (cost begins) | Hard cap (+1000) | Literature anchor |
|---|---|---|---|
| d3 (3'-anchored) | -4.0 | -6.0 (relax to -7.0 if infeasible) | NGS-PrimerPlex default -6; PrimerROC dimer-free cut-off -4.43 (own score units); SADDLE: 5 nt at 3' end gives dimers, which is about -6.3 (GC-rich) to -3.7 (mixed) at these conditions (4.2) |
| dany | -8.0 | -12.0 (untailed), -14 (tailed, tail-tail excluded) | NGS-PrimerPlex -10; varVAMP -9; PrimerPooler residual -6 at 60 °C |
| Pair-level end overlap | perfect 3'-terminal overlap of >= 6 nt is a hard fail regardless of dG | | varVAMP `has_end_overlap` (> 5 nt) |

Per-interaction cost `g(a,b) = max(0, -4 - d3)^2 + 0.25 * max(0, -8 - dany)^2 + 1000 * [d3 <= H3 or dany <= HA]`. For tailed primers compute d3 on tailed sequences but add a check that the gene-specific 3' end is the one paired (calc_end_stability already anchors the 3' end, which is always gene-specific).

Speed: P unique primers give P(P+1)/2 calls x 3. Measured 4.0 million interactions in 48 s on 16 cores (4.5). For more than about 3,000 primers, prefilter: only call primer3 for pairs where rc of one primer's last 4 nt occurs anywhere in the other (a hash lookup); all others get g = 0 for d3 and are scored with dany only if they share an 8-mer. Note that this prefilter is untested here.

### 5.3 Objective (one pool p, chosen pair s_i for each window i)

```
E(s) = sum_i U_i(s_i) + sum_{i<j, pool(i)=pool(j)} [ G(s_i, s_j) + w_off * OT(s_i, s_j) ]

U_i(x) = w_tm * [ (Tm_F - Tm0)^2 + (Tm_R - Tm0)^2 ]          # °C^2, w_tm = 1
       + w_fb * [ (fb_F - fb0)^2 + (fb_R - fb0)^2 ] / 0.05^2   # fraction bound at Ta, w_fb = 1 (when Ta is known)
       + w_len * ((L - 300) / 25)^2                         # w_len = 1
       + w_gc  * ((GC_amp - GC_ref) / 0.05)^2               # GC_ref = median amplicon GC of candidates, w_gc = 1
       + w_ss  * SS(x)                                      # amplicon/primer-site structure penalty from the structure module
       + 0.5 * primer3 pair penalty
       + g(F,R) + g(F,F) + g(R,R)                           # within-pair and self dimers
G(x, y) = sum over a in {F_x, R_x}, b in {F_y, R_y} of g(a, b)
OT(x, y) = 1000 if a primer of x and a primer of y form an off-target product < 2,000 bp on the reference, else 0
```

Weights were not tuned on wet-lab data. They put 1 °C of Tm deviation, 25 bp of length deviation, 5% of GC deviation and one kcal/mol beyond the soft dimer start on roughly the same scale. The user requirement "amplicon lengths within 150 bp of each other" is enforced as a hard filter (225 to 375 bp gives exactly a 150 bp band).

### 5.4 Optimiser

Exact search is only practical for tiny problems: the MILP took 158 s at 16 x 30, and enumeration is K^n. Use:

```
input: windows 1..n with candidate lists C_i (sorted by U_i), pools from 5.5
precompute U_i vectors and pairwise tables G between candidates of windows in the same pool
if prod |C_i| <= 1e6: enumerate exhaustively and return the optimum            # e.g. n <= 3 with K = 100
best = None
repeat R = 50 times:                                   # randomised greedy (multiply-style)
    order = random permutation of windows
    for i in order: s_i = argmin_k U_i(k) + sum_{j placed, same pool} G(k, s_j)
    keep best
for r in 1..R_sa (default 8, run in parallel):        # simulated annealing
    s = best greedy (r = 1) or random (r > 1); maintain local cost tables L_i(k)
    for t in 1..T_steps (default 30 * n * 10, at least 20,000):
        T = T0 * (T1 / T0)^(t / T_steps)               # T0 = 50, T1 = 0.05, geometric
        pick window i uniformly, candidate k uniformly from C_i
        dE = L_i(k) - L_i(s_i)
        if dE <= 0 or rand() < exp(-dE / T): update s_i and L_j for neighbours j
        track best
    descend: repeat { for each i: s_i = argmin_k L_i(k) } until no change   # ICM / 1-opt
if any hard violation remains:
    relax H3 from -6 to -7, rerun; if still violated, report the offending pairs and
    suggest moving one of the two windows to a second pool (graph colouring with the dimer edge added)
return best s, its E, worst d3 and dany in each pool, Tm / length / GC ranges
```

Measured: 16 windows: SA reaches the MILP optimum in all 8 restarts (0.44 s). 50 windows: SA 214.8 vs greedy 298.7 in 2.3 s. 100 windows: SA 816.6 (no hard violations) vs greedy 2292.1 (one violation) in 6.0 s. Optimisation time is negligible next to pair scoring, so spend it: 8 to 16 restarts in parallel, and step count about 300 x n.

### 5.5 Nearby targets: merge, tile, or separate pools

1. Merge positions into one window when the span between the first and last position is <= 250 bp (with primers of up to 25 nt and a 5 nt margin, a 375 bp amplicon can hold at most 315 bp of target; 302 bp failed in Primer3, 250 worked, section 4.3).
2. If the span is longer, tile it with overlapping sub-windows of <= 250 bp, overlap >= 30 bp, and give adjacent tiles different pools (PrimalScheme, Olivar, varVAMP all alternate pools for overlapping amplicons; PrimalScheme's stated reason is the short overlap product that would win in one tube).
3. Build a conflict graph over windows: edge if any candidate amplicons could overlap (window distance < 375 bp) or if the outer product (upstream F to downstream R) would be shorter than a spurious-product limit (1,500 bp used in exp3). Colour with DSATUR. If the graph has no edges, output one pool, which is what the user asked for.
4. Before accepting a second pool, try re-picking: for two windows 375 to 1,500 bp apart, restrict candidates so that the amplicons do not overlap and add an OT-style penalty for the outer product. This keeps one tube at the cost of a long competing product that a short extension time largely suppresses. This rule is a design judgement; I found no published validation of an outer-product length limit.
5. For mutation panels (as in the Mtb example), expect 0 or 2 pools: dense hotspots (rpoB RRDR, pncA, Rv0678, embB 306 to 406) need tiling.

## 6. Open questions

- No wet-lab data in this work tie the soft and hard dG thresholds to dimer read fractions under the user's polymerase and cycling; SADDLE's thresholds are k-mer based and NGS-PrimerPlex's -6 is a default, not a fitted value.
- Whether the 37 °C dimer scale should be replaced by a hot-start-aware temperature (for example the polymerase activation temperature) is untested; hot-start enzymes reduce set-up dimers, which may justify relaxing H3 to -7.
- The unary weights in 5.3 are uncalibrated.
- MFEprimer-3 `dimer` default cut-offs remain unverified (closed binary).
- PRISM details beyond the abstract were not read.
- The 4-mer prefilter in 5.2 and the outer-product limit in 5.5 have not been validated.
