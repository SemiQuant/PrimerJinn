# Multiplex primer design report

## PCR conditions

| Item | Value |
|---|---|
| genome | /Users/semiquant/Bioinformatics/Projects/Archived/PrimerJinn/example/ref.fasta |
| genome GC % | 65.61 |
| polymerase preset | kapa (KAPA2G Robust HotStart ReadyMix) |
| Tm scale | physical: SantaLucia 2004 + Owczarzy 2008 in the assumed buffer at the primer concentration |
| target primer Tm | 65.0 +/- 1.5 C |
| primer Tm range in pool | 64.1 to 67.8 C |
| same primers on the NEB Q5 calculator scale | 68.3 to 72.1 C |
| assumed buffer | 50 mM monovalent, 2 mM Mg2+, 0.8 mM dNTP total, 0% DMSO |
| primer concentration | 100 nM each |
| RECOMMENDED annealing temperature | 59 C |
| vendor rule Ta (KAPA: 55/60/65 C for primer Tm <55/55-65/>65) | 60.0 |
| equilibrium model Ta: all primers >= 90% bound on unfolded template, capped where structure leaves the weakest primer at 50% | 59 |
| temperature where every primer is >= 90% bound including structure (model) | 56 |
| temperature where every primer is >= 50% bound (model) | 59 |
| 90% binding temperature without template structure (model) | 60 |
| weakest-binding primer (model) | T12_pncA_F |
| dG_eff spread at 59 C (kcal/mol): full range / middle 50% | 2.28 / 0.59 |
| pools | 2 |
| amplicons | 15 |
| amplicon length range | 296 to 361 bp |
| adapter tails | none |
| optimiser | annealing, energy 206.5 |
| advice | Cycling (KAPA2G Robust data sheet): 95 C 3 min; then 95 C 15 s, anneal at the recommended Ta, 72 C 15 s (amplicons < 1 kb); final 72 C 1 min. KAPA advise a 15 s anneal because long anneals give non-specific products with this enzyme; multiplex pools often need longer, so compare 15 s with 60 s on the first pool. |
| advice | GC-rich amplicons (max GC 67%, max product Tm 95.1 C): KAPA recommend 5% DMSO for amplicons above 65% GC. DMSO lowers primer Tm by about 0.6 C per %; re-run with --dmso 5 so primers and Ta are scored with it. If GC-rich amplicons still fail or give extra bands, Musso et al. 2006 (Taq) needed all three of 1.3 M betaine, 5% DMSO and 50 uM 7-deaza-dGTP for 68 to 79% GC targets; single additives and betaine with DMSO were not enough. |
| advice | Annealing: 59 C, where template structure leaves the weakest primer at half binding (limiting primer T12_pncA_F). Run a gradient from 56 to 62 C. The KAPA data sheet rule for this primer Tm gives 60 C. |

## Coverage

All 24 requested positions lie between the primers of at least one amplicon.

Binding balance: these primers keep their Tm-matched core (uppercase) and carry extra template-matching bases at the 5' end (lowercase), because their template site is folded and they would otherwise bind more weakly than the rest of the pool at the annealing temperature. The 3' end and the covered positions are unchanged.

- T12_pncA_F: not extended (dG_eff -11.09, median -12.24): a strong dimer with T06_rpsL_K43+1_R (3' dG -6.0 kcal/mol); a strong dimer with T06_rpsL_K43+1_R (3' dG -6.5 kcal/mol)
- T13_pncA_R: +1 nt at the 5' end, dG_eff at 59 C -11.14 -> -12.74 kcal/mol (pool median -12.24)
- T10_fabG1_c-15+1_R: +3 nt at the 5' end, dG_eff at 59 C -11.14 -> -12.20 kcal/mol (pool median -12.24)
- T06_rpsL_K43+1_R: not extended (dG_eff -11.61, median -12.24): one more base overshoots the pool median
- T10_fabG1_c-15+1_F: not extended (dG_eff -11.68, median -12.24): one more base overshoots the pool median
- T09_rrs_1401+1_R: not extended (dG_eff -11.78, median -12.24): one more base overshoots the pool median
- T07_rplC_C154_R: not extended (dG_eff -11.85, median -12.24): one more base overshoots the pool median
- T04_Rv0678_R: not extended (dG_eff -11.85, median -12.24): one more base overshoots the pool median
- T11_katG_S315_F: not extended (dG_eff -11.90, median -12.24): one more base overshoots the pool median

## Primers

|   pool | primer             | sequence               |   length |   gc_percent |   tm_core_C |   tm_C |   dg_eff |   fraction_bound_model |   site_open_cost |   offtarget_primable_sites |
|-------:|:-------------------|:-----------------------|---------:|-------------:|------------:|-------:|---------:|-----------------------:|-----------------:|---------------------------:|
|      1 | T01_gyrB_N499_F    | CGGTCGCGATTCGATGTTCC   |       20 |         60   |        64.7 |   64.7 |   -12.99 |                  0.973 |             0.97 |                          1 |
|      1 | T01_gyrB_N499_R    | GCGGCGGTTGTGCCAAA      |       17 |         64.7 |        65.8 |   65.8 |   -12.34 |                  0.929 |             0.94 |                          0 |
|      1 | T03_rpoB_L430+4_F  | TCAAGGAGTTCTTCGGCACCAG |       22 |         54.5 |        65.6 |   65.6 |   -13.12 |                  0.977 |             0.83 |                          2 |
|      1 | T03_rpoB_L430+4_R  | CGGTCAGGTACACGATCTCGTC |       22 |         59.1 |        65.3 |   65.3 |   -12.33 |                  0.928 |             1.63 |                          0 |
|      1 | T04_Rv0678_F       | GACGCATACCGAACGTCACAG  |       21 |         57.1 |        64.3 |   64.3 |   -12.24 |                  0.918 |             1.91 |                          0 |
|      1 | T04_Rv0678_R       | GGCCGCAACCGGAAATAGG    |       19 |         63.2 |        65.1 |   65.1 |   -11.85 |                  0.863 |             1.78 |                          0 |
|      1 | T06_rpsL_K43+1_F   | GCCGGTAGATGCCAACCATC   |       20 |         60   |        64.5 |   64.5 |   -12.41 |                  0.936 |             1.32 |                          0 |
|      1 | T06_rpsL_K43+1_R   | CCGCGGATGATCTTGTAGCG   |       20 |         60   |        64.3 |   64.3 |   -11.61 |                  0.814 |             2.06 |                          2 |
|      1 | T07_rplC_C154_F    | CAGGTTGGGCAAGAGTTGACC  |       21 |         57.1 |        64.7 |   64.7 |   -12.16 |                  0.909 |             1.02 |                          0 |
|      1 | T07_rplC_C154_R    | CCGCACCCTTGATCAGCAG    |       19 |         63.2 |        64.7 |   64.7 |   -11.85 |                  0.862 |             1.8  |                          2 |
|      1 | T08_atpE_D28+1_F   | GGAGCTCGAAGAGGAACACCA  |       21 |         57.1 |        65.1 |   65.1 |   -12.19 |                  0.912 |             0.6  |                          1 |
|      1 | T08_atpE_D28+1_R   | GACGGGTGTAGCGAAGACGA   |       20 |         60   |        65.1 |   65.1 |   -12.81 |                  0.964 |             1.28 |                          2 |
|      1 | T09_rrs_1401+1_F   | TACGTTCCCGGGCCTTGTAC   |       20 |         60   |        65.4 |   65.4 |   -11.95 |                  0.88  |             1.53 |                          0 |
|      1 | T09_rrs_1401+1_R   | GTCCGAGTGTTGCCTCAGGA   |       20 |         60   |        65.3 |   65.3 |   -11.78 |                  0.849 |             1.15 |                          0 |
|      1 | T10_fabG1_c-15+1_F | CGCGGTCAGTTCCACACC     |       18 |         66.7 |        64.7 |   64.7 |   -11.68 |                  0.83  |             1.97 |                          2 |
|      1 | T10_fabG1_c-15+1_R | gttTGGCCCCTTCAGTGGCT   |       20 |         60   |        65.1 |   66.7 |   -12.2  |                  0.915 |             2.57 |                          1 |
|      1 | T11_katG_S315_F    | CGTCGGGGTGTTCGTCCATA   |       20 |         60   |        65.3 |   65.3 |   -11.9  |                  0.871 |             1.69 |                          1 |
|      1 | T11_katG_S315_R    | TGGGGCTGATCTACGTGAACC  |       21 |         57.1 |        65.2 |   65.2 |   -12.73 |                  0.96  |             0.79 |                          0 |
|      1 | T12_pncA_F         | GCCAACAGTTCATCCCGGTTC  |       21 |         57.1 |        65   |   65   |   -11.09 |                  0.665 |             2.32 |                          0 |
|      1 | T12_pncA_R         | GTCTGGACACGTCGGCAATC   |       20 |         60   |        64.5 |   64.5 |   -12.55 |                  0.947 |             0.97 |                         12 |
|      1 | T14_eis_c-10+1_F   | AGCACCACTTCACCAGGCA    |       19 |         57.9 |        65.4 |   65.4 |   -12.61 |                  0.952 |             0.51 |                          1 |
|      1 | T14_eis_c-10+1_R   | GCGTGATCCTTTGCCAGACAC  |       21 |         57.1 |        65.1 |   65.1 |   -12.32 |                  0.928 |             1.48 |                          2 |
|      1 | T15_embB_M306+1_F  | CGGACGACGGCTACATCCT    |       19 |         63.2 |        65   |   65   |   -13.37 |                  0.984 |             1.29 |                         11 |
|      1 | T15_embB_M306+1_R  | AGCACATAGGTGACCAGCGA   |       20 |         55   |        64.5 |   64.5 |   -12.24 |                  0.919 |             1.09 |                         21 |
|      2 | T02_gyrA_A90+1_F   | CCGGGTGCTCTATGCAATGTTC |       22 |         54.5 |        64.8 |   64.8 |   -12.41 |                  0.936 |             1.17 |                          0 |
|      2 | T02_gyrA_A90+1_R   | CGGCCGTCGTAGTTAGGGAT   |       20 |         60   |        64.8 |   64.8 |   -12.28 |                  0.923 |             1.77 |                          0 |
|      2 | T05_Rv0678_F       | CCGGATGCTGATCCAATTTGGG |       22 |         54.5 |        65   |   65   |   -12.06 |                  0.897 |             2.37 |                          0 |
|      2 | T05_Rv0678_R       | CGGTCAGATTGCGAGGTTGC   |       20 |         60   |        65.1 |   65.1 |   -12.28 |                  0.923 |             2.45 |                          0 |
|      2 | T13_pncA_F         | CCGGTGTAGGCACCCTTG     |       18 |         66.7 |        64.1 |   64.1 |   -11.96 |                  0.881 |             1.81 |                          3 |
|      2 | T13_pncA_R         | cGCGTCGGTAGGCAAACTGC   |       20 |         65   |        65.3 |   67.8 |   -12.74 |                  0.96  |             2.68 |                          1 |

## Amplicons

|   pool | target   | covers                                                                                                                                                         |   length |   gc_percent |   product_tm_C |   self_priming_top |   self_priming_bottom |   strongest_hairpin_dg_at_72C | flags                                                                                                                                      |
|-------:|:---------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------|---------:|-------------:|---------------:|-------------------:|----------------------:|------------------------------:|:-------------------------------------------------------------------------------------------------------------------------------------------|
|      1 | T01      | gyrB_N499 (NC_000962.3:6735)                                                                                                                                   |      296 |         57.8 |           90.9 |              0.05  |                 0.087 |                         -4.46 |                                                                                                                                            |
|      1 | T03      | rpoB_L430 (NC_000962.3:761095); rpoB_D435 (NC_000962.3:761110); rpoB_H445 (NC_000962.3:761140); rpoB_S450 (NC_000962.3:761155); rpoB_I491 (NC_000962.3:761277) |      318 |         65.4 |           94.1 |              0.029 |                 0.007 |                         -3.83 |                                                                                                                                            |
|      1 | T04      | Rv0678 (NC_000962.3:778990-779260, part of 778990-779487)                                                                                                      |      333 |         61.3 |           92.5 |              0.008 |                 0.127 |                         -1.91 |                                                                                                                                            |
|      1 | T06      | rpsL_K43 (NC_000962.3:781687); rpsL_K88 (NC_000962.3:781822)                                                                                                   |      307 |         64.2 |           93.6 |              0.253 |                 0.044 |                         -3.47 | top strand 3' end paired with p=0.25 (self-priming risk)                                                                                   |
|      1 | T07      | rplC_C154 (NC_000962.3:801268)                                                                                                                                 |      307 |         65.8 |           94.2 |              0.007 |                 0.018 |                         -3.05 | 50-bp window with 80% GC                                                                                                                   |
|      1 | T08      | atpE_D28 (NC_000962.3:1461127); atpE_A63 (NC_000962.3:1461232)                                                                                                 |      313 |         61.3 |           92.4 |              0.002 |                 0.014 |                         -3.27 |                                                                                                                                            |
|      1 | T09      | rrs_1401 (NC_000962.3:1473246); rrs_1484 (NC_000962.3:1473329)                                                                                                 |      305 |         59.7 |           91.7 |              0.085 |                 0.001 |                         -2.51 | fwd: site dG_open 2.3 > 1.5 kcal/mol at 59 C (penalise)                                                                                    |
|      1 | T10      | fabG1_c-15 (NC_000962.3:1673425); fabG1_c-8 (NC_000962.3:1673432)                                                                                              |      297 |         61.2 |           92.2 |              0.002 |                 0.045 |                         -2.49 |                                                                                                                                            |
|      1 | T11      | katG_S315 (NC_000962.3:2155168)                                                                                                                                |      305 |         66.6 |           94.5 |              0.008 |                 0.047 |                         -1.95 | fwd: site dG_open 1.5 > 1.5 kcal/mol at 59 C (penalise); 50-bp window with 80% GC                                                          |
|      1 | T12      | pncA (NC_000962.3:2288681-2288971, part of 2288681-2289241)                                                                                                    |      342 |         64   |           93.7 |              0.036 |                 0.019 |                         -2.58 |                                                                                                                                            |
|      1 | T14      | eis_c-10 (NC_000962.3:2715342); eis_c-14 (NC_000962.3:2715346)                                                                                                 |      306 |         63.4 |           93.2 |              0.012 |                 0.011 |                         -3.26 |                                                                                                                                            |
|      1 | T15      | embB_M306 (NC_000962.3:4247429-4247431); embB_G406 (NC_000962.3:4247730)                                                                                       |      361 |         67.3 |           95.1 |              0.001 |                 0.05  |                         -3.43 | 50-bp window with 80% GC                                                                                                                   |
|      2 | T02      | gyrA_A90 (NC_000962.3:7570); gyrA_D94 (NC_000962.3:7582)                                                                                                       |      318 |         66   |           94.4 |              0.065 |                 0.017 |                         -5.19 |                                                                                                                                            |
|      2 | T05      | Rv0678 (NC_000962.3:779224-779487, part of 778990-779487)                                                                                                      |      309 |         65.7 |           94.2 |              0.031 |                 0.039 |                         -1.03 | 50-bp window with 82% GC                                                                                                                   |
|      2 | T13      | pncA (NC_000962.3:2288958-2289241, part of 2288681-2289241)                                                                                                    |      353 |         65.1 |           94.2 |              0.015 |                 0.012 |                         -1.09 | fwd: site dG_open 1.7 > 1.5 kcal/mol at 59 C (penalise); rev: site dG_open 2.1 > 1.5 kcal/mol at 59 C (penalise); 50-bp window with 82% GC |

## Strongest dimers per pool (37 C)

|   pool | primer1          | primer2            |   d3_kcal |   dany_kcal | call   |
|-------:|:-----------------|:-------------------|----------:|------------:|:-------|
|      1 | T01_gyrB_N499_R  | T10_fabG1_c-15+1_R |     -5.47 |       -5.47 | soft   |
|      1 | T01_gyrB_N499_F  | T08_atpE_D28+1_F   |     -5.29 |       -5.29 | soft   |
|      2 | T05_Rv0678_R     | T13_pncA_R         |     -4.98 |       -4.98 | soft   |
|      1 | T11_katG_S315_R  | T12_pncA_R         |     -4.87 |       -6.15 | soft   |
|      1 | T06_rpsL_K43+1_F | T14_eis_c-10+1_F   |     -4.83 |       -4.83 | soft   |
|      1 | T08_atpE_D28+1_R | T10_fabG1_c-15+1_F |     -4.63 |       -4.63 | soft   |
|      2 | T02_gyrA_A90+1_F | T13_pncA_F         |     -4.56 |       -7.78 | soft   |
|      2 | T05_Rv0678_F     | T13_pncA_F         |     -3.46 |       -4.5  | ok     |
|      2 | T05_Rv0678_R     | T13_pncA_F         |     -3.42 |       -3.42 | ok     |
|      2 | T02_gyrA_A90+1_R | T05_Rv0678_F       |     -3.31 |       -4.11 | ok     |

## Other products (genome-wide, within pools)

kind: 'outer product' = forward primer of one designed amplicon with the reverse primer of a neighbouring one (normal in tiled pools, suppressed by a short extension); 'repeated locus' = the designed pair also amplifies another copy of the locus; 'unintended' = off-target priming.

|   pool | kind       | fwd               | rev               | contig      |   start |     end |   size |   fwd_mismatches |   rev_mismatches |   score | call     |
|-------:|:-----------|:------------------|:------------------|:------------|--------:|--------:|-------:|-----------------:|-----------------:|--------:|:---------|
|      1 | unintended | T01_gyrB_N499_R   | T07_rplC_C154_R   | NC_000962.3 |  801260 |  801400 |    141 |                4 |                0 |  0.0035 | possible |
|      1 | unintended | T07_rplC_C154_F   | T03_rpoB_L430+4_R | NC_000962.3 |  801094 |  802517 |   1424 |                0 |               10 |  0.003  | possible |
|      1 | unintended | T01_gyrB_N499_R   | T01_gyrB_N499_R   | NC_000962.3 |    6627 |    6941 |    315 |                3 |                0 |  0.0023 | possible |
|      1 | unintended | T03_rpoB_L430+4_R | T07_rplC_C154_R   | NC_000962.3 |  799853 |  801400 |   1548 |                7 |                0 |  0.0007 | possible |
|      1 | unintended | T06_rpsL_K43+1_F  | T14_eis_c-10+1_F  | NC_000962.3 |  781552 |  781846 |    295 |                0 |               10 |  0.0005 | possible |
|      1 | unintended | T01_gyrB_N499_F   | T12_pncA_R        | NC_000962.3 |    6646 |    6845 |    200 |                0 |               10 |  0.0005 | possible |
|      1 | unintended | T15_embB_M306+1_R | T03_rpoB_L430+4_R | NC_000962.3 |  760506 |  761385 |    880 |               11 |                0 |  0.0005 | possible |
|      1 | unintended | T12_pncA_F        | T15_embB_M306+1_R | NC_000962.3 | 4247123 | 4247766 |    644 |                8 |                0 |  0.0003 | possible |
|      1 | unintended | T15_embB_M306+1_R | T11_katG_S315_R   | NC_000962.3 | 2154808 | 2155438 |    631 |                8 |                0 |  0.0003 | possible |
|      1 | unintended | T03_rpoB_L430+4_F | T14_eis_c-10+1_R  | NC_000962.3 | 2715197 | 2715408 |    212 |               10 |                0 |  0.0002 | possible |
|      1 | unintended | T08_atpE_D28+1_R  | T04_Rv0678_R      | NC_000962.3 |  777428 |  779279 |   1852 |                7 |                0 |  0.0002 | possible |
|      1 | unintended | T14_eis_c-10+1_F  | T03_rpoB_L430+4_F | NC_000962.3 | 1010210 | 1011582 |   1373 |                7 |                7 |  0.0001 | possible |

## Notes

- tm_C is on the scale named in the conditions table; tm_neb_q5_scale_C reproduces the NEB Tm calculator for Q5 at 500 nM for comparison. dG_eff is the binding free energy at the recommended Ta (SantaLucia 2004 + Owczarzy 2008 in the assumed buffer) plus the cost of the primer's own hairpin and of opening template structure at its site; the pool is optimised for a narrow dG_eff spread.
- Model fraction bound ignores the polymerase. For Taq-type enzymes it is a direct estimate; Q5 and Phusion carry a DNA-binding domain and are run above the primer Tm by vendor rule, so for them use dG_eff to compare primers, not as an absolute occupancy. Buffer monovalent salt is assumed (not published by vendors), so absolute values carry about 1 to 3 C uncertainty.
- Structure, dimer and off-target thresholds are model-based starting points, not calibrated on sequencing data. See docs/DESIGN.md.
