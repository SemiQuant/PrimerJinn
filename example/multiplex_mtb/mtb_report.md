# Multiplex primer design report

## PCR conditions

| Item | Value |
|---|---|
| genome | /Users/semiquant/Bioinformatics/Projects/Archived/PrimerJinn/example/ref.fasta |
| genome GC % | 65.61 |
| polymerase preset | kapa (KAPA2G Robust HotStart ReadyMix) |
| Tm scale | physical: SantaLucia 2004 + Owczarzy 2008 in the assumed buffer at the primer concentration |
| target primer Tm | 65.0 +/- 1.5 C |
| primer Tm range in pool | 63.6 to 65.6 C |
| same primers on the NEB Q5 calculator scale | 68.0 to 70.2 C |
| assumed buffer | 50 mM monovalent, 2 mM Mg2+, 0.8 mM dNTP total, 0% DMSO |
| primer concentration | 100 nM each |
| RECOMMENDED annealing temperature | 60 C |
| vendor rule Ta (KAPA: 55/60/65 C for primer Tm <55/55-65/>65) | 60.0 |
| equilibrium model Ta: all primers >= 90% bound on unfolded template, capped where structure leaves the weakest primer at 50% | 60 |
| temperature where every primer is >= 90% bound including structure (model) | 55 |
| temperature where every primer is >= 50% bound (model) | 60 |
| 90% binding temperature without template structure (model) | 60 |
| weakest-binding primer (model) | T10_fabG1_c-15+1_F |
| dG_eff spread at 60 C (kcal/mol) | 1.88 |
| pools | 2 |
| amplicons | 15 |
| amplicon length range | 289 to 369 bp |
| adapter tails | none |
| optimiser | annealing, energy 117.3 |
| advice | Cycling (KAPA2G Robust data sheet): 95 C 3 min; then 95 C 15 s, anneal at the recommended Ta, 72 C 15 s (amplicons < 1 kb); final 72 C 1 min. KAPA advise a 15 s anneal because long anneals give non-specific products with this enzyme; multiplex pools often need longer, so compare 15 s with 60 s on the first pool. |
| advice | GC-rich amplicons (max GC 67%, max product Tm 95.0 C): KAPA recommend 5% DMSO for amplicons above 65% GC. DMSO lowers primer Tm by about 0.6 C per %; re-run with --dmso 5 so primers and Ta are scored with it. |
| advice | Annealing: 60 C, where every primer is at least 90% bound on unfolded template (limiting primer T10_fabG1_c-15+1_F). Run a gradient from 57 to 63 C. The KAPA data sheet rule for this primer Tm gives 60 C. |

## Coverage

All 24 requested positions lie between the primers of at least one amplicon.

## Primers

|   pool | primer             | sequence                 |   length |   gc_percent |   tm_C |   dg_eff |   fraction_bound_model |   site_open_cost |   offtarget_primable_sites |
|-------:|:-------------------|:-------------------------|---------:|-------------:|-------:|---------:|-----------------------:|-----------------:|---------------------------:|
|      1 | T01_gyrB_N499_F    | CGGTCGCGATTCGATGTTCC     |       20 |         60   |   64.7 |   -12.05 |                  0.889 |             0.92 |                          0 |
|      1 | T01_gyrB_N499_R    | GCGCTGCCACTTGAGTTTGT     |       20 |         55   |   65   |   -11.68 |                  0.822 |             1.36 |                          0 |
|      1 | T03_rpoB_L430+4_F  | CGCCGCGATCAAGGAGTTC      |       19 |         63.2 |   64.8 |   -11.46 |                  0.768 |             1.4  |                          2 |
|      1 | T03_rpoB_L430+4_R  | GGCGGTCAGGTACACGATCT     |       20 |         60   |   65   |   -11.39 |                  0.746 |             1.63 |                          1 |
|      1 | T04_Rv0678_F       | ACGCATACCGAACGTCACAGA    |       21 |         52.4 |   64.8 |   -11.32 |                  0.728 |             1.72 |                          0 |
|      1 | T04_Rv0678_R       | GCCGCAACCGGAAATAGGTG     |       20 |         60   |   65.2 |   -11.05 |                  0.639 |             1.68 |                          0 |
|      1 | T06_rpsL_K43+1_F   | CAAGGGTCGTCGGGACAAGA     |       20 |         60   |   65.3 |   -12.55 |                  0.945 |             0.56 |                          0 |
|      1 | T06_rpsL_K43+1_R   | CCTGTTTGCGGTTCTTGACACC   |       22 |         54.5 |   65.3 |   -11.63 |                  0.81  |             1.77 |                          0 |
|      1 | T07_rplC_C154_F    | CAGGTTGGGCAAGAGTTGACC    |       21 |         57.1 |   64.7 |   -11.52 |                  0.784 |             0.97 |                          0 |
|      1 | T07_rplC_C154_R    | GCACCCTTGATCAGCAGCAC     |       20 |         60   |   65.2 |   -11.32 |                  0.727 |             1.8  |                          1 |
|      1 | T08_atpE_D28+1_F   | GCGGTGGACTGATCATGGC      |       19 |         63.2 |   64.7 |   -11.44 |                  0.762 |             1.34 |                          1 |
|      1 | T08_atpE_D28+1_R   | GACTCGCCGCCTTCCTCT       |       18 |         66.7 |   65.4 |   -11.39 |                  0.749 |             1.55 |                          3 |
|      1 | T09_rrs_1401+1_F   | TACGTTCCCGGGCCTTGTAC     |       20 |         60   |   65.4 |   -11.7  |                  0.826 |             1.45 |                          0 |
|      1 | T09_rrs_1401+1_R   | TGTTGCCTCAGGACCCAACA     |       20 |         55   |   65.2 |   -10.99 |                  0.619 |             2    |                          0 |
|      1 | T10_fabG1_c-15+1_F | GTTACGCTCGTGGACATACCGA   |       22 |         54.5 |   65.1 |   -10.71 |                  0.515 |             2.62 |                          0 |
|      1 | T10_fabG1_c-15+1_R | CGGACCCTGGTGCTCTTCT      |       19 |         63.2 |   65.2 |   -11.19 |                  0.685 |             1.73 |                          0 |
|      1 | T11_katG_S315_F    | TCGGGGTGTTCGTCCATACG     |       20 |         60   |   65.3 |   -11.42 |                  0.756 |             1.71 |                          0 |
|      1 | T11_katG_S315_R    | TGGGGCTGATCTACGTGAACC    |       21 |         57.1 |   65.2 |   -12.41 |                  0.932 |             0.75 |                          0 |
|      1 | T12_pncA_F         | GCCAACAGTTCATCCCGGTTC    |       21 |         57.1 |   65   |   -10.8  |                  0.55  |             2.29 |                          0 |
|      1 | T12_pncA_R         | GGACACGTCGGCAATCGAG      |       19 |         63.2 |   64.4 |   -11.81 |                  0.849 |             0.88 |                          4 |
|      1 | T14_eis_c-10+1_F   | ACCGTCAACCGCAGATCCA      |       19 |         57.9 |   65.1 |   -12.59 |                  0.948 |             0.29 |                          0 |
|      1 | T14_eis_c-10+1_R   | CGCGTGATCCTTTGCCAGAC     |       20 |         60   |   65.1 |   -11.67 |                  0.819 |             1.43 |                          0 |
|      1 | T15_embB_M306+1_F  | GTCGGACGACGGCTACATC      |       19 |         63.2 |   63.6 |   -10.95 |                  0.603 |             1.36 |                          2 |
|      1 | T15_embB_M306+1_R  | TCGATCAGCACATAGGTGACCAG  |       23 |         52.2 |   64.8 |   -11.5  |                  0.779 |             1.44 |                          0 |
|      2 | T02_gyrA_A90+1_F   | CCGGGTGCTCTATGCAATGTTC   |       22 |         54.5 |   64.8 |   -11.72 |                  0.83  |             1.09 |                          0 |
|      2 | T02_gyrA_A90+1_R   | CGGCCGTCGTAGTTAGGGAT     |       20 |         60   |   64.8 |   -11.21 |                  0.694 |             1.71 |                          0 |
|      2 | T05_Rv0678_F       | CGGATGCTGATCCAATTTGGGTTC |       24 |         50   |   65.1 |   -11.27 |                  0.712 |             2.24 |                          0 |
|      2 | T05_Rv0678_R       | TCGCCACCTCGGTCAGATTG     |       20 |         60   |   65.6 |   -10.83 |                  0.561 |             2.45 |                          0 |
|      2 | T13_pncA_F         | TCCGGTGTAGGCACCCTTG      |       19 |         63.2 |   65.5 |   -10.95 |                  0.604 |             1.74 |                          1 |
|      2 | T13_pncA_R         | GCGTCGGTAGGCAAACTGC      |       19 |         63.2 |   65.3 |   -10.71 |                  0.515 |             2.39 |                          0 |

## Amplicons

|   pool | target   | covers                                                                                                                                                         |   length |   gc_percent |   product_tm_C |   self_priming_top |   self_priming_bottom |   strongest_hairpin_dg_at_72C | flags                                                                                                                                      |
|-------:|:---------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------|---------:|-------------:|---------------:|-------------------:|----------------------:|------------------------------:|:-------------------------------------------------------------------------------------------------------------------------------------------|
|      1 | T01      | gyrB_N499 (NC_000962.3:6735)                                                                                                                                   |      319 |         57.4 |           90.8 |              0.028 |                 0.076 |                         -4.46 |                                                                                                                                            |
|      1 | T03      | rpoB_L430 (NC_000962.3:761095); rpoB_D435 (NC_000962.3:761110); rpoB_H445 (NC_000962.3:761140); rpoB_S450 (NC_000962.3:761155); rpoB_I491 (NC_000962.3:761277) |      328 |         66.2 |           94.5 |              0.02  |                 0.075 |                         -3.83 |                                                                                                                                            |
|      1 | T04      | Rv0678 (NC_000962.3:778990-779258, part of 778990-779487)                                                                                                      |      331 |         61   |           92.4 |              0.011 |                 0.043 |                         -1.91 |                                                                                                                                            |
|      1 | T06      | rpsL_K43 (NC_000962.3:781687); rpsL_K88 (NC_000962.3:781822)                                                                                                   |      311 |         63   |           93.1 |              0.013 |                 0.013 |                         -3.47 | rev: site dG_open 1.9 > 1.5 kcal/mol at 60 C (penalise)                                                                                    |
|      1 | T07      | rplC_C154 (NC_000962.3:801268)                                                                                                                                 |      305 |         65.6 |           94.1 |              0.021 |                 0.016 |                         -3.05 | 50-bp window with 80% GC                                                                                                                   |
|      1 | T08      | atpE_D28 (NC_000962.3:1461127); atpE_A63 (NC_000962.3:1461232)                                                                                                 |      304 |         61.2 |           92.3 |              0.018 |                 0.005 |                         -3.27 |                                                                                                                                            |
|      1 | T09      | rrs_1401 (NC_000962.3:1473246); rrs_1484 (NC_000962.3:1473329)                                                                                                 |      298 |         59.4 |           91.5 |              0.015 |                 0.001 |                         -2.51 | fwd: site dG_open 2.2 > 1.5 kcal/mol at 60 C (penalise)                                                                                    |
|      1 | T10      | fabG1_c-15 (NC_000962.3:1673425); fabG1_c-8 (NC_000962.3:1673432)                                                                                              |      302 |         65.6 |           94.1 |              0.02  |                 0.016 |                         -4.06 |                                                                                                                                            |
|      1 | T11      | katG_S315 (NC_000962.3:2155168)                                                                                                                                |      303 |         66.3 |           94.4 |              0.008 |                 0.037 |                         -1.95 | fwd: site dG_open 1.5 > 1.5 kcal/mol at 60 C (penalise); 50-bp window with 80% GC                                                          |
|      1 | T12      | pncA (NC_000962.3:2288681-2288968, part of 2288681-2289241)                                                                                                    |      338 |         64.2 |           93.7 |              0.009 |                 0.017 |                         -2.58 |                                                                                                                                            |
|      1 | T14      | eis_c-10 (NC_000962.3:2715342); eis_c-14 (NC_000962.3:2715346)                                                                                                 |      289 |         63.7 |           93.2 |              0.028 |                 0.003 |                         -3.26 |                                                                                                                                            |
|      1 | T15      | embB_M306 (NC_000962.3:4247429-4247431); embB_G406 (NC_000962.3:4247730)                                                                                       |      369 |         66.9 |           95   |              0.016 |                 0.062 |                         -3.43 | 50-bp window with 80% GC                                                                                                                   |
|      2 | T02      | gyrA_A90 (NC_000962.3:7570); gyrA_D94 (NC_000962.3:7582)                                                                                                       |      318 |         66   |           94.4 |              0.055 |                 0.016 |                         -5.19 |                                                                                                                                            |
|      2 | T05      | Rv0678 (NC_000962.3:779227-779487, part of 778990-779487)                                                                                                      |      317 |         65.6 |           94.2 |              0.013 |                 0.008 |                         -1.03 | rev: site dG_open 2.6 > 1.5 kcal/mol at 60 C (penalise); 50-bp window with 82% GC                                                          |
|      2 | T13      | pncA (NC_000962.3:2288958-2289241, part of 2288681-2289241)                                                                                                    |      353 |         64.9 |           94.1 |              0.014 |                 0.002 |                         -1.09 | fwd: site dG_open 1.7 > 1.5 kcal/mol at 60 C (penalise); rev: site dG_open 1.9 > 1.5 kcal/mol at 60 C (penalise); 50-bp window with 82% GC |

## Strongest dimers per pool (37 C)

|   pool | primer1            | primer2           |   d3_kcal |   dany_kcal | call   |
|-------:|:-------------------|:------------------|----------:|------------:|:-------|
|      1 | T08_atpE_D28+1_F   | T12_pncA_F        |     -5.48 |       -5.48 | soft   |
|      1 | T11_katG_S315_F    | T15_embB_M306+1_F |     -5.31 |       -6.54 | soft   |
|      1 | T01_gyrB_N499_F    | T12_pncA_R        |     -5.24 |       -5.93 | soft   |
|      1 | T08_atpE_D28+1_F   | T11_katG_S315_F   |     -4.85 |       -4.85 | soft   |
|      1 | T10_fabG1_c-15+1_R | T15_embB_M306+1_R |     -4.8  |       -5.35 | soft   |
|      2 | T02_gyrA_A90+1_F   | T13_pncA_F        |     -4.56 |       -7.78 | soft   |
|      2 | T05_Rv0678_F       | T05_Rv0678_F      |     -4.22 |       -6.07 | soft   |
|      2 | T02_gyrA_A90+1_R   | T05_Rv0678_F      |     -3.31 |       -3.62 | ok     |
|      2 | T05_Rv0678_R       | T13_pncA_R        |     -2.94 |       -3.85 | ok     |
|      2 | T05_Rv0678_F       | T13_pncA_F        |     -2.64 |       -3.99 | ok     |

## Other products (genome-wide, within pools)

kind: 'outer product' = forward primer of one designed amplicon with the reverse primer of a neighbouring one (normal in tiled pools, suppressed by a short extension); 'repeated locus' = the designed pair also amplifies another copy of the locus; 'unintended' = off-target priming.

|   pool | kind       | fwd               | rev                | contig      |   start |     end |   size |   fwd_mismatches |   rev_mismatches |   score | call     |
|-------:|:-----------|:------------------|:-------------------|:------------|--------:|--------:|-------:|-----------------:|-----------------:|--------:|:---------|
|      1 | unintended | T14_eis_c-10+1_F  | T10_fabG1_c-15+1_R | NC_000962.3 | 2715121 | 2715224 |    104 |                0 |                4 |  0.0022 | possible |
|      1 | unintended | T12_pncA_F        | T15_embB_M306+1_R  | NC_000962.3 | 4247123 | 4247772 |    650 |                8 |                0 |  0.0002 | possible |
|      1 | unintended | T15_embB_M306+1_F | T15_embB_M306+1_R  | NC_000962.3 | 4247443 | 4247772 |    330 |                7 |                0 |  0.0001 | possible |

## Notes

- tm_C is on the scale named in the conditions table; tm_neb_q5_scale_C reproduces the NEB Tm calculator for Q5 at 500 nM for comparison. dG_eff is the binding free energy at the recommended Ta (SantaLucia 2004 + Owczarzy 2008 in the assumed buffer) plus the cost of the primer's own hairpin and of opening template structure at its site; the pool is optimised for a narrow dG_eff spread.
- Model fraction bound ignores the polymerase. For Taq-type enzymes it is a direct estimate; Q5 and Phusion carry a DNA-binding domain and are run above the primer Tm by vendor rule, so for them use dG_eff to compare primers, not as an absolute occupancy. Buffer monovalent salt is assumed (not published by vendors), so absolute values carry about 1 to 3 C uncertainty.
- Structure, dimer and off-target thresholds are model-based starting points, not calibrated on sequencing data. See docs/DESIGN.md.
