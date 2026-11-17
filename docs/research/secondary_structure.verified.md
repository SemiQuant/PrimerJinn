# Fact-check: secondary_structure.md

Checked 2026-11-11. Scripts and outputs are in `verify_ss2/` (this pass: `v2.py`, `v3.py`, `v4.py`, `p3tm.py`, `demo.out`, `pytest.out`, `dg.out`). A previous attempt left reruns in `verify_ss/` (`api.out`, `checks.out`, `bench/bench2_rerun.out`, `v1.out`, `gc_bias.out`), which I re-read and used only where this pass's own runs agree with them. Literature: full texts in `verify_ss/lit/` (Europe PMC XML) and PubMed metadata.

Units: kcal/mol, degrees C, mol/L.

## Verdicts on the key claims

| # | Claim | Verdict | Evidence / correction |
|---|---|---|---|
| 1 | Vienna Mathews2004 DNA set is non-physical at high T (rpoB MFE -45.5 / -64.8 / -129.0 at 55 / 80 / 99 C; NUPACK -36.2 / -8.9 / -1.2) | Confirmed | `ss2/t1.py` rerun gives identical numbers. Reproduced on two new sequences (`v2.out` section 2): M2004 MFE goes from -17.1 at 55 C to -119.5 at 99 C (Mtb 137013) while NUPACK goes to -1.8. 1x2 loop cost 3.3 / 1.59 / -0.02 / -1.72 reproduced (`verify_ss/checks.out`). |
| 2 | Vienna M1999 matches NUPACK G_ens with MAE 1.84 to 2.14 kcal/mol; M2004 error 9 to 35 | Corrected | The numbers reproduce on the report's 10 loci (56.7 to 68.7% GC), but they do **not** generalise. On new 300 nt sequences in PCR salt (`v2.out` section 1): 40 to 55% GC, Vienna is 3 to 9 kcal/mol **less** stable than NUPACK (Mtb 52% GC windows +6.5 to +8.8); about 80% GC Mtb windows, Vienna is 3 to 14 kcal/mol **more** stable. Only 60 to 70% GC sequences are within about 2 kcal/mol, and some (rnd60_2, rnd70_2) are off by 5 to 6. Cause, in part: M1999 helix energies differ from SantaLucia. A 20-bp AT-rich duplex at 65 C gives -3.0 (M1999), -4.6 (primer3 thal), -7.1 (NUPACK) (`v3.py`). The M2004 9 to 35 figure applies at 1 M salt on the bench loci only. |
| 3 | DNA param loaders return 1, set helical_rise 3.4 and backbone 6.76, leave noGU 0; rpoB 65 C fold has 10 G-T pairs | Confirmed | `v2.out` section 3. On a new sequence, 6 G-T pairs with noGU 0 and 0 with noGU 1. |
| 4 | md.salt is monovalent, in M, default 1.021; no Mg term; rpoB G_ens -33.9 at 1.021 M, -13.7 at 0.05 M | Confirmed | `model.h`: "Salt (monovalent) concentration (M)", `VRNA_MODEL_DEFAULT_SALT 1.021`. `salt.c` cites Chen and Znosko 2013 for duplex initiation. No md attribute mentions Mg. Salt scan reproduced. |
| 5 | NUPACK 4.0.2 uses Na_eq = Na + 3.3 sqrt(Mg) and enforces Na 0.05 to 1.1 M, Mg 0 to 0.2 M | Confirmed | On a new random sequence, at 55 C with Na 0.1 M, the fitted Na-only equivalent equals 0.1 + 3.3 sqrt(Mg) to 5 digits for Mg 1, 3 and 10 mM. Ranges match the assertion errors and docs.nupack.org/model. The docs do not state the formula itself; it is verified only numerically. |
| 6 | probs_window equals the hard-constraint dG_open; RNA.pfl_fold_up gives 0.0013 because it uses the global 37 C model | Confirmed, with a nuance | On a new sequence at 60 C, both methods give the same value at 3 sites (for example 0.7061 and 0.7061). pfl_fold_up returns 0.0003 with the default globals, but returns exactly the probs_window value (0.1158) after setting `RNA.cvar.temperature/noGU/salt`. So it does not ignore the model; it reads the globals. probs_window is still the right choice because it is thread-safe. |
| 7 | PCR salt, 65 C: G_ens -19.4 / -20.6, 31% / 27% paired, P22 0.28 / 0.17; 72 C 26% / 21% paired | Confirmed (for these 10 GC-rich loci only) | bench2 rerun identical apart from NUPACK sampling noise (±0.003). Does not transfer to lower-GC genomes (see #2). |
| 8 | Vienna-NUPACK site accessibility Spearman 0.83 to 0.74; single sites differ up to 5 to 20 fold | Corrected | Spearman reproduced (0.834 to 0.735). The per-site disagreement reaches at least 23-fold (7576 minus strand: 0.010 vs 0.233). The body text twice says "up to about 5" (Key decision 6, Open question 2), which is wrong. On a new sequence, site 0-20 gave 0.706 (Vienna) vs 0.217 (NUPACK). |
| 9 | Genomic site 1.4 to 6 times less accessible than at the amplicon 3' end; plfold W200/L150 reproduces full 900 nt fold with r = 0.98 | Corrected | On 12 new random Mtb sites (`v4.py`), amplicon/genomic ratio was 1.9 to 87, with r(log10) amplicon vs genomic only 0.50 (report: 0.96). plfold vs full fold r = 0.97 confirmed, but plfold overestimates P(open) by 1 to 10 fold (for example 0.0075 vs 0.0007). Also, a ±300 nt context is itself an approximation. |
| 10 | Gao 2006: 3 intramolecular bp halve k_on; k_on 12.0e5 to 2.0e5; >= 4 bp gives non-two-state kinetics | Confirmed, with a nuance | Table 2: 12.0, 7.2, 3.0, 2.0 (x10^5 M^-1 s^-1), in 0.5 M NaCl at **20 C**, 1 µM strands. Only P4T4 (4 bp in both strands) failed the two-state fit; P4T3 fitted two-state. Using k_on 1e6 at 65 C in PCR salt is an extrapolation. |
| 11 | Aird 2011: >70% GC depleted to <1/100; 3 min + 80 s denaturation + 2 M betaine; 94% GC windows > half coverage | Confirmed | Full text. Also relevant to Mtb: after 10 cycles, standard Phusion depleted >65% GC loci to about 1/100 (Fig 1e). This is library PCR, not target amplicon PCR. |
| 12 | Gimpel 2025: self-priming hairpins, not GC, main cause of low efficiency | Confirmed, with a nuance | The paper calls the mechanism **adapter-mediated** self-priming: a CGTG-type motif complementary to the common 5' adapter. GCall and GCfix pools behaved alike. AUROC 0.88. |
| 13 | Rees 1993: betaine isostabilising at about 5.2 M | Confirmed | PubMed abstract (doi:10.1021/bi00052a019). "Only partly reduces GC dependence at 1 to 2 M" is the report's inference, not from the abstract. |
| 14 | Huppert G3+N1-7 rule in 100 mM KCl; Woodford/Weitzmann K+-dependent arrest | Confirmed | Huppert full text: 100 mM KCl, 10 mM Tris pH 7.4. Woodford 1994 abstract: arrest only on the G-rich template, K+-dependent, removed by N7 blocking. Perrone 2017: 45 genes confirmed. |
| 15 | Product Tm formula and Na_eq; rpoB 94.5 C, katG 93.6 C; primer3-py about 2 C lower, cause not traced | Confirmed, cause found | oligotm.c `long_seq_tm` and `divalent_to_monovalent` (120 sqrt(Mg - dNTP)) match digit for digit, and I recomputed 94.53 and 93.57. **Cause of the 2 C gap:** primer3-py 2.0.3 `argdefaults.py` sets `formamide_conc = 0.8` mol/L, so the formamide term (0.453 GC - 2.88) x 0.8 is applied. With formamide_conc=0 the values agree exactly. The same default lowers **primer** Tm from `calc_tm` (22-mer: 70.21 vs 68.15 C). The primer3 manual default for PRIMER_FORMAMIDE_CONC is 0.0. |

Other citations checked against PubMed (all metadata correct): Henke 1997 NAR 25:3957 (doi:10.1093/nar/25.19.3957), Frey 2008 Nat Protoc 3:1312 (doi:10.1038/nprot.2008.112: >83% GC, 2.5 C/s ramp, 1.5 C/s cooling, 48 cycles), Chester and Marshak 1993 Anal Biochem 209:284 (doi:10.1006/abio.1993.1121), McConlogue 1988 NAR 16:9869 (PMID 3186456), Benita 2003 (84.3 / 94.8%, 61% GC, 21 bp), Musso 2006 (67 to 79% GC), Pythia (89% vs 51%, recall 97% vs 48% at 81% sensitivity). Primer3 manual DMSO factors (0.6 default, 0.5 Cullen, 0.675 Escara) match. I did not open the primary DMSO papers.

## Module check (`structure_reference.py`)

- The tests pass (8/8). The demo reproduces `ss3/demo.out` except for NUPACK sampling noise in the third decimal. f_bound recomputed by hand: dG_eff -12.34 at 65 C and 200 nM gives 0.950. Correct.
- The dG_eff sign convention is correct: adding -G_ens(primer) = RT ln Q is the cost of the fully open primer.
- Index conventions checked: p_seed3 and p_nuc windows sit on the 5'-most template nucleotides, as documented.
- Problems fixed in `structure_reference_fixed.py` (10/10 tests pass, demo unchanged):
  1. Stale docstring claimed agreement "to about 1 kcal/mol" and monotonic behaviour. It now states the GC-dependent bias.
  2. No `eval_structure` guard existed despite the recommendation. Added `safe_eval_structure` (segfault confirmed as SIGSEGV exit -11; a length mismatch silently returns 100000).
  3. Added `primer_tm()`, which passes formamide_conc=0.0 explicitly.
  4. `self_priming_prob` docstring said "all k paired" but the code returns "any of k paired". Docstring fixed.
  5. `nupack_crosscheck` did not clamp sodium to 1.1 M. Fixed.
  6. `load_params` cache can go stale. Added `force`.
  7. `local_melting_tm` documents that Vienna values above 90 C are unreliable. The demo's "Tm 101 C" hairpins should read "> 90 C".
- Not fixed, design-level: `recommend_ta` uses only the amplicon-strand context, not the genomic context that recommendation 6 asks for. dG_duplex (primer3 thal, SantaLucia) and dG_open (Vienna M1999) also come from different energy models. NUPACK gives perfect 22-mer duplexes 2.1 to 2.3 kcal/mol more stable than thal at 65 C in PCR salt (`dg.out`). That shifts K by about 25-fold and moves the recommended Ta.

## Recommendations: status

- Use M1999 with noGU=1, dangles=2 and salt=Na_eq: acceptable for GC-rich genomes such as Mtb, but energies are GC-biased (#2). For 40 to 55% GC genomes, use NUPACK as the reference where the licence allows, or calibrate.
- Never use M2004: supported.
- Na_eq = mono + 3.3 sqrt(Mg - dNTP): matches NUPACK. The dNTP subtraction is a primer3 or von Ahsen convention, not NUPACK's.
- probs_window over pfl_fold_up: supported (with nuance #6).
- Thresholds (1.5 / 3.0 kcal/mol, P_nuc 0.1, f 0.9, self-priming 0.2, hairpin -3): unverifiable. No primary source; the report says so itself.
- Additive advice: consistent with Aird, Musso and Frey.

## Missed points

1. primer3-py 2.0.3 defaults to 0.8 M formamide in `calc_tm`, which lowers all primer and product Tm values by about 2 C unless `formamide_conc=0` is passed. This affects the whole PrimerJinn Tm pipeline, not only this report.
2. Vienna-NUPACK agreement depends on GC content. The user's non-Mtb genomes (target Tm 55 to 60 C) are exactly where M1999 deviates most.
3. Gimpel's mechanism is adapter-mediated. In tailed amplicon-seq primers, also scan each amplicon for internal matches to the reverse complement of the universal tails, not only P(3' nt paired).
4. Mixing energy models (thal duplex plus Vienna folding) makes dG_eff internally inconsistent by about 2 kcal/mol.
5. Gao's k_on was measured at 20 C in 0.5 M NaCl. Rates at 60 to 72 C in about 0.16 M Na_eq are not established.
6. The report itself says Vienna is unreliable above 90 C, yet it quotes Vienna hairpin Tm values of 96 to 101 C.
7. Product reannealing and primer-primer dimers are not in f_bound. Pool-level dimers belong in a separate model.
