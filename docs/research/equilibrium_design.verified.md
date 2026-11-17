# Verification of equilibrium_design.md

Verifier work is in `research/verify_eq/` (fresh downloads in `dl/`, scripts in `code/`). Fixed module: `research/eqdesign/structcost_fixed.py` (output `structcost_fixed.out`). Python: /Users/semiquant/micromamba/bin/python3 (primer3-py 2.0.3, Biopython 1.85, NUPACK 4.0.2, ViennaRNA 2.7.0).

Primary sources fetched independently: ThermoPlex paper (Europe PMC full-text XML, PMC12598146) and a fresh clone of github.com/aagmata/ThermoPlex (commit 4e1efc6, 2025-09-09); Pythia (PMC2715258 XML); SADDLE main text (PMC9001684 XML; supplement S9 read from the earlier agent's local copy `eq/saddle/s1.txt`, not re-downloaded); DECIPHER 3.8.1 release tarball from Bioconductor; primer3 `src/oligotm.c` and `primer3_manual.htm` (2.6.1) from primer3-org GitHub main; primer3-py v2.0.3 `oligotm.c` from libnano GitHub tag; NC_000962.3:760900-761500 from NCBI efetch (identical to the file used).

## Verdicts on the key claims

| # | Claim | Verdict | Evidence / correction |
|---|---|---|---|
| 1 | ThermoPlex chi_t formula and psi | confirmed | Paper eq 5 and 6 read. Code `Xt` in spsprimer.m and SimulEq.m: b = 1 + Cr/Ct + exp(G/RT)/Ct, identical to 1 + alpha + (alpha/c0p)exp(dG/RT). Algebra checked: it is the smaller root of theta^2 - theta(1+alpha+1/(K T0)) + alpha = 0. fb.py reproduces K P0/(1+K P0) at T0 = 1 pM. |
| 2 | ThermoPlex defaults and sign typo | confirmed | Table 1: 55 C, 10 mmol/l MgCl2, 0.5 umol/l, psi 1000, chi_t 0.5, length 20, M 3, R 50. Paper eq 2-3 print Z = sum e^(dG/RT), dG = RT ln Z; ThermoDhyb.m lines 197-198 use exp(-E/RT) and -RT ln Z. |
| 3 | SiMulEq reduces to mass action; Newton; 40-70 C, 2 C steps; extension inhibition | confirmed | eqnsetLUM.m: ln(c_i/c0_i) + lambda_i = 0 and dG/RT + ln(c_ij 1e6/(c0_i c0_j)) + lambda_i + lambda_j = 0 (concentrations in uM, hence 1e6) eliminate to c_ij = K c_i c_j. Variables squared, finite-difference Jacobian, stop at max residual <= 1e-10 or 100 iterations (NewtonR.m; its message wrongly says 200). Temperature range is a function argument; 40-70 C in 2 C steps is stated in the paper text. Product correction prod(1 - chi_j) in SimulEq.m lines 118-121. |
| 4 | dG needed for theta at 200 nM, 60 C | confirmed | Recomputed: RT ln(P0(1-theta)/theta) = -10.21, -11.67, -12.51, -13.25 kcal/mol. |
| 5 | Width 0.5 to 0.99 is 6.5-7.3 C; T97 3.0-3.5 C below Tm | corrected (minor) | fb.py T_theta values agree with Biopython Tm_NN (DNA_NN3, saltcorr=7, concentration set to P0(1-theta)/theta) within 0.04 C on the report's primers and two fresh primers. T50-T99 = 6.5, 6.7, 7.3 C confirmed. Tm - T97 is 3.0, 3.0, 3.3 C, so "3.0 to 3.3 C" for these three primers; the report body says 6.5 to 7.5 C where its own data say 6.5 to 7.3 C. |
| 6 | DECIPHER equations and defaults | confirmed | DECIPHER 3.8.1: CalculateEfficiencyPCR.R lines 116-213 (K1, K2a, K2b, K3, K2eff branches, Kov, eff = P Kov/(1+P Kov) x eff_taq); DesignPrimers.R defaults annealingTemp 64, P 4e-7, monovalent 0.07, divalent 0.003, dNTPs 8e-4, minEfficiency 0.8, primerDimer 1e-7, length 17-26, product 75-1200; line 122 ions = monovalent + 3.3 sqrt(divalent - dNTPs) = 0.2248 M; Rd example ions=0.225; hybrid-min / hybrid-ss-min calls; the 99 % / 50 % documentation sentence exists. |
| 7 | SADDLE window and occupancy | confirmed, with a nuance | Main text: truncate until dG between -10.5 and -12.5 kcal/mol; -11.5 best trade-off; 0.25-0.75 GC; Vent (exo-), ThermoPol, 17 cycles 95 C 30 s / 60 C 3 min / 72 C 30 s; median 45 nM. Supplement S9: dS + 0.368 ln(0.18) per NN step, dG_init 2.09, fP1 -11.89 (I recomputed -11.90 from their Table S22). Nuance: Table S22 values are SantaLucia and Hicks 2004 (AA/TT -7.6/-21.3), not SantaLucia 1998, and no terminal AT penalties are applied; Fig. S24's written step list has 17 steps instead of 18 (summing it as printed gives -10.74), a typo in the supplement. Occupancy numbers 0.19/0.51/0.83 at 30 nM recomputed. |
| 8 | Pythia metric, 11 reactions, late stage, SVM thresholds | confirmed | All read in PMC2715258: min of left/right fraction bound, 11 reactions (Fig. 1), "late stage in an idealized PCR", Table 1 thresholds 0.8/0.85/0.9/0.95 with 642/1474/3056/10 498 examples, ~0.7 s per solve, 195 vs 106 pairs, 94/80 % vs 95/82 %. Wording fix for §3.1: the paper says "sensitivities of 81%" (not precision), and the 97 % vs 48 % recall is against Primer3 on RepeatMasked sequence, not against P316. |
| 9 | primer3 BOUND formula, defaults, Owczarzy error | confirmed | oligotm.c (2.6.1) lines 322-345: bound = 100/(1+sqrt(1/((DNA_nM/4e9) K))), 1e9 for symmetric; Owczarzy branch never sets ret.bound. Manual: PRIMER_OPT_BOUND 97.0, MIN -10, MAX 110, ANNEALING_TEMP -10, "95-98% should be bound", equal-concentration wording. primer3-py raises "Primer binding does not work with Owczarzy salt correction". Reimplementation matched primer3-py BOUND on two fresh primers (70.53 %, 91.94 % at 58 C) to 1e-12. |
| 10 | primer3-py 2.0.3 Owczarzy integer division; formamide default 0.8 | confirmed | libnano/primer3-py v2.0.3 oligotm.c line 530 `(1/(2 * (len - 1)))` (integer, = 0); primer3 main casts `((double) len)`. inspect.signature(calc_tm) shows formamide_conc=0.8 (also dv_conc 1.5, dntp 0.6, dna_conc 50). Own NN code agrees with Biopython (saltcorr=7) within 0.05 C on 5 sequences, so the 0.54-0.64 C offset is the primer3-py bug. |
| 11 | Hairpin shortcut K/(1+K_hp); homodimer -12 example | confirmed | Shortcut is exact algebra when T0 << P0; outputs reproduced. Note the -12 kcal/mol homodimer is hypothetical (the primer's real primer3 homodimer dG is -0.46). |
| 12 | Template hairpin costs 1.5/2.4/4.7/5.9; Vienna 1-2 kcal/mol higher | corrected | NUPACK numbers reproduce, but they were run at 2 mM total Mg; with free Mg 1.2 mM (the report's own advice) they are 1.45/2.30/4.50/5.65 at 60 C. The "Vienna 1 to 2 kcal/mol higher" statement is not general: in the report's own data the gap is +1.25 (stem 0), +1.87 (stem 8) and -0.11 (stem 12); on two fresh random constructs (free Mg 1.2 mM, Vienna 0.18 M) Vienna exceeded NUPACK by 0.7 to 3.7 kcal/mol. Costs are also strongly sequence-dependent (stem 8: 1.9 and 4.0 kcal/mol on fresh sites vs 2.4 in the report). Treat Vienna open-site cost as usually higher, not as a calibrated offset. |
| 13 | Mtb rpoB Tm-matched vs dG-matched spreads | confirmed | pool.py rerun gives identical output; sequence matches NCBI NC_000962.3:760900-761500 byte for byte. |
| 14 | NUPACK vs NN differ 1.7 at 55 C, 0.2 at 70 C | corrected | Those numbers compare NUPACK at 2 mM Mg with NN at 1.2 mM free Mg. At matched free Mg 1.2 mM the gaps are 2.01 (55 C), 1.46 (60 C), 0.97 (65 C), 0.54 (70 C) kcal/mol (structcost_fixed.out [C]). NUPACK's 1 M conversion (+RT ln rho_water, rho = 54.575 M at 60 C) was checked against a NUPACK tube_analysis: identical theta to 4 decimals at 60/66/70 C. |

## Other numbers and statements checked

- SantaLucia 1998 NN table and initiation terms in fb.py: identical to Biopython DNA_NN3, digit by digit.
- Owczarzy 2004 eq 22 and 2008 eq 16 constants, the 0.22 and 6.0 ratio branches and the a, d, g rescaling: identical to primer3 2.6.1 oligotm.c and Biopython MeltingTemp.py lines 556-580.
- Entropy form dS(salt) = dS(1 M) + dH corr: algebraically exact at Tm. Applying it at other temperatures is a modelling assumption (Owczarzy fits are on 1/Tm), not something the source establishes; the report calls it "consistent", which is fair, but it is untested away from Tm.
- primer3 "BOUND solves theta/(1-theta) = sqrt(KC/4), exact is theta/(1-theta)^2 = KC/2": algebra checked, both give 50 % at the C/4 Tm.
- Recommendation table, theta = 0.995 limit: -13.72 kcal/mol at 60 C and -13.92 at 65 C recomputed.
- Sensitivity 0.392 kcal/mol per C and 2.33 C per kcal/mol: reproduced.
- ViennaRNA 2.7.0: `RNA.params_load_DNA_Mathews2004()` exists; default `md.salt` is 1.021 M.

## Errors in the report's recommendations

1. **Hairpin threshold arithmetic is wrong.** §6 says "reject if dG_hp < -1 (costs > 0.3 kcal/mol)". With the report's own formula RT ln(1+exp(-dG_hp/RT)) at 60 C: dG_hp = +1 costs 0.13, 0 costs 0.46, -1 costs 1.13, -2 costs 2.03 kcal/mol. A 0.3 kcal/mol cost corresponds to dG_hp of about +0.37 kcal/mol. So a "-1 kcal/mol" hairpin rule already exceeds the "penalize above 1.0 kcal/mol" rule. Base the rule on the computed cost, not on dG_hp.
2. **Na-equivalent of 0.18 M** for ViennaRNA comes from primer3's 120 sqrt([Mg-dNTP] in mM) rule (0.181 M). DECIPHER's 3.3 sqrt rule gives 0.164 M for the same buffer. Neither was derived for ViennaRNA's salt model, which was fitted for RNA; the report does not say this.
3. **Pass free Mg to NUPACK** is advised but the report's own NUPACK runs used total Mg (2 mM). Use the fixed module.

## Errors in the ThermoPlex description (§2.1)

- "NN parameters SantaLucia 1998": the code comment says 1998, but `dnaenergies.m` uses SantaLucia and Hicks 2004 values (AA/TT -7.6/-21.3, initiation 0.2/-5.7, terminal AT 2.2/6.9; identical to Biopython DNA_NN4). The paper cites both.
- `dnaenergies.m` hard-codes the Owczarzy length term with N = 20 (`./(2.*(20-1))`) regardless of primer length, and uses total Mg (confirmed: no dNTP term). The report mentioned the 50 mM K assumption but not the N = 20 hard-coding.

## Module review

- `eqdesign/fb.py`: no bugs found. Independent check against Biopython on fresh inputs agrees within 0.05 C (Tm) and 0.04 C (T50, T97, T99). No fixed copy needed.
- `eqdesign/pool.py`: reproducible; no bugs found.
- `eqdesign/structcost.py`: bug. `vienna_open_cost` claims DNA Mathews 2004 parameters but never loads them, and sets no salt, so the "Vienna DNA site-opening cost" column in `structcost.out` (4.39 to 13.68) is actually RNA-parameter output. NUPACK gets total Mg. Fixed in `eqdesign/structcost_fixed.py` (loads DNA params, `md.salt` argument default 0.18 M, NUPACK free Mg 1.2 mM); output in `structcost_fixed.out`.

## Unverifiable

- SantaLucia 2007 chapter (Visual OMP) content: paywalled, not read; only the primer3 manual's reference to the "third myth" on page 18 was confirmed.
- Gimpel et al. 2025 and the Sci Adv 2026 paper numbers in §3.7: not rechecked.
- ThermoPlex "template-template interactions excluded" and "only one structure per nucleation site": not traced in code.
- Pythia's production threshold and simulation concentrations: in supplementary data, not obtained (the report says so too).

## Points the report missed

1. Template context: a primer on a long template gains dangling-end and coaxial-stacking terms that a two-state primer/complement NN model omits (NUPACK stem-0 cost is negative, about -0.1 to -0.4 kcal/mol, for this reason). Compute dG_dup with the 3' and 5' template flanks, or calibrate the threshold window against NUPACK with flanks.
2. Late-cycle competition from the complementary product strand (amplicon re-annealing) is the main competitor for the primer site once products reach nM levels; theta = K P0/(1+K P0) ignores it. Pythia evaluates at a late stage for this reason.
3. Tris contributes monovalent cations (Biopython counts Tris/2); typical buffers add 5 to 10 mM.
4. Equilibrium theta is not efficiency: polymerase binding, 3' terminal stability and extension kinetics are not in the model; DECIPHER multiplies by an elongation term for this.
5. The model and the threshold window must use the same thermodynamic set: SADDLE's -11.5 kcal/mol uses SantaLucia and Hicks 2004 at 0.18 M with no terminal penalties, so its window cannot be transplanted into an Owczarzy-corrected 1998 model without recomputation.
6. ThermoPlex's 10 mM MgCl2 simulation setting differs from its own wet-lab reaction composition; its thresholds are tied to that unusual setting.
