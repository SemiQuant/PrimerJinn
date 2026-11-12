# Fact-check of tm_models.md

Checker scripts: `research/tmv/xcheck.py`, `research/tmv/xcheck2.py` (run with `python3 -I <script> research/`). Vendor PDFs: `research/tmv/dl/`.

Bottom line: the report holds up. Every coefficient I could reach a primary or independent source for matches digit by digit. `thermo_reference.py` passes its self-tests and agrees with independent code on fresh inputs: Biopython 1.85 within 0.0033 C on 300 new primers in 9 buffers (all four salt regimes), primer3-py within 0.0033 C (Na-eq mode) and 0.0000 C (owczarzy mode with the length term removed), and NUPACK dna04 dG37 within 0.05 kcal/mol. I found no formula bug, so I did not write a `_fixed.py`. The corrections below are small: source statements, labels and caveats.

## Claim by claim

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| 1 | Owczarzy 2008 Mg equation, a-g, a/d/g(Mon) for 0.22 <= R < 6, R < 0.22 uses 2004 | confirmed | Patent OCR `dl_papers/us.txt`: Table IIIa (3.92e-5, -9.11e-6, 6.26e-5, 1.42e-5, -4.82e-4, 5.25e-4, 8.31e-5), Eq 18-20 written as (1-0.157 ...), (1+0.279 ...), (1-0.514 ...). The same values appear in primer3 oligotm.c lines 513-526 and Biopython method 7. The worked example gives 73.2 C in the patent, 73.15 C in Biopython (my run) and 73.2 C in the module. One labelling caution: in the patent the 1/Tm form is **Eq 13** with Table IIIa. Patent "Eq 16" is a different fit with Table IIIb coefficients (a' = -4.59 and so on). Never take coefficients from Table IIIb. primer3 calls the 1/Tm form "Eq 16" of the 2008 paper. I could not open the paper to check its numbering. |
| 2 | Owczarzy 2004 monovalent correction, 1.6 C average error | confirmed | The PubMed abstract (PMID 15035624) states 92 duplexes, 69 mM to 1.02 M and 1.6 C. The coefficients match primer3 lines 508-510 and Biopython method 6. |
| 3 | SL2004 differs from SL1998 only in AA/TT and initiation; within 0.14 C | confirmed | **SL1998 Table 2 read from PMC19045 (the PMC page is reachable; the report says it was blocked)**: all 10 stacks, init G-C 0.1/-2.8, init A-T 2.3/4.1 and symmetry 0/-1.4 match digit by digit. SL2004 values come from Biopython DNA_NN4 only, because the Annual Review is paywalled. An independent page (bioinf.fisica.ufmg.br comparetm, parameter list) says the 2004 set differs from 1998 "mainly in initiation and AA/TT". The NUPACK dna04 dG37 table matches the SL2004 dG37 (AA/TT -1.00). On 300 fresh primers, SL1998 minus SL2004 Tm was -0.14 to +0.08 C. |
| 4 | primer3 owczarzy drops the length term (integer division) | confirmed | oligotm.c line 323 declares `int len`, and line 530 computes `(1/(2 * (len - 1)))`, which is 0. My reimplementation without the term matches primer3 to 0.0000 C on 300 fresh primers. Effect at 50 mM K+, 2 mM Mg, 0.8 mM dNTP: dropping the term raises Tm by 0.1 to 0.9 C. |
| 5 | primer3-py 2.0.3 formamide default 0.8 M; about 2.1 C lower | confirmed | `inspect.signature` and `argdefaults.Primer3PyArguments.formamide_conc = 0.8`. ThermoAnalysis also uses 0.8. The primer3 2.6.1 manual gives `PRIMER_FORMAMIDE_CONC (float; default 0.0)`. Re-run: 58.83 vs 60.95 C (santalucia) and 58.91 vs 61.03 C (owczarzy). |
| 6 | primer3 santalucia and thal use Na_eq = Mon + 120 sqrt(Mg - dNTP); dS + 0.368(N-1)ln; C_eff = dna_conc/4 | confirmed, with one detail | oligotm.c lines 448-449 and 775-784; thal.c `saltCorrectS`. The /4 applies only to non-self-complementary oligos. Symmetric oligos use dna_conc/1 (Equation A). The 0.368 coefficient matches SL1998 Eq 7 (-0.114 kcal/mol per N ln[Na] divided by 0.310 = 0.368). Agreement with primer3 was 0.0033 C at most. |
| 7 | 342 measured Tm: MAE 1.25 (ref), about 2.8 (Na-eq), 2.6 (primer3 owczarzy) | confirmed (re-run) | I re-ran `patent_tableI.py` and got identical output. I spot-checked OCR rows against the patent header ("2 uM DNA, 2 mM / 10 mM Tris-HCl, no KCl") and against values (TACTTCCAGTGCTCAGCGTA at 1.5 mM is 337.8 K = 64.65 C). Caveats from the report still apply: these data were used to fit the salt model (in-sample), and 61 of 448 values were lost to OCR. |
| 8 | Table VIII: 0.6-0.8 C vs von Ahsen 1.5-2.9 C | confirmed | Patent OCR Table VIII: 0.6 (n = 484) and 0.8 (n = 456) with a, d, g varying; von Ahsen 1.5 and 2.9; Peyret 1.2 and 2.6; Mitsuhashi 2.1 and 3.8. |
| 9 | von Ahsen: 0.75 C per % DMSO, SE 1.76 C, n = 217 | confirmed | PubMed abstract PMID 11673362 (also SE 2.06 C for n = 475 and 1.87 C for the empirical formula). |
| 10 | Free Mg: patent says Mg - dNTP; Biopython Ka 3e4; 0.03 C difference | confirmed except Ka source | Patent Example 7 text: 3.0 mM Mg with 0.8 mM dNTP gives 2.2 mM, citing Sigel 1993. Biopython code uses `ka = 3e4`; its comment wrongly calls this a dissociation constant. On 300 fresh primers the difference was -0.05 to -0.02 C. **Whether Ka = 3e4 M^-1 appears in the 2008 paper is unverifiable**: the ACS PDF returned an HTML block page. |
| 11 | NEB Q5 API = SL1998 + Owczarzy 2004 at 150 mM, C_eff = primer conc | confirmed | 4 fresh API calls (3 new sequences, 2 concentrations) compared with Biopython `Tm_NN(DNA_NN3, Na=150, saltcorr=6, dnac1=C, dnac2=0)`: differences < 0.001 C. The API note "Ta = Tm_lower+1 C" is present. |
| 12 | R about 0.70 at 50 mM K+, 2 mM Mg, 0.8 mM dNTP; jumps < 1 C for Mon 10-100 mM | confirmed | sqrt(1.221e-3)/0.05 = 0.699. Independent check with Biopython (not the module), 20-mer: jumps from -0.78 to +0.72 C for Mon 10-100 mM; +15.56 C (R = 0.22, 2 mM), +3.39 C (R = 0.22, 5 mM), -3.44 C (R = 6, 2 mM). The abstract says PCR conditions (1.5-5 mM Mg, 20-100 mM monovalent) fall in the Mg-dominant range. |
| 13 | Vendor 1x Mg: Taq 1.5, Q5 2.0, Phusion HF 1.5, KAPA HiFi 2.0, SuperFi II 1.75 | confirmed (partly from primary PDFs) | I opened the PDFs: Thermo MAN0018859 says "1.75 mM MgCl2 in 1X", and the 5X buffer "includes 8.75 mM". MAN0018860 (master mix) also gives 1.75 mM. KAPA HiFi HotStart TDS: "Both buffers contain MgCl2 at a 1X concentration of 2 mM". NEB Taq (10 mM Tris-HCl, 50 mM KCl, 1.5 mM MgCl2, pH 8.3), Q5 2.0 mM and Phusion 1.5 mM come from search summaries of NEB pages only, because neb.com returned 403. An application note quotes 1.5 mM for SuperFi; the user guide (1.75 mM) is the authority. |
| 14 | Betaine isostabilizing about 5.2 M; no per-M coefficient | confirmed (first part); second part is an absence claim | PubMed abstract PMID 8418834: about 5.2 M. I also found no validated linear coefficient, but did not search exhaustively. |

## Other numbers checked

- thal parameters: `primer3_config/stack.dh` and `stack.ds`, entry 16 (AA/TT) is -7900 and -22.2, which are the 1998 values. thal.c has `AT_H = 2200`, `AT_S = 6.9`, `dplx_init 200 / -5.7`. The hybrid nature is confirmed. thal `calc_heterodimer(seq, revcomp)` Tm minus reference: +0.27 +/- 0.51 C (range -0.82 to +1.34), which is the Na-eq signature.
- primer3 with all defaults minus reference: -5.77 +/- 0.63 C (report: -5.73 +/- 0.66). The GC slope of primer3 santalucia (dna 800) is -3.99 C per unit fGC (report -4.14). The fresh sample is consistent.
- Formamide: 1 M is 45.04 g/L / 1.13 g/mL = 39.9 mL/L = 4.0 % v/v, as stated. The Blake and Delcourt form matches primer3 and Biopython (fmdmethod 2). The primary paper was not read (unverifiable).
- DMSO factors 0.6 (Musielski), 0.5 (Cullen), 0.675 (Escara) and 0.75 (von Ahsen) appear verbatim in primer3 manual 2.6.1. The Biopython docstring says the default is 0.65 while the code default is 0.75 (confirmed).
- Biopython Tm_NN uses R = 1.987 and the module uses 1.98720. This alone causes a constant offset of about 0.002 to 0.003 C, so unit-test tolerances must be at least 0.005 C.

## Module review (`thermo_reference.py`)

The self-tests pass. I found no formula error. Small robustness issues, none affecting the default PCR path, so I did not write a `_fixed.py`:
- `analyze(..., salt="owczarzy2004")` or `"santalucia1998"` with zero monovalent and Mg <= dNTP raises a math domain error (log 0). Raise a clear ValueError instead.
- `fraction_bound` for a self-complementary sequence ignores `target_nM`.
- The docstring calls the formula "Owczarzy 2008 Eq 16". Add "patent Eq 13, Table IIIa" so nobody uses Table IIIb.

## Recommendations: verdicts

All are supported, with these refinements:
- The Biopython agreement test should use a tolerance of 0.005 to 0.01 C, not exact equality, because of the R constant.
- "Warn when [Mon+] < 10 mM" will fire on harmless Mg-dominant buffers (for example the patent's 1 mM Tris+ data, where R >> 6). Warn only when Mon < 10 mM **and** R is within about 30 % of 0.22 or 6.
- Choosing Ta from theta >= 0.9 to 0.95 is a reasonable heuristic. It is a two-state equilibrium value: it ignores primer hairpins, template structure and annealing kinetics (see missed points).

## Missed points

1. **Template overhang stacking.** A primer on a long template is not an isolated duplex. The template continues past both primer ends, so the ends gain dangling or coaxial-stacking stabilisation. None of the Tm tools above model this. The direction is known (more stable), but I did not verify a size.
2. **Multistate competition.** Effective binding at Ta is reduced by primer hairpins, primer-dimers and template secondary structure. Primer-level theta should come from a coupled equilibrium (NUPACK complex concentrations), not from two-state theta alone. This links to the user's requirement that amplicon structure inform Ta.
3. **Kinetics.** Annealing steps of 10 to 30 s may not reach equilibrium, especially for structured templates. Equilibrium theta overstates binding.
4. **NUPACK convention.** `nupack.structure_energy` for a two-strand complex includes the -RT ln(55.14 M) water-standard-state term (-2.47 kcal/mol at 37 C). After removing it, NUPACK dna04 matches SL2004 dG37 (+0.03 +/- 0.05 kcal/mol). Anyone mixing NUPACK and NN dG must correct for this, or the comparison is off by about 2.4 kcal/mol.
5. **Heat capacity.** dH and dS are assumed temperature independent. This is fine near Tm but adds error when dG is evaluated at Ta 5 to 10 C below Tm. The size is unverified.
6. **Ready-mix Mg differs from kit buffers.** A search summary says KAPA HiFi HotStart ReadyMix contains 2.5 mM MgCl2 at 1x, versus 2 mM in the kit buffers. I did not open its data sheet, so this is unverified. Amplicon-seq protocols often use the ReadyMix, so the CLI should take Mg explicitly rather than infer it from the enzyme name.
7. **Mismatches.** The perfect-match Tm ignores lineage SNPs under primers (relevant for M. tuberculosis lineages). Mismatch NN parameters (Allawi and SantaLucia, Peyret) or thal/NUPACK handle this.
8. **Free Mg in large multiplex pools.** High total primer and template concentrations bind some Mg. The report's free Mg model counts only dNTP. This is not verified to matter at about 50 to 100 primers x 50 to 200 nM, but it is worth a sensitivity check.
