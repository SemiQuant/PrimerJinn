# Oligo duplex thermodynamics for PCR primers: models, tools, recommendation

All numbers below come from code I ran (Python 3.13.2, primer3-py 2.0.3, Biopython 1.85) or from primary text I read. Where I could not reach a primary source I say so.

Files (all under `research/`):

| file | what |
|---|---|
| `thermo_reference.py` | reference implementation, self-tests pass (`python -I thermo_reference.py`) |
| `tm3/scripts/compare200.py`, `tm3/compare200.out` | 200-primer comparison of implementations |
| `tm3/scripts/neb_fetch.py`, `tm3/neb/q5_tm.json` | 40 NEB Tm API calls (Q5) |
| `tm3/scripts/neb_fit.py`, `tm3/neb/fit.out` | fit of NEB Q5 to SantaLucia + salt models |
| `tm3/scripts/patent_tableI.py`, `tm3/patent_tableI.out` | accuracy vs 342 measured Tm values (US 8,055,451 Table I) |
| `tm3/p3v203/oligotm.c`, `thal.c` | primer3-py v2.0.3 C source from GitHub tag (md5 of oligotm.c identical to earlier cached copies) |
| `dl_papers/us.txt` | OCR text of US 8,055,451 B2 (Owczarzy et al., IDT), reused from earlier attempt |

## 1. Nearest-neighbour parameter sets

### 1.1 Table (1 M NaCl; dH kcal/mol, dS cal/(K mol))

| stack (5'->3'/3'->5') | SantaLucia 1998 unified dH | dS | SantaLucia and Hicks 2004 dH | dS |
|---|---|---|---|---|
| AA/TT | -7.9 | -22.2 | **-7.6** | **-21.3** |
| AT/TA | -7.2 | -20.4 | -7.2 | -20.4 |
| TA/AT | -7.2 | -21.3 | -7.2 | -21.3 |
| CA/GT | -8.5 | -22.7 | -8.5 | -22.7 |
| GT/CA | -8.4 | -22.4 | -8.4 | -22.4 |
| CT/GA | -7.8 | -21.0 | -7.8 | -21.0 |
| GA/CT | -8.2 | -22.2 | -8.2 | -22.2 |
| CG/GC | -10.6 | -27.2 | -10.6 | -27.2 |
| GC/CG | -9.8 | -24.4 | -9.8 | -24.4 |
| GG/CC | -8.0 | -19.9 | -8.0 | -19.9 |
| initiation | with terminal G-C: +0.1, -2.8; with terminal A-T: +2.3, +4.1 (one term per duplex end) | | +0.2, -5.7 (once) | |
| terminal A-T penalty | (folded into initiation above) | | +2.2, +6.9 per terminal A-T | |
| symmetry (self-complementary) | 0, -1.4 | | 0, -1.4 | |

How the two sets differ: only AA/TT changes (-7.9/-22.2 to -7.6/-21.3), and initiation is rewritten as one initiation term plus a per-end A-T penalty. Per-duplex sums are nearly identical: two G-C ends give 0.2/-5.6 (1998) vs 0.2/-5.7 (2004); two A-T ends 4.6/8.2 vs 4.6/8.1. On 200 random primers the two sets give the same Tm within -0.14 to +0.08 C (mean 0.00, SD 0.05).

Verification. I could not download the primary PDFs: PNAS / PMC19045 returned a Cloudflare challenge and the Annual Review (doi:10.1146/annurev.biophys.32.110601.141800) is paywalled. The values above were cross-checked in three independent code bases that cite those papers: primer3 `oligotm.c` (SL1998 table, read line by line, 1998 init and symmetry), Biopython 1.85 `DNA_NN3` (SL1998) and `DNA_NN4` (SantaLucia and Hicks 2004), and primer3 `primer3_config/stack.dh/.ds` read programmatically. All agree exactly.

What primer3 uses: `oligotm` (calc_tm) = SL1998 table with SL1998 initiation. `thal` (calc_hairpin/homodimer/heterodimer/end_stability) = 1998 stack values (AA -7.9/-22.2 in `stack.dh/ds`) combined with 2004-style initiation (`dplx_init_H = 200`, `dplx_init_S = -5.7`) and AT penalty (`AT_H = 2200`, `AT_S = 6.9`) in `thal.c` lines 125-131, 583-593. So thal is a hybrid of the two sets plus mismatch, dangling end and loop parameters.

### 1.2 Tm formula

Tm = dH x 1000 / (dS + R ln C_eff), R = 1.9872 cal/(K mol).

C_eff, the concentration term:
* self-complementary: C_T (total strand);
* non-self-complementary, equal strands: C_T/4;
* primer A in excess over target B: A0 - B0/2, which tends to A0 in PCR (primer 100-500 nM, template far lower). This is the physically correct choice for a PCR primer.

Conventions differ between tools, which shifts Tm by about 1 to 2 C for 200 nM primers:
* primer3: C_eff = dna_conc/4 (code: `log(DNA_nM/4000000000.0)`).
* Biopython Tm_NN: C_eff = dnac1 - dnac2/2 (default 25/25 nM gives 12.5 nM).
* NEB Tm API: C_eff = primer concentration (fit below).
* thermo_reference.py: C_eff = primer - target/2 (target default 0).

## 2. Salt corrections

### 2.1 SantaLucia 1998 entropy correction
dS(Na) = dS(1 M) + 0.368 x (N - 1) x ln[Na+], N = length, [Na+] in M. dH unchanged. (Implemented identically in primer3 `oligotm.c` line 449 and Biopython saltcorr=5.)

### 2.2 Owczarzy 2004 (monovalent), Biochemistry 43:3537-3554, doi:10.1021/bi034621r, PMID 15035624
1/Tm(Na) = 1/Tm(1 M) + (4.29 fGC - 3.95) x 10^-5 x ln[Na+] + 9.40 x 10^-6 x (ln[Na+])^2
Tm in K, [Na+] in M, fGC = fraction G-C pairs. Abstract (PubMed): 92 duplexes, 69 mM to 1.02 M Na+, average error 1.6 C on an independent literature validation set; correction independent of length and DNA concentration. Coefficients identical in primer3 (lines 508-510) and Biopython (saltcorr=6).

### 2.3 Owczarzy 2008 (Mg2+ and monovalent), Biochemistry 47:5336-5353, doi:10.1021/bi702363u, PMID 18422348
The paper is paywalled. Coefficients below were read from US patent 8,055,451 B2 by the same authors (Table IIIa, Eq 18-20, Example 6), and match primer3 `oligotm.c` lines 513-531 and Biopython `salt_correction(method=7)` exactly.

1/Tm(Mg) = 1/Tm(1 M Na+) + a + b ln[Mg2+] + fGC (c + d ln[Mg2+]) + [1 / (2 (Nbp - 1))] x (e + f ln[Mg2+] + g (ln[Mg2+])^2)

| coef | value (K^-1) | patent SE |
|---|---|---|
| a | 3.92 x 10^-5 | 0.2 x 10^-5 |
| b | -9.11 x 10^-6 | 0.5 x 10^-6 |
| c | 6.26 x 10^-5 | 0.4 x 10^-5 |
| d | 1.42 x 10^-5 | 0.08 x 10^-5 |
| e | -4.82 x 10^-4 | 0.7 x 10^-4 |
| f | 5.25 x 10^-4 | 0.2 x 10^-4 |
| g | 8.31 x 10^-5 | 0.2 x 10^-5 |

Regime switch with R = sqrt([Mg2+]free) / [Mon+] (both in M):
* R < 0.22: monovalent dominates, use Owczarzy 2004 (2.2) with [Mon+], ignore Mg.
* 0.22 <= R < 6.0: use the Mg equation with a, d, g replaced by
  * a = 3.92 x 10^-5 x (0.843 - 0.352 sqrt[Mon+] ln[Mon+])
  * d = 1.42 x 10^-5 x (1.279 - 4.03 x 10^-3 ln[Mon+] - 8.03 x 10^-3 (ln[Mon+])^2)
  * g = 8.31 x 10^-5 x (0.486 - 0.258 ln[Mon+] + 5.25 x 10^-3 (ln[Mon+])^3)
  (patent writes 0.843 as "1 - 0.157", 1.279 as "1 + 0.279", 0.486 as "1 - 0.514")
* R >= 6.0: Mg equation with the constants in the table.

Bookkeeping (patent Example 7 and abstract): [Mon+] = [Na+] + [K+] + [Tris+]; at pH 8.3, "10 mM Tris is equivalent to about 5 mM Na+". Na+ and K+ are equivalent. Free Mg: the patent says the Mg-dNTP binding constant is large enough that free Mg can be approximated as total Mg - total dNTP (3.0 mM Mg + 0.8 mM dNTP gives 2.2 mM). Biopython instead solves 1:1 binding with Ka = 3 x 10^4 M^-1 (2 mM Mg + 0.8 mM dNTP gives 1.221 mM free; 3.0/0.8 gives 2.212). I could not verify the 3 x 10^4 value against the 2008 paper itself. The two choices differ by only 0.03 C on 200 primers (mean -0.03, SD 0.01), so this does not matter in practice.

Typical PCR buffer (50 mM K+, 2 mM Mg, 0.8 mM dNTP): R = sqrt(0.00122)/0.05 = 0.70, so the mixed regime (a, d, g vary) applies. PubMed abstract: PCR and sequencing conditions (1.5-5 mM Mg, 20-100 mM monovalent) fall in the Mg-dominated range.

Reported accuracy (patent Table VIII, salt scaling from measured 1 M Tm): with a, d, g varying, mean |dTm| 0.6 C on 484 measurements in 50 mM KCl + Mg, 0.8 C on 456 with both K and Mg varying. Prior methods on the same data: Peyret 1.2 / 2.6 C, von Ahsen 1.5 / 2.9 C, Mitsuhashi 2.1 / 3.8 C.

The model is piecewise and not continuous. My measurement for a 20-mer: jumps at the boundaries are < 1 C for [Mon+] 10-100 mM, but 3.4 C (R = 0.22) and 15.6 C (R = 0.22, 2 mM Mon) at [Mon+] <= 5 mM, i.e. outside PCR use.

### 2.4 von Ahsen 2001 Na-equivalent, Clin Chem 47:1956-1961, PMID 11673362
[Na+eq] (mM) = [Mon+] + 120 sqrt([Mg2+] - [dNTP]), then apply SantaLucia 1998 entropy correction. With 0.75 C per % DMSO, SE 1.76 C for 217 matched duplexes (n = 475 incl. mismatches: 2.06 C) (PubMed abstract). primer3 manual (release 2.6.1) states it uses this formula for divalent cations and that dNTP > Mg removes the Mg effect.

### 2.5 Which is most accurate for PCR buffers
Owczarzy 2008 with the decision tree. Evidence:
* Patent Table VIII (above): 0.6-0.8 C vs 1.5-2.9 C for von Ahsen in K+ plus Mg buffers.
* My own test on 342 measured Tm values (patent Table I, 15-25-mers, 0.5-125 mM Mg, 1-5 mM Tris+, 2 uM total DNA; full stack NN prediction plus salt model; note these are the data Owczarzy 2008 was fitted on, so the salt part is in-sample):

| model | mean error (pred - meas), C | SD | MAE |
|---|---|---|---|
| SL2004 + Owczarzy 2008 (reference) | +0.50 | 1.51 | 1.25 |
| SL1998 + Owczarzy 2008 | +0.55 | 1.50 | 1.26 |
| SL2004 + von Ahsen Na-eq | +1.62 | 3.35 | 2.78 |
| primer3 `santalucia` | +1.69 | 3.32 | 2.77 |
| primer3 `owczarzy` (with length-term bug, see 4.3) | -1.97 | 2.35 | 2.61 |

Na-eq error grows with Mg (-1.1 C at 1.5 mM Mg, +6.8 C at 125 mM). primer3 `owczarzy` error grows the other way (-0.75 C at 1.5 mM to -3.6 C at 50 mM). Residual SD of about 1.5 C for the reference is dominated by NN error at 1 M Na (e.g. CAGTGGGCTCCTGGGCGTGCTGGTC is over-predicted by about 3 C).

## 3. Additives

| additive | correction | source | verified how | uncertainty |
|---|---|---|---|---|
| DMSO | Tm - 0.75 C per % v/v | von Ahsen 2001 (PMID 11673362) | PubMed abstract | fitted within an NN model with SE 1.76 C overall; no separate SE for the DMSO term in the abstract |
| DMSO | 0.6 C per % (primer3 default `dmso_fact`), alternatives 0.5 (Cullen 1976), 0.675 (Escara 1980) | Musielski 1981 via primer3 manual | primer3 manual text only; primary papers not read | published factors span 0.5-0.75 C per %, so 5% DMSO gives 2.5-3.75 C |
| formamide | dTm = (0.453 fGC - 2.88) x [formamide in mol/L] | Blake and Delcourt 1996 (NAR 24:2095) | primer3 `oligotm.c` line 441 and Biopython code; primary paper not read | 1 M formamide is about 4.0% v/v (45.04 g/mol, 1.13 g/mL); 0.6-0.72 C per % quoted by Biopython (McConaughy 1969) |
| betaine | no validated per-M linear coefficient found | Rees et al. 1993 (PMID 8418834): isostabilizing (AT = GC) at about 5.2 M | PubMed abstract | betaine reduces GC dependence, so a single linear shift is not correct; calibrate empirically if used |

Notes: Biopython `chem_correction` docstring says the DMSO default is 0.65, but the code default is `DMSOfactor=0.75` (doctest 70 - 3 x 0.75 = 67.75, run). primer3 manual mislabels the formamide tag as "concentration od DMSO in mol/l".

## 4. primer3-py 2.0.3: what it actually implements

### 4.1 Signatures (run with `inspect.signature`)
```
calc_tm(seq, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, dmso_conc=0.0, dmso_fact=0.6,
        formamide_conc=0.8, annealing_temp_c=-10.0, max_nn_length=60, tm_method='santalucia',
        salt_corrections_method='santalucia') -> float
calc_hairpin(seq, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0, max_loop=30, output_structure=False)
calc_homodimer(... same ...)
calc_heterodimer(seq1, seq2, ... same ...) -> ThermoResult
calc_end_stability(seq1, seq2, mv_conc, dv_conc, dntp_conc, dna_conc, temp_c, max_loop) -> ThermoResult
ThermoAnalysis(): attributes mv_conc, dv_conc, dntp_conc, dna_conc, temp_c, max_loop, tm_method,
  salt_correction_method, dmso_conc, dmso_fact, formamide_conc, annealing_temp_c; same defaults
```
Units: mv/dv/dntp mM, dna nM, dmso %, formamide mol/L, temps C. ThermoResult has tm (C), dg, dh (cal/mol), ds (cal/(K mol)).

### 4.2 Options mapped to models
| option | model |
|---|---|
| `tm_method='santalucia'` | SL1998 unified NN, SL1998 initiation and symmetry |
| `tm_method='breslauer'` | Breslauer 1986 |
| `salt_corrections_method='schildkraut'` | Tm + 16.6 log10[Na_eq] (Schildkraut and Lifson 1965), Na_eq = von Ahsen |
| `'santalucia'` (default) | von Ahsen Na_eq (Mon + 120 sqrt(Mg - dNTP)) then SL1998 entropy 0.368 (N-1) ln |
| `'owczarzy'` | Owczarzy 2004 for R < 0.22, Owczarzy 2008 Eq 16 for R >= 0.22 with a, d, g varying for R < 6; free Mg = Mg - dNTP; dNTP >= Mg sets Mg to 1e-11 M |
| `dmso_conc`, `dmso_fact` | Tm - dmso_fact x dmso_conc (all salt methods) |
| `formamide_conc` | Tm + (0.453 fGC - 2.88) x formamide (Blake and Delcourt) |
| `annealing_temp_c` | computes `bound` (percent) only in schildkraut / santalucia branches; `calc_tm` returns only Tm, so `bound` is not exposed by `calc_tm`; ignored in the owczarzy branch |
| thal functions | always von Ahsen Na_eq with SL1998 entropy correction (`saltCorrectS`, thal.c line 1039-1042); no Owczarzy option, no DMSO or formamide |

So the user's belief is partly wrong: primer3 CAN apply Owczarzy 2008 Mg correction for `calc_tm` (`salt_corrections_method='owczarzy'`), but (a) not for hairpins or dimers, and (b) its implementation has a bug.

### 4.3 Bug: length term dropped in primer3 Owczarzy 2008
`oligotm.c` line 530: `+ (1/(2 * (len - 1))) * (e + (f * log(free_divalent)) + g * ...)`. `len` is `int`, so `1/(2*(len-1))` is integer division and equals 0 for every oligo of 2 nt or more. The e, f, g length terms are never applied. Confirmed numerically: my Python re-implementation with that term removed reproduces primer3 `owczarzy` to 0.000 C (SD 0.000) on 200 primers. Effect at 50 mM K+, 2 mM Mg, 0.8 mM dNTP: about +0.1 to +0.7 C after correcting for the concentration convention; at higher Mg (10-50 mM) primer3 under-predicts by 3-3.6 C (Table I test).

### 4.4 Bug or trap: formamide default 0.8 M
primer3-py 2.0.3 `argdefaults.Primer3PyArguments.formamide_conc = 0.8`, while primer3 core manual says `PRIMER_FORMAMIDE_CONC` default 0.0. Measured: TACTTCCAGTGCTCAGCGTA (50 mM K+, 2 mM Mg, 0.8 mM dNTP, 200 nM) gives 58.83 C with the default and 60.95 C with `formamide_conc=0` (santalucia); 58.91 vs 61.03 (owczarzy). Over 200 primers `calc_tm` with all defaults is 5.73 C (SD 0.66) below the reference at our buffer (formamide plus lower Mg, lower DNA, C/4). Always pass `formamide_conc=0.0`.

## 5. Biopython 1.85 Bio.SeqUtils.MeltingTemp
* `DNA_NN3` = SantaLucia 1998 (table above); `DNA_NN4` = SantaLucia and Hicks 2004. Self-complementarity must be flagged by the caller (`selfcomp=True`).
* `Tm_NN`: R = 1.987; C term k = (dnac1 - dnac2/2) nM; selfcomp uses dnac1.
* `salt_correction(method)`: 1 = 16.6 log10[Na] (Schildkraut 1965); 2 = 16.6 log10([Na]/(1 + 0.7[Na])) (Wetmur 1991); 3 = 12.5 log10[Na] (SantaLucia 1996); 4 = 11.7 log10[Na] (SantaLucia 1998); 5 = dS + 0.368 (N-1) ln[Na] (SantaLucia 1998); 6 = Owczarzy 2004; 7 = Owczarzy 2008 with decision tree, Mon = Na + K + Tris/2, free Mg from Ka = 3 x 10^4 M^-1 quadratic (comment wrongly calls it a dissociation constant). Methods 1-6 use von Ahsen Na_eq when K, Tris, Mg or dNTP is non-zero; method 7 does not.
* method 7 keeps the length term (float division). Owczarzy worked example (83.4 C at 1 M, 0.5 mM Mg, 10 mM Tris): Biopython 73.15 C; published 73.2 C (measured 73.3 C).
* `chem_correction`: DMSO linear (default 0.75 C/%), formamide method 1 linear per % (0.65 default) or method 2 Blake and Delcourt per M.

## 6. Numerical comparison, 200 random primers
18-30 nt, GC 0.30-0.75, seed 20261111; buffer 50 mM K+, 2 mM Mg, 0.8 mM dNTP, 200 nM primer, template negligible. Reference = thermo_reference.py (SL2004 + Owczarzy 2008, Ka free Mg, C_eff = 200 nM). Reference Tm range 50.4-81.5 C, mean 68.46 C. Differences are method minus reference (C), output of `tm3/compare200.out`:

| method | mean | SD | min | max | Spearman |
|---|---|---|---|---|---|
| mine, free Mg by subtraction | -0.03 | 0.01 | -0.04 | -0.02 | 1.0000 |
| mine, SL1998 NN | 0.00 | 0.05 | -0.14 | 0.08 | 0.9999 |
| mine, SL2004 + von Ahsen Na_eq | +0.20 | 0.53 | -0.93 | 1.24 | 0.9985 |
| Biopython DNA_NN4 saltcorr=7, dnac1=200, dnac2=0 | 0.00 | 0.00 | 0.00 | 0.00 | 1.0000 |
| Biopython DNA_NN3 saltcorr=7 | 0.00 | 0.05 | -0.14 | 0.08 | 0.9999 |
| Biopython DNA_NN4 saltcorr=7, dnac1=dnac2=200 | -0.84 | 0.12 | -1.12 | -0.67 | 0.9997 |
| primer3 santalucia, dna_conc=200 (C/4 = 50 nM) | -1.44 | 0.57 | -2.83 | -0.18 | 0.9963 |
| primer3 owczarzy, dna_conc=200 | -1.23 | 0.16 | -1.62 | -0.96 | 0.9996 |
| primer3 santalucia, dna_conc=800 (C_eff 200 nM) | +0.24 | 0.53 | -0.87 | 1.37 | 0.9986 |
| primer3 owczarzy, dna_conc=800 | +0.45 | 0.10 | 0.21 | 0.68 | 0.9998 |
| primer3 calc_tm all defaults | -5.73 | 0.66 | -7.42 | -4.20 | 0.9936 |

Implementation checks: my re-implementation of primer3 owczarzy (length term dropped) vs primer3: 0.000 / 0.000 C; my SL1998 + Na_eq + C/4 vs primer3 santalucia: -0.002 / 0.000 C. Biopython DNA_NN4 + saltcorr 7 matches my independent code to < 0.005 C.

GC bias: primer3 santalucia (Na_eq) minus reference has slope -4.14 C per unit fGC (AT-rich primers over-predicted, GC-rich under-predicted, up to about 2.2 C apart across GC 0.3-0.75). This matters for a pool that must have uniform Tm across mixed GC.

## 7. Published test vectors (usable as unit tests)

1. Owczarzy 2008 / patent Example 6 (salt scaling only): CAGTGGGCTCCTGGGCGTGCTGGTC (25 bp, fGC 0.72), measured Tm 83.4 C in 1 M Na+ at 2 uM total strands; in 0.5 mM MgCl2 + 10 mM Tris-HCl ([Mon+] = 5 mM Tris+, R = 4.47): intermediate a = 3.82e-5, d = 1.53e-5, g = 8.91e-5; predicted 346.3 K = 73.2 C; measured 73.3 C. My code 73.2 (within 0.06), Biopython 73.15.
2. Patent Table I (2 uM total DNA, Tris-HCl, MgCl2, no KCl, Tm in K; 0.5-20 mM Mg with 1 mM Tris+, 50 and 125 mM Mg with 5 mM Tris+). Examples (measured C vs my full-stack prediction):

| seq | Mg mM | measured C | predicted C |
|---|---|---|---|
| TACTTCCAGTGCTCAGCGTA | 1.5 | 64.65 | 64.69 |
| TACTTCCAGTGCTCAGCGTA | 3.0 | 65.85 | 65.89 |
| TGATTCTACCTATGTGATTT | 1.5 | 54.55 | 54.72 |
| TGATTCTACCTATGTGATTT | 3.0 | 55.75 | 56.07 |
| ACGACAGGTTTACCA | 1.5 | 54.55 | 55.54 |
| CAGTGGGCTCCTGGGCGTGCTGGTC | 1.5 | 75.35 | 78.56 |

Caution: values were parsed from OCR text; I kept only cleanly parsed numbers (387 of 448). Use as tolerance tests (for example |error| < 2 C for 20-mers, set-level MAE < 1.5 C), not exact equality.
3. PCR-buffer data with 50 mM KCl (patent Table VI/VII, 80 duplexes) exist but the OCR lost the decimals; get them from the Owczarzy 2008 supporting information if a KCl + Mg test set is needed. Not verified here.

## 8. Polymerase buffers and the NEB calculator

| buffer | Mg2+ at 1x | monovalent | source |
|---|---|---|---|
| NEB Standard Taq | 1.5 mM | 50 mM KCl, 10 mM Tris-HCl pH 8.3 | NEB M0273 page and compatibility FAQ (via web search) |
| NEB Q5 | 2.0 mM | not disclosed | NEB protocol / protocols.io (via web search) |
| Phusion HF | 1.5 mM | not disclosed | NEB Phusion master mix protocol (via web search) |
| KAPA HiFi (Fidelity and GC) | 2.0 mM | not disclosed | Roche KAPA HiFi technical data sheets KR0368/KR0369 (via web search) |
| Platinum SuperFi II | 1.75 mM (5x buffer 8.75 mM) | not disclosed | Thermo Fisher MAN0018859 (via web search) |

Monovalent salt and other additives of the proprietary high-fidelity buffers are not public. I did not open the PDFs myself; the Mg values are from search-result summaries of vendor documents and should be checked against the lot's manual.

NEB Q5 calculator (tmapi.neb.com, version 0.19.0). 40 calls: 36 primers at 0.5 uM plus one primer at 0.1, 0.25, 1.0, 2.0 uM (`fmt=long` gives unrounded Tm). Grid fit over NN set, salt model, Mon 5-300 mM, Mg 0-3 mM, dNTP 0/0.8 mM, C_eff in {C, C/2, C/4}:
* best: SantaLucia 1998 NN + Owczarzy 2004 monovalent correction at [Mon+] = 150 mM, C_eff = primer concentration, Mg ignored. RMS 0.002 C, max 0.002 C over all 40 values (the tiny residual is consistent with rounding of R). Owczarzy 2008 fits only when R < 0.22, where it reduces to Owczarzy 2004.
* SL2004 at the same conditions: RMS 0.048 C. Best Na_eq fit: RMS 1.26 C.
* Example: TACTTCCAGTGCTCAGCGTA at 0.5 uM: NEB 65.36 C; the reference model at 50 mM K+, 2 mM Mg, 0.8 mM dNTP, 200 nM gives 62.47 C.
* API note returned for Q5: Ta rule "Ta = Tm_lower + 1 C".
So the NEB Q5 value is a vendor-calibrated number on a 150 mM monovalent scale, not a prediction from the true Q5 buffer composition. It is useful to match the Ta that NEB recommends, not as a physical Tm.

## 9. Recommended formula stack

1. dH, dS at 1 M Na+ from SantaLucia and Hicks 2004 NN (table 1.1), with init 0.2/-5.7, A-T end penalty 2.2/6.9, symmetry -1.4.
2. Tm(1 M) = 1000 dH / (dS + R ln C_eff), C_eff = [primer] - [target]/2 (about [primer] in PCR), R = 1.9872.
3. Cations: [Mon+] = [Na+] + [K+] + 0.5 [Tris]; free Mg from Ka = 3e4 M^-1 (or Mg - dNTP; difference 0.03 C).
4. Salt: Owczarzy 2008 Eq 16 with R decision tree and a, d, g of Mon (section 2.3); Owczarzy 2004 for R < 0.22 or Mg = 0.
5. Additives: DMSO -0.6 C per % (or 0.75 per von Ahsen; make the factor configurable), formamide Blake and Delcourt.
6. Effective dS = 1000 dH / Tm_final - R ln C_eff (dH unchanged); dG(T) = dH - T dS_eff.
7. Fraction bound: K = exp(-dG(T)/RT); primer excess: theta = K A0 / (1 + K A0); finite target: smaller root of x^2 - (A0 + B0 + 1/K) x + A0 B0 = 0, theta = x / B0.

`thermo_reference.py` implements exactly this. Self-tests (all pass): Owczarzy worked example (73.2 C), intermediate a, d, g values, 1998 vs 2004 dG37 agreement, hand-summed self-complementary duplex, free Mg values, theta = 0.5 at Tm for both concentration modes, regime selection, and boundary jumps < 1 C for Mon 10-100 mM. Output for TACTTCCAGTGCTCAGCGTA at the default buffer: dH -153.9 kcal/mol, dS(1 M) -415.4, dS_eff -427.9 cal/(K mol), Tm(1 M) 71.88 C, Tm 62.47 C (mixed regime), dG(60 C) -11.34 kcal/mol, theta(60 C) 0.847.

For hairpins and dimers keep primer3 thal (it has the loop and mismatch parameters) but be aware it uses Na_eq salt handling; or compute those with ViennaRNA/NUPACK DNA parameters. Do not mix Tm scales: rank and filter primers on one scale.

## 10. Open issues
* Owczarzy 2008 Ka = 3e4 M^-1 for Mg-dNTP: verified in the paper (Eq 17, p5346), for dNTP within about 20% of Mg; the default there is simple subtraction (p5346, p5351).
* SantaLucia 1998 / 2004 tables not read in the original PDFs (blocked/paywalled); verified via three code bases.
* Monovalent content of Q5, Phusion, KAPA, SuperFi II buffers unknown; any Tm for them is an assumption. Measure Tm (qPCR melt of a primer/complement pair in the real buffer) to calibrate a single offset if absolute Tm matters.
* Betaine has no validated Tm coefficient for primers.
* Should primer3-py bugs be reported upstream (integer division in oligotm.c Owczarzy branch; formamide default 0.8 in primer3-py)? Not done.

## Sources
* SantaLucia 1998 PNAS 95:1460, doi:10.1073/pnas.95.4.1460 (PMID 9465037)
* SantaLucia and Hicks 2004 Annu Rev Biophys Biomol Struct 33:415, doi:10.1146/annurev.biophys.32.110601.141800 (PMID 15139820)
* Owczarzy et al. 2004 Biochemistry 43:3537, doi:10.1021/bi034621r (PMID 15035624)
* Owczarzy et al. 2008 Biochemistry 47:5336, doi:10.1021/bi702363u (PMID 18422348)
* US 8,055,451 B2 (Owczarzy et al., 2011), local OCR `dl_papers/us.txt`
* von Ahsen et al. 2001 Clin Chem 47:1956 (PMID 11673362)
* Rees et al. 1993 Biochemistry 32:137, doi:10.1021/bi00052a019 (PMID 8418834)
* primer3-py v2.0.3 source: https://github.com/libnano/primer3-py/tree/v2.0.3/primer3/src/libprimer3
* primer3 manual release 2.6.1 (local `eq_downloads/primer3_src/primer3_manual.htm`)
* NEB Tm API: https://tmapi.neb.com (OpenAPI `tm_dl/reference.json`)
* NEB M0273: https://www.neb.com/en/products/m0273-taq-dna-polymerase-with-standard-taq-buffer ; NEB PCR optimization guide: https://www.neb.com/en/tools-and-resources/usage-guidelines/guidelines-for-pcr-optimization-with-thermophilic-dna-polymerases ; Roche KAPA HiFi kits: https://sequencing.roche.com/us/en/products/group/kapa-hifi-kits.html ; Thermo SuperFi II UG: https://assets.thermofisher.com/TFS-Assets/LSG/manuals/MAN0018860_Platinum_SuperFi_II_PCR_MM_UG.pdf
