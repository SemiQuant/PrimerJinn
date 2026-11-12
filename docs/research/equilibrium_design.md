# Equilibrium (fraction bound) primer design versus Tm matching

Scope: what an equilibrium-based multiplex designer should compute for each primer at the annealing temperature Ta, what the literature tools actually do (exact equations, defaults, thresholds), and concrete targets for a 20-plex at 200 nM per primer. Units: kcal/mol (1 kcal = 4.184 kJ), K or °C, mol/L. R = 1.98720 × 10⁻³ kcal mol⁻¹ K⁻¹.

All scripts and outputs are in `research/eqdesign/` (`fb.py`/`fb.out`, `structcost.py`/`structcost.out`, `structcost_B_dna.out`, `structcost_B_salt.out`, `pool.py`/`pool.out`, `intdiv.out`). Run with `/Users/semiquant/micromamba/bin/python3` (primer3-py 2.0.3, NUPACK 4.0.2.0, ViennaRNA 2.7.0). Earlier attempts' material that was reused after reading and re-testing: ThermoPlex repo clone (`dl/ThermoPlex`, commit 4e1efc6, 2025-09-09), ThermoPlex paper and supplement text, Pythia PDF text, DECIPHER R source, SADDLE XML and supplement, primer3 2.6.1 manual and `oligotm.c`, primer3-py's `oligotm.c` (`dl_primer3/`), and two 2025/2026 papers (`eq2/dl/`). The earlier `eq/code/eq_lib.py` was read; its formulas are correct, but I rewrote a self-contained version (`fb.py`) so results do not depend on unreviewed code.

---

## 1. Core math for one primer at Ta

### 1.1 Binding constant
For primer + site ⇌ duplex, with ΔG°(Ta) the standard-state (1 M) free energy at the reaction's salt:

K = exp(−ΔG°/RT)  (units M⁻¹)

Two-state NN: ΔG°(T) = ΔH° − T·ΔS°, with ΔH°, ΔS° from SantaLucia 1998 unified NN (1 M Na⁺) and the salt correction applied to ΔS° only. With the Owczarzy 2008 Mg²⁺ correction, which is expressed on 1/Tm, the consistent entropy form is

ΔS°(salt) = ΔS°(1 M Na⁺) + ΔH° · corr, corr = a + b·ln[Mg²⁺free] + f_GC·(c + d·ln[Mg²⁺free]) + (e + f·ln[Mg²⁺free] + g·ln²[Mg²⁺free]) / (2(N_bp − 1))

a = 3.92e-5, b = −9.11e-6, c = 6.26e-5, d = 1.42e-5, e = −4.82e-4, f = 5.25e-4, g = 8.31e-5; when 0.22 ≤ √[Mg²⁺free]/[Mon⁺] < 6, a, d, g are rescaled for the monovalent ion (a × (0.843 − 0.352·√[Mon⁺]·ln[Mon⁺]), d × (1.279 − 4.03e-3·ln[Mon⁺] − 8.03e-3·ln²[Mon⁺]), g × (0.486 − 0.258·ln[Mon⁺] + 5.25e-3·ln³[Mon⁺])); below a ratio of 0.22 the monovalent-only Owczarzy 2004 eq 22 is used. [Mg²⁺free] = [Mg²⁺] − [dNTP]. Source: primer3 2.6.1 `oligotm.c` (code read), Owczarzy et al. 2008 Biochemistry 47:5336, PMID 18422348.

### 1.2 Fraction of template sites occupied (primer in excess)
Mass action with P₀ ≫ T₀ (true in early cycles: 200 nM primer vs sub-pM target sites):

θ = K·P₀ / (1 + K·P₀)

Exact two-species solution (needed only late in PCR, when product approaches primer concentration), with α = P₀/T₀:

θ = ½[(1 + α + 1/(K·T₀)) − √((1 + α + 1/(K·T₀))² − 4α)]

This is ThermoPlex eq 5 (paper) / eq 20 (supplement): there 1/(K·T₀) is written (α/c₀,p)·exp(ΔG/RT). Verified numerically: at T₀ = 1 pM and P₀ = 200 nM it equals K·P₀/(1+K·P₀) to 3 decimals across 50 to 78 °C (`fb.out` section 3).

Inverting for a design target θ*:

ΔG°(Ta) ≤ RT·ln(P₀·(1 − θ*)/θ*)  (1)

T_θ* = ΔH° / (ΔS° + R·ln(P₀(1 − θ*)/θ*))  (2)

Equation 2 is the temperature at which the primer reaches θ*; the pool Ta must not exceed min over primers of T_θ*.

Note on "Tm": the usual Tm = ΔH°/(ΔS° + R ln(C_T/4)) is the 50% point for two strands each at C_T/2. With primer in excess, θ = 0.5 occurs at ΔG = RT ln P₀, i.e. 2.0 to 2.3 °C above the C_T/4 Tm for 18 to 20-mers at 200 nM (`fb.out` section 4). So "Ta = Tm − 5 °C" is not a fixed fraction bound; it depends on ΔH.

### 1.3 Coupled equilibria: primer hairpin, primer dimer, template structure
Species: free primer p, primer hairpin H = K_hp·p (K_hp = exp(−ΔG_hp/RT), dimensionless), homodimer D = K_d·p², free site s, duplex C = K·p·s. Balances: P₀ = p(1 + K_hp) + 2K_d·p² + C; T₀ = s + C.

- Primer hairpin only, primer in excess: K_eff = K / (1 + K_hp), i.e. ΔG_eff = ΔG_dup + RT·ln(1 + exp(−ΔG_hp/RT)). If ΔG_hp ≪ 0 this tends to ΔG_dup − ΔG_hp (the full cost of opening the hairpin).
- Homodimer/heterodimers: they lower free p. Solve the quadratic (or the full multi-strand system) for p, then θ = K·p/(1 + K·p).
- Template site structure: K_eff = K·p_open, where p_open is the equilibrium probability that the whole site is single-stranded in the template ensemble; ΔG_eff = ΔG_dup + ΔG_open, ΔG_open = −RT ln p_open = G(site constrained unpaired) − G(unconstrained). This treats binding as all-or-none on an open site, so it is an upper bound on the true cost (the complex ensemble also allows partial invasion).
- Combined (DECIPHER form, see §3): K_ov = K / ((1 + K₂,eff)(1 + K₃)), θ = P·K_ov/(1 + P·K_ov).

Verified (`structcost.out` [A], 70% GC 20-mer at 62 °C, P₀ = 200 nM): the full mass-action solver and the shortcut K/(1+K_hp) agree to 4 decimals for hairpin ΔG from +1 to −2 kcal/mol. A homodimer with ΔG = −12 kcal/mol leaves only 17.6% of primer free yet occupancy falls only from 0.9994 to 0.9967, because K·P₀ is so large for this primer. Occupancy loss from sequestration matters only for primers near the edge (θ ≈ 0.9 to 0.97).

---

## 2. ThermoPlex (Agmata et al. 2025, Biol Methods Protoc 10:bpaf074, doi 10.1093/biomethods/bpaf074, PMC12598146; repo github.com/aagmata/ThermoPlex)

Read in full: paper text (PMC HTML and JATS XML), supplementary PDF text, and the MATLAB source (`ThermoDhyb.m`, `SimulEq.m`, `eqnsetLUM.m`, `eqnsetKeq.m`, `prodsolKeq.m`, `spsprimer.m`, `dnaenergies.m`, `NewtonR.m`, `sens_fxn.m`).

### 2.1 ThermoDHyb (hybridization ΔG)
- Smith–Waterman-like DP over a doublet matrix P (and a reverse matrix P′), cell recurrence ΔG(m,n) = min{bulge from (m−bu, n−1), bulge from (m−1, n−bl), NN stack + loop from (m−1, n−1), 0} (eq 1). NN parameters SantaLucia 1998; bulge init ΔG table `[4.0 2.9 3.1 3.2 3.3 3.5 3.7 3.9 4.1 4.3, then +0.2/nt]`, internal-loop init `[0 3.6 0.8 0.4 0.1 0.3 0.2 …]` kcal/mol, both converted to entropy at 310.15 K (code lines 69 to 73).
- Each local minimum of the bidirectional DP is one "microstate"; ensemble ΔG from a partition function over those minima. Code: `Z = sum(exp(-(E*1000)/(1.986*T)))`, `dGensemble = -1.986*T*log(Z)/1000` (ThermoDhyb.m lines 197 to 198). The paper's eq 2 and 3 (and supplement) print Z = Σ exp(ΔG/RT) and ΔG = RT ln Z without minus signs: a typesetting error; the code uses the standard signs.
- Only one structure per nucleation site is kept (stated limitation). Intra-strand structure is not included in ΔG; instead a primer is rejected if MATLAB `rnafold` predicts any base pair (`spsprimer.m`: reject if any '(' in the fold), regardless of stability. Note: that is an RNA folding routine applied to DNA.
- Salt: Owczarzy 2008 Mg correction applied to NN entropies with a hard-coded 50 mM monovalent (KCl) assumption (`dnaenergies.m` line 4 onward). No dNTP chelation term in that file.
- Benchmark vs Visual OMP at 37 °C, 1 M NaCl, 99 random pairs 20 to 40 nt: slope 0.9959, R² 0.9770, intercept 0.0192 kcal/mol; structure concordance 84%. NUPACK showed a ~3.41 kcal/mol offset. Runtime O(M·N).

### 2.2 Fractional conversion χ_t and specificity ψ
- χ_t = ½[1 + α + (α/c₀,p)·exp(ΔG/RT) − √((1 + α + (α/c₀,p)·exp(ΔG/RT))² − 4α)], α = c₀,p/c₀,t (eq 5; supplement eq 15 to 21 derives it from K = c_tp/(c_t·c_p) = exp(−ΔG/RT)). Code (`Xt` in `spsprimer.m`, `SimulEq.m`) uses R = 1.987 cal/mol/K and ΔG in kcal/mol.
- ψ = χ_t(ΔG_target) / χ_t(ΔG_non-target) (eq 6). Multiplex version ψ_i(T) = c_i,t / Σ_{j≠i} c_j,t (eq 7).
- Sensitivity (Monte Carlo n = 1000): χ is most sensitive to T, then [MgCl₂], then α; α matters only below ~30. Typical PCR α ≥ 10⁸ (citing Cha and Thilly 1993).

### 2.3 SiMulEq (multi-reaction equilibrium)
- Gibbs-energy minimisation with Lagrange multipliers, giving (eq 4a to 4c): c₀,i − c_i − c_ii − Σ_k c_ij = 0; ln(c_i/c₀,i) + λ_i/RT = 0; ΔG_ij/RT + ln(c_ij/(c₀,i c₀,j)) + (λ_i + λ_j)/RT = 0. Eliminating λ gives exactly c_ij = K_ij·c_i·c_j, so the system is ordinary mass action (an implementer can solve it with any mass-action solver, e.g. NUPACK `tube` concentrations or a Newton solver on log-concentrations).
- Solver: Newton with finite-difference Jacobian, variables squared to keep concentrations positive, stop at max residual ≤ 1e-10 or 100 iterations; initial guess from a simpler mass-action system (`eqnsetKeq.m`). All pairwise species including homodimers are generated (n(n+1)/2 products); template–template interactions excluded.
- Extension inhibition: a primer binding inside another primer's extension path blocks that product: χ_i = χ₀,i·Π_j(1 − χ_j) over primers binding downstream (supplement eq 23; `SimulEq.m` applies it by sorted binding position).
- Temperature sweep 40 to 70 °C in 2 °C steps, output equilibrium product distribution (EPD) curves.

### 2.4 Defaults used in the paper (Table 1)
T_simulation 55 °C; MgCl₂ 10 mmol/L; primer 0.5 µmol/L each; ψ threshold 1000; χ_t threshold 0.5; initial primer length 20; mismatch parameter M = 3; amplicon resolution R = 50 bp. Wet lab: 10 µL GoTaq reactions; 4-plex and 5-plex sardine COI assays.

Discrepancy: Fig. 5 caption says "# of mismatch < M", but `spsprimer.m` keeps a primer only if mismatches to the non-target are ≥ M (`if numbermismatch >= mismatchparameter`). The code logic (require mismatches for specificity) is the sensible one.

### 2.5 Stated limitations
Needs aligned, equal-length, gap-free input; fails when polymorphisms are sparse; no quantitative (qPCR) validation; no ranking/scoring of candidate sets; no primer–dimer penalty in selection ("does not yet incorporate explicit penalties for primer–primer hybridizations"); models only the first annealing event; tested on 4 to 5-plex only.

Relevance for us: ThermoPlex is a specificity tool for species discrimination, with χ_t ≥ 0.5 at 55 °C as the efficiency gate. It does not model template secondary structure, uses a binary hairpin filter, and a χ_t of 0.5 is a weak efficiency criterion. We should adopt its χ formula and the mass-action multi-species idea, not its thresholds.

---

## 3. Related equilibrium approaches

### 3.1 Pythia (Mann et al. 2009 NAR 37:e95, doi 10.1093/nar/gkp443)
- Binding and folding free energies from partition-function DP (pseudoknot-free, O(n⁴) folding), same NN parameter set for both.
- Chemical-reaction equilibrium by Gibbs-energy minimisation (gradient descent), chemical potentials μ_i = ΔG_i + RT ln(n_i/(n_A n_B)) for dimers and μ_i = ΔG_i + RT ln(n_i/n_A) for folded species.
- 11 competing reactions (Fig. 1): primer folding, primer homo/heterodimers, primers binding the template outside the priming region, and the desired priming-site binding.
- Evaluated "at a late stage of an idealized PCR" so that problems at high template concentration show up.
- Metric: equilibrium efficiency = min(fraction of left primer on its site, fraction of right primer on its site).
- SVM (quadratic kernel on the 11 ΔG values) as a fast prefilter; trained at thresholds 0.8, 0.85, 0.9, 0.95 (training sets 642 to 10 498; ROC > 0.99). Full equilibrium ~0.7 s per pair.
- Results: tiling repeats near TSSs, Pythia 195 pairs vs P316 (Primer3 + 16-nt 3′ check) 106; success 94% permissive / 80% stringent vs 95% / 82%; Pythia recall 97% vs 48%, precision 81% for both.
- Not verified: the threshold and the concentrations/temperature actually used in the tiling runs are in the supplementary data, which I did not obtain.

### 3.2 DECIPHER DesignPrimers / CalculateEfficiencyPCR (Wright et al. 2014, PMID 24750536)
Read from the R source (`eq/dec/CalculateEfficiencyPCR.R`, `dl/DP.R`; package version not recorded):
- Sodium equivalent: ions = monovalent + 3.3·√(divalent − dNTPs) (M). With defaults 0.07 M, 0.003 M, 0.0008 M this gives 0.225 M, matching the documented example `ions=0.225`.
- ΔG values from OligoArrayAux (`hybrid-min` for duplexes, `hybrid-ss-min` for folding) at the annealing temperature.
- K₁ = exp(−ΔG_dup/RT) (primer–target), K₂a primer folding, K₂b primer–primer, K₃ target folding over the site region.
- K₂,eff: 0 if K₂b = 0; K₂a if K₂b·P/K₂a < 0.01; K₂b if 8K₂b·P > (1 + K₂a)²; otherwise 4K₂b·P / (−1 − K₂a + √((1 + K₂a)² − 8K₂b·P)) − 1; floored at 0.
- K_ov = K₁ / ((1 + K₂,eff)(1 + K₃)); efficiency = P·K_ov/(1 + P·K_ov) × elongation efficiency (Taq 3′-mismatch model, `terminalMismatch`; set to 1 for proofreading polymerases).
- Defaults: annealingTemp 64 °C, P = 4e-7 M (400 nM), monovalent 0.07 M, divalent 0.003 M, dNTPs 8e-4 M, minEfficiency 0.8, primerDimer 1e-7 (max acceptable dimer amplification efficiency), length 17 to 26, product 75 to 1200 bp.
- Documentation: "an efficiency of 99% (0.99) will greatly lower predicted specificity ... an efficiency of 50% (0.5) may be too low in actuality".
- Caveat: K₃ is folding of the target substring only (local), not the full template.

### 3.3 SADDLE (Xie et al. 2022 Nat Commun 13:1881, doi 10.1038/s41467-022-29500-4, PMC9001684)
- Candidate primers truncated from the 3′ end until ΔG° is between −10.5 and −12.5 kcal/mol, target ≈ −11.5 kcal/mol ("best tradeoff between amplification efficiency/uniformity and nonspecific hybridization").
- ΔG° at 60 °C and "salinity" 0.18 (M Na⁺ equivalent), SantaLucia NN, ΔS corrected by 0.368·ln(S) per NN step, initiation 2.09 kcal/mol (Supplementary S9, example fP1 = −11.89 kcal/mol). GC 0.25 to 0.75.
- Wet lab: Vent (exo−), ThermoPol buffer, 30 nM each primer in Table S1 (main text quotes a median of 45 nM for the 96-plex), 17 cycles of 95 °C 30 s / 60 °C 3 min / 72 °C 30 s.
- Translated to occupancy (`fb.out` section 5), primer excess, 60 °C: at 30 nM, −10.5 / −11.5 / −12.5 kcal/mol give θ = 0.19 / 0.51 / 0.83; at 45 nM 0.26 / 0.61 / 0.88; at 200 nM 0.61 / 0.88 / 0.97. So SADDLE's working point is roughly θ ≈ 0.5 at its own concentrations and a long (3 min) anneal. That is a lower target than primer3's 97%, and it was chosen for uniformity and specificity in a 96-plex, not for maximum yield.

### 3.4 Visual OMP / SantaLucia (2007) Methods Mol Biol 402:3–34, doi 10.1007/978-1-59745-528-2_1
- Not verified from the chapter text (paywalled). Secondary sources: the primer3 manual cites the chapter's "third myth" (primers should be matched on fraction bound at Ta, not Tm). A 2010 GEN trade article by the vendor (genengnews.com/insights/design-tools-for-use-in-assay-optimization/) says Visual OMP uses "a multistate coupled equilibrium model" and ranks primers by sensitivity, defined as percent bound to target; its examples report 94 to 96% for good primers, 77 to 82% for weaker Primer3 primers, and > 99.9% bound with < 0.01% to off-targets for an 8-plex. These are example results, not stated thresholds.

### 3.5 primer3 PRIMER_ANNEALING_TEMP "bound" (primer3 2.6.1 manual and `oligotm.c`; primer3-py 2.0.3)
- Output PRIMER_{LEFT,RIGHT}_n_BOUND: "fraction of primers bound to a template at the PRIMER_ANNEALING_TEMP in percent if the primer and template would be at equal concentration".
- Formula (code and manual): ΔG = ΔH − T·ΔS, K = exp(−ΔG/RT), bound = 100 / (1 + √(1/((PRIMER_DNA_CONC/4e9)·K))) (non-self-complementary; 1e9 for symmetric). This is exactly 50% at the C_T/4 Tm but is an approximation elsewhere: it solves θ/(1−θ) = √(K·C/4), not the exact equal-concentration θ/(1−θ)² = K·C/2.
- Works only with SantaLucia or Schildkraut salt; with Owczarzy salt (PRIMER_SALT_CORRECTIONS=2) primer3-py raises "Primer binding does not work with Owczarzy salt correction". The Owczarzy branch in `oligotm.c` never sets `ret.bound`.
- Defaults: PRIMER_ANNEALING_TEMP −10 (off), PRIMER_OPT_BOUND 97.0, PRIMER_MIN_BOUND −10, PRIMER_MAX_BOUND 110, weights 0. Manual: "at the reduced annealing temperature 95-98% should be bound"; feature described as experimental.
- Reproduced exactly (`fb.out` section 6): primer GACCACCCAGGACGTGGAGG, 50 mM K⁺, 2 mM Mg²⁺, 0.8 mM dNTP, 200 nM, SantaLucia salt: primer3 BOUND 98.54 / 91.65 / 65.36% at 55 / 60 / 65 °C; my reimplementation identical; Tm 66.84 °C identical.
- Because it assumes equal concentrations, primer3's BOUND is much lower than the primer-excess occupancy relevant to early PCR cycles (`fb.out` section 3: 70% GC 20-mer at 66 °C: 0.835 primer3 vs 0.990 primer-excess).

### 3.6 Two primer3-py pitfalls found while checking (both verified by running code)
1. `primer3.bindings.calc_tm` default `formamide_conc=0.8`, which lowers Tm by (0.453·f_GC − 2.88)·0.8 ≈ 2.0 °C for a 70% GC 20-mer. Always pass `formamide_conc=0`.
2. primer3-py 2.0.3's bundled `oligotm.c` Owczarzy branch computes `(1/(2 * (len - 1)))` in integer arithmetic, which is 0, so the length-dependent term of Owczarzy 2008 eq 16 is dropped. Tm comes out 0.54 to 0.64 °C too high (`intdiv.out`: 55.87 vs 56.41, 70.72 vs 71.30, 65.33 vs 65.97 °C). Primer3 2.6.1's `oligotm.c` casts to double and is correct. Our own NN code reproduces primer3-py exactly when that term is dropped, and agrees with primer3 at 1 M Na⁺ (65.26 °C both).

### 3.7 Recent work (2020 to 2026)
- Gimpel et al. 2025 Nat Commun 16:9187, doi 10.1038/s41467-025-64221-4 (PMC12533003): 1D-CNN predicts sequence-specific efficiency in multi-template PCR (AUROC 0.88, AUPRC 0.44); about 2% of random sequences amplify poorly (as low as 80% of mean efficiency), driven mainly by adapter-mediated self-priming motifs next to priming sites. The "state of the art" filter it beats is GC 40 to 60%, homopolymer ≤ 4, mfold ΔG > −15 kcal/mol. Usable rule: check amplicon-internal sequences that can prime from the primer/adapter region (self-priming), not just primer–primer dimers.
- Sci Adv 2026, doi 10.1126/sciadv.aec9175 (PMC13015887): programmable amplification via primer tags; 2483 measurements; a ΔG°-only model gave R² = 0.62 and an ML model with sequence features R² = 0.86; "primers with identical ΔG° under the same reaction conditions exhibited substantial variation in amplification efficiency". Their simulation: a primer with equivalent efficiency 57.8% at 60 °C gives 95.6% at 55 °C and 15.2% at 65 °C, illustrating the steepness of θ(T).
- No 2020 to 2026 source I found gives a validated universal θ threshold; thresholds remain tool defaults (0.5 ThermoPlex, 0.8 DECIPHER, 0.95 to 0.98 primer3 manual).

---

## 4. Numeric sanity checks (2 mM Mg²⁺, 50 mM K⁺, 200 nM primer)

### 4.1 ΔG°(Ta) needed for a target occupancy (eq 1), kcal/mol (`fb.out` [2])

| P₀ | Ta | θ=0.5 | θ=0.9 | θ=0.97 | θ=0.99 |
|---|---|---|---|---|---|
| 200 nM | 55 °C | −10.06 | −11.49 | −12.33 | −13.06 |
| 200 nM | 60 °C | −10.21 | −11.67 | −12.51 | −13.25 |
| 200 nM | 65 °C | −10.37 | −11.84 | −12.70 | −13.45 |
| 45 nM | 60 °C | −11.20 | −12.65 | −13.50 | −14.24 |
| 500 nM | 55 °C | −9.46 | −10.89 | −11.73 | −12.46 |

Each tenfold drop in primer concentration needs 2.303·RT ≈ 1.5 kcal/mol more binding energy (at 60 °C).

### 4.2 θ versus temperature (primer excess; free Mg²⁺ 1.2 mM = 2 mM − 0.8 mM dNTP) (`fb.out` [3])

| T (°C) | 50 | 54 | 58 | 60 | 62 | 64 | 66 | 68 | 70 | 72 | 74 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 20-mer 45% GC ACTAGCATACTGGACATCTG | 0.996 | 0.937 | 0.475 | 0.186 | 0.055 | 0.015 | 0.004 | 0.001 | 0 | 0 | 0 |
| 20-mer 70% GC GGCACCAGCCAGCTGAGCCA | 1.000 | 1.000 | 1.000 | 1.000 | 0.999 | 0.998 | 0.990 | 0.962 | 0.867 | 0.627 | 0.307 |
| 18-mer 67% GC CACCAGCCAGCTGAGCCA | 1.000 | 1.000 | 0.998 | 0.991 | 0.970 | 0.902 | 0.725 | 0.436 | 0.187 | 0.065 | 0.021 |

Temperatures where θ reaches each target (`fb.out` [4]):

| primer | ΔH (kcal/mol) | ΔS (cal/mol/K) | Tm (C_T/4) | T50 | T90 | T97 | T99 | T50 − T99 |
|---|---|---|---|---|---|---|---|---|
| 45% GC 20-mer | −151.0 | −425.5 | 55.9 | 57.9 | 54.7 | 52.9 | 51.4 | 6.5 |
| 70% GC 20-mer | −159.1 | −429.3 | 70.7 | 72.8 | 69.5 | 67.7 | 66.0 | 6.7 |
| 67% GC 18-mer | −141.3 | −384.0 | 65.3 | 67.6 | 64.0 | 62.0 | 60.2 | 7.3 |

So θ goes from 0.5 to 0.99 over ~6.5 to 7.5 °C, and T97 is 3 to 3.5 °C below the C_T/4 Tm. Free Mg²⁺ matters: with 0.2 mM dNTP (1.8 mM free) the 45% GC primer at 58 °C rises from 0.475 to 0.622.

### 4.3 Tm matching versus ΔG(Ta) matching on real Mtb sequence (`pool.out`)
Sequence: Mtb H37Rv NC_000962.3:760900–761500 (rpoB, 601 nt, GC 0.667); all 16 to 30-mers on both strands (17 370).
- Tm-matched set (Tm 64.5 to 65.5 °C, n = 556, ΔH −202 to −125 kcal/mol): T97 spans 60.93 to 63.06 °C (2.1 °C spread); θ at 63 °C spans 0.908 to 0.971.
- ΔG(63 °C)-matched set (±0.25 kcal/mol around the θ = 0.97 value, −12.63, n = 688): θ at 63 °C spans 0.957 to 0.979; T97 spans 62.39 to 63.57 °C; Tm spans 65.0 to 67.0 °C.
So matching ΔG at Ta halves the spread in the quantity that matters, while Tm differs by up to 2 °C among equally good primers.

Sensitivity: dΔG/dT = −ΔS ≈ +0.39 kcal/mol per °C for an 18-mer; at fixed θ = 0.97 a 1 kcal/mol structure penalty lowers T_θ by 1/|ΔS + R ln(P₀(1−θ)/θ)| ≈ 2.3 °C.

### 4.4 Pool Ta choice with structure (`pool.out` [3])
Toy pool of 40 ΔG-matched primers with illustrative structure penalties (0, 0.5, 1, 2 kcal/mol): highest Ta at which all primers reach θ ≥ 0.9 / 0.97 / 0.99 is 64.4 / 62.4 / 60.6 °C without structure and 59.7 / 57.7 / 55.9 °C with structure. One primer with a 2 kcal/mol penalty pulls the whole pool down ~4.7 °C. This is why the designer should screen out (or redesign) the worst-structured sites rather than lower Ta for them.

### 4.5 Template structure cost (`structcost.out` [B], `structcost_B_dna.out`, `structcost_B_salt.out`)
Engineered strand: 150 nt random + blocker (complement of the site's 5′ k nt) + TTTTT + 20-nt site at the 3′ end; primer = rc(site). NUPACK 4 (dna04, 50 mM Na⁺, 2 mM Mg²⁺), ΔG_bind = G(complex) − G(primer) − G(strand) + RT ln ρ_water (ρ = 54.57 M at 60 °C) to convert to the 1 M standard state.

| stem k | NUPACK cost 60 °C | θ 60 °C (iso → full) | NUPACK cost 65 °C | θ 65 °C (iso → full) | Vienna open-site cost 60 °C (0.18 M Na) |
|---|---|---|---|---|---|
| 0 | −0.17 | 0.999 → 0.999 | −0.13 | 0.970 → 0.975 | 1.08 |
| 6 | 1.54 | 0.999 → 0.987 | 1.06 | 0.970 → 0.868 | n/a |
| 8 | 2.44 | 0.999 → 0.951 | 1.72 | 0.970 → 0.710 | 4.31 |
| 10 | 4.68 | 0.999 → 0.396 | 3.71 | 0.970 → 0.114 | n/a |
| 12 | 5.85 | 0.999 → 0.100 | 4.72 | 0.970 → 0.027 | 5.74 |

Findings: (a) a ≥ 10 bp stem that sequesters the 5′ half of the site makes the primer fail (θ < 0.4) even though the isolated duplex predicts θ = 0.999; (b) the cost is temperature-dependent and smaller at higher T; (c) the ViennaRNA "force site unpaired" cost (DNA Mathews 2004 parameters, `RNA.params_load_DNA_Mathews2004()`, md.salt) overestimates the NUPACK complex cost by ~1 to 2 kcal/mol, so it is a conservative screen; (d) salt matters greatly for Vienna: the same open-site cost is 3.00 / 1.08 / 0.22 kcal/mol at 1.021 / 0.18 / 0.05 M Na⁺ (stem 0). Without loading DNA parameters (RNA defaults), costs were 4 to 14 kcal/mol, so the parameter file must be loaded explicitly. An earlier scan of 30 real 20-nt sites on a 300-nt rpoB amplicon (`eq_code/eq_check2.out`) found median cost 0.29 to 0.40 kcal/mol and maximum 2.3 to 3.3 kcal/mol (55 to 65 °C), with the worst sites in the GC-rich core; NUPACK took ~24 ms per primer + 300-nt evaluation.

### 4.6 Model disagreement on ΔG itself (`structcost.out` [C], 70% GC 20-mer with its complement)

| T | NUPACK dna04 (50 mM Na, 2 mM Mg) | NN + Owczarzy, free Mg 1.2 mM | NN + Owczarzy, free Mg 2 mM | primer3 thal |
|---|---|---|---|---|
| 55 | −16.54 | −18.23 | −18.63 | −17.98 |
| 60 | −14.97 | −16.09 | −16.49 | −15.83 |
| 65 | −13.33 | −13.94 | −14.35 | −13.68 |
| 70 | −11.63 | −11.80 | −12.21 | −11.53 |

The models agree near 70 °C and diverge by up to 1.7 kcal/mol at 55 °C (NUPACK has a smaller effective ΔH). NUPACK does not model dNTP chelation, so pass free Mg²⁺ to it. Use one model consistently for thresholds; the NN + Owczarzy two-state model is the best-validated for perfect duplexes in Mg²⁺ buffers, NUPACK for structure costs (as differences).

---

## 5. How to choose the pool annealing temperature

1. For each primer i compute ΔH_i, ΔS_i (NN + Owczarzy, free Mg²⁺ = Mg − dNTP) and ΔG_struct,i(T) = primer hairpin term RT ln(1 + e^{−ΔG_hp/RT}) + template site cost (NUPACK complex vs isolated site, or Vienna open-site cost as a conservative bound), evaluated at a trial Ta.
2. T_θ*,i is the root of ΔH_i − T ΔS_i + ΔG_struct,i(T) = RT ln(P₀(1 − θ*)/θ*).
3. Ta_pool = min_i T_θ*,i, then iterate once or twice (structure costs depend on T). Report Ta_pool and the per-primer θ at Ta_pool.
4. Off-target bound: at Ta_pool require that the strongest off-target 3′-anchored duplex for each primer has θ_off ≤ 0.01 (ψ ≥ ~100 at θ_on ≈ 0.97 to 0.99). ThermoPlex used ψ ≥ 1000 for species-level discrimination; for genome-wide off-target sites in a single-genome panel, 100 is a reasonable start (not validated).
5. Design the pool so that the spread of T_θ* is small (≤ 2 °C); then Ta_pool sits close to every primer's T97, so primers are both efficient (θ ≥ 0.97) and not overly stable (θ ≤ ~0.995, which limits mismatch tolerance and off-target priming).

Target fractions in the literature: 0.5 at Tm by definition (and ThermoPlex χ_t gate); 0.5 to 0.6 effective for SADDLE's −11.5 kcal/mol at 30 to 45 nM, 60 °C; 0.8 DECIPHER default; 0.8 to 0.95 Pythia SVM training range; 0.95 to 0.98 primer3 manual (default PRIMER_OPT_BOUND 97). Practical for a 20-plex: θ ≥ 0.95 for every primer at Ta_pool, aim 0.97 to 0.99.

---

## 6. Recommended numeric targets for a 20-plex, 200 nM per primer
Conditions assumed: 50 mM K⁺, 2 mM Mg²⁺, 0.8 mM total dNTP (1.2 mM free Mg²⁺); expose all as parameters.

| quantity | target |
|---|---|
| θ at Ta_pool, every primer, including structure | ≥ 0.95 hard, 0.97 to 0.99 preferred |
| equivalent ΔG°(Ta) for θ = 0.97 / 0.99 | −12.5 / −13.25 kcal/mol at 60 °C; −12.7 / −13.45 at 65 °C (scale with eq 1) |
| upper stability (avoid over-stable primers) | θ_on ≤ 0.995, i.e. ΔG°(Ta) ≥ −13.7 kcal/mol at 60 °C, ≥ −13.9 kcal/mol at 65 °C (eq 1 with θ = 0.995) |
| spread of ΔG_eff(Ta) across pool | ≤ 1.0 kcal/mol (≈ ±1 °C in T97) |
| primer hairpin | ΔG_hp(Ta) ≥ 0 kcal/mol preferred; reject if ΔG_hp < −1 (costs > 0.3 kcal/mol in ΔG_eff) |
| template site structure cost (NUPACK) | ≤ 1.0 kcal/mol preferred; reject > 2.0 |
| any primer–primer dimer (all 40 × 40 pairs) at Ta | free primer fraction ≥ 0.9; and no 3′-extensible dimer with ΔG ≤ −6 kcal/mol (heuristic, see open questions) |
| off-target site occupancy | θ_off ≤ 0.01 at Ta_pool |
| Ta_pool for Mtb (65% GC) | aim 63 to 66 °C: choose primer lengths so T97 ≈ Ta_target + 1 °C |
| Ta_pool for ~50% GC genomes | aim 58 to 60 °C |

---

## 7. Code snippet (core functions, from `eqdesign/fb.py`)

```python
R = 1.98720425864083e-3  # kcal/(mol K)
def K(dg, Tc): return math.exp(-dg/(R*(Tc+273.15)))
def theta_excess(dg, Tc, P): x = K(dg, Tc)*P; return x/(1+x)
def dG_needed(th, Tc, P): return R*(Tc+273.15)*math.log(P*(1-th)/th)
def T_at_theta(s, th, P):           # eq 2
    dH, dS = thermo(s); return dH/(dS + R*math.log(P*(1-th)/th)) - 273.15
def chi_thermoplex(dg, Tc, c0p, c0t):  # ThermoPlex suppl. eq 20
    a = c0p/c0t; b = 1 + a + (a/c0p)*math.exp(dg/(R*(Tc+273.15)))
    return 0.5*(b - math.sqrt(b*b - 4*a))
```

---

## 8. Open questions
- Pythia's actual threshold and simulation conditions (supplementary data not obtained).
- Visual OMP's exact multi-state model and default % bound targets (chapter not read).
- Which ΔG model (NN + Owczarzy vs NUPACK) better predicts real primer θ at 55 to 65 °C; they differ by up to 1.7 kcal/mol.
- Kinetics: annealing times of 15 to 60 s may not reach equilibrium for structured templates; equilibrium θ is an upper bound.
- How θ thresholds translate into amplicon read-depth uniformity: no quantitative calibration found; needs a pilot run.
- Dimer thresholds in ΔG terms for 3′-extensible dimers are heuristics, not derived from θ.
