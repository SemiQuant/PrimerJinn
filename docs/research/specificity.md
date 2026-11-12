# Genome-wide specificity for a multiplex primer pool

Scope: off-target priming sites and unintended products (including cross-pair products) for one multiplex pool, fast on 1 to 10 Mb bacterial genomes and usable on larger ones. All numbers below come from code run on this machine (macOS arm64, 16 cores, Python 3.13.2, primer3-py 2.0.3, numba, BLAST+ 2.16.0) against `/Users/semiquant/Bioinformatics/Projects/Archived/PrimerJinn/example/ref.fasta` (M. tuberculosis H37Rv, NC_000962.3, 4,411,532 bp), unless a literature source is cited.

Deliverables
- Module: `/private/tmp/claude-501/-Users-semiquant--claude-spares-spare-n9Plyg/5e124d61-7974-4f46-b4a7-ff4fbdb28e2f/scratchpad/research/specificity_reference.py`
- Tests: `/private/tmp/claude-501/-Users-semiquant--claude-spares-spare-n9Plyg/5e124d61-7974-4f46-b4a7-ff4fbdb28e2f/scratchpad/research/test_specificity_reference.py` (copy of `spec_work/t/test_specificity.py`; "ALL TESTS PASSED")
- Experiment scripts and outputs: `spec_work/t/` (t1 to t12, `pool200.tsv`, `truth_broad.npy`, `truth_broad_scored.pkl`, `t11_out.txt`)

Units: temperature C, dG kcal/mol, salts mM, oligo nM, lengths nt/bp.

---

## 1. What existing tools do

### 1.1 NCBI Primer-BLAST
Source: Ye J et al. 2012, BMC Bioinformatics 13:134, doi:10.1186/1471-2105-13-134 (PMC3412702), and the live form downloaded 2026-11-10 (`dl_primerblast/index.html`).

- Search: BLAST, then Needleman-Wunsch re-alignment of each hit to the full primer length, because BLAST local alignments truncate the primer ends.
- Sensitivity: "can detect a target that contains up to 35% mismatches" (7 mismatches for a 20-mer).
- BLAST defaults (paper and form): E value 30,000 (primer-only case), word size 7, max database sequences 50,000, reward:penalty 1:1. Low-complexity filter on (form).
- Products: checks F-R, F-F and R-R combinations; max target amplicon size `defval="4000"` on the form.
- Stringency (form, `selected` options): a primer pair counts as specific only if at least one primer has at least 2 total mismatches to the unintended target, including at least 2 mismatches within the last 5 bp at the 3' end; targets with 6 or more mismatches to a primer are ignored. Note: the HTML carries `defVal="1"` for the two mismatch selects while the `selected` option is 2; the user-visible default is 2. I could not resolve which the server uses when the field is omitted.
- Paper's rule of thumb: "a two base mismatch at the 3' end generally prevents amplification"; single mismatches, even at the 3' terminus, still allow amplification at reduced efficiency.

### 1.2 MFEprimer-3.0
Source: Wang K et al. 2019, NAR 47:W610, doi:10.1093/nar/gkz351 (PMC6602485); CLI docs for v3.1 (`spec_work/dl/mfe31.html`).

- k-mer index of the genome; the primer's 3'-terminal k-mer (k = 9) is expanded into all variants with 1 mismatch except at the 3'-terminal base (9*3 + 1 = 28 sequences), and these are looked up exactly. Then binding stability (Tm/dG) of the full primer is computed.
- CLI defaults: `-k 9`, `--misStart 1`, `--misEnd 9`, `--misMatch 1`, `--maxSize 2000`, `-t/--tm 30` (Tm cutoff for amplicons), `--mono 50`, `--diva 1.5`, `--dntp 0.25`, `--oligo 50`.
- Limitation: a site needs a near-exact 3' 9-mer; terminal mismatches are never found (fine for Taq, wrong for proofreading enzymes, see 1.4).

### 1.3 ThermoAlign
Source: Francis F et al. 2017, Sci Rep 7:44437, doi:10.1038/srep44437 (PMC5353602).

- Off-target search: blastn, `e_value 30000; word_size 7; gapopen 2; gapextend 2; reward 1; penalty -1; perc_identity 70; max_target_seqs 13; max_hsps 20`.
- Each HSP is edited into a full-length ungapped "thermoalignment" (end-filled) before computing Tm with primer3 thal. Native truncated BLAST alignments had Tm on average 7 C higher than their thermoalignment, range -14 to +272 C, so raw HSPs are unusable for thermodynamics.
- Mismatch count "is not a suitable proxy for the potential of mispriming".
- Pair selection: lower primer Tm minus the best off-target Tm of either primer must be at least 10 C.
- Admitted limitation: BLAST "does not allow for mismatches with the sequence used to seed an alignment".

### 1.4 Polymerase extension from mismatched 3' ends

| Source | Finding used here |
|---|---|
| Huang MM, Arnheim N, Goodman MF 1992, NAR 20:4567, doi:10.1093/nar/20.17.4567 (PMC334186) | Taq, 3'-terminal mispair, relative extension efficiency (primer base first): transitions A.C, C.A, G.T, T.G 1e-3 to 1e-4; T.C, T.T 1e-4 to 1e-5; A.A about 1e-6; G.A, A.G, G.G, C.C below 1e-6; C.T about 1e-2. Independent of temperature (45, 55, 70 C). |
| Stadhouders R et al. 2010, J Mol Diagn 12:109, doi:10.2353/jmoldx.2010.090035 (PMC2797725) | Real-time PCR: single 3'-region mismatches cost < 1.5 Ct (A-C, C-A, T-G, G-T) up to > 7 Ct (A-A, G-A, A-G, C-C). Up to sevenfold differences between master mixes. |
| Lefever S et al. 2013, Clin Chem 59:1470, doi:10.1373/clinchem.2013.203653 | >= 4 mismatches in one primer, or 3 + 2 across both primers, needed to block qPCR completely; single mismatches > 5 bp from the 3' end have a moderate effect. |
| Rejali NA, Moric E, Wittwer CT 2018, Clin Chem 64:801, doi:10.1373/clinchem.2017.282285 | Taq (deletion mutant) extension at 65 C: 177 to 205 nt/s matched vs 0.008 to 4.55 nt/s with 3'-terminal mismatch; a C.C mismatch 6 nt from the 3' end still cut the rate by up to 84 %. KCl inhibits mismatch extension. |
| Gohl DM et al. 2021, NAR 49:e87, doi:10.1093/nar/gkab471 (PMC8421222) | Proofreading enzymes edit primer 3' ends: near-complete editing of up to 6 terminal bases for KAPA HiFi, only the last 2 for Phusion; Q5, Pfu, Vent also edit; no editing by three Taq products. A single phosphorothioate bond blocks editing beyond it (KAPA HiFi, Q5, Phusion). With Taq, primer variants within the last 5 bases were recovered at 0.922, 0.774, 0.404, 0.318, 0.189 of expected (the paper lists the values without restating the position order; from context the strongest penalty is nearest the 3' end, but I did not verify the mapping in the figure). |

Abstracts only (not full text) were read for Stadhouders, Lefever and Rejali; the Huang numbers are from the abstract. The Kwok 1990 paper was not retrieved and is not relied on.

Consequences for the model
- Taq: the 3'-terminal base and, less strongly, positions 2 to 5 control extension; 2 mismatches in the last 5 nt nearly blocks priming. A mismatch-only rule (Primer-BLAST) is a coarse version of this.
- Proofreading (Q5, Phusion, KAPA HiFi): a terminal mismatch is not protective because the exonuclease removes it; off-target risk then depends on duplex stability of the rest of the primer. Treat 3' mismatches as protective only for Taq-type enzymes, or when primers carry 3' phosphorothioate bonds.

---

## 2. How primer3-py thal treats mismatches and dangling ends (run here)

Script `spec_work/exp2_thal.py`, output `spec_work/exp2_out.txt`. Conditions: 50 mM monovalent, 2 mM Mg2+, 0.8 mM dNTP, 250 nM oligo, 60 C. Primer `GAGGAACACCACTAGTACCG`.

```
[1] perfect                                dg= -10.87 kcal/mol tm= 59.76
[2] perfect + 5nt template flanks          dg= -11.67 kcal/mol tm= 61.39
[2b] perfect + 1nt flank at primer 3' side dg= -11.47
[2c] perfect + 1nt flank at primer 5' side dg= -11.27
mm at  1 from 3': het dg= -9.64 | end dg= -3.38 | flanked het dg= -9.85 | flanked end dg= -6.87
mm at  2 from 3': het dg= -8.82 | end dg= -6.94 | flanked het dg= -9.02 | flanked end dg= -7.74
mm at  3 from 3': het dg= -7.43 | end dg= -7.43
mm at 10 from 3': het dg= -9.70 ; identity at pos 10: template A -8.34, C -7.97, G -9.70 (G.T wobble), T -10.87 (match)
[7] 1-nt deletion in template (bulge)      dg= -6.28
12 kb template -> error "Target sequence length > maximum allowed (10000)"
calc_heterodimer 30-nt target: 68.3 us/call; calc_end_stability: 67.2 us/call
```

Findings
- `calc_heterodimer` (thal ANY) finds the best gapped alignment, including internal mismatches (type-specific, G.T is mild), bulges and dangling ends. Template flanks add 0.4 to 0.8 kcal/mol (dangling ends / terminal mismatches), so the template passed to thal should include 2 nt on each side.
- `calc_end_stability` (thal END1) forces the primer 3' end to be paired. With a 3'-terminal mismatch it returns a much weaker duplex (-3.38 vs -9.64) and, with flanks, a slipped alignment (-6.87, identical for all three wrong bases). It therefore cannot rank terminal mismatches; use ANY binding plus a polymerase factor.
- thal source (`src_p3/thal.c`): salt term `0.368 * ln((mv + 120*sqrt(max(0, dv - dntp)))/1000)` per (paired bases - 1) (line 1042, 3027); duplex initiation dH 200 cal/mol, dS -5.7 cal/K/mol (lines 588-589); Tm uses `R ln(C/4)` for non-self-complementary pairs (line 593); dG is reported without the concentration term. So `frac_bound` must add concentration itself (done in the module).
- primer3-py 2.0.3 thal is not thread-safe: 8 threads calling `calc_heterodimer` died with exit 138 (SIGBUS) in 2 of 2 runs; an 8-process pool gave identical sums and a 3.75x speedup (`spec_work/t/gil2.py`).

Fraction bound and dG thresholds (two-state, primer in excess at C = 250 nM, theta = 1/(1 + exp(dG/RT)/C)):

| Ta | theta 0.5 | theta 0.01 | theta 1e-4 |
|---|---|---|---|
| 55 C | -9.91 | -6.92 | -3.91 |
| 60 C | -10.06 | -7.02 | -3.97 |
| 65 C | -10.22 | -7.13 | -4.03 |

A 10-nt GC-rich perfect 3' match is already near theta 0.3 at 60 C: thal gave dG -9.39 kcal/mol for `...GGCCACCCGC` (10 nt paired, 11 nt unpaired), and a hand nearest-neighbour sum (SantaLucia 2004 stacks + thal salt term) gave -9.3. In a 65 % GC genome such short perfect 3' matches are common, which is why mismatch counting alone underestimates risk.

---

## 3. Implementation (specificity_reference.py)

Pipeline
1. `find_sites_seed()`: seed table of three 8-mers taken 0, 5 and 10 nt from the primer 3' end, each with all 1-mismatch variants (25 per seed, both strands), stored in a CSR table indexed by the 2-bit 8-mer code. One numba pass per genome chunk (16 threads, numba `nogil`). Each hit is extended to the full ungapped primer footprint; for each footprint the scan computes mismatches and an optimistic dG estimate at Ta: initiation + every Watson-Crick stack whose two pairs both match (SantaLucia & Hicks 2004 values) + thal salt term; mismatch loops cost nothing. Keep if `dG_est <= est_max` (default -5.0).
2. `score_sites()`: thal ANY on primer vs template footprint +- 2 nt; thal END1 also when the 3' base matches and the enzyme is Taq. Fraction bound at Ta (`conc_divisor = 1`, primer in excess). Extension factor for Taq (below). `score = theta * ext`. Process pool (fork).
3. `enumerate_products()`: sites with score >= `site_min_score` are merged per locus (same primer and strand, 3' ends within 20 nt, keep best) to remove shifted duplicates, then every + site is paired with every - site whose footprint lies downstream within 40 to 2000 bp, for all primers in the pool (F-R, F-F, R-R, cross-pair). `product score = score_plus * score_minus`. Optional circular contigs.

Extension factor (Taq-type only; `polymerase="proofreading"` sets it to 1)
- 3'-terminal mismatch: 0.01 if primer:template is A:A, A:G, G:A, G:G or C:C (Huang < 1e-6 or about 1e-6; Stadhouders > 7 Ct), else 0.3.
- Single mismatch at position 2, 3, 4, 5 from the 3' end: 0.3, 0.4, 0.8, 0.9.
- Two or more mismatches in the last 5: capped at 0.01.
These factors are heuristic: the kinetic ratios (1e-2 to 1e-6) are per extension event, and in PCR a product needs one successful extension per cycle among many binding events, so effective penalties are much smaller (Stadhouders 1.5 to 7 Ct, i.e. about 3- to 130-fold). I chose values in that Ct range; they are not fitted to data.

Changes made in this session to the earlier, unreviewed prototype (backup at `spec_work/specificity_reference_v1_backup.py`)
- G:G added to the severe terminal set (Huang 1992 lists G.G below 1e-6; the prototype omitted it).
- Hamming caps (7 mismatches, 35 %) removed; they dropped primable sites such as an 11-mismatch footprint whose 3' 10 nt match (dG -9.39 at 60 C). Replaced by the dG prefilter.
- Two seeds (offsets 0, 8) changed to three (0, 5, 10) after the sensitivity test below.
- Thread pool for thal replaced by a fork process pool (thread crash above).
- Linear contigs by default (`circular=False`); the prototype treated every contig as circular.
- Locus merging before product enumeration (a primer 6 nt off its true site re-found the true duplex with score 0.93 and created a duplicate "unintended" product).
- `keep_min_score` (1e-6) to cap memory.

---

## 4. Sensitivity and speed (200 primers, 4.4 Mb)

Test pool: 100 primer3 pairs designed in random 700-bp windows of H37Rv, Tm 63 to 67 C, 18 to 26 nt, products 250 to 350 bp (`spec_work/t/mk_primers.py`, `pool200.tsv`).

Ground truth (`t4_truth.py`): exhaustive numba scan of every position, primer and strand, keeping a footprint if total mismatches <= 7, or <= 2 mismatches in the 3'-terminal 12 nt, or <= 2 mismatches in positions 7 to 18 (what remains after up to 6 nt of primer editing), the last two with total <= 0.6 L. 518,677 sites in 11.3 s; all scored with thal (6.5 to 11.8 s per condition on 16 processes). Primable = score >= 0.01.

| Condition | primable sites (score >= 0.01) | score >= 1e-4 |
|---|---|---|
| Taq, Ta 60 | 6,248 | 62,943 |
| Taq, Ta 65 | 1,459 | 28,049 |
| Proofreading, Ta 60 | 27,758 | 198,408 |
| Proofreading, Ta 65 | 5,912 | 105,538 |

Recall of primable truth sites, locus level (a hit within +- 5 nt counts; `t9_locus.py`). Search time is candidate generation only, 16 threads.

| Method | Ta | candidates | time | Taq >= 0.01 | Taq >= 1e-4 | Proof >= 0.01 | Proof >= 1e-4 |
|---|---|---|---|---|---|---|---|
| blastn-short r1/p-1 e1000 | 60 | 763,125 | 11.0 s | 71.0 % | 58.8 % | 79.8 % | 58.4 % |
| blastn-short r1/p-1 e30000 | 60 | 1,990,394 | 29.2 s | 71.1 % | 59.5 % | 80.2 % | 59.4 % |
| seed 2 x 8-mer (0, 8), est <= -3 | 60 | 1,438,402 | 1.5 s | 95.0 % | 82.4 % | 79.7 % | 62.0 % |
| seed 3 x 8-mer (0, 5, 10), est <= -3 | 60 | 2,063,560 | 2.2 s | 96.1 % | 87.5 % | 93.3 % | 78.6 % |
| blastn-short r1/p-1 e1000 | 65 | 763,125 | 11.0 s | 83.8 % | 62.7 % | 88.1 % | 66.2 % |
| seed 3 x 8-mer, est <= -3 | 65 | 1,374,077 | 1.4 s | 98.2 % | 91.4 % | 96.6 % | 84.1 % |

Seed of 9-mers (0, 6, 12) was faster (0.9 s) but recall fell to 76 % (Taq, 60 C), so 8-mers are kept.

Effect of the dG prefilter cut (`t10_cut.py`; fraction of primable truth sites whose estimate passes, and seed candidates):

| Ta | est_max | Taq >= 0.01 kept | Taq >= 1e-4 kept | Proof >= 0.01 kept | candidates |
|---|---|---|---|---|---|
| 60 | -3 | 100.0 % | 99.7 % | 100.0 % | 2,063,560 |
| 60 | -4 | 99.9 % | 96.4 % | 99.8 % | 1,178,854 |
| 60 | -5 | 99.5 % | 83.2 % | 98.5 % | 531,938 |
| 60 | -6 | 93.4 % | 58.3 % | 86.5 % | 185,436 |
| 65 | -4 | 99.8 % | 91.5 % | 99.8 % | 595,258 |
| 65 | -5 | 99.3 % | 66.5 % | 97.4 % | 195,097 |

The estimate equals thal within 0.27 kcal/mol for perfect duplexes (50 primers, mean thal - est = -0.06 at 60 C), but thal can be up to 7.3 kcal/mol more negative for mismatched footprints because it finds gapped or shifted alignments, so the estimate is a fast filter, not a bound. Default `est_max = -5.0`; use -4.0 when the "possible" tier matters.

BLAST notes
- blastn-short defaults (verified with `blastn -help` and an outfmt 0 footer): E 1000, dust off, reward 1, penalty -3, gap 5/2, word size >= 4 allowed, max_target_seqs 500. With these defaults (t3, older caps) 3-mismatch sites were found 125/135 and 7-mismatch 30,997/258,709.
- Best practical settings for this task: `-task blastn-short -word_size 7 -evalue 1000 -reward 1 -penalty -1 -gapopen 2 -gapextend 2 -dust no -soft_masking false -strand both -max_target_seqs 100000 -max_hsps 100000`, outfmt 6 with qstart/qend/sstart/send, then anchor on the 3' end and end-fill each HSP to the full ungapped footprint (ThermoAlign style). Raising E to 30,000 cost 2.7x time for < 1 % more recall; word size 6 took 112 s; perc_identity is not useful after end-filling (filter thermodynamically instead).
- `-num_threads 16` gave no speedup (10.52 s vs 10.39 s with 1 thread) for a one-sequence database, so BLAST is about 5x (e1000) to 13x (e30000) slower than the seed scan here and less sensitive for Taq-primable sites. BLAST's advantage is scale: an indexed database and no per-run table; it remains a reasonable option for > 100 Mb genomes.

End-to-end runtime of `check_pool()` with defaults (seed 3 x 8-mer, est <= -5, 16 workers; `t11_runtime.py`, output `spec_work/t/t11_out.txt`):

| Polymerase | Ta | candidates | find_sites | score_sites | products | total | primable sites | likely products (intended) |
|---|---|---|---|---|---|---|---|---|
| Taq | 60 | 531,938 | 0.57 s | 13.9 s | 0.04 s | 14.6 s | 7,689 | 933 (103) |
| Taq | 65 | 195,097 | 0.28 s | 4.9 s | 0.01 s | 5.2 s | 1,475 | 175 (103) |
| Proofreading | 60 | 531,938 | 0.69 s | 11.2 s | 0.50 s | 12.4 s | 54,247 | 19,251 (103) |
| Proofreading | 65 | 195,097 | 0.28 s | 3.7 s | 0.03 s | 4.2 s | 8,133 | 608 (103) |

Peak RSS 1.2 GB for the four runs in one process. thal is about 95 % of runtime. For a 10 Mb genome expect about 2.3x these times; for a 3 Gb genome candidate counts grow about 700x, so thal must only see candidates passing a stricter prefilter (for example est <= -6 and a primable-only mode), or BLAST with a database should be used.

What the "likely" products are (Taq, `t12_cat.py`)
- 103 intended products for 100 pairs: pairs A023 and A071 amplify a second perfect copy (400201-400550 and 606560-606909; 1095866-1096159 and 2944400-2944693), and A026 gives a 1,700 bp product between perfect sites (429837-431537). These are real repeats in H37Rv, correctly flagged.
- Ta 65: 72 unintended likely; 61 have one perfectly matched primer and one off-target site, 7 have two off-target sites, 2 are cross-pair products between perfectly matched sites of different amplicons (e.g. A048_F/A093_R, 1,520 bp). A typical off-target site: A070_R with mismatches only at positions 17 to 21 (3' 16 nt perfect), dG -10.0 kcal/mol, score 0.42.
- Ta 60: 830 unintended likely. Dropping Ta from 65 to 60 C multiplies primable sites by 5.2 (Taq) and likely products by 11.

---

## 5. Recommended thresholds

| Quantity | Default | Reason |
|---|---|---|
| Ta for scoring | the recommended PCR Ta of the pool (about 65 C for Mtb pools) | risk changes 5 to 10x between 60 and 65 C |
| Oligo concentration | per-primer working concentration (default 250 nM), primer in excess (divisor 1) | thal dG excludes the concentration term |
| Site "primable" | score = theta(Ta) * ext >= 0.01 (theta 0.01 is dG <= -7.0 kcal/mol at 60 C, 250 nM) | 1 % occupancy with an intact 3' end can seed product over 30 cycles |
| Product "likely" | score_plus * score_minus >= 0.01, size 40 to 2000 bp | same scale as a single primable site; both sites needed |
| Product "possible" | product score >= 1e-4 | report, do not reject |
| Taq terminal mismatch | factor 0.01 for A:A, A:G, G:A, G:G, C:C; 0.3 otherwise | Huang 1992, Stadhouders 2010 |
| Taq 2+ mismatches in last 5 | cap 0.01 | Primer-BLAST rule of thumb, Lefever 2013 |
| Proofreading enzymes | no 3' mismatch protection (ext 1); thal ANY only | Gohl 2021 (editing up to 6 nt) |
| Candidate search | seeds 8-mers at 0, 5, 10 nt from 3' end, 1 mismatch each; est_max -5.0 | 96 to 98 % locus recall for Taq-primable sites, 0.3 to 0.7 s |
| Rejection rule in design | reject or replace a primer if it creates any unintended "likely" product shorter than about 1.5x the longest intended amplicon, or more than N (e.g. 3) primable off-target loci | short off-target products out-compete intended ones; not validated here |

These thresholds are calibrated against a thermodynamic model, not against wet-lab outcomes. They should be checked against at least one sequenced pool (off-target read fraction per primer).

---

## 6. Bugs and limitations of the current PrimerJinn/PCRinSilico.py

File: `/Users/semiquant/Bioinformatics/Projects/Archived/PrimerJinn/PrimerJinn/PCRinSilico.py` (286 lines). Verified items were run in `spec_work/t` (see notes).

Correctness bugs
1. L88: off-target Tm is computed from `sseq` alone as a perfect duplex, so mismatches never lower Tm. Run: `Tm_NN` gives 60.8 C for the perfect site and also 60.8 C for the same site with 2 3'-terminal mismatches.
2. L85: `req_five` only requires the alignment to start at query position 1 (5' end); nothing requires the 3' end (`qend == len(primer)`), which is the end that matters. HSPs truncated at the 3' end are accepted.
3. L20: `--req_five` is `type=str, default=True`; any value, including "False", is truthy, so it cannot be turned off.
4. L82-83 and L99: forward sites are positioned at their 5' end (`sstart`) but reverse sites at `send`, the 3'-end side; product sizes are short by about one primer length and `ref_amplicon` omits the reverse primer site.
5. L237: `full_amplicon = pf + ref_amplicon + pr` duplicates the forward primer, because `ref_amplicon` already starts at the forward primer 5' end.
6. L224: `pr_seq = str(pf_seq.reverse_complement())` uses `pf_seq`; the reverse primer is replaced by the reverse complement of the forward primer (and raises AttributeError when `pf_seq` is already a str).
7. L105-108: in the `elif` branch the keys are crossed: `qseq1` is row2 but `qseq1_input` and `mismatch1` come from row1, so primer names, sequences and mismatch counts are mismatched in output.
8. L125: primer-dimer Tm is `Tm_NN(seq1 + seq2)`, the Tm of the concatenated 40-mer as a perfect duplex. Run: 73.0 C for a pair whose thal heterodimer Tm is -14.1 C, so nearly every pair is flagged.
9. L64: Q5 branch with `seq2` never assigns `tm` (UnboundLocalError), hidden by the bare `except: continue` at L266, so the interactions sheet is empty in Q5 mode.
10. L66/L260: non-Q5 branch calls `Tm_NN(amplicon1, c_seq=amplicon2)` with unrelated sequences of different length. Run: raises `ValueError: no thermodynamic data for neighbors 'CG/TG'`, again swallowed; the amplicon-interaction output is meaningless or empty.
11. L213: `SeqIO.read` accepts only single-record FASTA; multi-contig assemblies crash after BLAST (which did search all contigs), and the contig ID of each product is ignored when slicing.
12. L18: salt help says nM; the value is passed to Biopython as mM Na+. Mg2+ and dNTP are ignored outside the Q5 branch.
13. L40: annealing temperature is silently `target_tm - 4`; L17 declares `type=int` with a float default.

Method limitations
14. L72: blastn-short with default scoring (reward 1, penalty -3, E 1000, max_target_seqs 500); no end-filling, so local alignments truncate mismatched ends (ThermoAlign: truncated HSP Tm on average 7 C too high). In my test the default scoring found 3-mismatch sites 125/135 and 7-mismatch sites 12 %.
15. No polymerase or 3' mismatch model; no Taq vs proofreading distinction.
16. Criterion is Tm >= Ta, a hard cut on a single-duplex Tm, instead of occupancy at Ta with primer concentration.
17. L93-111: O(n^2) pandas `.loc` over all hit pairs; slow with thousands of hits. Products between hits on different contigs are excluded correctly, but circular replicons are not handled.
18. L28-34: argparse and the blastn check run at import, so the module cannot be imported or tested; L12 imports `requests` for an unused NEB API helper (L49-55).
19. L142, L156-158: only ACGTN accepted (no IUPAC); an invalid line `break`s the loop and the run continues with a truncated primer list. L146-150: the name counter `i` is incremented once, outside the loop.
20. Uses deprecated primer3-py camelCase names (`calcHeterodimer`, `bindings.calcTm`); they still exist in 2.0.3 (checked) but are deprecated.

---

## 7. Open issues

- Extension factors for Taq positions 2 to 5 and the 0.01 / 0.3 terminal factors are heuristics in the range of published Ct shifts, not fitted.
- Two-state equilibrium at Ta ignores kinetics (annealing time, polymerase binding of transient duplexes) and competition between products; product score = product of site scores is a ranking, not a yield.
- Per-enzyme editing depth (Phusion about 2 nt, KAPA HiFi about 6 nt) is not modelled; proofreading mode assumes full trimming of mismatched 3' bases.
- Recall was measured against a truth set that itself has inclusion rules (section 4); sites with > 7 mismatches outside both 12-nt windows are not in it.
- For > 100 Mb genomes thal on every candidate is too slow; a tighter prefilter or a learned dG surrogate is needed.
- The Primer-BLAST default stringency (form shows 2/2/5 selected while defVal attributes say 1) was not confirmed against server behaviour.
