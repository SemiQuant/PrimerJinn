import math, primer3, pytest
from Bio import SeqIO
import structure_reference as sr
G = str(next(SeqIO.parse('/Users/semiquant/Bioinformatics/Projects/Archived/PrimerJinn/example/ref.fasta', 'fasta')).seq).upper()
AMP = G[761005:761305]

def test_probs_window_equals_constraint():
    n = len(AMP)
    up = sr.unpaired_windows(AMP, 65.0, 22)
    for a, b in ((n - 22, n), (100, 122), (0, 8)):
        d1 = sr.dg_open(up[b][b - a], 65.0)
        d2 = sr.dg_open_constrained(AMP, a, b, 65.0)
        assert abs(d1 - d2) < 0.02

def test_ensemble_energy_monotonic_in_T():
    gs = [sr.fold_stats(AMP, t).g_ensemble for t in (55, 60, 65, 68, 72, 80, 90)]
    assert all(x < y for x, y in zip(gs, gs[1:])), gs

def test_no_gt_pairs():
    s = sr.fold_stats(AMP, 60.0).structure
    import RNA
    pt = RNA.ptable(s)
    assert not any(pt[i] > i and {AMP[i-1], AMP[pt[i]-1]} == {'G', 'T'} for i in range(1, len(AMP)+1))

def test_product_tm_matches_primer3():
    c = sr.Conditions()
    # formula of primer3 oligotm.c long_seq_tm(); note primer3-py 2.0.3 calc_tm() returns about 2 C lower
    gc = (AMP.count('G') + AMP.count('C')) / len(AMP)
    na = c.mono_mM + 120 * math.sqrt(c.mg_mM - c.dntp_mM)
    ref = 81.5 + 16.6 * math.log10(na / 1000) + 41 * gc - 600 / len(AMP)
    assert abs(sr.duplex_tm_long(AMP, c) - ref) < 1e-9

def test_na_eq_matches_nupack():
    nupack = pytest.importorskip('nupack')
    c = sr.Conditions()
    a = nupack.pfunc([AMP], model=nupack.Model(material='dna', celsius=65, sodium=0.05, magnesium=c.mg_free_M))[1]
    b = nupack.pfunc([AMP], model=nupack.Model(material='dna', celsius=65, sodium=c.na_eq_M))[1]
    assert abs(a - b) < 1e-3

def test_fraction_bound_limits():
    fwd = AMP[:22]
    m = sr.site_metrics(sr.revcomp(AMP), 278, 300, fwd, 65.0)
    assert 0 <= m.f_bound <= m.f_bound_nostruct <= 1
    assert m.dg_eff >= m.dg_duplex

def test_g4_finder():
    s = 'AAAA' + 'GGGTTGGGTTGGGTTGGG' + 'AAAA'
    assert sr.g4_motifs(s) == [(4, 22, '+')]
    assert sr.g4_motifs(sr.revcomp(s)) == [(4, 22, '-')]

def test_local_tm_salt_dependence():
    hp = 'GCATGCTTTTGCATGC'
    t_hi = sr.local_melting_tm(hp, sr.Conditions(mono_mM=1000, mg_mM=0, dntp_mM=0))
    t_lo = sr.local_melting_tm(hp, sr.Conditions())
    assert t_lo < t_hi
