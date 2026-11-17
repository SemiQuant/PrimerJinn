"""Unit tests for specificity_reference.py on a synthetic genome with planted sites.
Run: python3 test_specificity.py"""
import os, random, sys, tempfile
sys.path.insert(0, "/private/tmp/claude-501/-Users-semiquant--claude-spares-spare-n9Plyg/5e124d61-7974-4f46-b4a7-ff4fbdb28e2f/scratchpad/research")
import numpy as np
import specificity_reference as S

random.seed(7)
F = "GTTTCTGTACGGCCCACTCGA"      # 21 nt
R = "GAATCCCGTCGTCGCAAACAG"      # 21 nt
X = "ACATCAATTACGTGCTGGCGC"      # 21 nt
rc = S.revcomp


def mut(seq, pos_from_3p, to=None):
    """Change the base `pos_from_3p` (1 = 3' terminal) of a primer-sense sequence."""
    i = len(seq) - pos_from_3p
    b = to or {"A": "C", "C": "A", "G": "T", "T": "G"}[seq[i]]
    assert b != seq[i]
    return seq[:i] + b + seq[i + 1:]


def plant(g, pos, s):
    g[pos:pos + len(s)] = list(s)


N = 200_000
c1 = [random.choice("ACGT") for _ in range(N)]
c2 = [random.choice("ACGT") for _ in range(50_000)]
# 1. intended product on contig 1: F at 1000 (+), R at 1300 (-)  -> size 1300 + 21 - 1000 = 321
plant(c1, 1000, F); plant(c1, 1300, rc(R))
# 2. cross-pair off-target on contig 1: X with 1 mismatch at position 10 (+) at 50000, F perfect (-) at 50500
plant(c1, 50000, mut(X, 10)); plant(c1, 50500, rc(F))
# 3. same-primer product: R (+) at 90000 and R (-) at 90200, both perfect
plant(c1, 90000, R); plant(c1, 90200, rc(R))
# 4. Taq-blocked site: F with C:C terminal mismatch. F ends in A, so use X (ends in C): template base must be C
#    -> genome base in primer sense = G (complement of template C)
plant(c1, 120000, mut(X, 1, "G")); plant(c1, 120300, rc(R))
# 5. two mismatches in the last 5 (positions 2 and 4) for F at 150000, R (-) at 150300
plant(c1, 150000, mut(mut(F, 2), 4)); plant(c1, 150300, rc(R))
# 6. wrap-around product on contig 2 (circular): F (+) 100 nt before the end, R (-) at 150
plant(c2, 50_000 - 100, F); plant(c2, 150, rc(R))
# 7. site with 3 mismatches all in the 5' half (positions 14, 17, 20) for X (+) at 170000
plant(c1, 170000, mut(mut(mut(X, 14), 17), 20))

tmp = tempfile.mkdtemp()
fa = os.path.join(tmp, "g.fa")
with open(fa, "w") as fh:
    fh.write(">c1\n" + "".join(c1) + "\n>c2\n" + "".join(c2) + "\n")
g = S.load_genome(fa)
primers = [F, R, X]; names = ["F", "R", "X"]
P = S.Params(ta_c=60.0, circular=True)

seed = S.find_sites_seed(primers, g, P)
brute = S.find_sites_brute(primers, g, 10)
blast = S.find_sites_blast(primers, fa, g, P)
ks = lambda a: set(map(tuple, a[:, :4].tolist()))
expected = {  # (contig, primer, strand, start): mismatches
    (0, 0, 1, 1000): 0, (0, 1, -1, 1300): 0, (0, 2, 1, 50000): 1, (0, 0, -1, 50500): 0,
    (0, 1, 1, 90000): 0, (0, 1, -1, 90200): 0, (0, 2, 1, 120000): 1, (0, 1, -1, 120300): 0,
    (0, 0, 1, 150000): 2, (0, 1, -1, 150300): 0, (1, 0, 1, 49900): 0, (1, 1, -1, 150): 0, (0, 2, 1, 170000): 3}
for name, arr in (("seed", seed), ("brute", brute), ("blast", blast)):
    d = {tuple(r[:4]): r[4] for r in arr.tolist()}
    for k, mm in expected.items():
        assert k in d, (name, "missing", k)
        assert d[k] == mm, (name, k, d[k], mm)
    print(f"{name}: all {len(expected)} planted sites found with exact coordinates and mismatch counts ({len(arr)} sites total)")
assert {k for k in ks(seed) if dict((tuple(r[:4]), r[4]) for r in seed.tolist())[k] <= 10} <= ks(brute)
b3 = {tuple(r[:4]) for r in brute.tolist() if r[4] <= 3}
assert b3 <= ks(seed), "seed missed a <=3 mismatch site"
print("seed is a subset of brute and contains every brute site with <= 3 mismatches")

sc = S.score_sites(seed, primers, g, P)
by = {(s["contig"], s["primer"], s["strand"], s["start"]): s for s in sc}
assert by[(0, 0, 1, 1000)]["score"] > 0.5 and by[(0, 0, 1, 1000)]["ext"] == 1.0
s4 = by[(0, 2, 1, 120000)]; assert s4["mm_pos"] == "1" and s4["ext"] == P.taq_term_severe, s4
s5 = by[(0, 0, 1, 150000)]; assert s5["mm_pos"] == "2,4" and s5["ext"] == P.taq_two_in_last5, s5
s2 = by[(0, 2, 1, 50000)]; assert s2["mm_pos"] == "10" and s2["ext"] == 1.0
print("scores: perfect %.3f | 1 mm pos 10 %.3f (dG %.2f) | C:C terminal %.4f | 2 mm in last 5 %.5f | 3 mm 5' half %.4f"
      % (by[(0, 0, 1, 1000)]["score"], s2["score"], s2["dg_any"], s4["score"], s5["score"], by[(0, 2, 1, 170000)]["score"]))

prods = S.enumerate_products(sc, primers, g, P, names=names, pairs=[(0, 1)], only_primable=False)
key = {}
for p in prods:  # prods are sorted by score, keep the best product per (contig, start, fwd, rev)
    key.setdefault((p["contig"], p["start"], p["fwd"], p["rev"]), p)
p1 = key[("c1", 1000, "F", "R")]; assert p1["size"] == 321 and p1["intended"] and p1["call"] == "likely", p1
p2 = key[("c1", 50000, "X", "F")]; assert p2["size"] == 521 and not p2["intended"], p2
p3 = key[("c1", 90000, "R", "R")]; assert p3["size"] == 221 and p3["same_primer"] and p3["call"] == "likely"
p4 = key[("c1", 120000, "X", "R")]; assert p4["call"] != "likely", p4
p5 = key[("c1", 150000, "F", "R")]; assert p5["call"] != "likely", p5
p6 = key[("c2", 49900, "F", "R")]; assert p6["wraps_origin"] and p6["size"] == 100 + 150 + 21, p6
print("products: intended 321 bp likely; cross-pair X/F 521 bp call=%s score=%.3f; same-primer R/R 221 bp likely; "
      "C:C terminal call=%s; 2-in-last-5 call=%s; wrap-around %d bp call=%s"
      % (p2["call"], p2["score"], p4["call"], p5["call"], p6["size"], p6["call"]))

# proofreading polymerase: 3' mismatches no longer protect
Pq = S.Params(ta_c=60.0, polymerase="proofreading", circular=True)
scq = S.score_sites(seed, primers, g, Pq)
byq = {(s["contig"], s["primer"], s["strand"], s["start"]): s for s in scq}
assert byq[(0, 2, 1, 120000)]["ext"] == 1.0 and byq[(0, 2, 1, 120000)]["score"] > s4["score"] * 10
print("proofreading: C:C terminal site score %.3f (Taq %.4f); 2-in-last-5 site score %.4f (Taq %.5f)"
      % (byq[(0, 2, 1, 120000)]["score"], s4["score"], byq[(0, 0, 1, 150000)]["score"], s5["score"]))
# G:G terminal mismatch is severe (Huang et al. 1992), C:T is mild
assert S.extension_factor([1], "G", "G", P) == P.taq_term_severe
assert S.extension_factor([1], "C", "T", P) == P.taq_term_mild
# linear contigs: no wrap-around product by default
prods_lin = S.enumerate_products(sc, primers, g, S.Params(ta_c=60.0), names=names, only_primable=False)
assert not any(p["wraps_origin"] for p in prods_lin)
# process-pool scoring gives the same numbers as serial scoring
big = np.repeat(seed, 150, axis=0)
a = S.score_sites(big, primers, g, P, threads=1); b = S.score_sites(big, primers, g, P, threads=8)
assert [x["score"] for x in a] == [x["score"] for x in b]
print("G:G severe, linear default, process pool == serial (%d sites)" % len(big))
print("ALL TESTS PASSED")
