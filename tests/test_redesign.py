"""Redesign with kept primers: fixed pairs, excluded 3' ends, preserved pools."""
import os

import pandas as pd
import pytest

from primerJinn.multiplex import design, redesign, report
from primerJinn.multiplex.candidates import design_target, fixed_target
from primerJinn.multiplex.config import Config
from primerJinn.multiplex.optimize import dsatur
from primerJinn.multiplex.redesign import FixedPair, FixedPrimer
from primerJinn.multiplex.targets import Site, Target

REF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "example", "ref.fasta")
META = {"fasta": REF, "gc_percent": 65.61}


def test_dsatur_preset():
    edges = {0: {1}, 1: {0, 2}, 2: {1}, 3: set()}
    c = dsatur(edges, preset={0: 1, 2: 1, 3: 0})
    assert c[0] == 1 and c[2] == 1 and c[1] != 1 and c[3] == 0
    c2 = dsatur(edges, preset={0: 0, 1: 0})        # a conflicting preset is dropped for the second node
    assert c2[0] == 0 and c2[1] != 0


def test_fixed_target_scores_like_design_target():
    store = design.SeqStore(REF)
    cfg = Config(specificity=False, threads=1).resolve(0.656)
    t = Target("T01", "NC_000962.3", 2155168, 2155168, [Site("NC_000962.3", 2155168, 2155168)])
    region, rs = store.region(t.seqid, t.start - 625, t.end + 625)
    res = design_target(t, region, rs, store.lengths[t.seqid], cfg, keep_pairs=5)
    pr = res["pairs"][0]
    fp = FixedPair("T01", FixedPrimer(pr.fwd.seq, "F", t.seqid, pr.fwd.start, pr.fwd.end),
                   FixedPrimer(pr.rev.seq, "R", t.seqid, pr.rev.start, pr.rev.end), 1)
    fx = fixed_target(t, region, rs, store.lengths[t.seqid], cfg, fp)
    assert fx["fail"] is None and len(fx["pairs"]) == 1
    q = fx["pairs"][0]
    assert q.fixed and q.fwd.fixed and q.rev.fixed
    for a, b in ((pr.fwd, q.fwd), (pr.rev, q.rev)):
        assert b.dg_eff == pytest.approx(a.dg_eff, abs=1e-6)
        assert b.site_open_genomic == pytest.approx(a.site_open_genomic, abs=1e-6)
    assert q.pair_d3 == pytest.approx(pr.pair_d3, abs=1e-6)
    # a primer that does not match the genome is refused, not redesigned silently
    bad = FixedPair("T01", FixedPrimer("ACGT" * 5, "F", t.seqid, pr.fwd.start, pr.fwd.start + 19), fp.rev, 1)
    assert "does not match the genome" in fixed_target(t, region, rs, store.lengths[t.seqid], cfg, bad)["fail"]


def test_read_previous_and_tails(tmp_path):
    rows = []
    for t, nm, d, seq, s, e, pool in (("T01", "T01_F", "forward", "ACGTACGTACGTACGTAC", 100, 117, 1),
                                      ("T01", "T01_R", "reverse", "TTGCATGCATGCATGCAA", 400, 417, 1),
                                      ("T02", "T02_F", "forward", "aaCGTACGTACGTACGTAC", 1100, 1118, 2),
                                      ("T02", "T02_R", "reverse", "TTGCATGCATGCATGCAA", 1400, 1417, 2)):
        rows.append(dict(target=t, primer=nm, direction=d, sequence=seq, order_sequence=seq, added_5prime_nt=seq.count("a"),
                         seqid="chr", start=s, end=e, three_prime_end=e if d == "forward" else s, pool=pool))
    p = tmp_path / "x_primers.tsv"
    pd.DataFrame(rows).to_csv(p, sep="\t", index=False)
    prev = redesign.read_previous(str(p), replace=["T01"])
    assert set(prev.fixed) == {"T02"} and prev.prev_pools == {"T02": 2}
    assert prev.fixed["T02"].fwd.seq == "AACGTACGTACGTACGTAC" and prev.fixed["T02"].fwd.ext5 == 2
    assert prev.avoid_3prime == {("chr", "F", 117), ("chr", "R", 400)}
    assert prev.avoid_sites == [("chr", "F", 100, 117), ("chr", "R", 400, 417)]
    with pytest.raises(ValueError, match="not in"):
        redesign.read_previous(str(p), replace=["T09"])
    with pytest.raises(ValueError, match="both kept and replaced"):
        redesign.read_previous(str(p), replace=["T01"], keep=["T01", "T02"])
    with pytest.raises(ValueError, match="same --tails"):
        redesign.read_previous(str(p), replace=["T01"], tail_fwd="GGGG")


def test_redesign_keeps_primers_and_pools(tmp_path):
    """Design two targets, replace one: the kept pair is byte-identical, the new one avoids the old 3' ends."""
    sites = [Site("NC_000962.3", 761155, 761155, "rpoB_S450"), Site("NC_000962.3", 2155168, 2155168, "katG_S315")]
    cfg = Config(specificity=False, threads=2, candidates_per_target=12, sa_restarts=2)
    d1 = design.design_pool(sites, META, cfg, verbose=False)
    t1 = report.tables(d1)["primers"]
    prefix = str(tmp_path / "d1")
    t1.to_csv(prefix + "_primers.tsv", sep="\t", index=False)
    files = report.write_order_sheets(t1, prefix)
    assert os.path.exists(prefix + "_order_sheet.tsv") and len(files) == 1
    full = report.write(d1, str(tmp_path / "full" / "d1"))      # the whole output set, map included
    assert any(f.endswith("_map.png") for f in full) and os.path.getsize(tmp_path / "full" / "d1_map.png") > 5000
    assert "![Amplicon map](d1_map.png)" in (tmp_path / "full" / "d1_report.md").read_text()
    prev = redesign.read_previous(prefix + "_primers.tsv", replace=["T01"])
    cfg2 = Config(specificity=False, threads=2, candidates_per_target=12, sa_restarts=2)
    redesign.apply(cfg2, prev)
    d2 = design.design_pool(sites, META, cfg2, verbose=False)
    t2 = report.tables(d2)["primers"]
    cols = ["primer", "sequence", "order_sequence", "start", "end", "three_prime_end", "pool", "added_5prime_nt"]
    k1 = t1[t1.target == "T02"][cols].reset_index(drop=True)
    k2 = t2[t2.target == "T02"][cols].reset_index(drop=True)
    assert k1.equals(k2), (k1, k2)
    assert t2[t2.target == "T02"]["kept_from_previous"].all()
    assert not t2[t2.target == "T01"]["kept_from_previous"].any()
    old3 = set(t1[t1.target == "T01"].three_prime_end)
    new3 = set(t2[t2.target == "T01"].three_prime_end)
    assert old3.isdisjoint(new3)
    assert not d2.coverage_problems and any("keep their primers" in n for n in d2.redesign_notes)
    files = report.write_order_sheets(t2, str(tmp_path / "d2"))
    new = pd.read_csv(tmp_path / "d2_order_sheet_new.tsv", sep="\t")
    assert list(new.fwd_name) == ["T01_rpoB_S450_F"] and len(files) == 2
    # the kept pair still appears in the final dimer table against the new primers
    dm = pd.DataFrame(d2.dimers)
    assert ((dm.primer1.str.startswith("T02")) & (dm.primer2.str.startswith("T01"))).any() or \
           ((dm.primer1.str.startswith("T01")) & (dm.primer2.str.startswith("T02"))).any()
