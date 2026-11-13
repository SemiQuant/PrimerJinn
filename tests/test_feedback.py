"""Sequencing feedback: reading the counts table, matching, decisions and the rebalance rule."""
import os

import numpy as np
import pandas as pd
import pytest

from primerJinn.multiplex import feedback as fb

HERE = os.path.dirname(os.path.abspath(__file__))
COUNTS = os.path.join(HERE, "..", "example", "sequencing_run", "counts.tsv")

PRIMER_COLS = ["pool", "target", "target_name", "covers", "primer", "direction", "sequence", "order_sequence",
               "added_5prime_nt", "seqid", "start", "end", "three_prime_end", "length"]


def _prev(tmp_path, rows, conc=None):
    """A previous-design primers table: rows of (target, primer, direction, seq, pool)."""
    recs = []
    for t, nm, d, seq, pool in rows:
        recs.append(dict(pool=pool, target=t, target_name=nm.rsplit("_", 1)[0].split("_", 1)[-1], covers="",
                         primer=nm, direction=d, sequence=seq, order_sequence=seq, added_5prime_nt=0, seqid="chr",
                         start=1, end=len(seq), three_prime_end=len(seq), length=len(seq)))
    p = tmp_path / "prev_primers.tsv"
    pd.DataFrame(recs, columns=PRIMER_COLS).to_csv(p, sep="\t", index=False)
    cpath = None
    if conc:
        cpath = tmp_path / "conc.tsv"
        pd.DataFrame(list(conc.items()), columns=["primer", "conc_nM"]).to_csv(cpath, sep="\t", index=False)
    return fb.load_previous(str(p), str(cpath) if cpath else None, default_nM=100.0)


def _run(tmp_path, rows, sep="\t", name="counts.tsv"):
    p = tmp_path / name
    pd.DataFrame(rows).to_csv(p, sep=sep, index=False)
    return fb.load_run(str(p))


ROWS = [("T01", "T01_a_F", "forward", "ACGTACGTACGTACGTAC", 1), ("T01", "T01_a_R", "reverse", "TTGCATGCATGCATGCAA", 1),
        ("T02", "T02_b_F", "forward", "GGGTACGTACGTACGTAC", 1), ("T02", "T02_b_R", "reverse", "TTGCATGCATGCATGAAA", 1),
        ("T03", "T03_c_F", "forward", "CCGTACGTACGTACGTAC", 2), ("T03", "T03_c_R", "reverse", "TTGCATGCATGCATGCCC", 2)]


def test_factor_rule():
    th = fb.Thresholds()
    assert fb.factor(4.0, 100, th) == pytest.approx(0.5)
    assert fb.factor(0.25, 100, th) == pytest.approx(2.0)
    assert fb.factor(100.0, 100, th) == pytest.approx(0.25)       # capped
    assert fb.factor(0.0, 0, th) == 4.0                          # dropout: the cap


def test_status_rule():
    assert fb.status_of(1.0, 10) == "good" and fb.status_of(3.0, 10) == "warn"
    assert fb.status_of(0.1, 10) == "bad" and fb.status_of(0.0, 0) == "drop"


def test_load_counts_aliases_other_rows_and_csv(tmp_path):
    rows = [dict(target="T01", count=1000), dict(target="T02", count=250), dict(target="T03", count=0),
            dict(target="Unaligned", count=500), dict(target="primer dimer", count=50), dict(target="other", count=10)]
    run = _run(tmp_path, rows, sep=",", name="counts.csv")
    u = run.units.set_index("amplicon")
    assert list(u.index) == ["T01", "T02", "T03"]
    assert run.other == {"unaligned": 500, "dimer": 50, "off_target": 10}
    # median over amplicons with reads (1000, 250) = 625; the other rows are not amplicons
    assert u.loc["T01", "fold"] == pytest.approx((1000 + 0.5) / (625 + 0.5), rel=1e-3)
    assert u.loc["T03", "status"] == "drop" and u.loc["T01", "status"] == "good"
    assert not run.notes
    run2 = _run(tmp_path, [dict(amplicon="T01", reads=5)], name="c2.tsv")
    assert any("no unaligned row" in n for n in run2.notes)
    with pytest.raises(ValueError, match="lacks column"):
        _run(tmp_path, [dict(foo="T01", bar=5)], name="c3.tsv")
    with pytest.raises(ValueError, match="not a number"):
        _run(tmp_path, [dict(amplicon="T01", reads="many")], name="c4.tsv")


def test_pooling_over_samples(tmp_path):
    rows = [dict(sample="a", amplicon="T01", reads=1000), dict(sample="a", amplicon="T02", reads=100),
            dict(sample="b", amplicon="T01", reads=100), dict(sample="b", amplicon="T02", reads=1000)]
    run = _run(tmp_path, rows)
    u = run.units.set_index("amplicon")
    assert u.loc["T01", "reads"] == 1100 and u.loc["T01", "samples"] == 2
    assert u.loc["T01", "fold"] == pytest.approx(u.loc["T02", "fold"], rel=1e-6)   # symmetric: geometric mean


def test_decisions_and_rebalance(tmp_path):
    prev = _prev(tmp_path, ROWS)
    run = _run(tmp_path, [dict(amplicon="T01", reads=5000), dict(amplicon="T02", reads=300),
                          dict(amplicon="T03", reads=0), dict(amplicon="unaligned", reads=1500)])
    out = fb.analyse(run, prev)
    t = out.targets.set_index("target")
    assert t.loc["T01", "action"] == "rebalance"
    assert t.loc["T02", "action"] == "replace" and "above 500" in t.loc["T02", "reason"]   # 0.11x at 100 nM
    assert t.loc["T03", "action"] == "replace" and t.loc["T03", "status"] == "drop"
    assert out.replace == ["T02", "T03"]
    s = dict(out.summary)
    assert s["total reads"] == "6,800"
    assert s["designed reads %"] == f"{100 * 5300 / 6800:.1f}"
    assert s["unaligned or off-target reads %"] == f"{100 * 1500 / 6800:.1f}"
    assert not any("not from a designed amplicon" in w for w in out.warnings)      # 22 % < 30 %
    # overrides
    out2 = fb.analyse(run, prev, keep=["T02"], replace=["T01"])
    t2 = out2.targets.set_index("target")
    assert t2.loc["T02", "action"] == "rebalance" and t2.loc["T01", "reason"] == "replaced on request"
    p = out2.primers.set_index("primer")
    med = np.median([5000, 300])
    fac = min(4.0, ((med + 0.5) / (300 + 0.5)) ** 0.5)
    assert p.loc["T02_b_F", "action"] == "raise"
    assert p.loc["T02_b_F", "suggested_nM"] == pytest.approx(100 * fac, abs=0.1)
    assert p.loc["T02_b_F", "ul_stock_per_rxn"] == pytest.approx(100 * fac * 25 / 10000, abs=0.01)
    conc = fb.conc_table(out2, new_primers=["X_F"], default_nM=80)
    assert set(conc[conc.primer == "X_F"].conc_nM) == {80}
    assert "T01_a_F" not in set(conc.primer)       # replaced primers are not in the concentration table
    with pytest.raises(ValueError, match="not in"):
        fb.analyse(run, prev, replace=["T09"])


def test_concentration_decides_raise_or_redesign(tmp_path):
    """A weak amplicon is raised when its primers were dilute, redesigned when it would need > conc_max."""
    prev = _prev(tmp_path, ROWS)
    # T02 at 0.1x: at 100 nM it would need ~1000 nM -> replace; at 25 nM it needs ~250 nM -> raise
    run = _run(tmp_path, [dict(amplicon="T01", reads=3000, conc_nM=100), dict(amplicon="T02", reads=300, conc_nM=25),
                          dict(amplicon="T03", reads=3000, conc_nM=100)])
    out = fb.analyse(run, prev)
    t = out.targets.set_index("target")
    assert t.loc["T02", "action"] == "rebalance" and "raise the concentration" in t.loc["T02", "reason"]
    assert t.loc["T02", "current_nM"] == 25
    p = out.primers.set_index("primer")
    assert p.loc["T02_b_F", "action"] == "raise" and p.loc["T02_b_F", "conc_source"] == "counts table"
    assert p.loc["T02_b_F", "suggested_nM"] == pytest.approx(25 * ((3000.5 / 300.5) ** 0.5), abs=0.5)
    run2 = _run(tmp_path, [dict(amplicon="T01", reads=3000, conc_nM=100), dict(amplicon="T02", reads=300, conc_nM=100),
                           dict(amplicon="T03", reads=3000, conc_nM=100)], name="c2.tsv")
    t2 = fb.analyse(run2, prev).targets.set_index("target")
    assert t2.loc["T02", "action"] == "replace" and "above 500" in t2.loc["T02", "reason"]
    # --conc per primer wins over the column; a pair already at conc_max and under 0.5x is replaced
    prev3 = _prev(tmp_path, ROWS, conc={"T02_b_F": 500, "T02_b_R": 500})
    run3 = _run(tmp_path, [dict(amplicon="T01", reads=1000), dict(amplicon="T02", reads=400),
                           dict(amplicon="T03", reads=1000)], name="c3.tsv")
    out3 = fb.analyse(run3, prev3)
    t3 = out3.targets.set_index("target")
    assert t3.loc["T02", "action"] == "replace" and "already at 500" in t3.loc["T02", "reason"]
    # a dilute weak pair is raised but never above conc_max
    th = fb.Thresholds(conc_max=60)
    out4 = fb.analyse(run, prev, th)
    p4 = out4.primers.set_index("primer")
    assert out4.targets.set_index("target").loc["T02", "action"] == "replace"      # 250 needed > 60
    # not weak (0.3x) but already at 400 nM: raised, capped at conc_max
    run5 = _run(tmp_path, [dict(amplicon="T01", reads=3000, conc_nM=100), dict(amplicon="T02", reads=900, conc_nM=400),
                           dict(amplicon="T03", reads=3000, conc_nM=100)], name="c5.tsv")
    p5 = fb.analyse(run5, prev).primers.set_index("primer")
    assert p5.loc["T02_b_F", "action"] == "raise" and p5.loc["T02_b_F", "suggested_nM"] == 500


def test_mostly_off_target_warns(tmp_path):
    prev = _prev(tmp_path, ROWS)
    run = _run(tmp_path, [dict(amplicon="T01", reads=100), dict(amplicon="T02", reads=100),
                          dict(amplicon="T03", reads=100), dict(amplicon="unaligned", reads=900)])
    out = fb.analyse(run, prev)
    assert any("75% of reads are not from a designed amplicon" in w for w in out.warnings)


def test_write_files(tmp_path):
    prev = _prev(tmp_path, ROWS, conc={"T01_a_F": 50})
    run = _run(tmp_path, [dict(amplicon="T01", reads=100), dict(amplicon="T02", reads=100),
                          dict(amplicon="T03", reads=100)])
    out = fb.analyse(run, prev)
    assert set(out.primers.action) == {"keep"}
    p = out.primers.set_index("primer")
    assert p.loc["T01_a_F", "current_nM"] == 50 and p.loc["T01_a_F", "conc_source"] == "given"
    files = fb.write(out, str(tmp_path / "o" / "x"))
    names = {os.path.basename(f) for f in files}
    assert names == {"x_rebalance.tsv", "x_rebalance_targets.tsv", "x_conc.tsv", "x_feedback.tsv", "x_feedback.md"}
    fbt = pd.read_csv(tmp_path / "o" / "x_feedback.tsv", sep="\t")
    assert "amplicons within 2x of the median" in set(fbt["item"])


def test_matching_by_primer_name_target_name_and_merged(tmp_path):
    """Rows may use the primer name without _F/_R, the target name, or 'A/B' for merged amplicons."""
    prev = _prev(tmp_path, ROWS)       # target names a, b, c
    run = _run(tmp_path, [dict(amplicon="T01_a", reads=800), dict(amplicon="b/T03_c_F", reads=200),
                          dict(amplicon="rpoB", reads=5)])
    unit_of, warn = fb.match(run, prev)
    assert unit_of == {"T01": "T01_a", "T02": "b/T03_c_F", "T03": "b/T03_c_F"}
    assert len(warn) == 1 and "rpoB" in warn[0]
    out = fb.analyse(run, prev)
    t = out.targets.set_index("target")
    assert t.loc["T02", "reads"] == 200 and t.loc["T03", "amplicon"] == "b/T03_c_F"


@pytest.mark.skipif(not os.path.exists(COUNTS), reason="example counts missing")
def test_example_counts_table(tmp_path):
    run = fb.load_run(COUNTS)
    assert run.other == {"unaligned": 578, "off_target": 13301}
    u = run.units.set_index("amplicon")
    assert len(u) == 23 and u.loc["ddn", "status"] == "drop"
    live = u[u["reads"] > 0]["reads"]
    assert u.loc["katG", "fold"] == pytest.approx((3697 + 0.5) / (np.median(live) + 0.5), rel=1e-3)
    # a design named after the run's amplicons: dropouts and weak amplicons are replaced
    recs = []
    for i, nm in enumerate(u.index):
        recs += [(f"T{i + 1:02d}", f"T{i + 1:02d}_{nm.replace('/', '_')}_F", "forward", "ACGTACGTACGTACGTAC", 1),
                 (f"T{i + 1:02d}", f"T{i + 1:02d}_{nm.replace('/', '_')}_R", "reverse", "TTGCATGCATGCATGCAA", 1)]
    prev = _prev(tmp_path, recs)
    out = fb.analyse(run, prev)
    t = out.targets.set_index("target_name")
    assert t.loc["ddn", "action"] == "replace" and t.loc["atpE", "action"] == "rebalance"
    assert out.primers.set_index("primer").loc["T01_atpE_F", "action"] == "lower"
