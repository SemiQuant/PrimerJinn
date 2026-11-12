import os

import pytest

from primerJinn.multiplex.targets import Site, coverage_gaps, group_sites, span_limits

HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLE = os.path.join(HERE, "..", "example")


META = {"sequences": [{"seqid": "NC_000962.3", "length": 4411532, "aliases": ["NC_000962"]}],
        "n_sequences": 1, "default_seqid": "NC_000962.3", "default_unambiguous": True}


def read_positions(src):
    from primerJinn.multiplex.cli import sites_from_positions
    return sites_from_positions(src, META)


def test_read_inputs(tmp_path):
    sites = read_positions("7570,7582 761155;NC_000962:5-9")
    assert [(s.start, s.end) for s in sites] == [(7570, 7570), (7582, 7582), (761155, 761155), (5, 9)]
    assert all(s.seqid == "NC_000962.3" for s in sites)
    legacy = read_positions(os.path.join(EXAMPLE, "primer_regions.tsv"))
    assert (legacy[0].start, legacy[0].end) == (490900, 491416) and len(legacy) == 8
    mtb = read_positions(os.path.join(EXAMPLE, "mtb_positions.tsv"))
    assert Site("NC_000962.3", 761155, 761155, "rpoB_S450") in mtb
    assert Site("NC_000962.3", 2288681, 2289241, "pncA") in mtb
    f = tmp_path / "r.bed"
    f.write_text("NC_000962.3\t99\t200\tx\n")
    assert read_positions(str(f)) == [Site("NC_000962.3", 100, 200, "x")]
    with pytest.raises(ValueError):
        read_positions("5000000")          # beyond the end of the genome
    with pytest.raises(ValueError):
        read_positions("chrZ:100")         # unknown sequence


def _check_partition(sites, targets, span_hard):
    inserts = [(t.seqid, t.start, t.end) for t in targets]
    assert coverage_gaps(sites, inserts) == []
    for t in targets:
        assert t.span <= span_hard
        assert t.sites, "every target holds at least one site"
    for a, b in zip(targets, targets[1:]):
        if a.seqid == b.seqid:
            assert a.end < b.start


def test_group_single_and_close_positions():
    comfort, hard, gap = span_limits(300, 375)
    sites = [Site("", p, p) for p in (100000, 100050, 100180, 300000)]
    targets = group_sites(sites, comfort, hard, gap)
    assert [(t.start, t.end) for t in targets] == [(100000, 100180), (300000, 300000)]
    _check_partition(sites, targets, hard)
    assert targets[0].id == "T01"


def test_group_splits_at_largest_gap():
    comfort, hard, gap = span_limits(300, 375)
    # 0, 40, 260, 300 cannot be one target (span 301 > hard). Two targets; the
    # cut must fall in the 220 bp gap, not between close neighbours.
    sites = [Site("", p, p) for p in (10000, 10040, 10260, 10300)]
    targets = group_sites(sites, comfort, hard, gap)
    assert [(t.start, t.end) for t in targets] == [(10000, 10040), (10260, 10300)]


def test_group_tiles_long_range():
    comfort, hard, gap = span_limits(300, 375)
    sites = [Site("", 2288681, 2289241, "pncA")]
    targets = group_sites(sites, comfort, hard, gap)
    assert len(targets) == 2  # 561 bp cannot fit one target, fits two
    _check_partition(sites, targets, hard)
    sites = [Site("", 1000, 2999, "long")]
    targets = group_sites(sites, comfort, hard, gap)
    # one extra amplicon is allowed when it keeps every tile near 300 bp
    assert len(targets) in (-(-2000 // hard), -(-2000 // hard) + 1)
    _check_partition(sites, targets, hard)


def test_group_mtb_example_and_multiseq():
    comfort, hard, gap = span_limits(300, 375)
    sites = read_positions(os.path.join(EXAMPLE, "mtb_positions.tsv"))
    targets = group_sites(sites, comfort, hard, gap)
    _check_partition(sites, targets, hard)
    members = [{s.name for s in t.sites} for t in targets]
    assert any({"rpoB_L430", "rpoB_S450", "rpoB_I491"} <= m for m in members)  # one amplicon for the rpoB sites
    assert any({"gyrA_A90", "gyrA_D94"} <= m for m in members)

    multi = [Site("chrA", 100, 100), Site("chrB", 120, 120), Site("chrA", 150, 150)]
    targets = group_sites(multi, comfort, hard, gap)
    assert sorted((t.seqid, t.start, t.end) for t in targets) == [("chrA", 100, 150), ("chrB", 120, 120)]


def test_coverage_gaps_reports_missing_bases():
    sites = [Site("c", 100, 200, "r"), Site("c", 500, 500, "p")]
    gaps = coverage_gaps(sites, [("c", 90, 150), ("c", 160, 220)])
    assert gaps == [(sites[0], [(151, 159)]), (sites[1], [(500, 500)])]
    assert coverage_gaps(sites, [("c", 90, 150), ("c", 151, 220), ("c", 400, 600)]) == []
