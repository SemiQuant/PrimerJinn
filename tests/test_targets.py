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


# ----------------------------------------------------------------------------
# junction amplicons for long deletions (2026-11-12)
# ----------------------------------------------------------------------------
def test_vcf_kind_and_snps_only(tmp_path):
    from primerJinn.multiplex.genome import vcf_kind, parse_positions

    r = lambda pos, ref, alt, end=None, info=".": dict(line=1, chrom="c", pos=pos, id=".", ref=ref, alt=alt,
                                                       filter=".", info=info, end=end or pos + len(ref) - 1,
                                                       ref_block=False)
    assert vcf_kind(r(10, "A", "G")) == "snp"
    assert vcf_kind(r(10, "A", "G,T")) == "snp"
    assert vcf_kind(r(10, "AC", "A")) == "indel"
    assert vcf_kind(r(10, "A" * 200, "A"), 100) == "long_deletion"
    assert vcf_kind(r(10, "A" * 200, "A"), 0) == "indel"           # threshold off
    assert vcf_kind(r(10, "A" * 200, "A"), None) == "indel"
    assert vcf_kind(r(10, "N", "<DEL>", end=500), 100) == "long_deletion"
    assert vcf_kind(r(10, "N", "<INV>", end=500), 100) == "sv"
    assert vcf_kind(r(10, "A", "<INS>"), 100) == "sv"
    vcf = tmp_path / "x.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                   "c\t100\ts1\tA\tG\t.\t.\t.\n"
                   "c\t300\td1\t" + "A" * 151 + "\tA\t.\t.\t.\n"
                   "c\t600\ti1\tAC\tA\t.\t.\t.\n")
    rows = parse_positions(str(vcf), None, long_variant=100)
    assert [(r["start"], r["end"], r.get("gap")) for r in rows] == [(100, 100, None), (300, 451, (301, 450)),
                                                                     (600, 601, None)]
    rows = parse_positions(str(vcf), None, snps_only=True)
    assert [r["name"] for r in rows] == ["s1"]
    rows = parse_positions(str(vcf), None, long_variant=0)
    assert rows[1]["end"] == 450 and "gap" not in rows[1]


def test_junction_site_and_targets():
    from primerJinn.multiplex.config import Config
    from primerJinn.multiplex.targets import Site, build_targets, junction_targets

    with pytest.raises(ValueError):
        Site("c", 100, 110, gap=(100, 105))          # the gap must lie strictly inside
    s = Site("c", 300, 451, "del", gap=(301, 450))
    assert s.gap_length == 150 and "deletion c:301-450, 150 bp" in s.label()
    sites = [Site("c", 1000, 1000, "snp"), s, Site("c", 5000, 5000, "snp2"), Site("c", 9000, 9601, "d2", gap=(9001, 9600))]
    cfg = Config(specificity=False).resolve(0.5)
    targets, windows, notes = build_targets(sites, cfg)
    ids = [t.id for t in targets]
    assert ids[:2] == ["T01", "T02"] and ids[2:] == ["J01", "J02"]
    j = targets[2]
    assert j.gap == (301, 450) and j.start == 300 and j.end == 451 and j.sites == [s]
    assert any("junction" in n for n in notes)
    assert junction_targets([Site("c", 10, 10)]) == []


def test_conflict_graph_sees_targets_inside_a_long_junction_span():
    """A junction target spans its deletion; targets starting inside it must conflict with it."""
    from primerJinn.multiplex.config import Config
    from primerJinn.multiplex.optimize import conflict_graph
    from primerJinn.multiplex.targets import Target

    cfg = Config().resolve(0.65)
    j = Target("J01", "c", 10000, 12401, gap=(10001, 12400))
    inside = Target("T01", "c", 12100, 12100)
    just_after = Target("T02", "c", 12410, 12410)
    far = Target("T03", "c", 20000, 20000)
    e = conflict_graph([j, inside, just_after, far], cfg)
    assert 1 in e[0] and 2 in e[0] and 3 not in e[0]


def test_vcf_max_indel_filter(tmp_path):
    from primerJinn.multiplex.genome import parse_positions, vcf_indel_length

    vcf = tmp_path / "x.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                   "c\t100\ts1\tA\tG\t.\t.\t.\n"
                   "c\t200\ti5\tACGTAC\tA\t.\t.\t.\n"
                   "c\t300\td60\t" + "A" * 61 + "\tA\t.\t.\t.\n"
                   "c\t400\tsv\tN\t<DEL>\t.\t.\tEND=900\n"
                   "c\t950\tins\tA\tACCCCCCCCCCC\t.\t.\t.\n")
    names = lambda rows: [r["name"] for r in rows]
    assert names(parse_positions(str(vcf), None)) == ["s1", "i5", "d60", "sv", "ins"]
    assert names(parse_positions(str(vcf), None, max_indel=10)) == ["s1", "i5"]
    assert names(parse_positions(str(vcf), None, max_indel=60)) == ["s1", "i5", "d60", "ins"]
    assert names(parse_positions(str(vcf), None, max_indel=20, long_variant=50)) == ["s1", "i5", "ins"]
    assert vcf_indel_length(dict(pos=1, end=1, ref="A", alt="G,T")) == 0
