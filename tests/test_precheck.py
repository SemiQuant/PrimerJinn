"""Precheck and VCF/BED position input. Offline: small synthetic FASTA files."""
import os
import random

import pytest

from primerJinn.multiplex import cli
from primerJinn.multiplex import genome as gmod
from primerJinn.multiplex.precheck import run as precheck

HDR = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


@pytest.fixture()
def fasta(tmp_path):
    """A 20 kb sequence 'seq1' (random, 60% GC) with an N gap at 9001-9100."""
    rnd = random.Random(7)
    s = "".join(rnd.choices("ACGT", weights=[2, 3, 3, 2], k=20000))
    s = s[:9000] + "N" * 100 + s[9100:]
    p = tmp_path / "g.fa"
    p.write_text(">seq1 test\n" + "\n".join(s[i:i + 80] for i in range(0, len(s), 80)) + "\n")
    return str(p), s


def _levels(pc, level, topic=None):
    return [i.msg for i in pc.items if i.level == level and (topic is None or i.topic == topic)]


def test_vcf_positions_names_refs_and_deletion(tmp_path, fasta):
    path, s = fasta
    v = tmp_path / "t.vcf"
    v.write_text(HDR + f"seq1\t5000\trs1\t{s[4999]}\tA\t.\tPASS\t.\n"
                       f"seq1\t6000\t.\t{s[5999:6002]}\t{s[5999]}\t.\tPASS\t.\n")
    meta = gmod.resolve_genome(path, str(tmp_path / "cache"))
    rows = gmod.parse_positions(str(v), meta)
    assert [(r["name"], r["start"], r["end"]) for r in rows] == [("rs1", 5000, 5000), (f"{s[5999:6002]}6000{s[5999]}", 6000, 6002)]
    assert rows[0]["ref"] == s[4999] and rows[0]["line"] == 3
    pc = precheck(path, str(v), genome_dir=str(tmp_path / "cache"), check_env=False)
    assert pc.ok, pc.text()
    assert any("REF alleles match" in m for m in _levels(pc, "OK", "positions"))
    assert len(pc.targets) == 2


def test_vcf_snv_without_id_is_named_by_change(tmp_path, fasta):
    path, s = fasta
    v = tmp_path / "t.vcf"
    alt = "A" if s[4999] != "A" else "C"
    v.write_text(HDR + f"seq1\t5000\t.\t{s[4999]}\t{alt}\t.\tPASS\t.\n")
    meta = gmod.resolve_genome(path, str(tmp_path / "cache"))
    assert gmod.parse_positions(str(v), meta)[0]["name"] == f"{s[4999]}5000{alt}"


def test_ref_mismatch_majority_is_error_minority_warn(tmp_path, fasta):
    path, s = fasta
    flip = {"A": "C", "C": "G", "G": "T", "T": "A"}
    bad = HDR + "".join(f"seq1\t{p}\t.\t{flip[s[p - 1]]}\tA\t.\tPASS\t.\n" for p in (2000, 4000, 6000))
    v = tmp_path / "bad.vcf"
    v.write_text(bad)
    pc = precheck(path, str(v), genome_dir=str(tmp_path / "c"), check_env=False)
    assert not pc.ok and any("different reference" in m for m in _levels(pc, "ERROR"))
    mixed = HDR + "".join(f"seq1\t{p}\t.\t{s[p - 1]}\tA\t.\tPASS\t.\n" for p in (2000, 4000, 6000)) \
        + f"seq1\t8000\t.\t{flip[s[7999]]}\tA\t.\tPASS\t.\n"
    v.write_text(mixed)
    pc = precheck(path, str(v), genome_dir=str(tmp_path / "c"), check_env=False)
    assert pc.ok and any("1 of 4" in m for m in _levels(pc, "WARN"))


def test_every_bad_row_is_reported_and_unknown_names_grouped(tmp_path, fasta):
    path, s = fasta
    v = tmp_path / "t.vcf"
    v.write_text(HDR + "chrX\t100\t.\tA\tC\t.\t.\t.\n" + "chrX\t200\t.\tA\tC\t.\t.\t.\n"
                 + "seq1\tabc\t.\tA\tC\n" + "seq1\t25000\t.\tA\tC\t.\t.\t.\n"
                 + f"seq1\t5000\t.\t{s[4999]}\tA\t.\t.\t.\n")
    pc = precheck(path, str(v), genome_dir=str(tmp_path / "c"), check_env=False)
    errs = _levels(pc, "ERROR", "positions")
    assert sum("unknown sequence 'chrX' (2 row(s)" in e for e in errs) == 1
    assert any("VCF line 5" in e for e in errs)          # malformed POS
    assert any("beyond the end" in e for e in errs)
    assert any("1 position(s) read" in m for m in _levels(pc, "OK"))   # the good row is still read


def test_vcf_header_contig_length_and_version(tmp_path, fasta):
    path, s = fasta
    v = tmp_path / "t.vcf"
    v.write_text("##fileformat=VCFv4.2\n##contig=<ID=seq1,length=19999>\n" + HDR.split("\n", 1)[1]
                 + f"seq1\t5000\t.\t{s[4999]}\tA\t.\t.\t.\n")
    pc = precheck(path, str(v), genome_dir=str(tmp_path / "c"), check_env=False)
    assert any("19,999 bp" in e for e in _levels(pc, "ERROR"))


def test_bed_and_two_bp_table_hint(tmp_path, fasta):
    path, _ = fasta
    b = tmp_path / "t.bed"
    b.write_text("track name=x\nseq1\t4999\t5000\tsnp\nseq1\t6000\t6200\tregion\n")
    pc = precheck(path, str(b), genome_dir=str(tmp_path / "c"), check_env=False)
    assert pc.ok and [(r["start"], r["end"]) for r in pc.rows] == [(5000, 5000), (6001, 6200)]
    t = tmp_path / "t.tsv"
    t.write_text("name\tstart\tend\na\t4999\t5000\nb\t6999\t7000\nc\t7999\t8000\n")
    pc = precheck(path, str(t), genome_dir=str(tmp_path / "c"), check_env=False)
    assert any("exactly 2 bp" in w for w in _levels(pc, "WARN"))


def test_n_gap_and_sequence_end_warnings(tmp_path, fasta):
    path, _ = fasta
    pc = precheck(path, "9050,50,15000", genome_dir=str(tmp_path / "c"), check_env=False)
    warns = " ".join(_levels(pc, "WARN"))
    assert "N gap" in warns and "end of seq1" in warns
    assert pc.ok


def test_regions_mode(tmp_path, fasta):
    path, s = fasta
    b = tmp_path / "r.bed"
    b.write_text("seq1\t1000\t1500\ta\nseq1\t1400\t1700\tb\nseq1\t3000\t4500\tbig\n")
    pc = precheck(path, str(b), regions_mode=True, check_env=False, product_max=800,
                  genome_dir=str(tmp_path / "c"))
    assert any("overlaps" in e for e in _levels(pc, "ERROR", "regions"))
    assert any("more than --product_size_max 800" in w for w in _levels(pc, "WARN"))
    two = tmp_path / "two.fa"
    two.write_text(f">a\n{s[:5000]}\n>b\n{s[5000:10000]}\n")
    pc = precheck(str(two), "a:100", regions_mode=True, check_env=False, genome_dir=str(tmp_path / "c"))
    assert any("one sequence" in e for e in _levels(pc, "ERROR", "genome"))


def test_output_folder_checks(tmp_path, fasta):
    path, _ = fasta
    ro = tmp_path / "ro"
    ro.mkdir()
    os.chmod(ro, 0o500)
    try:
        if os.access(ro, os.W_OK):
            pytest.skip("running as a user who can write anywhere")
        pc = precheck(path, "5000", out=str(ro / "sub" / "x"), genome_dir=str(tmp_path / "c"), check_env=False)
        assert any("no write permission" in e for e in _levels(pc, "ERROR", "output"))
    finally:
        os.chmod(ro, 0o700)
    pc = precheck(path, "5000", out=str(tmp_path / "nodir" / "x"), regions_mode=True,
                  genome_dir=str(tmp_path / "c"), check_env=False)
    assert any("does not exist" in e for e in _levels(pc, "ERROR", "output"))


def test_cli_check_exit_codes_and_no_design(tmp_path, fasta, capsys):
    path, _ = fasta
    out = tmp_path / "o" / "d"
    base = ["--genome", path, "--genome_dir", str(tmp_path / "c"), "--out", str(out), "--check"]
    assert cli.main(base + ["--positions", "5000"]) == 0
    assert "all checks passed" in capsys.readouterr().out
    assert not list(tmp_path.glob("o/d_*"))               # --check designs nothing
    assert cli.main(base + ["--positions", "seq9:5000,25000"]) == 2
    err = capsys.readouterr().err
    assert "unknown sequence 'seq9'" in err and "beyond the end" in err


def test_read_vcf_strict_and_lenient(tmp_path):
    v = tmp_path / "t.vcf"
    v.write_text(HDR + "s\tx\t.\tA\tC\n" + "s\t5\t.\tA\t<DEL>\t.\t.\tSVTYPE=DEL;END=20\n")
    with pytest.raises(ValueError, match="VCF line 3"):
        gmod.read_vcf(str(v))
    errs = []
    recs, _ = gmod.read_vcf(str(v), errs)
    assert len(errs) == 1 and recs[0]["end"] == 20
    assert gmod.is_vcf(str(v))
    plain = tmp_path / "noext"
    plain.write_text(HDR)
    assert gmod.is_vcf(str(plain))


# ---- regressions found by the adversarial review ---------------------------------------
def test_regions_mode_keeps_old_getmultiprimerset_inputs(tmp_path, fasta):
    path, s = fasta
    uneven = tmp_path / "uneven.fa"           # SeqIO read this; faidx cannot
    uneven.write_text(">seq1\n" + s[:70] + "\n" + s[70:150] + "\n\n" + s[150:] + "\n")
    t = tmp_path / "r.tsv"                    # touching regions, float coordinates, a seqid column
    t.write_text("chrom\tstart\tend\nchr1\t1000.0\t1500.0\nchr1\t1500\t1800\n")
    pc = precheck(str(uneven), str(t), regions_mode=True, check_env=False, product_max=800,
                  genome_dir=str(tmp_path / "c"))
    assert pc.ok, pc.text()
    assert [(r["start"], r["end"]) for r in pc.rows] == [(1000, 1500), (1500, 1800)]
    assert any("rewrapped" in m for m in _levels(pc, "INFO", "genome"))
    assert any("taken as its only sequence" in w for w in _levels(pc, "WARN"))
    assert gmod.fetch(pc.meta, "seq1", 1, 200) == s[:200]


def test_corrupt_gzip_and_bcf_are_reported_not_raised(tmp_path, fasta):
    import gzip
    path, s = fasta
    v = tmp_path / "t.vcf.gz"
    with gzip.open(v, "wt") as fh:
        fh.write(HDR + "".join(f"seq1\t{p}\t.\t{s[p - 1]}\tA\t.\t.\t.\n" for p in range(1000, 9000, 50)))
    trunc = tmp_path / "trunc.vcf.gz"
    trunc.write_bytes(v.read_bytes()[:60])
    pc = precheck(path, str(trunc), variants=str(trunc), genome_dir=str(tmp_path / "c"), check_env=False)
    errs = " ".join(_levels(pc, "ERROR"))
    assert "truncated or corrupt" in errs and any(i.topic == "variants" for i in pc.errors)
    bcf = tmp_path / "x.bcf"
    bcf.write_bytes(b"BCF\x02\x02\x00\x00")
    pc = precheck(path, str(bcf), genome_dir=str(tmp_path / "c"), check_env=False)
    assert any("bcftools view" in e for e in _levels(pc, "ERROR"))


def test_gvcf_blocks_svlen_af_and_long_deletions(tmp_path, fasta):
    from primerJinn.multiplex.variants import load_variants
    path, s = fasta
    meta = gmod.resolve_genome(path, str(tmp_path / "c"))
    g = tmp_path / "s.g.vcf"
    g.write_text(HDR + "seq1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=1099\n"
                 f"seq1\t1100\t.\t{s[1099]}\tG,<NON_REF>\t.\t.\tAF=.\n"
                 "seq1\t1101\t.\tA\t<*>\t.\t.\tEND=4000\n"
                 "seq1\t7000\t.\tA\t<DEL>\t.\t.\tSVLEN=-40\n"
                 "seq1\t7500\t.\tA\tC,G\t.\t.\tAF=0.01,0.4\n"
                 f"seq1\t10000\t.\t{s[9999:16000]}\t{s[9999]}\t.\t.\t.\n"      # 6001-bp deletion, plain REF
                 "seq1\t13000\t.\tA\t<DEL>\t.\t.\tEND=19500\n"
                 "seq1\tbad\n")
    m, notes = load_variants(str(g), meta, min_af=0.1)
    assert not m.within("seq1", 1, 1099) and not m.within("seq1", 1101, 4000)   # reference blocks
    assert m.within("seq1", 1100, 1100) and m.within("seq1", 7500, 7500)         # AF '.' kept; max AF used
    assert m.within("seq1", 7000, 7040) == list(range(7000, 7041))               # SVLEN
    assert len(m.within("seq1", 10000, 16000)) == 6001                           # long plain deletion masked
    assert not m.within("seq1", 17000, 19000)                                    # long symbolic SV skipped...
    assert any("ignored" in n and "longer than" in n for n in notes)              # ...and reported
    assert any("malformed" in n for n in notes) and any("reference block" in n for n in notes)
    rows = gmod.parse_positions(str(g), meta, errors=[])
    assert all(r["start"] >= 1100 for r in rows) and not any(r["start"] in (1, 1101) for r in rows)


def test_header_only_contig_does_not_block(tmp_path, fasta):
    path, s = fasta
    v = tmp_path / "t.vcf"
    v.write_text("##fileformat=VCFv4.2\n##contig=<ID=seq1,length=20000>\n##contig=<ID=chrUn_7,length=5>\n"
                 "##contig=<ID=seq1.2,length=10>\n" + HDR.split("\n", 1)[1] + f"seq1\t5000\t.\t{s[4999]}\tA\t.\t.\t.\n")
    pc = precheck(path, str(v), genome_dir=str(tmp_path / "c"), check_env=False)
    assert pc.ok, pc.text()


def test_quiet_design_still_prints_warnings(tmp_path, fasta, capsys, monkeypatch):
    import primerJinn.multiplex.design as design

    def stop(*a, **k):
        raise ValueError("stopped after the check")
    monkeypatch.setattr(design, "design_pool", stop)
    path, _ = fasta
    base = ["--genome", path, "--genome_dir", str(tmp_path / "c"), "--out", str(tmp_path / "d"), "--quiet"]
    assert cli.main(base + ["--positions", "9050"]) == 2
    out = capsys.readouterr().out
    assert "N gap" in out and "OK " not in out


def test_second_review_regressions(tmp_path, fasta):
    from primerJinn.multiplex.variants import load_variants
    path, s = fasta
    meta = gmod.resolve_genome(path, str(tmp_path / "c"))
    # a numeric-looking name is not a coordinate; integral floats are read in coordinate columns
    t = tmp_path / "p.tsv"
    t.write_text("12.0\t7570\n")
    assert [(r["name"], r["start"]) for r in gmod.parse_positions(str(t), meta)] == [("12.0", 7570)]
    t.write_text("name\tpos\n12.0\t7570.0\nr\t1000.0-1100.0\n")
    assert [(r["name"], r["start"], r["end"]) for r in gmod.parse_positions(str(t), meta)] == \
        [("12.0", 7570, 7570), ("r", 1000, 1100)]
    # <INS> SVLEN is inserted length; a symbolic ALT among several still counts as an SV
    v = tmp_path / "v.vcf"
    v.write_text(HDR + f"seq1\t9000\t.\t{s[8999]}\t<INS>\t.\t.\tSVTYPE=INS;SVLEN=3000\n"
                 f"seq1\t1000\t.\t{s[999]}\tC,<DEL>\t.\t.\tEND=19000\n"
                 f"seq1\t5000\t.\t{s[4999]}\tA\t.\t.\t.\n" + "seq1\tXX\t.\tA\tC\n")
    m, notes = load_variants(str(v), meta)
    assert m.within("seq1", 9000, 12000) == [9000] and not m.within("seq1", 1001, 4000)
    # one malformed line in --variants does not block the design
    pc = precheck(path, "5000", variants=str(v), genome_dir=str(tmp_path / "c"), check_env=False)
    assert pc.ok, pc.text()
    # every record skipped: the error says why
    g = tmp_path / "g.vcf"
    g.write_text(HDR + "seq1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=20000\n")
    with pytest.raises(ValueError, match="reference block"):
        load_variants(str(g), meta)


def test_exclude_positions_by_name_locus_and_target(tmp_path):
    """--exclude drops positions by name, locus or target id before grouping."""
    import os
    from primerJinn.multiplex.config import Config
    from primerJinn.multiplex.precheck import run

    ref = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "example", "ref.fasta")
    pos = tmp_path / "p.tsv"
    pos.write_text("gyrA_A90\t7570\ngyrA_D94\t7582\nrpoB_S450\t761155\nkatG_S315\t2155168\n")
    pc = run(ref, str(pos), Config(), out=str(tmp_path / "x"), exclude=["rpoB_S450", "NC_000962.3:2155168", "T09"])
    assert [s.name for s in pc.sites] == ["gyrA_A90", "gyrA_D94"]
    msgs = " ".join(i.msg for i in pc.items)
    assert "2 position(s) left out" in msgs and "no position or target matches t09" in msgs
    pc2 = run(ref, str(pos), Config(), out=str(tmp_path / "y"), exclude=["T01"])   # gyrA pair is T01
    assert [s.name for s in pc2.sites] == ["rpoB_S450", "katG_S315"]
