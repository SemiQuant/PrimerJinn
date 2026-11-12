"""Regression tests for defects found in the 2026-11-11 validation (docs/research is the evidence base)."""
import pytest

from primerJinn.multiplex import structure, thermo
from primerJinn.multiplex.candidates import Pair, Primer
from primerJinn.multiplex.config import Config
from primerJinn.multiplex.design import classify_product


def _pair(seqid="chr", f=("ACGTACGTACGTACGTAC", 1254, 1271), r=("TTGCATGCATGCATGCAA", 1533, 1550)):
    fp = Primer(f[0], "F", seqid, f[1], f[2], 60.0, 0.5)
    rp = Primer(r[0], "R", seqid, r[1], r[2], 60.0, 0.5)
    return Pair("T01", fp, rp, rp.end - fp.start + 1, 0.5)


def test_designed_product_is_contig_aware():
    pr = _pair()
    seqs = [pr.fwd.seq, pr.rev.seq]
    d = dict(contig="chr", start=1253, end=1550, size=297)
    assert classify_product(d, 0, 1, pr, pr, seqs, "chr") == "designed"
    # identical second copy on another contig is a repeated locus, not the designed product
    d2 = dict(contig="plasmid_copy", start=1253, end=1550, size=297)
    assert classify_product(d2, 0, 1, pr, pr, seqs, "plasmid_copy") == "repeat"


def test_outer_product_between_designed_amplicons():
    a = _pair()
    b = _pair(f=("GGGTACGTACGTACGTAC", 1900, 1917), r=("CCCCATGCATGCATGCAA", 2180, 2197))
    seqs = [a.fwd.seq, b.rev.seq]
    d = dict(contig="chr", start=a.fwd.start - 1, end=b.rev.end, size=b.rev.end - a.fwd.start + 1)
    assert classify_product(d, 0, 1, a, b, seqs, "chr") == "outer"


def test_duplex_over_60_nt_does_not_crash():
    tail = "GACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
    dh, ds = structure.duplex_dh_ds(tail + "ACCATCATAGTTATCGATAAAATCTATTG")
    assert dh < 0 and ds < 0


def test_q5_vendor_rule_uses_reaction_concentration():
    # NEB API q5/0.1 for this pair gives Tm 57/57, Ta 58.2 (validation, E. coli run)
    tm = min(thermo.tm_neb(s, 100.0) for s in ("GCCTGTAAATTCACGAACA", "GTTAGTCTGTGCGTACAC"))
    assert abs(round(tm + 1.0, 1) - 58.2) <= 0.3
    cfg = Config(polymerase="q5").resolve(0.51)
    assert cfg.ta < cfg.tm   # design Ta follows the 100 nM rule, not the 500 nM scale


def test_tail_validation():
    with pytest.raises(ValueError):
        Config(tails="custom").resolve(0.5)
    with pytest.raises(ValueError):
        Config(tails="illumina", tail_fwd="ACGT").resolve(0.5)
    with pytest.raises(ValueError):
        Config(tails="custom", tail_fwd="ACGN").resolve(0.5)


def test_at_rich_defaults():
    cfg = Config().resolve(0.193)
    assert cfg.max_poly == 5 and cfg.primer_gc_min < 0.15
    assert Config().resolve(0.656).primer_gc_min == 0.20


def test_plasmodb_chromosome_names():
    from primerJinn.multiplex import genome as G
    meta = {"sequences": [{"seqid": f"NC_0043{i:02d}.1", "length": 3_000_000, "aliases": [], "type": "Chromosome",
                           "name": str(i)} for i in range(1, 15)],
            "n_sequences": 14, "default_seqid": "NC_004301.1", "default_unambiguous": False}
    r = G.parse_positions(["Pf3D7_13_v3:1725259"], meta)
    assert r[0]["seqid"] == "NC_004313.1"
    with pytest.raises(ValueError, match="Chromosome names"):
        G.parse_positions(["Pf3D7_99_v3:100"], meta)
