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


def test_template_dangling_ends_match_primer3_thal():
    """Our SantaLucia-Hicks 2004 Table 3 dangling-end term should explain what thal
    adds when a primer binds inside a longer template (orientation check)."""
    import os
    import random
    import statistics

    import primer3
    from Bio import SeqIO

    from primerJinn.multiplex.candidates import Primer, revcomp, template_dangles

    ref = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "example", "ref.fasta")
    g = str(SeqIO.read(ref, "fasta").seq)

    def sub(a, b):
        return g[a - 1:b] if a >= 1 else ""

    kw = dict(mv_conc=1000, dv_conc=0, dntp_conc=0, dna_conc=50, temp_c=37)
    rng = random.Random(3)
    res = []
    for _ in range(80):
        s = rng.randint(1000, 4_000_000)
        e = s + 19
        for side in "FR":
            if side == "F":
                seq, flank, bare = sub(s, e), revcomp(sub(s - 1, e + 1)), revcomp(sub(s, e))
            else:
                seq, flank, bare = revcomp(sub(s, e)), sub(s - 1, e + 1), sub(s, e)
            dh, ds = template_dangles(Primer(seq, side, "x", s, e, 0, 0), sub)
            thal = (primer3.calc_heterodimer(seq, flank, **kw).dg - primer3.calc_heterodimer(seq, bare, **kw).dg) / 1000
            res.append(thal - (dh - 310.15 * ds / 1000))
    assert abs(statistics.mean(res)) < 0.1
    assert statistics.pstdev(res) < 0.3


def test_balance_extension_geometry():
    from primerJinn.multiplex.balance import _extended, display
    from primerJinn.multiplex.candidates import Primer, revcomp

    class Store:
        seq = "AAACCCGGGTTTACGTACGTACGTACGTAGCTAGCTAGCTAGGATCCAAATTT"
        lengths = {"c": len(seq)}

        def region(self, sid, a, b):
            return self.seq[a - 1:b], a

    st = Store()
    f = Primer(st.seq[10:30], "F", "c", 11, 30, 60.0, 0.5)
    seq, s, e = _extended(f, 3, st)
    assert (s, e) == (8, 30) and seq == st.seq[7:30] and seq.endswith(f.seq)
    r = Primer(revcomp(st.seq[20:40]), "R", "c", 21, 40, 60.0, 0.5)
    seq, s, e = _extended(r, 4, st)
    assert (s, e) == (21, 44) and seq == revcomp(st.seq[20:44]) and seq.endswith(r.seq)
    r.seq, r.ext5 = seq, 4
    assert display(r) == seq[:4].lower() + seq[4:]
    assert _extended(f, 11, st) is None   # cannot run off the sequence start


def test_resolve_is_idempotent_with_preset_tails():
    c = Config(tails="illumina").resolve(0.5)
    tail = c.tail_fwd
    assert c.resolve(0.5) is c and c.tail_fwd == tail


def test_rescore_loop_copies_config_and_survives_tails(monkeypatch):
    """The re-score pass used to crash with preset tails (second resolve) and
    mutated the caller's Config."""
    from types import SimpleNamespace

    from primerJinn.multiplex import design

    seen = []

    def fake_once(sites, meta, cfg, verbose):
        cfg.resolve(0.5)
        seen.append(cfg.ta)
        rec = 52.0 if len(seen) == 1 else cfg.ta   # first pass recommends a different Ta
        return SimpleNamespace(ta={"recommended": rec}, timings={"total": 1.0})

    monkeypatch.setattr(design, "_design_once", fake_once)
    user = Config(tails="illumina")
    d = design.design_pool([], {}, user, verbose=False)
    assert seen[1] == 52.0 and len(seen) == 2 and d.passes == 2
    assert user.ta is None and not user._resolved      # caller's Config untouched


def test_no_amplicon_crosses_an_N_gap():
    import random

    from primerJinn.multiplex.candidates import design_target
    from primerJinn.multiplex.targets import Site, Target

    rng = random.Random(5)
    seq = "".join(rng.choice("ACGT") for _ in range(3000))
    seq = seq[:1500] + "N" * 40 + seq[1540:]
    cfg = Config(tm=60.0, tm_scale="physical", specificity=False, threads=1).resolve(0.5)
    t = Target("T01", "c", 1480, 1480, [Site("c", 1480, 1480)])
    res = design_target(t, seq, 1, len(seq), cfg, keep_pairs=20)
    for pr in res["pairs"]:
        assert "N" not in seq[pr.fwd.start - 1:pr.rev.end]


def test_bad_options_fail_early():
    for kw in (dict(primer_nM=0), dict(primer_min=30, primer_max=20), dict(mg_mM=0, mono_mM=0), dict(dmso_pct=-1)):
        with pytest.raises(ValueError):
            Config(**kw).resolve(0.5)
