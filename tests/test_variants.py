import gzip

import pytest

from primerJinn.multiplex.config import Config
from primerJinn.multiplex.variants import VariantMask, load_variants

META = {"sequences": [{"seqid": "NC_000962.3", "length": 4411532, "aliases": ["NC_000962", "chr1"]}],
        "n_sequences": 1, "default_seqid": "NC_000962.3", "default_unambiguous": True}


def test_vcf_tsv_bed(tmp_path):
    vcf = tmp_path / "v.vcf.gz"
    with gzip.open(vcf, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        fh.write("NC_000962\t7570\t.\tC\tT\t.\tPASS\tAF=0.9\n")
        fh.write("NC_000962.3\t7600\t.\tACG\tA\t.\tPASS\tAF=0.5\n")      # deletion: 7600-7602
        fh.write("chrX\t10\t.\tA\tG\t.\tPASS\tAF=1\n")                   # unknown sequence
        fh.write("NC_000962.3\t9000\t.\tA\tG\t.\tPASS\tAF=0.01\n")       # rare
    m, notes = load_variants(str(vcf), META, min_af=0.05)
    assert m.within("NC_000962.3", 1, 10000) == [7570, 7600, 7601, 7602]
    assert any("chrX" in n for n in notes)
    tsv = tmp_path / "v.tsv"
    tsv.write_text("seqid\tpos\nchr1\t100\nNC_000962.3\t200-202\n")
    m, _ = load_variants(str(tsv), META)
    assert m.within("NC_000962.3", 1, 1000) == [100, 200, 201, 202]
    bed = tmp_path / "v.bed"
    bed.write_text("NC_000962.3\t99\t100\n")
    m, _ = load_variants(str(bed), META)
    assert m.within("NC_000962.3", 1, 1000) == [100]
    with pytest.raises(ValueError):
        load_variants(str(tsv), {**META, "sequences": [{"seqid": "other", "length": 10, "aliases": []}],
                                 "default_seqid": "other"})


def test_primer_hit_distances():
    m = VariantMask({"c": [100, 110]})
    assert sorted(m.primer_hits("F", "c", 91, 110)) == [0, 10]   # F 3' end at 110
    assert sorted(m.primer_hits("R", "c", 100, 119)) == [0, 10]  # R 3' end at 100


def test_candidates_avoid_variants_near_3prime():
    import os

    from primerJinn.multiplex import design
    from primerJinn.multiplex.candidates import design_target
    from primerJinn.multiplex.targets import Site, Target

    ref = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "example", "ref.fasta")
    store = design.SeqStore(ref)
    t = Target("T01", "NC_000962.3", 2155168, 2155168, [Site("NC_000962.3", 2155168, 2155168)])
    region, rs = store.region(t.seqid, t.start - 625, t.end + 625)
    # mask every 7th base around the target: many primers lose their 3' end
    mask = VariantMask({"NC_000962.3": list(range(t.start - 500, t.end + 500, 7))})
    cfg = Config(specificity=False, threads=1, variants=mask).resolve(0.656)
    res = design_target(t, region, rs, store.lengths[t.seqid], cfg, keep_pairs=10)
    assert res["pairs"]
    for pr in res["pairs"]:
        for p in (pr.fwd, pr.rev):
            hits = mask.primer_hits(p.side, p.seqid, p.start, p.end)
            assert all(h >= cfg.variant_3prime_nt for h in hits)
            assert p.n_variants == len(hits)
