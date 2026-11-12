import pytest

pytest.importorskip("streamlit")

from primerJinn.multiplex import gui  # noqa: E402
from primerJinn.multiplex.cli import build_parser  # noqa: E402


def test_build_args_parse_with_cli():
    opts = dict(polymerase="q5", primer_nM=50, tails="custom", tail_fwd="ACGT", tail_rev="TTGA", amp_min=200,
                amp_opt=300, amp_max=380, tm=62.0, tm_window=1.5, mono_mM=60, mg_mM=None, dntp_mM=0.8, dmso=2,
                ta=None, outer_min=900, balance=False, max_ext=6, specificity=False, threads=4, variants="v.vcf")
    a = build_parser().parse_args(gui.build_args("g.fa", "p.txt", "out/x", opts))
    assert (a.polymerase, a.primer_nM, a.tails, a.tail_fwd, a.amp_max, a.tm, a.mg_mM, a.dmso) == \
        ("q5", 50.0, "custom", "ACGT", 380, 62.0, None, 2.0)
    assert a.no_balance and a.no_specificity and a.variants == "v.vcf" and a.threads == 4 and a.ta is None


def test_sequence_html_and_positions():
    h = gui.sequence_html("gtcACGT<")
    assert "fef3c7" in h and "gtc" in h and "&lt;" in h
    assert gui.requested_positions("a (NC_1:100); b (NC_1:200-210, part of 150-300)") == [(100, 100), (200, 210)]


def test_report_sections():
    s = gui.report_sections("# T\n\n## Coverage\n\nAll 3 covered.\n## Warnings\n\n- x\n")
    assert s["Coverage"] == "All 3 covered." and s["Warnings"] == "- x"
