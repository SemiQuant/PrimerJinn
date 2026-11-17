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


def test_build_args_sequencing_round_trip():
    opts = dict(sequencing="counts.tsv", previous="old_primers.tsv", conc="conc.tsv", replace="T01,T05", keep="T02",
                redesign_below=0.25, reaction_ul=20.0, stock_uM=100.0, no_redesign=True)
    a = gui.build_args("g.fa", "p.txt", "out/x", opts)
    pairs = dict(zip(a[::2], a[1::2]))
    assert pairs["--sequencing"] == "counts.tsv" and pairs["--previous"] == "old_primers.tsv"
    assert pairs["--conc"] == "conc.tsv" and pairs["--replace"] == "T01,T05" and pairs["--keep"] == "T02"
    assert pairs["--redesign_below"] == "0.25" and "--offtarget_max" not in pairs
    assert pairs["--reaction_ul"] == "20.0" and pairs["--stock_uM"] == "100.0" and "--no_redesign" in a
    # nothing of it without the keys, and no --genome without a genome
    b = gui.build_args(None, None, "out/x", {})
    assert not any(x.startswith("--seq") or x in ("--no_redesign", "--genome", "--positions") for x in b)
    parser = build_parser()
    if any("--sequencing" in act.option_strings for act in parser._actions):
        p = parser.parse_args(gui.build_args("g.fa", "p.txt", "out/x", opts))
        assert p.sequencing == "counts.tsv" and p.no_redesign and p.redesign_below == 0.25


def test_sequence_html_and_positions():
    h = gui.sequence_html("gtcACGT<")
    assert 'class="sq-ext"' in h and "gtc" in h and "&lt;" in h and "fef3c7" not in h
    assert gui.requested_positions("a (NC_1:100); b (NC_1:200-210, part of 150-300)") == [(100, 100), (200, 210)]


def test_report_sections():
    s = gui.report_sections("# T\n\n## Coverage\n\nAll 3 covered.\n## Warnings\n\n- x\n")
    assert s["Coverage"] == "All 3 covered." and s["Warnings"] == "- x"


def test_header_html_escapes():
    h = gui.header_html("<b>x</b> & y", "sub <i>", logo="data:image/png;base64,AAAA")
    assert "&lt;b&gt;x&lt;/b&gt; &amp; y" in h and "sub &lt;i&gt;" in h and "<b>x</b>" not in h
    assert 'class="sq-header"' in h and "sq-badge-offline" in h and "All tools" in h and "https://drdx.me/" in h
    assert 'src="data:image/png;base64,AAAA"' in h
    assert gui.header_html("t", "s").count("data:image/png;base64,") == 1  # default logo from static/


def test_style_css_substitutes_placeholders_and_has_tokens():
    css = gui.style_css()
    import re
    assert not re.search(r"\{\{[a-z0-9-]+\}\}", css)
    assert css.count("data:font/woff2;base64,") == 5
    for tok in ("--sq-body: rgb(39, 43, 48)", "--sq-header: rgb(70, 80, 90)", "--sq-accent: rgb(198, 253, 168)",
                "--sq-heading: #a9dff9", "--sq-warning: rgb(240, 80, 210)", "--sq-focus:", "--sq-box-radius: 5px"):
        assert tok in css
    assert 'data-testid="stSidebar"' in css and "prefers-reduced-motion" in css and ":focus-visible" in css


def test_pill_and_table_html():
    assert gui.pill_html("replace") == '<span class="sq-pill sq-pill-danger">replace</span>'
    assert "sq-pill-info" in gui.pill_html("Raise") and "sq-pill-warning" in gui.pill_html("lower")
    assert "sq-pill-muted" in gui.pill_html("keep") and "sq-pill-muted" in gui.pill_html("<odd>")
    assert "&lt;odd&gt;" in gui.pill_html("<odd>")
    pd = pytest.importorskip("pandas")
    df = pd.DataFrame({"name": ["a<b", "c"], "x": [1.5, float("nan")], "n": [3, 4], "act": ["<i>", "<j>"]})
    t = gui.table_html(df, raw_cols=("act",))
    assert "a&lt;b" in t and "<i>" in t and 'class="sq-right">1.50<' in t and "sq-na" in t and ">3<" in t
    assert t.startswith('<div class="sq-table-wrap"><table class="sq-table">')


def test_stat_tiles_and_footer():
    h = gui.stat_tiles_html([("Total reads", "12 345", ""), ("Bad", "<x>", "bad")])
    assert "sq-stat-bad" in h and "&lt;x&gt;" in h and h.count("sq-stat-label") == 2
    assert gui.footer_html().startswith('<footer class="sq-footer">DrDxMe offline tools')
