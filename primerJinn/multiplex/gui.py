"""primerJinn multiplex designer: a small local web GUI in the DrDx.Me dashboard style.

Start with `primerjinn-gui` (after `pip install .[gui]`) or
`streamlit run primerJinn/multiplex/gui.py`. The design itself runs the
command-line tool in a separate process, so results are identical to the CLI.
"""
from __future__ import annotations

import base64
import functools
import html
import io
import os
import re
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

# Streamlit runs this file as a script: make the package importable either way
_ROOT = Path(__file__).resolve().parents[2]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

STATIC = Path(__file__).resolve().parent / "static"
HOME = "https://drdx.me/"
OFFLINE_TITLE = ("Everything runs on this computer. No data leaves it and no internet connection is needed "
                 "except for genome lookup by name or accession.")
THEME = {"base": "dark", "primaryColor": "#c6fda8", "backgroundColor": "#272b30",
         "secondaryBackgroundColor": "#343e48", "textColor": "#cdcdcd"}
PLOT_RANGE = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#008300", "#9085e9", "#e66767"]

POLY = {"kapa": "Generic Taq master mix (default)", "taq": "Standard Taq (1.5 mM MgCl2)",
        "q5": "NEB Q5", "hifi": "Proofreading high-fidelity (model Ta rule)"}
TAILS = {"none": "None", "illumina": "Illumina TruSeq partial adapters", "nextera": "Illumina Nextera overhangs",
         "custom": "Custom"}
PILL_CLASS = {"keep": "muted", "raise": "info", "lower": "warning", "replace": "danger"}
# feedback.tsv rows shown as stat tiles: (item name starts with, tile label, unit)
FEEDBACK_TILES = (("total reads", "Total reads", ""), ("designed reads", "Designed reads", " %"),
                  ("unaligned or off-target", "Unaligned or off-target", " %"),
                  ("amplicons within 2x", "Within 2x of median", ""),
                  ("targets replaced", "Targets replaced", ""))
SEQ_KEYS = ("sequencing", "previous", "conc", "replace", "keep", "redesign_below", "conc_max", "reaction_ul",
            "stock_uM", "no_redesign")


# ----------------------------------------------------------------------------
# helpers (no Streamlit needed; tested in tests/test_gui.py)
# ----------------------------------------------------------------------------
def build_args(genome: str, positions: str, out_prefix: str, opts: dict) -> list:
    """Command-line arguments for primerJinn.multiplex.cli from GUI settings."""
    a = ["--out", out_prefix]
    if genome:
        a += ["--genome", genome]
    if positions:
        a += ["--positions", positions]
    a += ["--polymerase", opts.get("polymerase", "kapa"), "--primer_nM", str(opts.get("primer_nM", 100.0))]
    for key in ("amp_min", "amp_opt", "amp_max", "mono_mM", "dntp_mM", "dmso", "tm_window", "max_ext", "outer_min"):
        if opts.get(key) is not None:
            a += [f"--{key}", str(opts[key])]
    for key in ("tm", "mg_mM", "ta"):
        if opts.get(key) not in (None, "", 0, 0.0):
            a += [f"--{key}", str(opts[key])]
    tails = opts.get("tails", "none")
    a += ["--tails", tails]
    if tails == "custom":
        if opts.get("tail_fwd"):
            a += ["--tail_fwd", opts["tail_fwd"]]
        if opts.get("tail_rev"):
            a += ["--tail_rev", opts["tail_rev"]]
    if opts.get("variants"):
        a += ["--variants", opts["variants"]]
    if not opts.get("specificity", True):
        a.append("--no_specificity")
    if not opts.get("balance", True):
        a.append("--no_balance")
    if opts.get("threads"):
        a += ["--threads", str(opts["threads"])]
    if opts.get("one_pool"):
        a += ["--one_pool", "--one_pool_max", str(opts.get("one_pool_max", 1000))]
    if opts.get("vcf_snps_only"):
        a.append("--vcf_snps_only")
    if opts.get("vcf_max_indel"):
        a += ["--vcf_max_indel", str(int(opts["vcf_max_indel"]))]
    if opts.get("one_pool_drop"):
        a.append("--one_pool_drop")
    if opts.get("exclude"):
        ex = opts["exclude"]
        a += ["--exclude", ",".join(str(x) for x in ex) if isinstance(ex, (list, tuple)) else str(ex)]
    if opts.get("long_variant") is not None:
        a += ["--long_variant", str(int(opts["long_variant"]))]
    # from sequencing (rebalance and redesign)
    if opts.get("sequencing"):
        a += ["--sequencing", str(opts["sequencing"])]
    for key in ("previous", "conc", "replace", "keep"):
        if opts.get(key):
            a += [f"--{key}", str(opts[key])]
    for key in ("redesign_below", "conc_max", "reaction_ul", "stock_uM"):
        if opts.get(key) is not None:
            a += [f"--{key}", str(opts[key])]
    if opts.get("no_redesign"):
        a.append("--no_redesign")
    return a


def sequence_html(seq: str) -> str:
    """Primer sequence with lowercase (5' extension) bases highlighted."""
    out = []
    for run in re.findall(r"[a-z]+|[^a-z]+", seq):
        if run.islower():
            out.append(f'<span class="sq-ext">{html.escape(run)}</span>')
        else:
            out.append(html.escape(run))
    return f'<span class="sq-mono">{"".join(out)}</span>'


def pill_html(action: str) -> str:
    """A coloured status pill for a rebalance action (keep, raise, lower, replace)."""
    key = str(action or "").strip().lower()
    return f'<span class="sq-pill sq-pill-{PILL_CLASS.get(key, "muted")}">{html.escape(key or "n/a")}</span>'


def table_html(df, raw_cols=(), float_fmt="{:.2f}") -> str:
    """A pandas frame as an .sq-table (numbers right-aligned, NaN muted). raw_cols hold HTML already."""
    import pandas as pd

    numeric = {c for c in df.columns if pd.api.types.is_numeric_dtype(df[c]) and c not in raw_cols}
    head = "".join(f'<th{" class=sq-right" if c in numeric else ""}>{html.escape(str(c))}</th>' for c in df.columns)
    rows = []
    for _, r in df.iterrows():
        cells = []
        for c in df.columns:
            v = r[c]
            if c in raw_cols:
                cells.append(f"<td>{v}</td>")
            elif v is None or (isinstance(v, float) and v != v):
                cells.append('<td class="sq-na">n/a</td>')
            elif c in numeric:
                s = float_fmt.format(v) if isinstance(v, float) and not float(v).is_integer() else str(int(v)) \
                    if isinstance(v, float) else str(v)
                cells.append(f'<td class="sq-right">{html.escape(s)}</td>')
            else:
                cells.append(f"<td>{html.escape(str(v))}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    return (f'<div class="sq-table-wrap"><table class="sq-table"><thead><tr>{head}</tr></thead>'
            f'<tbody>{"".join(rows)}</tbody></table></div>')


def stat_tiles_html(tiles) -> str:
    """Stat tiles: an iterable of (label, value, status) with status in good, warn, bad, muted or ''."""
    out = []
    for label, value, status in tiles:
        cls = f" sq-stat-{status}" if status else ""
        out.append(f'<div class="sq-stat{cls}"><div class="sq-stat-label">{html.escape(str(label))}</div>'
                   f'<div class="sq-stat-value">{html.escape(str(value))}</div></div>')
    return f'<div class="sq-stats">{"".join(out)}</div>'


def requested_positions(covers: str) -> list:
    """(start, end) pairs from an amplicons-table 'covers' cell."""
    out = []
    for m in re.finditer(r":(\d+)(?:-(\d+))?(?:,|\))", covers or ""):
        a = int(m.group(1))
        out.append((a, int(m.group(2) or a)))
    return out


def report_sections(md: str) -> dict:
    """Split the markdown report into {heading: text}."""
    parts, cur, buf = {}, "_top", []
    for line in md.splitlines():
        if line.startswith("## "):
            parts[cur] = "\n".join(buf).strip()
            cur, buf = line[3:].strip(), []
        else:
            buf.append(line)
    parts[cur] = "\n".join(buf).strip()
    return parts


def zip_outputs(prefix: str) -> bytes:
    bio = io.BytesIO()
    base = Path(prefix)
    with zipfile.ZipFile(bio, "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(base.parent.glob(base.name + "*")):
            if f.is_file():
                z.write(f, f.name)
    return bio.getvalue()


def version() -> str:
    try:
        from importlib.metadata import version as v
        return v("primerJinn")
    except Exception:
        return ""


# ----------------------------------------------------------------------------
# look: fonts, logo and stylesheet (static/), header and footer markup
# ----------------------------------------------------------------------------
@functools.lru_cache(maxsize=None)
def data_uri(name: str) -> str:
    mime = {"woff2": "font/woff2", "png": "image/png"}[name.rsplit(".", 1)[-1]]
    return f"data:{mime};base64," + base64.b64encode((STATIC / name).read_bytes()).decode("ascii")


FONT_FILES = {"montserrat-600": "montserrat-latin-600-normal.woff2", "montserrat-700": "montserrat-latin-700-normal.woff2",
              "roboto-400": "roboto-latin-400-normal.woff2", "roboto-500": "roboto-latin-500-normal.woff2",
              "roboto-mono-400": "roboto-mono-latin-400-normal.woff2"}


@functools.lru_cache(maxsize=None)
def style_css() -> str:
    """static/style.css with every {{placeholder}} replaced by a data: URI (works offline)."""
    css = (STATIC / "style.css").read_text()
    for key, fname in FONT_FILES.items():
        css = css.replace("{{%s}}" % key, data_uri(fname))
    return css.replace("{{logo}}", data_uri("DrDxMe_white.png"))


def inject_style(st) -> None:
    st.markdown(f"<style>{style_css()}</style>", unsafe_allow_html=True)


_ICON_HOME = ('<svg class="sq-icon" width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" '
              'stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">'
              '<path d="M3 11l9-8 9 8"/><path d="M5 10v10h5v-6h4v6h5V10"/></svg>')
_ICON_CHECK = ('<svg class="sq-icon" width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" '
               'stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">'
               '<path d="M5 12l5 5L20 7"/></svg>')


def header_html(title: str, subtitle: str, logo: str | None = None) -> str:
    """The DrDx.Me .sq-header bar: logo, tool name and subtitle, All tools link, OFFLINE pill."""
    logo = data_uri("DrDxMe_white.png") if logo is None else logo
    return (
        '<header class="sq-header">'
        f'<a class="sq-brand" href="{HOME}" target="_blank" rel="noopener" title="DrDx.Me">'
        f'<img class="sq-logo" src="{logo}" alt="DrDxMe" width="128" height="36"></a>'
        f'<div class="sq-titles"><span class="sq-tool-name">{html.escape(title)}</span>'
        f'<span class="sq-tool-sub">{html.escape(subtitle)}</span></div>'
        '<div class="sq-header-actions">'
        f'<a class="sq-hbtn" href="{HOME}" target="_blank" rel="noopener" title="More DrDx.Me tools">{_ICON_HOME}'
        '<span>All tools</span></a>'
        f'<span class="sq-badge-offline" title="{html.escape(OFFLINE_TITLE)}">{_ICON_CHECK}<span>Offline</span></span>'
        '</div></header>')


def footer_html() -> str:
    v = version()
    return (f'<footer class="sq-footer">DrDxMe offline tools · works without internet · primerJinn'
            f'{" " + html.escape(v) if v else ""}</footer>')


def enable_plot_theme(alt) -> None:
    """The dashboard's dark plot theme for Altair (registered once per process)."""
    config = {"config": {
        "background": "rgb(52,62,72)", "font": "Roboto, Helvetica, Arial, sans-serif",
        "view": {"stroke": "transparent"},
        "axis": {"labelColor": "rgb(215,215,215)", "titleColor": "rgb(215,215,215)", "gridColor": "rgba(255,255,255,.08)",
                 "domainColor": "rgba(255,255,255,.25)", "tickColor": "rgba(255,255,255,.25)", "labelFontSize": 12,
                 "titleFontSize": 13, "titleFontWeight": 500},
        "legend": {"labelColor": "rgb(215,215,215)", "titleColor": "rgb(215,215,215)", "labelFontSize": 12,
                   "orient": "bottom", "direction": "horizontal", "columns": 4, "labelLimit": 260},
        "title": {"color": "rgb(240,240,240)", "fontSize": 15, "fontWeight": 600},
        "header": {"labelColor": "rgb(215,215,215)", "titleColor": "rgb(215,215,215)"},
        "range": {"category": PLOT_RANGE, "ordinal": {"scheme": "blues"}},
    }}
    try:  # Altair 5.5+
        if "sq_dark" not in alt.theme.names():
            alt.theme.register("sq_dark", enable=True)(lambda: config)
        alt.theme.enable("sq_dark")
    except AttributeError:  # Altair 5.0 to 5.4
        alt.themes.register("sq_dark", lambda: config)
        alt.themes.enable("sq_dark")


# ----------------------------------------------------------------------------
# the app
# ----------------------------------------------------------------------------
def app():
    import altair as alt
    import pandas as pd
    import streamlit as st

    st.set_page_config(page_title="primerJinn multiplex", page_icon="🧬", layout="wide",
                       initial_sidebar_state="expanded")
    inject_style(st)
    enable_plot_theme(alt)
    st.markdown(header_html("primerJinn", "Multiplex amplicon-sequencing pool designer"), unsafe_allow_html=True)
    ss = st.session_state
    if "workdir" not in ss:
        ss.workdir = tempfile.mkdtemp(prefix="primerjinn_gui_")
        ss.run = 0
        preload = os.environ.get("PRIMERJINN_GUI_RESULT")  # screenshots: show an existing result prefix
        if preload and (Path(preload + "_report.md").exists() or Path(preload + "_rebalance.tsv").exists()):
            ss.result = preload
            if Path(preload + "_primers.tsv").exists():
                ss.last_primers = preload + "_primers.tsv"
    work = Path(ss.workdir)

    st.title("Multiplex pool design")
    st.markdown('<p class="sq-intro">Design a multiplex amplicon-sequencing primer pool from a genome and a list of '
                'positions. Every position is covered; primers are matched on Tm and on binding at the annealing '
                'temperature, including template structure. After a sequencing run, rebalance the pool or redesign '
                'the primers that failed.</p>', unsafe_allow_html=True)

    # ------------------------------------------------------------------ inputs
    with st.sidebar:
        st.markdown('<div class="sq-sidebar-title">Inputs</div>', unsafe_allow_html=True)
        st.header("Genome")
        gmode = st.radio("Source", ["Name or accession", "Upload FASTA"], horizontal=True, label_visibility="collapsed")
        genome = None
        if gmode == "Name or accession":
            genome = st.text_input("Organism name, assembly (GCF_/GCA_) or nucleotide accession",
                                   value="Mycobacterium tuberculosis H37Rv").strip() or None
        else:
            up = st.file_uploader("FASTA (.fa, .fasta, .fna, optionally .gz)",
                                  type=["fa", "fasta", "fna", "gz"])
            if up is not None:
                path = work / up.name
                path.write_bytes(up.getvalue())
                genome = str(path)

        st.header("Positions")
        st.caption("1-based. One per line or comma separated: 7570, NC_000962.3:7570, 761095-761277, "
                   "or a table, BED or VCF file.")
        ptext = st.text_area("Positions", height=150, label_visibility="collapsed",
                             placeholder="gyrA_A90\t7570\nrpoB_S450\t761155\npncA\t2288681-2289241")
        pfile = st.file_uploader("or upload a positions file (table, BED or VCF)",
                                 type=["tsv", "csv", "txt", "xlsx", "bed", "vcf", "gz"])
        vcf_snps_only = st.checkbox("VCF: SNPs only (drop indels and structural variants)", value=False)
        # options are keyed by the target's positions (seqid:start-end|...), not by its id: ids renumber
        # once something is left out, positions do not
        ex_sel = st.multiselect("Leave out targets (run Check positions first, then pick)",
                                options=[t[0] for t in ss.get("check_targets", [])],
                                format_func=lambda k: dict(ss.get("check_targets", [])).get(k, k))
        ex_text = st.text_input("Leave out positions by name or locus (comma separated)", value="").strip()
        exclude = [loc for key in ex_sel for loc in key.split("|")] \
            + [x.strip() for x in ex_text.split(",") if x.strip()]
        long_variant = st.number_input("VCF: deletions longer than this (bp) get a junction amplicon, 0 = tile them",
                                       0, 100000, 100, 10)

        st.header("Reaction")
        polymerase = st.selectbox("Polymerase", list(POLY), format_func=POLY.get)
        primer_nM = st.number_input("Each primer (nM)", 5.0, 2000.0, 100.0, 5.0)
        tails = st.selectbox("Adapter tails", list(TAILS), format_func=TAILS.get)
        tail_fwd = tail_rev = ""
        if tails == "custom":
            tail_fwd = st.text_input("Forward tail (5'->3')").strip().upper()
            tail_rev = st.text_input("Reverse tail (5'->3')").strip().upper()
        vfile = st.file_uploader("Known strain variants (optional): VCF, BED or positions",
                                 type=["vcf", "gz", "bed", "tsv", "csv", "txt"])

        st.header("From sequencing")
        st.caption("Count the reads per amplicon with your own pipeline and upload one table: columns "
                   "`amplicon` (the target id or primer name without _F/_R) and `reads`, plus a row named "
                   "`unaligned` with the reads that match no amplicon. Add a `conc_nM` column with the "
                   "concentration each pair was used at, so a weak amplicon is raised rather than redesigned "
                   "when its primers were dilute. See the README.")
        sfile = st.file_uploader("Read counts per amplicon (TSV or CSV: amplicon, reads)",
                                 type=["tsv", "txt", "csv"])
        last_primers = ss.get("last_primers")
        pprev = st.file_uploader("Previous design (_primers.tsv)", type=["tsv", "txt"])
        if last_primers and pprev is None:
            st.caption(f"Using the primers designed in this session ({Path(last_primers).name}) unless you upload "
                       "another table.")
        cfile = st.file_uploader("Concentrations used per primer (optional, instead of a conc_nM column): "
                                 "primer, conc_nM", type=["tsv", "txt", "csv"])
        replace = st.text_input("Replace these targets anyway (comma separated, optional)").strip()
        keep = st.text_input("Keep these targets anyway (comma separated, optional)").strip()

        with st.expander("Advanced"):
            c1, c2, c3 = st.columns(3)
            amp_min = c1.number_input("Amplicon min", 100, 2000, 225, 5)
            amp_opt = c2.number_input("optimum", 100, 2000, 300, 5)
            amp_max = c3.number_input("max", 100, 2000, 375, 5)
            tm = st.number_input("Target primer Tm (0 = from genome GC)", 0.0, 80.0, 0.0, 0.5)
            tm_window = st.number_input("Tm window (+/- C)", 0.5, 5.0, 1.5, 0.5)
            mono = st.number_input("Monovalent cations (mM: K+ + Na+ + half the Tris)", 0.0, 500.0, 50.0, 5.0)
            mg = st.number_input("Mg2+ (mM, 0 = polymerase preset)", 0.0, 20.0, 0.0, 0.25)
            dntp = st.number_input("Total dNTP (mM)", 0.0, 5.0, 0.8, 0.1)
            dmso = st.number_input("DMSO (%)", 0.0, 10.0, 0.0, 0.5)
            ta = st.number_input("Force annealing temperature (0 = derive)", 0.0, 80.0, 0.0, 0.5)
            outer_min = st.number_input("Separate pools for amplicons closer than (bp)", 0, 5000, 1000, 100)
            one_pool = st.checkbox("Force one pool (close targets merge into larger amplicons; the rest stay short)",
                                   value=False)
            one_pool_max = st.number_input("Largest amplicon with one pool (bp)", 200, 5000, 1000, 25)
            one_pool_drop = st.checkbox("One pool: drop targets that would need a second pool (instead of adding it)",
                                        value=False)
            vcf_max_indel = st.number_input("VCF: leave out indels and structural variants longer than (bp, 0 = keep all)",
                                            0, 100000, 0, 10)
            balance = st.checkbox("Binding balance (5' extensions, lowercase)", value=True)
            max_ext = st.number_input("Most bases added per primer", 1, 15, 8, 1)
            specificity = st.checkbox("Genome-wide specificity check", value=True)
            threads = st.number_input("Threads (0 = all but one)", 0, 256, 0, 1)
            st.markdown('<div class="sq-group-title">From sequencing</div>', unsafe_allow_html=True)
            redesign_below = st.number_input("Weak amplicon: below this fold of the median", 0.0, 1.0, 0.2, 0.05)
            conc_max = st.number_input("Highest primer concentration to suggest (nM); a weak amplicon needing more "
                                       "is redesigned", 50.0, 5000.0, 500.0, 50.0)
            reaction_ul = st.number_input("Reaction volume (uL)", 1.0, 200.0, 25.0, 1.0)
            stock_uM = st.number_input("Primer stock (uM)", 0.5, 1000.0, 10.0, 0.5)

    # positions source
    pos_arg = None
    if pfile is not None:
        p = work / pfile.name
        p.write_bytes(pfile.getvalue())
        pos_arg = str(p)
    elif ptext.strip():
        p = work / "positions.txt"
        p.write_text(ptext.strip() + "\n")
        pos_arg = str(p)
    var_arg = None
    if vfile is not None:
        v = work / vfile.name
        v.write_bytes(vfile.getvalue())
        var_arg = str(v)
    seq_arg = None
    if sfile is not None:
        p = work / sfile.name
        p.write_bytes(sfile.getvalue())
        seq_arg = str(p)
    prev_arg = last_primers
    if pprev is not None:
        p = work / "previous_primers.tsv"
        p.write_bytes(pprev.getvalue())
        prev_arg = str(p)
    conc_arg = None
    if cfile is not None:
        p = work / "conc.tsv"
        p.write_bytes(cfile.getvalue())
        conc_arg = str(p)

    opts = dict(polymerase=polymerase, primer_nM=primer_nM, tails=tails, tail_fwd=tail_fwd, tail_rev=tail_rev,
                amp_min=amp_min, amp_opt=amp_opt, amp_max=amp_max, tm=tm or None, tm_window=tm_window,
                mono_mM=mono, mg_mM=mg or None, dntp_mM=dntp, dmso=dmso, ta=ta or None, outer_min=outer_min,
                balance=balance, max_ext=max_ext, specificity=specificity, threads=threads or None,
                one_pool=one_pool, one_pool_max=one_pool_max, one_pool_drop=one_pool_drop,
                vcf_snps_only=vcf_snps_only, long_variant=long_variant, vcf_max_indel=vcf_max_indel,
                exclude=exclude or None,
                variants=var_arg)
    seq_opts = dict(sequencing=seq_arg, previous=prev_arg, conc=conc_arg, replace=replace or None, keep=keep or None,
                    redesign_below=redesign_below, conc_max=conc_max, reaction_ul=reaction_ul, stock_uM=stock_uM)

    # ------------------------------------------------------------------ actions
    can_design = bool(genome and pos_arg)
    can_rebalance = bool(seq_arg and prev_arg)
    with st.container(horizontal=True, gap="small"):
        check = st.button("Check positions", disabled=not can_design)
        run = st.button("Design pool", type="primary", disabled=not can_design)
        rebalance = st.button("Rebalance", disabled=not can_rebalance,
                              help="Suggest new concentrations from the sequencing run; no new primers.")
        redesign = st.button("Rebalance and redesign", disabled=not (can_rebalance and can_design),
                             help="Rebalance, and design new primers for the targets that failed.")
    if not can_design:
        st.info("Give a genome and at least one position in the sidebar. To rebalance a pool after sequencing, "
                "add the read counts per amplicon and the previous primer table under From sequencing.")
    elif not can_rebalance and seq_arg:
        st.info("Add the previous design's _primers.tsv to rebalance.")

    if check:
        with st.spinner("Checking the genome, positions, variants and settings"):
            try:
                from primerJinn.multiplex.cli import _id_list, build_parser, config_from_args
                from primerJinn.multiplex.precheck import run as precheck

                a = build_parser().parse_args(build_args(genome, pos_arg, str(work / "check" / "design"), opts))
                pc = precheck(a.genome, a.positions, config_from_args(a), variants=a.variants,
                              variant_min_af=a.variant_min_af, out=a.out, exclude=_id_list(a.exclude))
                ss.pc = pc
                ss.pc_args = build_args(genome, pos_arg, "x", opts)   # to notice settings changed since
                # the sidebar's "leave out" list offers the targets of this check; excluded ones stay listed
                known = dict(ss.get("check_targets", []))
                for t in pc.targets:
                    key = "|".join(f"{s.seqid}:{s.start}-{s.end}" for s in t.sites)
                    known.setdefault(key, f"{t.id}: " + ", ".join(s.label() for s in t.sites[:3])
                                     + (f" +{len(t.sites) - 3}" if len(t.sites) > 3 else ""))
                new_targets = sorted(known.items(), key=lambda kv: kv[1])
                if new_targets != ss.get("check_targets"):
                    ss.check_targets = new_targets
                    st.rerun()          # so the sidebar's "leave out" list offers these targets
            except Exception as exc:  # show the message, keep the app alive
                ss.pc = None
                st.error(str(exc))
    pc = ss.get("pc")
    if pc is not None and not (run or rebalance or redesign):
        with st.container():
            try:
                show = {"OK": st.success, "INFO": st.info, "WARN": st.warning, "ERROR": st.error}
                ne, nw = len(pc.errors), len(pc.warnings)
                if ss.get("pc_args") != build_args(genome, pos_arg, "x", opts):
                    st.warning("The inputs or settings have changed since this check: run Check positions again "
                               "before designing.")
                (st.error if ne else st.warning if nw else st.success)(
                    f"{ne} error(s), {nw} warning(s)" + (": fix the errors before designing." if ne else "."))
                for it in pc.items:
                    if it.level in ("WARN", "ERROR"):
                        show[it.level](f"**{it.topic}**: {it.msg}")
                with st.expander("All checks", expanded=not (ne or nw)):
                    st.markdown(table_html(pd.DataFrame([{"level": i.level, "topic": i.topic, "finding": i.msg}
                                                         for i in pc.items])), unsafe_allow_html=True)
                if pc.targets:
                    st.markdown("**Amplicon targets**")
                    st.markdown(table_html(pd.DataFrame([{"target": t.id, "sequence": t.seqid, "start": t.start,
                                                          "end": t.end, "span": t.span,
                                                          "positions": ", ".join(s.label() for s in t.sites)}
                                                         for t in pc.targets])), unsafe_allow_html=True)
                    st.caption("To leave a target out, pick it under Positions in the sidebar (Leave out targets) "
                               "and check again; its positions are then not designed.")
            except Exception as exc:  # show the message, keep the app alive
                st.error(str(exc))

    if run or rebalance or redesign:
        ss.run += 1
        prefix = str(work / f"run{ss.run}" / "design")
        Path(prefix).parent.mkdir(parents=True, exist_ok=True)
        all_opts = dict(opts)
        if rebalance or redesign:
            all_opts.update(seq_opts, no_redesign=bool(rebalance))
        cmd = [sys.executable, "-m", "primerJinn.multiplex.cli"] + build_args(genome, pos_arg, prefix, all_opts)
        label = ("Rebalancing the pool from the sequencing run" if rebalance else
                 "Rebalancing and redesigning (1 to 4 minutes for a bacterial genome)" if redesign else
                 "Designing the pool (1 to 4 minutes for a bacterial genome)")
        log_lines = []
        with st.status(label, expanded=True) as status:
            box = st.empty()
            env = dict(os.environ, PYTHONUNBUFFERED="1")
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env,
                                    cwd=str(_ROOT))
            for line in proc.stdout:
                log_lines.append(line.rstrip())
                box.code("\n".join(log_lines[-25:]), language=None)
            code = proc.wait()
            made = Path(prefix + "_report.md").exists() or Path(prefix + "_rebalance.tsv").exists()
            ok = code in (0, 1) and made
            status.update(label=("Done" if code == 0 else "Done, with uncovered positions" if ok else
                                 "The run failed"), state="complete" if ok else "error", expanded=not ok)
        ss.result = prefix if ok else None
        if ok and Path(prefix + "_primers.tsv").exists() and Path(prefix + "_report.md").exists():
            ss.last_primers = prefix + "_primers.tsv"
        ss.log = "\n".join(log_lines)

    # ------------------------------------------------------------------ results
    prefix = ss.get("result")
    if not prefix:
        if ss.get("log"):
            st.error("The run did not finish. Log:")
            st.code(ss.log[-4000:], language=None)
    else:
        show_results(prefix, st, pd, alt)
    st.markdown(footer_html(), unsafe_allow_html=True)


DOWNLOADS = (("_primers.tsv", "Primer table (TSV)", "text/tab-separated-values"),
             ("_order_sheet.tsv", "Order sheet (TSV)", "text/tab-separated-values"),
             ("_order_sheet_new.tsv", "Order sheet, new primers only (TSV)", "text/tab-separated-values"),
             ("_conc.tsv", "Concentration table (TSV)", "text/tab-separated-values"),
             ("_rebalance.tsv", "Rebalance table (TSV)", "text/tab-separated-values"),
             (".xlsx", "All tables (XLSX)", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"),
             ("_amplicons.fasta", "Amplicons (FASTA)", "text/plain"),
             (".bed", "Amplicons (BED)", "text/plain"),
             ("_report.md", "Report (Markdown)", "text/markdown"))


def download_buttons(prefix, st, suffixes=None, zip_all=True):
    for suffix, label, mime in DOWNLOADS:
        if suffixes is not None and suffix not in suffixes:
            continue
        f = Path(prefix + suffix)
        if f.exists():
            st.download_button(label, f.read_bytes(), file_name="primerjinn" + suffix, mime=mime,
                               key=f"dl{suffix}{'' if suffixes is None else '_fb'}")
    if zip_all:
        st.download_button("Everything (ZIP)", zip_outputs(prefix), file_name="primerjinn_design.zip",
                           mime="application/zip", key="dl_zip" + ("" if suffixes is None else "_fb"))


def show_feedback(prefix, st, pd):
    """The Sequencing feedback tab: stat tiles, warnings, the rebalance table with action pills, downloads."""
    fb = Path(prefix + "_feedback.tsv")
    items = []
    if fb.exists():
        t = pd.read_csv(fb, sep="\t", dtype=str, keep_default_na=False)
        items = list(zip(t.iloc[:, 0], t.iloc[:, 1])) if len(t.columns) >= 2 else []
    tiles = []
    for keys, label, unit in FEEDBACK_TILES:
        keys = (keys,) if isinstance(keys, str) else keys
        v = next((v for k, v in items if k.strip().lower().startswith(keys)), None)
        if v is None:
            continue
        v, status = str(v).strip(), ""
        try:
            x = float(v.replace("%", "").replace(",", "").split("/")[0].split()[0])
            if label == "Unaligned or off-target":
                status = "good" if x < 10 else "warn" if x < 30 else "bad"
            elif label == "Designed reads":
                status = "good" if x >= 80 else "warn" if x >= 50 else "bad"
            elif label == "Targets replaced":
                status = "muted" if x == 0 else "warn"
        except (ValueError, IndexError):
            pass
        tiles.append((label, v if not unit or v.endswith(unit.strip()) else v + unit, status))
    if tiles:
        st.markdown(stat_tiles_html(tiles), unsafe_allow_html=True)
    for k, v in items:
        if k.strip().lower() == "warning":
            st.warning(v)
    notes = [(k, v) for k, v in items if k.strip().lower() in ("rule", "concentration source", "source",
                                                                 "primers raised / lowered / kept")]
    if notes:
        st.caption("  \n".join(f"{k}: {v}" for k, v in notes))
    fmd = Path(prefix + "_feedback.md")
    if fmd.exists():
        with st.expander("Feedback report"):
            st.markdown(fmd.read_text())
    rb = pd.read_csv(prefix + "_rebalance.tsv", sep="\t").copy()
    if "action" in rb:
        rb["action"] = rb["action"].map(pill_html)
    if "unit" in rb and "target" in rb and (rb["unit"].astype(str) == rb["target"].astype(str)).all():
        rb = rb.drop(columns="unit")
    rb = rb.drop(columns=[c for c in ("conc_source",) if c in rb])
    why = None
    if "reason" in rb:  # the long text goes in its own table, one row per target
        why = rb[["target", "action", "reason"]].drop_duplicates("target") if "target" in rb else None
        rb = rb.drop(columns="reason")
    st.markdown("**Rebalance table**: suggested concentration per primer and what to do with each target.")
    st.markdown(table_html(rb, raw_cols=("action",)), unsafe_allow_html=True)
    if why is not None and why["reason"].notna().any():
        st.markdown("**Why**")
        st.markdown(table_html(why, raw_cols=("action",)), unsafe_allow_html=True)
    st.markdown("**Downloads**")
    download_buttons(prefix, st, suffixes=("_rebalance.tsv", "_conc.tsv", "_order_sheet.tsv", "_order_sheet_new.tsv"),
                     zip_all=False)


def show_results(prefix, st, pd, alt):
    has_design = Path(prefix + "_report.md").exists()
    has_feedback = Path(prefix + "_rebalance.tsv").exists()
    if not has_design:  # rebalance only: no report, only the feedback
        st.subheader("Sequencing feedback")
        show_feedback(prefix, st, pd)
        return

    rd = lambda name: pd.read_csv(f"{prefix}_{name}.tsv", sep="\t")
    primers, amps = rd("primers"), rd("amplicons")
    conds = rd("conditions")
    cmap = dict(zip(conds["item"], conds["value"]))
    md = Path(prefix + "_report.md").read_text()
    sec = report_sections(md)

    c = st.columns(4)
    c[0].metric("Amplicons", len(amps))
    c[1].metric("Pools", int(amps["pool"].max()) if len(amps) else 0)
    c[2].metric("Annealing temperature", cmap.get("RECOMMENDED annealing temperature", ""))
    c[3].metric("Amplicon length", f"{amps['length'].min()}-{amps['length'].max()} bp" if len(amps) else "")
    cov = sec.get("Coverage", "").splitlines()[0] if sec.get("Coverage") else ""
    (st.success if cov.startswith("All") else st.warning)(cov or "Coverage not reported")
    adv = [v for k, v in zip(conds["item"], conds["value"]) if k == "advice"]
    if adv:
        with st.expander("PCR conditions and advice", expanded=True):
            for a in adv:
                st.markdown(f"- {a}")
    if sec.get("Warnings"):
        with st.expander("Warnings", expanded=False):
            st.markdown(sec["Warnings"])

    names = ["Primers", "Amplicon map", "Binding vs temperature", "Dimers and products", "Conditions", "Report"]
    if has_feedback:
        names.append("Sequencing feedback")
    names.append("Downloads")
    tabs = dict(zip(names, st.tabs(names)))
    with tabs["Primers"]:
        st.caption("Uppercase: Tm-matched core. Highlighted lowercase: 5' bases added to balance binding "
                   "(they match the template; order the sequence as written).")
        show = primers[["pool", "primer", "order_sequence", "length", "gc_percent", "tm_core_C", "tm_C", "dg_eff",
                        "fraction_bound_model", "site_open_cost", "offtarget_primable_sites"]
                       + (["variants_under_primer"] if "variants_under_primer" in primers
                          and primers["variants_under_primer"].notna().any() else [])].copy()
        show["order_sequence"] = show["order_sequence"].map(sequence_html)
        st.markdown(table_html(show, raw_cols=("order_sequence",)), unsafe_allow_html=True)
    with tabs["Amplicon map"]:
        if len(amps):
            amps = amps.copy()
            amps["mid"] = (amps["amplicon_start"] + amps["amplicon_end"]) / 2
            amps["pool"] = amps["pool"].astype(str)
            st.markdown("**Genome overview** (each mark is an amplicon; hover for details)")
            overview = alt.Chart(amps).mark_circle(size=120).encode(
                x=alt.X("mid:Q", title="genome position"), y=alt.Y("pool:N", title="pool"),
                color="pool:N", tooltip=["target", "covers", "length", "gc_percent", "seqid"])
            if amps["seqid"].nunique() > 1:
                overview = overview.facet(row="seqid:N")
            st.altair_chart(overview, width="stretch", theme=None)
            tgt = st.selectbox("Zoom on an amplicon", amps["target"] + "  " + amps["covers"].str.slice(0, 80))
            row = amps.iloc[[i for i, v in enumerate(amps["target"] + "  " + amps["covers"].str.slice(0, 80))
                             if v == tgt][0]]
            pr = primers[primers["target"] == row["target"]]
            parts = [dict(track="amplicon", start=row["amplicon_start"], end=row["amplicon_end"], label="amplicon"),
                     dict(track="reads cover", start=row["insert_start"], end=row["insert_end"], label="insert")]
            for _, p in pr.iterrows():
                parts.append(dict(track="primers", start=p["start"], end=p["end"],
                                  label=f"{p['primer']} {p['direction']}"))
            for a, b in requested_positions(row["covers"]):
                a, b = max(a, row["amplicon_start"]), min(b, row["amplicon_end"])
                parts.append(dict(track="requested", start=a, end=b + 1, label=f"{a}-{b}" if b > a else str(a)))
            dd = pd.DataFrame(parts)
            zoom = alt.Chart(dd).mark_bar(height=14).encode(
                x=alt.X("start:Q", scale=alt.Scale(zero=False), title=f"{row['seqid']} position"), x2="end:Q",
                y=alt.Y("track:N", sort=["amplicon", "primers", "reads cover", "requested"], title=None),
                color=alt.Color("track:N", legend=None), tooltip=["label", "start", "end"])
            st.altair_chart(zoom, width="stretch", theme=None)
            st.caption(f"{row['length']} bp, GC {row['gc_percent']} %, product Tm {row['product_tm_C']} C. "
                       f"{row.get('flags') if isinstance(row.get('flags'), str) else ''}")
    with tabs["Binding vs temperature"]:
        bt = rd("binding_vs_temperature")
        rec = re.findall(r"[\d.]+", str(cmap.get("RECOMMENDED annealing temperature", "")))
        weakest = primers.sort_values("dg_eff", ascending=False)["primer"].head(4).tolist()
        sel = st.multiselect("Primers", sorted(bt["primer"].unique()), default=weakest)
        sub = bt[bt["primer"].isin(sel)]
        if len(sub):
            base = alt.Chart(sub).encode(x=alt.X("temp_C:Q", title="temperature (C)"), color="primer:N")
            lines = base.mark_line().encode(y=alt.Y("theta:Q", title="fraction bound (with structure)"),
                                            tooltip=["primer", "temp_C", "theta", "theta_nostruct", "site_open"])
            dashed = base.mark_line(strokeDash=[4, 3], opacity=0.5).encode(y="theta_nostruct:Q")
            chart = lines + dashed
            if rec:
                chart = chart + alt.Chart(pd.DataFrame({"Ta": [float(rec[0])]})).mark_rule(
                    color="#e66767").encode(x="Ta:Q")
            st.altair_chart(chart, width="stretch", theme=None)
            st.caption("Solid: with template and primer structure. Dashed: unfolded template. Red: recommended "
                       "annealing temperature. Model of the assumed buffer; for Q5 use it to compare primers.")
    with tabs["Dimers and products"]:
        st.markdown("**Strongest primer dimers (37 C)**")
        st.markdown(table_html(rd("dimers").head(30)), unsafe_allow_html=True)
        st.markdown("**Other products** (outer products, repeated loci, unintended)")
        st.markdown(sec.get("Other products (genome-wide, within pools)", "")[-3000:])
    with tabs["Conditions"]:
        st.markdown(table_html(conds[conds["item"] != "advice"]), unsafe_allow_html=True)
    with tabs["Report"]:
        st.markdown(md)
    if has_feedback:
        with tabs["Sequencing feedback"]:
            show_feedback(prefix, st, pd)
    with tabs["Downloads"]:
        download_buttons(prefix, st)


def main():
    """Entry point: start the Streamlit server on this file."""
    from streamlit.web import cli as stcli

    theme = [x for k, v in THEME.items() for x in (f"--theme.{k}", v)]
    sys.argv = ["streamlit", "run", str(Path(__file__).resolve()), "--browser.gatherUsageStats", "false"] \
        + theme + sys.argv[1:]
    sys.exit(stcli.main())


if __name__ == "__main__":
    app()
