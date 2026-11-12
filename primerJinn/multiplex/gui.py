"""primerJinn multiplex designer: a small local web GUI.

Start with `primerjinn-gui` (after `pip install .[gui]`) or
`streamlit run primerJinn/multiplex/gui.py`. The design itself runs the
command-line tool in a separate process, so results are identical to the CLI.
"""
from __future__ import annotations

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

POLY = {"kapa": "KAPA2G Robust HotStart ReadyMix (default)", "taq": "Standard Taq (1.5 mM MgCl2)",
        "q5": "NEB Q5", "hifi": "Proofreading high-fidelity (model Ta rule)"}
TAILS = {"none": "None", "illumina": "Illumina TruSeq partial adapters", "nextera": "Illumina Nextera overhangs",
         "custom": "Custom"}


# ----------------------------------------------------------------------------
# helpers (no Streamlit needed; tested in tests/test_gui.py)
# ----------------------------------------------------------------------------
def build_args(genome: str, positions: str, out_prefix: str, opts: dict) -> list:
    """Command-line arguments for primerJinn.multiplex.cli from GUI settings."""
    a = ["--genome", genome, "--positions", positions, "--out", out_prefix,
         "--polymerase", opts.get("polymerase", "kapa"), "--primer_nM", str(opts.get("primer_nM", 100.0))]
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
    return a


def sequence_html(seq: str) -> str:
    """Primer sequence with lowercase (5' extension) bases highlighted."""
    out = []
    for run in re.findall(r"[a-z]+|[^a-z]+", seq):
        if run.islower():
            out.append(f'<span style="color:#b45309;background:#fef3c7;font-weight:600">{html.escape(run)}</span>')
        else:
            out.append(html.escape(run))
    return f'<span style="font-family:monospace">{"".join(out)}</span>'


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


# ----------------------------------------------------------------------------
# the app
# ----------------------------------------------------------------------------
def app():
    import altair as alt
    import pandas as pd
    import streamlit as st

    st.set_page_config(page_title="primerJinn multiplex", page_icon="🧬", layout="wide")
    ss = st.session_state
    if "workdir" not in ss:
        ss.workdir = tempfile.mkdtemp(prefix="primerjinn_gui_")
        ss.run = 0
    work = Path(ss.workdir)

    st.title("primerJinn multiplex designer")
    st.caption("Design a multiplex amplicon-sequencing primer pool from a genome and a list of positions. "
               "Every position is covered; primers are matched on Tm and on binding at the annealing "
               "temperature, including template structure.")

    # ------------------------------------------------------------------ inputs
    with st.sidebar:
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
                   "or a table / BED file.")
        ptext = st.text_area("Positions", height=150, label_visibility="collapsed",
                             placeholder="gyrA_A90\t7570\nrpoB_S450\t761155\npncA\t2288681-2289241")
        pfile = st.file_uploader("or upload a positions file", type=["tsv", "csv", "txt", "xlsx", "bed"])

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
            balance = st.checkbox("Binding balance (5' extensions, lowercase)", value=True)
            max_ext = st.number_input("Most bases added per primer", 1, 15, 8, 1)
            specificity = st.checkbox("Genome-wide specificity check", value=True)
            threads = st.number_input("Threads (0 = all but one)", 0, 256, 0, 1)

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

    opts = dict(polymerase=polymerase, primer_nM=primer_nM, tails=tails, tail_fwd=tail_fwd, tail_rev=tail_rev,
                amp_min=amp_min, amp_opt=amp_opt, amp_max=amp_max, tm=tm or None, tm_window=tm_window,
                mono_mM=mono, mg_mM=mg or None, dntp_mM=dntp, dmso=dmso, ta=ta or None, outer_min=outer_min,
                balance=balance, max_ext=max_ext, specificity=specificity, threads=threads or None,
                variants=var_arg)

    # ------------------------------------------------------------------ actions
    b1, b2, _ = st.columns([1, 1, 4])
    check = b1.button("Check positions", disabled=not (genome and pos_arg))
    run = b2.button("Design pool", type="primary", disabled=not (genome and pos_arg))
    if not (genome and pos_arg):
        st.info("Give a genome and at least one position in the sidebar.")

    if check:
        with st.spinner("Resolving the genome and grouping positions"):
            try:
                from primerJinn.multiplex import genome as gmod
                from primerJinn.multiplex.cli import sites_from_positions
                from primerJinn.multiplex.config import Config, primer_length_range
                from primerJinn.multiplex.targets import group_sites, span_limits

                meta = gmod.resolve_genome(genome)
                sites = sites_from_positions(pos_arg, meta)
                gc = (meta.get("gc_percent") or 50) / 100
                lo, _hi = primer_length_range(gc)
                comfort = span_limits(amp_opt, amp_max, primer_len=lo + 4, margin=3, flex=20)[0]
                _, hard, gap = span_limits(amp_opt, amp_max, primer_len=lo, margin=3, flex=20)
                targets = group_sites(sites, comfort, max(hard, comfort), gap)
                cfg = Config(**({"tm": tm} if tm else {})).resolve(gc)
                st.success(f"{meta.get('organism') or genome}: {meta['total_length']:,} bp, "
                           f"{meta['n_sequences']} sequence(s), GC {meta['gc_percent']}%. "
                           f"{len(sites)} position(s) -> {len(targets)} amplicon target(s); target primer Tm "
                           f"{cfg.tm:.1f} C.")
                st.dataframe(pd.DataFrame([{"target": t.id, "sequence": t.seqid, "start": t.start, "end": t.end,
                                            "span": t.span, "positions": ", ".join(s.label() for s in t.sites)}
                                           for t in targets]), hide_index=True, use_container_width=True)
                for n in meta.get("notes", []) + meta.get("warnings", []):
                    st.caption(n)
            except Exception as exc:  # show the message, keep the app alive
                st.error(str(exc))

    if run:
        ss.run += 1
        prefix = str(work / f"run{ss.run}" / "design")
        Path(prefix).parent.mkdir(parents=True, exist_ok=True)
        cmd = [sys.executable, "-m", "primerJinn.multiplex.cli"] + build_args(genome, pos_arg, prefix, opts)
        log_lines = []
        with st.status("Designing the pool (1 to 4 minutes for a bacterial genome)", expanded=True) as status:
            box = st.empty()
            env = dict(os.environ, PYTHONUNBUFFERED="1")
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env,
                                    cwd=str(_ROOT))
            for line in proc.stdout:
                log_lines.append(line.rstrip())
                box.code("\n".join(log_lines[-25:]), language=None)
            code = proc.wait()
            ok = code in (0, 1) and Path(prefix + "_report.md").exists()
            status.update(label=("Done" if code == 0 else "Done, with uncovered positions" if ok else
                                 "Design failed"), state="complete" if ok else "error", expanded=not ok)
        ss.result = prefix if ok else None
        ss.log = "\n".join(log_lines)

    # ------------------------------------------------------------------ results
    prefix = ss.get("result")
    if not prefix:
        if ss.get("log"):
            st.error("The design did not finish. Log:")
            st.code(ss.log[-4000:], language=None)
        return
    show_results(prefix, st, pd, alt)


def show_results(prefix, st, pd, alt):
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

    tabs = st.tabs(["Primers", "Amplicon map", "Binding vs temperature", "Dimers and products", "Conditions",
                    "Report", "Downloads"])
    with tabs[0]:
        st.caption("Uppercase: Tm-matched core. Highlighted lowercase: 5' bases added to balance binding "
                   "(they match the template; order the sequence as written).")
        show = primers[["pool", "primer", "order_sequence", "length", "gc_percent", "tm_core_C", "tm_C", "dg_eff",
                        "fraction_bound_model", "site_open_cost", "offtarget_primable_sites"]
                       + (["variants_under_primer"] if "variants_under_primer" in primers
                          and primers["variants_under_primer"].notna().any() else [])].copy()
        show["order_sequence"] = show["order_sequence"].map(sequence_html)
        st.markdown(show.to_html(escape=False, index=False, float_format=lambda x: f"{x:.2f}"),
                    unsafe_allow_html=True)
    with tabs[1]:
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
            st.altair_chart(overview, use_container_width=True)
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
            st.altair_chart(zoom, use_container_width=True)
            st.caption(f"{row['length']} bp, GC {row['gc_percent']}%, product Tm {row['product_tm_C']} C. "
                       f"{row.get('flags') if isinstance(row.get('flags'), str) else ''}")
    with tabs[2]:
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
                chart = chart + alt.Chart(pd.DataFrame({"Ta": [float(rec[0])]})).mark_rule(color="red").encode(x="Ta:Q")
            st.altair_chart(chart, use_container_width=True)
            st.caption("Solid: with template and primer structure. Dashed: unfolded template. Red: recommended "
                       "annealing temperature. Model of the assumed buffer; for Q5 use it to compare primers.")
    with tabs[3]:
        st.markdown("**Strongest primer dimers (37 C)**")
        st.dataframe(rd("dimers").head(30), hide_index=True, use_container_width=True)
        st.markdown("**Other products** (outer products, repeated loci, unintended)")
        st.markdown(sec.get("Other products (genome-wide, within pools)", "")[-3000:])
    with tabs[4]:
        st.dataframe(conds[conds["item"] != "advice"], hide_index=True, use_container_width=True)
    with tabs[5]:
        st.markdown(md)
    with tabs[6]:
        for suffix, label, mime in (("_primers.tsv", "Primer order sheet (TSV)", "text/tab-separated-values"),
                                    (".xlsx", "All tables (XLSX)",
                                     "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"),
                                    ("_amplicons.fasta", "Amplicons (FASTA)", "text/plain"),
                                    (".bed", "Amplicons (BED)", "text/plain"),
                                    ("_report.md", "Report (Markdown)", "text/markdown")):
            f = Path(prefix + suffix)
            if f.exists():
                st.download_button(label, f.read_bytes(), file_name="primerjinn" + suffix, mime=mime)
        st.download_button("Everything (ZIP)", zip_outputs(prefix), file_name="primerjinn_design.zip",
                           mime="application/zip")


def main():
    """Entry point: start the Streamlit server on this file."""
    from streamlit.web import cli as stcli

    sys.argv = ["streamlit", "run", str(Path(__file__).resolve()), "--browser.gatherUsageStats", "false"] \
        + sys.argv[1:]
    sys.exit(stcli.main())


if __name__ == "__main__":
    app()
