"""Write a Design to TSV, XLSX, FASTA, BED and a markdown report."""
from __future__ import annotations

import os
from typing import Dict, List

import pandas as pd

from .config import POLYMERASES
from .balance import display
from .design import Design


def _r(x, n=2):
    return None if x is None else round(float(x), n)


def tables(d: Design) -> Dict[str, pd.DataFrame]:
    cfg, ta = d.cfg, d.ta
    tmap = {t.id: t for t in d.targets}
    at = ta["at_recommended"]
    prim_rows: List[dict] = []
    amp_rows: List[dict] = []
    order = sorted(d.chosen, key=lambda tid: (d.pools[tid], tmap[tid].seqid, tmap[tid].start))
    for tid in order:
        pr = d.chosen[tid]
        t = tmap[tid]
        ins = pr.insert
        parts = []
        for s in t.sites:
            a, b = max(s.start, ins[0]), min(s.end, ins[1])
            if a > b:
                continue
            loc = str(a) if a == b else f"{a}-{b}"
            tag = "" if (a, b) == (s.start, s.end) else f", part of {s.start}-{s.end}"
            parts.append(f"{s.name} ({s.seqid}:{loc}{tag})" if s.name else f"{s.seqid}:{loc}{tag}")
        covered = "; ".join(parts)
        for side, p in (("F", pr.fwd), ("R", pr.rev)):
            nm = d.names[tid][side]
            a = at.get(nm, {})
            prim_rows.append({
                "pool": d.pools[tid], "target": tid, "target_name": t.name, "covers": covered,
                "primer": nm, "direction": "forward" if side == "F" else "reverse",
                "sequence": display(p), "order_sequence": (cfg.tail_fwd if side == "F" else cfg.tail_rev) + display(p),
                "added_5prime_nt": p.ext5,
                "seqid": p.seqid, "start": p.start, "end": p.end, "three_prime_end": p.three_prime,
                "length": p.length, "gc_percent": _r(100 * p.gc, 1),
                "tm_C": _r(p.tm, 1),
                "tm_core_C": _r(p.tm_core if p.tm_core is not None else p.tm, 1),
                "tm_neb_q5_scale_C": _r(p.tm_neb, 1),
                "tm_physical_C": _r(p.tm_phys, 1),
                f"dg_duplex_at_{ta['rec_grid']:g}C": _r(a.get("dg_dup")),
                "hairpin_cost": _r(a.get("hairpin")),
                "site_open_cost": _r(a.get("site_open")),
                "dg_eff": _r(a.get("dg_eff")),
                "fraction_bound_model": _r(a.get("theta"), 3),
                "self_dimer_3prime_dg37": _r(p.self_d3),
                "offtarget_primable_sites": len(d.offtarget_sites.get(nm, [])) if cfg.specificity else None,
                "variants_under_primer": p.n_variants if cfg.variants is not None else None,
                "kept_from_previous": bool(p.fixed),
                "deletion": f"{pr.gap[0]}-{pr.gap[1]}" if pr.gap else "",
                "notes": "; ".join(p.notes),
            })
        rep = d.amplicon_reports.get(tid)
        strands = {(s.strand, s.celsius): s for s in rep.strands} if rep else {}
        top = strands.get(("top", float(ta["recommended"])))
        bot = strands.get(("bottom", float(ta["recommended"])))
        ins = pr.insert
        amp_rows.append({
            "pool": d.pools[tid], "target": tid, "target_name": t.name, "covers": covered,
            "forward": d.names[tid]["F"], "reverse": d.names[tid]["R"], "seqid": pr.fwd.seqid,
            "amplicon_start": pr.start, "amplicon_end": pr.end, "length": pr.length,
            "kind": "junction" if pr.gap else "standard",
            "deletion": f"{pr.fwd.seqid}:{pr.gap[0]}-{pr.gap[1]} ({pr.gap_length} bp)" if pr.gap else "",
            "wildtype_product_bp": pr.genomic_length if pr.gap else "",
            "size_class_bp": "{}-{}".format(cfg.amp(tid)[0], cfg.amp(tid)[2]),
            "length_with_tails": pr.length + len(cfg.tail_fwd) + len(cfg.tail_rev),
            "insert_start": ins[0], "insert_end": ins[1],
            "gc_percent": _r(100 * pr.gc, 1),
            "max_gc_50bp_percent": _r(100 * rep.max_gc_50, 0) if rep else None,
            "product_tm_C": _r(rep.duplex_tm, 1) if rep else None,
            "ensemble_dg_top": _r(top.g_ensemble, 1) if top else None,
            "ensemble_dg_bottom": _r(bot.g_ensemble, 1) if bot else None,
            "fraction_paired_top": _r(top.frac_paired, 2) if top else None,
            "self_priming_top": _r(pr.self_prime_top, 3), "self_priming_bottom": _r(pr.self_prime_bottom, 3),
            "fwd_site_open_amplicon": _r(pr.site_open_amp_f), "rev_site_open_amplicon": _r(pr.site_open_amp_r),
            "strongest_hairpin_dg_at_ta": _r(min([h.dg for h in (rep.hairpin_top, rep.hairpin_bottom) if h] or [0])) if rep else None,
            "strongest_hairpin_dg_at_72C": _r(min([h.dg_72 for h in (rep.hairpin_top, rep.hairpin_bottom) if h] or [0])) if rep else None,
            "strongest_hairpin_tm_C": (lambda tm: ">95" if tm > 95 else round(tm))(
                max([h.tm for h in (rep.hairpin_top, rep.hairpin_bottom) if h] or [0])) if rep else None,
            "g_quadruplex_motifs": len(rep.g4) if rep else None,
            "pair_cost": _r(pr.cost),
            "relaxed": d.relaxed.get(tid, ""),
            "flags": "; ".join(pr.flags + (rep.flags if rep else [])),
        })
    out = {
        "primers": pd.DataFrame(prim_rows),
        "amplicons": pd.DataFrame(amp_rows),
        "dimers": pd.DataFrame(d.dimers),
        "products": pd.DataFrame(d.products, columns=["pool", "kind", "fwd", "rev", "contig", "start", "end", "size",
                                                      "fwd_mismatches", "rev_mismatches", "score", "call"]),
        "offtarget_sites": pd.DataFrame([r for rows in d.offtarget_sites.values() for r in rows],
                                        columns=["primer", "contig", "start", "strand", "mismatches",
                                                 "mismatch_positions_from_3prime", "dg_kcal", "score"]),
    }
    # binding curves
    rows = []
    for nm, row in ta["table"].items():
        for T, v in row.items():
            rows.append({"primer": nm, "temp_C": T, **{k: _r(x, 4) for k, x in v.items()}})
    out["binding_vs_temperature"] = pd.DataFrame(rows)
    out["conditions"] = pd.DataFrame(conditions(d), columns=["item", "value"])
    return out


def conditions(d: Design) -> List[tuple]:
    cfg, ta = d.cfg, d.ta
    g = d.genome
    rows = [
        ("genome", f"{g.get('organism') or ''} {g.get('accession') or g.get('fasta')}".strip()),
        ("genome GC %", g.get("gc_percent")),
        ("polymerase preset", f"{cfg.polymerase} ({POLYMERASES[cfg.polymerase]['label']})"),
        ("Tm scale", "physical: SantaLucia 2004 + Owczarzy 2008 in the assumed buffer at the primer "
                     "concentration" if cfg.tm_scale == "physical" else "NEB Tm calculator (Q5), 500 nM"),
        ("target primer Tm", f"{cfg.tm:.1f} +/- {cfg.tm_window} C"),
        ("primer Tm range in pool", f"{ta['tm_min']:.1f} to {ta['tm_max']:.1f} C"),
        ("same primers on the NEB Q5 calculator scale", f"{ta['tm_neb_min']:.1f} to {ta['tm_neb_max']:.1f} C"),
        ("assumed buffer", f"{cfg.mono_mM:g} mM monovalent, {cfg.mg_mM:g} mM Mg2+, {cfg.dntp_mM:g} mM dNTP total, "
                           f"{cfg.dmso_pct:g}% DMSO"),
        ("primer concentration", f"{cfg.primer_nM:g} nM each"),
        ("RECOMMENDED annealing temperature", f"{ta['recommended']:g} C"),
        ("temperature primers, structure and off-target sites were scored at",
         f"{cfg.ta:.1f} C" + (f" (after {d.passes} passes)" if d.passes > 1 else "")),
        ("genome-wide specificity", "checked" if cfg.specificity else "NOT checked (--no_specificity)"),
        ("known strain variants kept out of primer 3' ends",
         f"{len(cfg.variants)} positions from {cfg.variants.source} (rejected within {cfg.variant_3prime_nt} nt of the "
         f"3' end, cost {cfg.variant_penalty:g} each elsewhere)" if cfg.variants is not None else "none given"),
        ("vendor rule Ta" + (f" (NEB Q5: lowest Tm at {cfg.primer_nM:g} nM + 1, max 72)" if cfg.polymerase == "q5"
                             else " (generic Taq master mix: 55/60/65 C for primer Tm <55/55-65/>65)" if cfg.polymerase == "kapa"
                             else ""), ta["vendor_ta"] if ta["vendor_ta"] is not None else "n/a (model rule)"),
    ]
    if cfg.polymerase == "q5":
        rows.append(("physical-model temperatures (no polymerase effect, for comparison only)",
                     f"90% bound without structure {ta['t90_without_structure']} C, 50% bound with structure "
                     f"{ta['t50']} C"))
    else:
        rows += [
            ("equilibrium model Ta: all primers >= 90% bound on unfolded template, capped where structure leaves "
             "the weakest primer at 50%", ta["model_ta"]),
            ("temperature where every primer is >= 90% bound including structure (model)", ta["t90"]),
            ("temperature where every primer is >= 50% bound (model)", ta["t50"]),
            ("90% binding temperature without template structure (model)", ta["t90_without_structure"]),
        ]
    rows += [
        ("primer that limits the model Ta (first below the binding criterion)", ta["limiting_primer"]),
        ("weakest primer at the recommended Ta, with template structure",
         f"{ta['weakest']} ({100 * ta['weakest_theta']:.0f}% bound)"),
        (f"dG_eff spread at {ta['rec_grid']:g} C (kcal/mol): full range / middle 50%",
         f"{ta['dg_eff_spread']:.2f} / {ta['dg_eff_iqr']:.2f}"),
        ("pools", (max(d.pools.values()) if d.pools else 0)),
        *[("pool layout", n) for n in d.layout_notes],
        *[("redesign", n) for n in d.redesign_notes],
        ("amplicons", len(d.chosen)),
        ("amplicon length range", f"{min(p.length for p in d.chosen.values())} to "
                                  f"{max(p.length for p in d.chosen.values())} bp"),
        ("adapter tails", cfg.tails),
        ("optimiser", f"{d.method}, energy {d.energy:.1f}"),
    ]
    for a in ta.get("advice", []):
        rows.append(("advice", a))
    return rows


def _reasons(d: Design) -> List[str]:
    r = list(getattr(d, "dropped_reasons", []) or [])
    return r + ["N (unknown) in the reference"] * (len(d.dropped_sites) - len(r))


def uncovered(d: Design) -> pd.DataFrame:
    """Requested positions that no amplicon covers, with the reason."""
    rows = [dict(name=s.name, seqid=s.seqid, start=s.start, end=s.end,
                 reason=f"not between the primers of any amplicon {gaps}") for s, gaps in d.coverage_problems]
    rows += [dict(name=s.name, seqid=s.seqid, start=s.start, end=s.end, reason=why)
             for s, why in zip(d.dropped_sites, _reasons(d))]
    return pd.DataFrame(rows, columns=["name", "seqid", "start", "end", "reason"])


def write(d: Design, prefix: str, fasta_with_primers: bool = True) -> List[str]:
    os.makedirs(os.path.dirname(os.path.abspath(prefix)), exist_ok=True)
    tabs = tables(d)
    written = []
    for name, df in tabs.items():
        path = f"{prefix}_{name}.tsv"
        df.to_csv(path, sep="\t", index=False)
        written.append(path)
    unc = uncovered(d)
    if len(unc):
        path = f"{prefix}_uncovered.tsv"
        unc.to_csv(path, sep="\t", index=False)
        written.append(path)
    xlsx = f"{prefix}.xlsx"
    def _safe(df):
        # a cell starting with = + - @ would be read as a formula by spreadsheet programs
        def fix(v):
            return "'" + v if isinstance(v, str) and v[:1] in ("=", "+", "-", "@") else v
        return df.apply(lambda col: col.map(fix) if col.dtype == object else col)

    with pd.ExcelWriter(xlsx, engine="openpyxl") as xw:
        for name in ("conditions", "primers", "amplicons", "dimers", "products", "offtarget_sites",
                     "binding_vs_temperature"):
            _safe(tabs[name]).to_excel(xw, sheet_name=name[:31], index=False)
    written.append(xlsx)
    # FASTA and BED
    tmap = {t.id: t for t in d.targets}
    from .design import SeqStore

    store = SeqStore(d.genome["fasta"], d.genome.get("fai"))
    fa = f"{prefix}_amplicons.fasta"
    bed = f"{prefix}.bed"
    with open(fa, "w") as fh, open(bed, "w") as bh:
        for tid in sorted(d.chosen, key=lambda x: (d.pools[x], tmap[x].seqid, tmap[x].start)):
            pr = d.chosen[tid]
            if fasta_with_primers:
                a, b = pr.start, pr.end
            else:
                a, b = pr.insert
            fetch = lambda x, y, s=pr.fwd.seqid: store.region(s, x, y)[0]
            seq = fetch(a, pr.gap[0] - 1) + fetch(pr.gap[1] + 1, b) if pr.gap else fetch(a, b)
            note = f" junction deletion={pr.gap[0]}-{pr.gap[1]} (deleted allele shown)" if pr.gap else ""
            fh.write(f">{tid} pool={d.pools[tid]} {pr.fwd.seqid}:{a}-{b} len={len(seq)} "
                     f"covers={tmap[tid].name}{note}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")
            bh.write(f"{pr.fwd.seqid}\t{pr.start - 1}\t{pr.end}\t{tid}_pool{d.pools[tid]}"
                     f"{'_junction' if pr.gap else ''}\t0\t+\t{pr.insert[0] - 1}\t{pr.insert[1]}\n")
    written += [fa, bed]
    written += write_order_sheets(tabs["primers"], prefix)
    from . import mapfig

    png = mapfig.draw(d, f"{prefix}_map.png")
    if png:
        written.append(png)
    md = f"{prefix}_report.md"
    with open(md, "w") as fh:
        fh.write(markdown(d, tabs, prefix))
    written.append(md)
    return written


def order_sheet(primers: pd.DataFrame) -> pd.DataFrame:
    """The 4-column order sheet (forward name, forward sequence, reverse name,
    reverse sequence; one amplicon per row) that most oligo vendors accept."""
    rows = []
    for tid, grp in primers.groupby("target", sort=False):
        f = grp[grp["direction"] == "forward"].iloc[0]
        r = grp[grp["direction"] == "reverse"].iloc[0]
        rows.append({"fwd_name": f["primer"], "fwd_seq": str(f["order_sequence"]).upper(),
                     "rev_name": r["primer"], "rev_seq": str(r["order_sequence"]).upper()})
    return pd.DataFrame(rows, columns=["fwd_name", "fwd_seq", "rev_name", "rev_seq"])


def write_order_sheets(primers: pd.DataFrame, prefix: str) -> List[str]:
    """`<prefix>_order_sheet.tsv` for every design; `<prefix>_order_sheet_new.tsv`
    (only primers that are not kept from a previous design) after a redesign."""
    out = []
    path = f"{prefix}_order_sheet.tsv"
    order_sheet(primers).to_csv(path, sep="\t", index=False)
    out.append(path)
    if "kept_from_previous" in primers and primers["kept_from_previous"].any():
        new = primers[~primers["kept_from_previous"].astype(bool)]
        path = f"{prefix}_order_sheet_new.tsv"
        order_sheet(new).to_csv(path, sep="\t", index=False)
        out.append(path)
    return out


def markdown(d: Design, tabs, prefix: str = "") -> str:
    L = ["# Multiplex primer design report", ""]
    L.append("## PCR conditions")
    L.append("")
    L.append("| Item | Value |")
    L.append("|---|---|")
    for k, v in conditions(d):
        L.append(f"| {k} | {v} |")
    L.append("")
    if d.warnings:
        L.append("## Warnings")
        L.append("")
        for w in d.warnings:
            L.append(f"- {w}")
        L.append("")
    L.append("## Coverage")
    L.append("")
    n_req = len(d.sites) + len(d.dropped_sites)
    if d.coverage_problems or d.dropped_sites:
        L.append(f"{len(d.sites) - len(d.coverage_problems)} of {n_req} requested positions are covered.")
        for s, gaps in d.coverage_problems:
            L.append(f"- NOT covered: {s.label()} {gaps}")
        for s, why in zip(d.dropped_sites, _reasons(d)):
            L.append(f"- NOT designed: {s.label()}: {why}")
        L.append("")
        L.append("The same list is in `_uncovered.tsv`.")
    else:
        L.append(f"All {n_req} requested positions lie between the primers of at least one amplicon.")
    L.append("")
    ext = [n_ for n_ in d.balance_notes if "not extended" not in n_]
    blocked = [n_ for n_ in d.balance_notes if "not extended" in n_ and "overshoots" not in n_]
    routine = [n_ for n_ in d.balance_notes if "overshoots" in n_]
    if ext:
        L.append("Binding balance: these primers keep their Tm-matched core (uppercase) and carry extra "
                 "template-matching bases at the 5' end (lowercase), because they would otherwise bind more "
                 "weakly than the rest of the pool at the annealing temperature (usually a folded template site). "
                 "The 3' end and the covered positions are unchanged.")
        L.append("")
        for n_ in ext:
            L.append(f"- {n_}")
        L.append("")
    if blocked:
        L.append("Weak primers that could not be extended:")
        L.append("")
        for n_ in blocked:
            L.append(f"- {n_}")
        L.append("")
    if routine:
        L.append(f"{len(routine)} slightly weak primer(s) were left as they are because one more base would "
                 f"overshoot the pool median.")
        L.append("")
    if d.redesign_notes:
        L.append("Redesign:")
        L.append("")
        for n_ in d.redesign_notes:
            L.append(f"- {n_}")
        L.append("")
    if d.relaxed:
        L.append("Targets designed with fallbacks:")
        L.append("")
        for k, v in d.relaxed.items():
            L.append(f"- {k}: {v}")
        L.append("")
    L.append("## Amplicon map")
    L.append("")
    L.append(f"![Amplicon map]({os.path.basename(prefix)}_map.png)" if prefix else "(map not written)")
    L.append("")
    L.append("Top: amplicons along the genome by pool. Below: each amplicon with its primers (arrows), the insert "
             "the reads cover (coloured), the requested positions (red ticks) and, for a junction amplicon, the "
             "spanned deletion (dashed).")
    L.append("")
    L.append("## Primers")
    L.append("")
    cols = ["pool", "primer", "sequence", "length", "gc_percent", "tm_core_C", "tm_C", "dg_eff", "fraction_bound_model",
            "site_open_cost",
            "offtarget_primable_sites"]
    L.append(tabs["primers"][cols].to_markdown(index=False))
    L.append("")
    L.append("## Amplicons")
    L.append("")
    cols = ["pool", "target", "covers", "length", "gc_percent", "product_tm_C", "self_priming_top",
            "self_priming_bottom", "strongest_hairpin_dg_at_72C", "flags"]
    if (tabs["amplicons"]["kind"] == "junction").any():
        cols = cols[:3] + ["kind", "deletion", "wildtype_product_bp"] + cols[3:]
    L.append(tabs["amplicons"][cols].to_markdown(index=False))
    L.append("")
    L.append("## Strongest dimers per pool (37 C)")
    L.append("")
    dm = tabs["dimers"]
    if len(dm):
        L.append(dm.groupby("pool").head(5).to_markdown(index=False))
    L.append("")
    pr = tabs["products"]
    L.append("## Other products (genome-wide, within pools)")
    L.append("")
    L.append("kind: 'outer product' = forward primer of one designed amplicon with the reverse primer of a "
             "neighbouring one (normal in tiled pools, suppressed by a short extension); 'repeated locus' = the "
             "designed pair also amplifies another copy of the locus; 'unintended' = off-target priming.")
    L.append("")
    if not d.cfg.specificity:
        L.append("Not checked (--no_specificity): off-target sites and unintended products were not searched.")
    else:
        L.append(pr.to_markdown(index=False) if len(pr) else "None predicted.")
    L.append("")
    L.append("## Notes")
    L.append("")
    L.append("- tm_C is on the scale named in the conditions table; tm_neb_q5_scale_C reproduces the NEB Tm "
             "calculator for Q5 at 500 nM for comparison. dG_eff is the binding free energy "
             "at the recommended Ta (SantaLucia 2004 + Owczarzy 2008 in the assumed buffer) plus the cost of the "
             "primer's own hairpin and of opening template structure at its site; the pool is optimised for a "
             "narrow dG_eff spread.")
    L.append("- Model fraction bound ignores the polymerase. For Taq-type enzymes it is a direct estimate; Q5 and "
             "Phusion carry a DNA-binding domain and are run above the primer Tm by vendor rule, so for them use "
             "dG_eff to compare primers, not as an absolute occupancy. Buffer monovalent salt is assumed (not "
             "published by vendors), so absolute values carry about 1 to 3 C uncertainty.")
    L.append("- Structure, dimer and off-target thresholds are model-based starting points, not calibrated on "
             "sequencing data. See docs/DESIGN.md.")
    return "\n".join(L) + "\n"
