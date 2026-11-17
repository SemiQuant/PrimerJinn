"""Sequencing feedback: rebalance primer concentrations and pick targets to redesign.

Input is one small table of read counts per amplicon that any pipeline can produce:

    amplicon            reads
    T01_gyrB_N499       12400
    T02_gyrA_A90        3100
    ...
    unaligned           5200      reads that match no designed amplicon (optional, but recommended)

`amplicon` is the primerJinn target id (T01), the primer name without its _F/_R, or the target
name; `reads` is the number of reads that are the designed product. Optional columns: `conc_nM`,
the concentration each primer of that amplicon was used at (so a weak amplicon whose primers were
dilute is raised rather than redesigned), and `sample` for one row per sample and amplicon (folds
are pooled over samples). Rows named unaligned
(or unmapped, off_target, other, dimer, background) are the reads that are not from a designed
amplicon; they are not amplicons but are needed to see whether the library is mostly on target.

Every function here is pure (no Streamlit); the CLI and the GUI call them.
"""
from __future__ import annotations

import io
import math
import os
import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

ALIASES = {
    "amplicon": ("amplicon", "unit", "target", "amplicon_id", "name", "region", "primer_pair", "pair"),
    "reads": ("reads", "read_count", "n_reads", "count", "counts", "n", "depth"),
    "sample": ("sample", "barcode", "sample_id", "run"),
    "conc_nM": ("conc_nm", "conc", "concentration", "nm", "primer_nm", "conc_used", "used_nm"),
    "primer": ("primer", "name", "oligo"),
}
# rows of the counts table that are not amplicons: reads that match no designed product
OTHER_ROWS = {
    "unaligned": "unaligned", "unmapped": "unaligned", "no_hit": "unaligned", "nohit": "unaligned",
    "off_target": "off_target", "off-target": "off_target", "offtarget": "off_target", "off": "off_target",
    "other": "off_target", "background": "off_target", "chimera": "off_target", "mispriming": "off_target",
    "dimer": "dimer", "dimers": "dimer", "primer_dimer": "dimer", "primer_dimers": "dimer",
}


@dataclass
class Thresholds:
    redesign_below: float = 0.2     # fold vs the median below which a target is weak
    conc_max: float = 500.0         # nM: a weak amplicon whose primers would need more than this is redesigned
    alpha: float = 0.5              # damping: factor = (1 / fold) ** alpha
    cap: float = 4.0                # largest change of one primer's concentration
    reaction_ul: float = 25.0
    stock_uM: float = 10.0
    min_reads: int = 200            # fewer reads per amplicon on average: rebalancing unreliable
    other_warn: float = 0.3         # warn when more than this fraction of reads is not a designed amplicon


@dataclass
class RunData:
    amplicons: pd.DataFrame                 # normalised: sample, amplicon, reads[, conc_nM]
    units: pd.DataFrame                     # pooled over samples: amplicon, reads, fold, status, samples, conc_nM
    other: Dict[str, int] = field(default_factory=dict)   # unaligned / off_target / dimer -> reads
    source: str = ""
    notes: List[str] = field(default_factory=list)


@dataclass
class PrevDesign:
    primers: pd.DataFrame                   # the previous `_primers.tsv`
    conc: Dict[str, float]                  # primer -> nM actually used
    default_nM: float
    source: str = ""
    conc_source: str = "assumed"


@dataclass
class Feedback:
    primers: pd.DataFrame                   # one row per primer (the rebalance table)
    targets: pd.DataFrame                   # one row per target
    replace: List[str]
    summary: List[Tuple[str, str]]
    warnings: List[str]
    thresholds: Thresholds


# ----------------------------------------------------------------------------
# reading
# ----------------------------------------------------------------------------
def _norm(df: pd.DataFrame, wanted: Iterable[str]) -> pd.DataFrame:
    """Rename alias columns to their canonical names (case-insensitive)."""
    low = {c.lower().strip(): c for c in df.columns}
    ren = {}
    for canon in wanted:
        if canon in df.columns:
            continue
        for alias in ALIASES.get(canon, (canon,)):
            if alias in low:
                ren[low[alias]] = canon
                break
    return df.rename(columns=ren)


def read_table(path: str) -> pd.DataFrame:
    """A TSV or CSV with one header row (a UTF-8 BOM and '#' comment lines are fine)."""
    with open(path, "rb") as fh:
        text = fh.read().decode("utf-8-sig", errors="replace")
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    if not lines:
        raise ValueError(f"{path} is empty")
    sep = "\t" if "\t" in lines[0] else "," if "," in lines[0] else r"\s+"
    return pd.read_csv(io.StringIO("\n".join(lines)), sep=sep, engine="python")


def status_of(fold: float, reads: float) -> str:
    """good 0.5-2x of the median, warn 0.2-5x, bad beyond, drop = no reads."""
    if reads == 0 or not np.isfinite(fold):
        return "drop"
    if 0.5 <= fold <= 2:
        return "good"
    if 0.2 <= fold <= 5:
        return "warn"
    return "bad"


def pool_units(amp: pd.DataFrame) -> pd.DataFrame:
    """One row per amplicon over all samples. fold per sample is (reads + 0.5) / (median over
    the amplicons with reads + 0.5), pooled over samples as the geometric mean; reads are summed."""
    amp = amp.copy()
    amp["reads"] = pd.to_numeric(amp["reads"], errors="coerce").fillna(0).astype(float)
    rows = []
    for unit, g in amp.groupby("amplicon", sort=False):
        lf = []
        for s, gs in g.groupby("sample"):
            live = amp[(amp["sample"] == s) & (amp["reads"] > 0)]
            if len(live):
                med = float(np.median(live["reads"]))
                lf.append(math.log(max((float(gs["reads"].sum()) + 0.5) / (med + 0.5), 1e-6)))
        reads = float(g["reads"].sum())
        fold = math.exp(float(np.mean(lf))) if lf else float("nan")
        if reads == 0:
            fold = 0.0
        conc = float("nan")
        if "conc_nM" in g.columns:
            c = pd.to_numeric(g["conc_nM"], errors="coerce").dropna()
            if len(c):
                conc = float(c.iloc[0])
        rows.append(dict(amplicon=unit, reads=int(reads), fold=round(fold, 4), status=status_of(fold, reads),
                         samples=int(g["sample"].nunique()), conc_nM=conc))
    return pd.DataFrame(rows, columns=["amplicon", "reads", "fold", "status", "samples", "conc_nM"])


def load_run(path) -> RunData:
    """Read the counts table (amplicon, reads[, sample]; rows named unaligned etc. set aside)."""
    if isinstance(path, (list, tuple)):
        if len(path) != 1:
            raise ValueError("give one counts table (amplicon, reads)")
        path = path[0]
    path = str(path)
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    df = _norm(read_table(path), ("amplicon", "reads", "sample", "conc_nM"))
    missing = [c for c in ("amplicon", "reads") if c not in df.columns]
    if missing:
        raise ValueError(f"{os.path.basename(path)} lacks column(s) {', '.join(missing)}; expected a header with "
                         "amplicon and reads (accepted names: "
                         + "; ".join(f"{k}: {', '.join(v)}" for k, v in ALIASES.items() if k in missing) + ")")
    if "sample" not in df.columns:
        df["sample"] = "run"
    df["amplicon"] = df["amplicon"].astype(str).str.strip()
    df["reads"] = pd.to_numeric(df["reads"], errors="coerce")
    bad = df[df["reads"].isna()]
    if len(bad):
        raise ValueError(f"{os.path.basename(path)}: reads is not a number for {', '.join(bad['amplicon'].head(5))}")
    key = df["amplicon"].str.lower().str.replace(" ", "_")
    is_other = key.isin(OTHER_ROWS)
    other: Dict[str, int] = {}
    for k, r in zip(key[is_other], df["reads"][is_other]):
        other[OTHER_ROWS[k]] = other.get(OTHER_ROWS[k], 0) + int(r)
    cols = ["sample", "amplicon", "reads"] + (["conc_nM"] if "conc_nM" in df.columns else [])
    amp = df[~is_other][cols].reset_index(drop=True)
    if amp.empty:
        raise ValueError(f"{os.path.basename(path)} has no amplicon rows")
    notes: List[str] = []
    if not other:
        notes.append("no unaligned row in the counts table: the fraction of reads from the designed amplicons "
                     "is unknown (add a row named unaligned with the reads that match no amplicon)")
    return RunData(amplicons=amp, units=pool_units(amp), other=other, source=os.path.basename(path), notes=notes)


def load_previous(primers_tsv: str, conc_tsv: Optional[str] = None, default_nM: Optional[float] = None) -> PrevDesign:
    """The previous design's primers table and the concentrations actually used."""
    from .redesign import read_primers_table

    df = read_primers_table(primers_tsv)
    df["primer"] = df["primer"].astype(str)
    if default_nM is None:
        default_nM = 100.0
    cond = re.sub(r"_primers\.(tsv|xlsx)$", "_conditions.tsv", str(primers_tsv))
    if cond != str(primers_tsv) and os.path.exists(cond):   # the previous design's own setting
        try:
            c = pd.read_csv(cond, sep="\t")
            v = c[c["item"] == "primer concentration"]["value"]
            if len(v):
                default_nM = float(re.findall(r"[\d.]+", str(v.iloc[0]))[0])
        except Exception:
            pass
    conc: Dict[str, float] = {}
    source = "assumed"
    if conc_tsv:
        c = _norm(read_table(conc_tsv), ("primer",))
        vcol = next((x for x in c.columns if x != "primer"), None)
        if "primer" not in c.columns or vcol is None:
            raise ValueError(f"{conc_tsv}: expected columns primer and conc_nM")
        for _, r in c.iterrows():
            try:
                conc[str(r["primer"]).strip()] = float(r[vcol])
            except (TypeError, ValueError):
                continue
        source = "given"
    return PrevDesign(primers=df, conc=conc, default_nM=float(default_nM), source=str(primers_tsv),
                      conc_source=source)


# ----------------------------------------------------------------------------
# matching the run to the design
# ----------------------------------------------------------------------------
def _strip_side(name: str) -> str:
    return re.sub(r"_(F|R|fwd|rev|FWD|REV|L)$", "", name)


def match(run: RunData, prev: PrevDesign) -> Tuple[Dict[str, str], List[str]]:
    """target id -> amplicon name in the run, and warnings. A run amplicon is recognised by the
    target id, the primer name without _F/_R, or the target name (case-insensitive); 'A/B'
    names two merged amplicons."""
    warnings: List[str] = []
    targets = list(dict.fromkeys(prev.primers["target"].astype(str)))
    keys: Dict[str, str] = {}
    for _, r in prev.primers.iterrows():
        tid = str(r["target"])
        for k in (tid, _strip_side(str(r["primer"])), str(r["primer"]),
                  str(r["target_name"]) if "target_name" in r and pd.notna(r["target_name"]) else ""):
            if k:
                keys.setdefault(k.lower(), tid)
    unit_of: Dict[str, str] = {}
    unmatched: List[str] = []
    for unit in run.units["amplicon"].astype(str):
        hit = set()
        for nm in unit.split("/"):
            nm = nm.strip().lower()
            tid = keys.get(nm) or keys.get(_strip_side(nm))
            if tid:
                hit.add(tid)
        if not hit:
            unmatched.append(unit)
        for t in hit:
            unit_of[t] = unit
    if unmatched:
        warnings.append(f"{len(unmatched)} amplicon(s) in the counts table match no target of the previous design "
                        f"(use the target id, the primer name without _F/_R, or the target name): "
                        + ", ".join(unmatched[:8]) + (" ..." if len(unmatched) > 8 else ""))
    missing = [t for t in targets if t not in unit_of]
    if missing:
        warnings.append(f"{len(missing)} target(s) of the previous design have no row in the counts table "
                        f"(treated as no reads): {', '.join(missing[:8])}" + (" ..." if len(missing) > 8 else ""))
    return unit_of, warnings


# ----------------------------------------------------------------------------
# decisions
# ----------------------------------------------------------------------------
def factor(fold: float, reads: float, th: Thresholds) -> float:
    """Damped concentration change: (1 / fold) ^ alpha, within 1/cap .. cap."""
    if reads == 0 or not np.isfinite(fold) or fold <= 0:
        return th.cap
    return float(np.clip((1.0 / fold) ** th.alpha, 1.0 / th.cap, th.cap))


def analyse(run: RunData, prev: PrevDesign, th: Optional[Thresholds] = None, replace: Iterable[str] = (),
            keep: Iterable[str] = ()) -> Feedback:
    th = th or Thresholds()
    unit_of, warnings = match(run, prev)
    warnings = list(run.notes) + warnings
    units = {str(r["amplicon"]): r for _, r in run.units.iterrows()}
    force_replace = {str(x) for x in replace}
    force_keep = {str(x) for x in keep}
    targets = list(dict.fromkeys(prev.primers["target"].astype(str)))
    unknown = sorted((force_replace | force_keep) - set(targets))
    if unknown:
        raise ValueError(f"target id(s) not in {prev.source}: {', '.join(unknown)}")

    designed = int(run.units["reads"].sum())
    other = int(sum(run.other.values()))
    total = designed + other
    if run.other and total and other / total > th.other_warn:
        warnings.append(f"{100 * other / total:.0f}% of reads are not from a designed amplicon "
                        f"({', '.join(f'{k} {v:,}' for k, v in run.other.items())}): check the library "
                        f"(primer dimers, annealing temperature) before trusting the balance")

    trows, prows = [], []
    replaced: List[str] = []
    n_within = n_live = 0
    for tid in targets:
        grp = prev.primers[prev.primers["target"].astype(str) == tid]
        unit = unit_of.get(tid)
        u = units.get(unit) if unit else None
        reads = int(u["reads"]) if u is not None else 0
        fold = float(u["fold"]) if u is not None else 0.0
        st = str(u["status"]) if u is not None else "missing"
        if u is not None:
            n_live += 1
            n_within += st == "good"
        # the concentration the pair was used at: --conc per primer, else the counts table, else the default
        names = [str(x) for x in grp["primer"]]
        unit_conc = float(u["conc_nM"]) if u is not None and np.isfinite(float(u["conc_nM"])) else None
        cur_of = {nm: prev.conc.get(nm, unit_conc if unit_conc is not None else prev.default_nM) for nm in names}
        cur_pair = max(cur_of.values())
        known = any(nm in prev.conc for nm in names) or unit_conc is not None
        reasons: List[str] = []
        notes: List[str] = []
        if st == "drop":
            reasons.append("no reads")
        elif st == "missing":
            reasons.append("no row in the counts table")
        elif fold < th.redesign_below:
            need = cur_pair / fold       # the concentration that would bring it to the median, undamped
            if need > th.conc_max:
                reasons.append(f"{fold:.2f}x the median amplicon at {cur_pair:g} nM: it would need about "
                               f"{need:.0f} nM, above {th.conc_max:g}")
            else:
                notes.append(f"{fold:.2f}x the median at {cur_pair:g} nM" + ("" if known else " (assumed)")
                             + ": raise the concentration before redesigning")
        elif fold < 0.5 and cur_pair >= th.conc_max:
            reasons.append(f"{fold:.2f}x the median although already at {cur_pair:g} nM")
        action = "replace" if reasons else "rebalance"
        if tid in force_keep:
            action = "rebalance"
            reasons = [r + " (kept on request)" for r in reasons]
        elif tid in force_replace:
            action = "replace"
            reasons = reasons or ["replaced on request"]
        if action == "replace":
            replaced.append(tid)
        fac = factor(fold, reads, th)
        pool = int(grp["pool"].iloc[0]) if "pool" in grp else 1
        trows.append(dict(target=tid, target_name=str(grp["target_name"].iloc[0]) if "target_name" in grp else "",
                          amplicon=unit or "", pool=pool, reads=reads, fold=round(fold, 4), status=st,
                          current_nM=round(cur_pair, 2), action=action, factor=round(fac, 3),
                          reason="; ".join(reasons or notes)))
        for _, r in grp.iterrows():
            nm = str(r["primer"])
            cur = cur_of[nm]
            if action == "replace":
                sug, pact, ul = float("nan"), "replace", float("nan")
            else:
                sug = min(cur * fac, th.conc_max) if fac > 1 else cur * fac
                pact = "raise" if fac > 1.05 else "lower" if fac < 0.95 else "keep"
                ul = sug * th.reaction_ul / (th.stock_uM * 1000.0)
            prows.append(dict(target=tid, amplicon=unit or "", pool=pool, primer=nm, direction=str(r["direction"]),
                              reads=reads, fold=round(fold, 3), status=st,
                              action=pact, reason="; ".join(reasons) if action == "replace" else
                              (f"{fold:.2f}x the median, concentration x{fac:.2f}" if u is not None else ""),
                              current_nM=round(cur, 2), suggested_nM=round(sug, 1) if np.isfinite(sug) else None,
                              ul_stock_per_rxn=round(ul, 3) if np.isfinite(ul) else None,
                              conc_source="given" if nm in prev.conc else "counts table" if unit_conc is not None
                              else "assumed"))
    primers = pd.DataFrame(prows)
    targets_df = pd.DataFrame(trows)

    # sanity on depth and the total oligo load
    if n_live and designed / max(n_live, 1) < th.min_reads:
        warnings.append(f"only {designed / max(n_live, 1):.0f} reads per amplicon on average "
                        f"(fewer than {th.min_reads}): fold changes are noisy, rebalance with care")
    kept = primers[primers["action"] != "replace"]
    if len(kept):
        tot_ul = float(kept["ul_stock_per_rxn"].sum())
        if tot_ul > 0.5 * th.reaction_ul:
            warnings.append(f"the suggested primers add up to {tot_ul:.1f} uL of {th.stock_uM:g} uM stock per "
                            f"{th.reaction_ul:g} uL reaction: make a concentrated pool or lower every "
                            f"concentration by the same factor")
    summary: List[Tuple[str, str]] = []
    if run.other:
        summary.append(("total reads", f"{total:,}"))
        summary.append(("designed reads %", f"{100 * designed / total:.1f}" if total else "0"))
        summary.append(("unaligned or off-target reads %", f"{100 * other / total:.1f}" if total else "0"))
        for k, v in run.other.items():
            summary.append((f"{k} reads", f"{v:,}"))
    else:
        summary.append(("total reads", f"{designed:,} (designed amplicons only; no unaligned row)"))
    summary.append(("amplicons within 2x of the median", f"{n_within}/{n_live}"))
    summary.append(("targets replaced", f"{len(replaced)}" + (f" ({', '.join(replaced)})" if replaced else "")))
    summary.append(("primers raised / lowered / kept",
                    " / ".join(str(int((primers['action'] == a).sum())) for a in ("raise", "lower", "keep"))))
    src = sorted(set(primers["conc_source"])) if len(primers) else []
    summary.append(("concentration source", ", ".join(src).replace("assumed", f"assumed {prev.default_nM:g} nM")))
    summary.append(("rule", f"new = current x (median / reads) ^ {th.alpha:g}, at most x{th.cap:g} either way and "
                            f"at most {th.conc_max:g} nM; replace with no reads, or below {th.redesign_below:g}x "
                            f"the median when bringing it to the median would take more than {th.conc_max:g} nM"))
    summary.append(("source", run.source))
    for w in warnings:
        summary.append(("warning", w))
    return Feedback(primers=primers, targets=targets_df, replace=replaced, summary=summary, warnings=warnings,
                    thresholds=th)


# ----------------------------------------------------------------------------
# writing
# ----------------------------------------------------------------------------
def conc_table(fb: Feedback, new_primers: Iterable[str] = (), default_nM: float = 100.0) -> pd.DataFrame:
    """`primer conc_nM`: suggested values for kept primers, the design default for new ones."""
    rows = [(r["primer"], r["suggested_nM"]) for _, r in fb.primers.iterrows() if r["action"] != "replace"]
    rows += [(p, default_nM) for p in new_primers]
    return pd.DataFrame(rows, columns=["primer", "conc_nM"])


def markdown(fb: Feedback) -> str:
    L = ["# Sequencing feedback", "", "| Item | Value |", "|---|---|"]
    for k, v in fb.summary:
        L.append(f"| {k} | {v} |")
    L += ["", "## Targets", "", fb.targets.to_markdown(index=False), "", "## Primers", "",
          fb.primers.drop(columns=["reason"]).to_markdown(index=False), "",
          "Actions: keep, raise or lower the concentration (suggested_nM, ul_stock_per_rxn of the stock per "
          "reaction); replace = redesign this amplicon's primers (see the reason column).", ""]
    return "\n".join(L)


def write(fb: Feedback, prefix: str, new_primers: Iterable[str] = (), default_nM: float = 100.0) -> List[str]:
    os.makedirs(os.path.dirname(os.path.abspath(prefix)), exist_ok=True)
    out = []
    p = f"{prefix}_rebalance.tsv"
    fb.primers.to_csv(p, sep="\t", index=False)
    out.append(p)
    p = f"{prefix}_rebalance_targets.tsv"
    fb.targets.to_csv(p, sep="\t", index=False)
    out.append(p)
    p = f"{prefix}_conc.tsv"
    conc_table(fb, new_primers, default_nM).to_csv(p, sep="\t", index=False)
    out.append(p)
    p = f"{prefix}_feedback.tsv"
    pd.DataFrame(fb.summary, columns=["item", "value"]).to_csv(p, sep="\t", index=False)
    out.append(p)
    p = f"{prefix}_feedback.md"
    with open(p, "w") as fh:
        fh.write(markdown(fb))
    out.append(p)
    return out
