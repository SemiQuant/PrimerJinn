"""Redesign part of a pool while keeping the primers that work.

A previous design's `<out>_primers.tsv` names every primer with its target id,
coordinates and sequence. Kept targets become fixed pairs: the designer scores
them with the normal machinery (binding, structure, dimers, products,
specificity) but never replaces or extends them, so new primers are chosen to
fit the existing pool. Replaced targets are designed afresh, and the 3' ends of
the primers that failed are excluded for them.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Set, Tuple

import pandas as pd

REQUIRED = ["target", "direction", "sequence", "order_sequence", "added_5prime_nt", "seqid", "start", "end",
            "three_prime_end", "pool"]


@dataclass
class FixedPrimer:
    seq: str            # full sequence, upper case, extension included
    side: str           # 'F' or 'R'
    seqid: str
    start: int
    end: int
    ext5: int = 0
    name: str = ""


@dataclass
class FixedPair:
    target_id: str
    fwd: FixedPrimer
    rev: FixedPrimer
    pool: int
    gap: Optional[Tuple[int, int]] = None    # junction amplicon: the deletion spanned


@dataclass
class Previous:
    fixed: Dict[str, FixedPair] = field(default_factory=dict)
    avoid_3prime: Set[Tuple[str, str, int]] = field(default_factory=set)      # (seqid, side, three_prime)
    avoid_sites: List[Tuple[str, str, int, int]] = field(default_factory=list)  # (seqid, side, start, end)
    prev_pools: Dict[str, int] = field(default_factory=dict)
    replaced: List[str] = field(default_factory=list)
    source: str = ""


def read_primers_table(path: str) -> pd.DataFrame:
    """The previous `_primers.tsv` (or the primers sheet of the XLSX)."""
    p = str(path)
    if p.lower().endswith((".xlsx", ".xls")):
        df = pd.read_excel(p, sheet_name="primers")
    else:
        df = pd.read_csv(p, sep="\t")
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(f"{path} is not a primerJinn primers table: missing column(s) {', '.join(missing)}")
    return df


def read_previous(path: str, replace: Iterable[str], keep: Optional[Iterable[str]] = None,
                  tail_fwd: str = "", tail_rev: str = "") -> Previous:
    """Fixed pairs for kept targets and the sites to avoid for replaced ones."""
    df = read_primers_table(path)
    ids = list(dict.fromkeys(df["target"].astype(str)))
    replace = {str(x) for x in replace}
    keep = set(ids) - replace if keep is None else {str(x) for x in keep}
    unknown = sorted((replace | keep) - set(ids))
    if unknown:
        raise ValueError(f"target id(s) not in {path}: {', '.join(unknown)} (ids there: {', '.join(ids)})")
    both = sorted(replace & keep)
    if both:
        raise ValueError(f"target(s) both kept and replaced: {', '.join(both)}")
    prev = Previous(source=str(path), replaced=sorted(replace, key=ids.index))
    for tid in ids:
        rows = df[df["target"].astype(str) == tid]
        prim = {}
        for _, r in rows.iterrows():
            side = "F" if str(r["direction"]).lower().startswith("f") else "R"
            seq = str(r["sequence"]).strip().upper()
            prim[side] = FixedPrimer(seq, side, str(r["seqid"]), int(r["start"]), int(r["end"]),
                                     int(r["added_5prime_nt"]) if pd.notna(r["added_5prime_nt"]) else 0,
                                     str(r["primer"]) if "primer" in r else "")
            if tid in keep:
                tail = tail_fwd if side == "F" else tail_rev
                ordered = str(r["order_sequence"]).strip().upper()
                if ordered != tail.upper() + seq:
                    raise ValueError(f"primer {r.get('primer', tid)}: order_sequence {ordered} is not the chosen tails "
                                     f"({tail or 'none'}) plus the sequence; run the redesign with the same --tails "
                                     f"as the previous design")
        if "F" not in prim or "R" not in prim:
            raise ValueError(f"target {tid} in {path} does not have one forward and one reverse primer")
        if tid in keep:
            gap = None
            if "deletion" in rows.columns and pd.notna(rows["deletion"].iloc[0]) and str(rows["deletion"].iloc[0]).strip():
                a, b = str(rows["deletion"].iloc[0]).split("-")
                gap = (int(a), int(b))
            prev.fixed[tid] = FixedPair(tid, prim["F"], prim["R"], int(rows["pool"].iloc[0]), gap)
            prev.prev_pools[tid] = int(rows["pool"].iloc[0])
        else:
            for side, p in prim.items():
                prev.avoid_3prime.add((p.seqid, side, p.end if side == "F" else p.start))
                prev.avoid_sites.append((p.seqid, side, p.start, p.end))
    return prev


def apply(cfg, prev: Previous) -> None:
    """Put a Previous into a Config (the fields the designer reads)."""
    cfg.fixed = dict(prev.fixed)
    cfg.avoid_3prime = set(prev.avoid_3prime)
    cfg.avoid_sites = list(prev.avoid_sites)
    cfg.prev_pools = dict(prev.prev_pools)
