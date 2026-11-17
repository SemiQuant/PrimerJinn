"""Known strain variants to keep out from under primers.

Read from a VCF (optionally gzipped), a BED file, or a TSV/CSV/TXT of positions
(any format genome.parse_positions accepts: '7570', 'NC_000962.3:7570', ranges,
tables with seqid/pos/start/end columns). Positions are 1-based; a VCF deletion
covers every reference base of its REF allele.

Rule (config): a primer with a known variant in its 3'-terminal
`variant_3prime_nt` bases is rejected (a mismatch there blocks or biases
extension); a variant elsewhere under the primer costs `variant_penalty` each.
"""
from __future__ import annotations

import bisect
import re
from typing import Dict, List, Tuple

MAX_RANGE = 5000   # longest single range expanded into masked positions


class VariantMask:
    """Sorted masked positions per sequence id."""

    def __init__(self, positions: Dict[str, List[int]], source: str = ""):
        self.pos = {k: sorted(set(v)) for k, v in positions.items()}
        self.source = source

    def __len__(self) -> int:
        return sum(len(v) for v in self.pos.values())

    def within(self, seqid: str, start: int, end: int) -> List[int]:
        """Masked positions in [start, end] (1-based, inclusive)."""
        v = self.pos.get(seqid)
        if not v:
            return []
        i = bisect.bisect_left(v, start)
        j = bisect.bisect_right(v, end)
        return v[i:j]

    def primer_hits(self, side: str, seqid: str, start: int, end: int) -> List[int]:
        """Distances from the primer 3' end (0 = 3'-terminal base) of masked
        positions under a primer occupying [start, end] on the top strand."""
        hits = self.within(seqid, start, end)
        if side == "F":
            return [end - x for x in hits]
        return [x - start for x in hits]


def _af_ok(info: str, min_af: float) -> bool:
    """Keep a record unless every numeric AF value (one per ALT) is below min_af; AF '.' or absent keeps it."""
    m = re.search(r"(?:^|;)AF=([^;]+)", info)
    if not m:
        return True
    vals = []
    for v in m.group(1).split(","):
        try:
            vals.append(float(v))
        except ValueError:
            pass
    return not vals or max(vals) >= min_af


def _read_vcf(path, lookup, min_af: float) -> Tuple[Dict[str, List[int]], List[str]]:
    """Masked positions per sequence, and notes. Each note about records not used contains
    'ignored', which the input check shows as a warning."""
    from .genome import read_vcf

    out: Dict[str, List[int]] = {}
    skipped: List[str] = []
    bad: List[str] = []
    recs, _ = read_vcf(path, bad)
    long_sv, ref_blocks = [], 0
    for r in recs:
        if r["ref_block"]:          # gVCF reference block or no-call: not a variant
            ref_blocks += 1
            continue
        s = lookup.get(r["chrom"].lower())
        if s is None:
            skipped.append(r["chrom"])
            continue
        if min_af > 0 and not _af_ok(r["info"], min_af):
            continue
        symbolic = any(a.startswith("<") for a in r["alt"].split(","))
        if symbolic and r["end"] - r["pos"] + 1 > MAX_RANGE:
            long_sv.append(f"VCF line {r['line']} ({r['pos']}-{r['end']})")
            continue
        out.setdefault(s["seqid"], []).extend(range(r["pos"], r["end"] + 1))
    notes = []
    if skipped:
        u = sorted(set(skipped))
        notes.append(f"variants on unknown sequences ignored: {', '.join(u[:5])}{' ...' if len(u) > 5 else ''}")
    if bad:
        notes.append(f"{len(bad)} malformed VCF line(s) ignored, first: {bad[0]}")
    if long_sv:
        notes.append(f"{len(long_sv)} symbolic structural variant(s) longer than {MAX_RANGE} bp ignored: "
                     f"{', '.join(long_sv[:3])}")
    if ref_blocks:
        notes.append(f"{ref_blocks} gVCF reference block(s) or no-call record(s) skipped (no alternative allele)")
    return out, notes


def load_variants(path: str, genome_meta: dict, min_af: float = 0.0) -> Tuple[VariantMask, List[str]]:
    """Returns (mask, notes). Raises ValueError with a clear message on bad input."""
    from . import genome as gmod

    lookup = gmod._seq_lookup(genome_meta)
    notes: List[str] = []
    if gmod.is_vcf(path):
        pos, vnotes = _read_vcf(path, lookup, min_af)
        notes.extend(vnotes)
    else:
        rows = gmod.parse_positions(path, genome_meta, allow_ambiguous_default=False, max_span=MAX_RANGE)
        pos = {}
        for r in rows:
            pos.setdefault(r["seqid"], []).extend(range(r["start"], r["end"] + 1))
    mask = VariantMask(pos, str(path))
    if not len(mask):
        why = "; ".join(notes)
        raise ValueError(f"no variant positions read from {path}" + (f" ({why})" if why else ""))
    notes.append(f"{len(mask)} known variant position(s) masked from {path}")
    return mask, notes
