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
import gzip
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


def _open(path):
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)


def _read_vcf(path, lookup, min_af: float) -> Tuple[Dict[str, List[int]], List[str]]:
    out: Dict[str, List[int]] = {}
    skipped: List[str] = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            chrom, pos, ref = f[0], int(f[1]), f[3]
            s = lookup.get(chrom.lower())
            if s is None:
                skipped.append(chrom)
                continue
            if min_af > 0 and len(f) > 7:
                m = re.search(r"(?:^|;)AF=([0-9.eE+-]+)", f[7])
                if m and float(m.group(1)) < min_af:
                    continue
            for x in range(pos, pos + max(1, len(ref))):
                out.setdefault(s["seqid"], []).append(x)
    return out, sorted(set(skipped))


def load_variants(path: str, genome_meta: dict, min_af: float = 0.0) -> Tuple[VariantMask, List[str]]:
    """Returns (mask, notes). Raises ValueError with a clear message on bad input."""
    from . import genome as gmod

    lookup = gmod._seq_lookup(genome_meta)
    notes: List[str] = []
    low = str(path).lower()
    if low.endswith((".vcf", ".vcf.gz", ".vcf.bgz")):
        pos, skipped = _read_vcf(path, lookup, min_af)
        if skipped:
            notes.append(f"variants on unknown sequences ignored: {', '.join(skipped[:5])}")
    else:
        rows = gmod.parse_positions(path, genome_meta, allow_ambiguous_default=False, max_span=MAX_RANGE)
        pos = {}
        for r in rows:
            pos.setdefault(r["seqid"], []).extend(range(r["start"], r["end"] + 1))
    mask = VariantMask(pos, str(path))
    if not len(mask):
        raise ValueError(f"no variant positions read from {path}")
    notes.append(f"{len(mask)} known variant position(s) masked from {path}")
    return mask, notes
