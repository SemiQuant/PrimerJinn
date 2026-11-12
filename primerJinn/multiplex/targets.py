"""Group requested positions into amplicon targets and check coverage.

Positions are read by genome.parse_positions (validated against the genome).
Coordinates are 1-based and inclusive everywhere in this module. A *site* is something the user wants sequenced: a single
position or a range. A *target* is a group of sites that one amplicon must
cover between its primers.

Grouping is an exact dynamic programme over the sorted required bases. It
trades the number of amplicons against how tightly each one is packed (tight
groups leave little room to place primers and push amplicons past the
optimum length), and prefers to cut where neighbouring amplicons will not
have to overlap (overlapping amplicons need separate pools).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Sequence

__all__ = [
    "Site",
    "Target",
    "group_sites",
    "coverage_gaps",
    "span_limits",
]


@dataclass(frozen=True)
class Site:
    """A required genome interval, 1-based inclusive. ``seqid`` is '' until
    it has been resolved against the genome (see ``Genome.resolve_sites``)."""

    seqid: str
    start: int
    end: int
    name: str = ""

    def __post_init__(self):
        if self.start < 1 or self.end < self.start:
            raise ValueError(f"bad interval {self.seqid}:{self.start}-{self.end}")

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def label(self) -> str:
        loc = str(self.start) if self.start == self.end else f"{self.start}-{self.end}"
        loc = f"{self.seqid}:{loc}" if self.seqid else loc
        return f"{self.name} ({loc})" if self.name else loc


@dataclass
class Target:
    """Bases ``start..end`` on ``seqid`` that one amplicon must contain
    between its primers."""

    id: str
    seqid: str
    start: int
    end: int
    sites: list = field(default_factory=list)
    origins: frozenset = field(default_factory=frozenset)   # requested sites this target came from

    @property
    def span(self) -> int:
        return self.end - self.start + 1

    @property
    def name(self) -> str:
        names = []
        for s in self.sites:
            if s.name and s.name not in names:
                names.append(s.name)
        if not names:
            return self.id
        shown = "+".join(names[:3])
        return shown + (f"+{len(names) - 3}more" if len(names) > 3 else "")


# --------------------------------------------------------------------------
# Grouping
# --------------------------------------------------------------------------


def span_limits(opt_len: int, max_len: int, primer_len: int = 22, margin: int = 5, flex: int = 30):
    """Comfortable and hard limits for the required span of one target.

    An amplicon of length L with two primers of ``primer_len`` and ``margin``
    clear bases next to each primer can cover L - 2*(primer_len + margin)
    required bases. ``flex`` bases are held back so the primer search has
    room to move. Also returns the smallest gap between two targets that
    still lets their amplicons sit side by side without overlapping.
    """
    comfort = max(1, opt_len - 2 * (primer_len + margin) - flex)
    hard = max(comfort, max_len - 2 * (primer_len + margin) - flex)
    min_gap = 2 * (primer_len + margin) + 10
    return comfort, hard, min_gap


def _merge(intervals):
    out = []
    for a, b in sorted(intervals):
        if out and a <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [(a, b) for a, b in out]


def group_sites(
    sites: Sequence[Site],
    span_comfort: int = 216,
    span_hard: int = 291,
    min_gap: int = 64,
    atom: int = 10,
    overlap_penalty: float = 0.4,
) -> list:
    """Partition the required bases into targets.

    Exact dynamic programme per sequence. Cost of a solution is
    (number of targets) + (stretch cost for targets wider than
    ``span_comfort``, at most 0.5 each) + ``overlap_penalty`` for every pair
    of neighbouring targets closer than ``min_gap`` (their amplicons will
    overlap and need separate pools). No target is wider than ``span_hard``.
    Ranges longer than ``span_hard`` are tiled.
    """
    if span_hard < span_comfort:
        raise ValueError("span_hard must be >= span_comfort")
    by_seq: dict = {}
    for s in sites:
        by_seq.setdefault(s.seqid, []).append(s)

    targets: list = []
    for seqid, seq_sites in by_seq.items():
        merged = _merge((s.start, s.end) for s in seq_sites)
        atoms = []
        for a, b in merged:
            pos = a
            while pos <= b:
                atoms.append((pos, min(b, pos + atom - 1)))
                pos += atom
        n = len(atoms)
        INF = float("inf")
        best = [INF] * (n + 1)
        back = [0] * (n + 1)
        best[0] = 0.0
        for j in range(1, n + 1):
            end = atoms[j - 1][1]
            i = j
            while i >= 1:
                span = end - atoms[i - 1][0] + 1
                if span > span_hard and i < j:
                    break
                stretch = 0.0
                if span > span_comfort and span_hard > span_comfort:
                    stretch = 0.5 * ((span - span_comfort) / (span_hard - span_comfort)) ** 2
                link = 0.0
                if i > 1 and atoms[i - 1][0] - atoms[i - 2][1] - 1 < min_gap:
                    link = overlap_penalty
                cost = best[i - 1] + 1.0 + stretch + link
                if cost < best[j] - 1e-12:
                    best[j] = cost
                    back[j] = i
                i -= 1
        groups = []
        j = n
        while j > 0:
            i = back[j]
            groups.append((atoms[i - 1][0], atoms[j - 1][1]))
            j = i - 1
        groups.reverse()
        for a, b in groups:
            members = [s for s in seq_sites if s.start <= b and s.end >= a]
            members = sorted(members, key=lambda s: (s.start, s.end))
            targets.append(Target("", seqid, a, b, members,
                                  frozenset((s.seqid, s.start, s.end, s.name) for s in members)))

    width = max(2, len(str(len(targets))))
    for k, t in enumerate(targets, 1):
        t.id = f"T{k:0{width}d}"
    return targets


def coverage_gaps(sites: Sequence[Site], inserts: Iterable) -> list:
    """Required intervals not covered by any insert.

    ``inserts`` is an iterable of (seqid, start, end), 1-based inclusive,
    giving the bases between the primers of each amplicon. Returns a list of
    (site, [(gap_start, gap_end), ...]) for every site that is not fully
    covered. An empty list means every position is covered.
    """
    cover: dict = {}
    for seqid, a, b in inserts:
        if b >= a:
            cover.setdefault(seqid, []).append((a, b))
    cover = {k: _merge(v) for k, v in cover.items()}
    problems = []
    for s in sites:
        gaps = []
        pos = s.start
        for a, b in cover.get(s.seqid, []):
            if b < pos:
                continue
            if a > s.end:
                break
            if a > pos:
                gaps.append((pos, min(a - 1, s.end)))
            pos = max(pos, b + 1)
            if pos > s.end:
                break
        if pos <= s.end:
            gaps.append((pos, s.end))
        if gaps:
            problems.append((s, gaps))
    return problems


def group_for_config(sites: Sequence[Site], cfg) -> list:
    """Targets as the designer builds them for cfg's amplicon and primer lengths: a
    comfortable span from a realistic primer length, the hard limit from the shortest."""
    comfort = span_limits(cfg.amp_opt, cfg.amp_max, primer_len=min(cfg.primer_max, cfg.primer_min + 4),
                          margin=cfg.margin, flex=20)[0]
    _, hard, gap = span_limits(cfg.amp_opt, cfg.amp_max, primer_len=cfg.primer_min, margin=cfg.margin, flex=20)
    return group_sites(sites, comfort, max(hard, comfort), gap)


class OnePoolError(ValueError):
    """The targets cannot share one pool within the largest allowed amplicon."""


def _merge_conflicts(targets: list, cfg, hard: int) -> tuple:
    """Merge targets that may not share a pool (overlapping or too-short outer product)
    into single targets while the merged span fits `hard`. Returns (targets, blocking
    pairs as (target a, target b, merged span))."""
    from .optimize import conflict_graph

    ts = sorted(targets, key=lambda t: (t.seqid, t.start))
    while True:
        edges = conflict_graph(ts, cfg)
        blocking, merged = [], False
        for i in range(len(ts)):
            js = sorted(j for j in edges[i] if j > i)
            if not js:
                continue
            j = js[-1]                                   # merge everything from i to its farthest conflict
            span = ts[j].end - ts[i].start + 1
            if span <= hard:
                group = ts[i:j + 1]
                sites = sorted({id(x): x for t in group for x in t.sites}.values(), key=lambda x: (x.start, x.end))
                origins = frozenset().union(*(t.origins for t in group))
                ts = ts[:i] + [Target("", ts[i].seqid, ts[i].start, max(t.end for t in group), sites, origins)] \
                    + ts[j + 1:]
                merged = True
                break
            blocking.append((ts[i], ts[j], span))
        if not merged:
            return ts, blocking


def one_pool_layout(sites: Sequence[Site], cfg, step: int = 25) -> tuple:
    """Targets for a single-pool design, in two size classes.

    Targets that may not share a tube (overlapping amplicons, or an outer product
    shorter than cfg.outer_min) are merged into one target. Targets that still fit
    the normal window keep it (amp_min..amp_max). Merged targets that need more
    room form a large class with its own window of the same spread, ending at the
    longest length needed (rounded up to `step` bp, at most cfg.one_pool_max), so
    the large amplicons are balanced among themselves and the rest stay short.
    Returns (targets, {target id: (min, opt, max)} for the large class, note).
    Raises OnePoolError naming the targets that cannot be merged."""
    pad = 2 * (cfg.primer_min + cfg.margin) + 20         # primers, margins and search room around a span
    hard_cap = cfg.one_pool_max - pad
    targets, blocking = _merge_conflicts(group_for_config(sites, cfg), cfg, hard_cap)

    def describe(pairs):
        out = []
        for a, b, span in pairs[:10]:
            da = ", ".join(x.label() for x in a.sites[:2]) or f"{a.seqid}:{a.start}-{a.end}"
            db = ", ".join(x.label() for x in b.sites[:2]) or f"{b.seqid}:{b.start}-{b.end}"
            out.append(f"{da} and {db} would need one amplicon of about {span + pad} bp")
        return out

    if blocking:
        raise OnePoolError(
            f"the targets cannot share one pool with amplicons up to {cfg.one_pool_max} bp: "
            + "; ".join(describe(blocking))
            + ". Raise --one_pool_max, lower --outer_min (shortest outer product tolerated in one tube), or drop "
              "--one_pool to use several pools.")
    width = max(2, len(str(len(targets))))
    for k, t in enumerate(targets, 1):
        t.id = f"T{k:0{width}d}"
    need = {t.id: t.span + pad for t in targets}
    large = [t for t in targets if need[t.id] > cfg.amp_max]
    if not large:
        return targets, {}, f"one pool: all amplicons {cfg.amp_min}-{cfg.amp_max} bp (no target needed merging)"
    # size classes: from the smallest need up, one class per spread of lengths, each ending at
    # the longest need in it (rounded up to `step`), so no target is made longer than its class needs
    spread, opt_off = cfg.amp_max - cfg.amp_min, cfg.amp_max - cfg.amp_opt
    windows, classes = {}, []
    todo = sorted(large, key=lambda t: need[t.id])
    while todo:
        first = need[todo[0].id]
        members = [t for t in todo if need[t.id] <= first + spread]
        hi = min(cfg.one_pool_max, -(-max(need[t.id] for t in members) // step) * step)
        win = (max(cfg.amp_min, hi - spread), max(cfg.amp_opt, hi - opt_off), hi)
        for t in members:
            windows[t.id] = win
        classes.append((win, members))
        todo = todo[len(members):]

    def names(ts):
        return "; ".join(f"{t.id}: " + ", ".join(x.label() for x in t.sites[:3])
                         + (f" +{len(t.sites) - 3}" if len(t.sites) > 3 else "") for t in ts[:6])

    parts = [f"{len(targets) - len(large)} amplicon(s) of {cfg.amp_min}-{cfg.amp_max} bp"]
    parts += [f"{len(m)} of {w[0]}-{w[2]} bp ({names(m)})" for w, m in classes]
    note = (f"one pool, {len(classes) + 1} size classes: " + ", ".join(parts)
            + f". Targets were merged so that no two amplicons in the tube overlap or give an outer product shorter "
              f"than {cfg.outer_min} bp. Longer amplicons usually amplify less well in a mixed pool: compare read "
              f"depth after a first run and raise those primers' concentration if needed.")
    return targets, windows, note
