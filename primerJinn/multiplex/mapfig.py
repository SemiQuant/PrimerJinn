"""A static amplicon map for the report: `<out>_map.png`.

Top: the genome overview, one row per pool, each amplicon a mark at its position with
its id. Below: one lane per amplicon in its own coordinates showing the primers
(arrows), the insert the reads cover, the requested positions (ticks) and, for a
junction amplicon, the deletion it spans.
"""
from __future__ import annotations

from typing import List, Optional, Tuple

POOL_COLOURS = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#9085e9", "#008300", "#e66767"]


def _requested(target) -> List[Tuple[int, int]]:
    return [(s.start, s.end) for s in target.sites]


def draw(d, path: str, dpi: int = 130) -> Optional[str]:
    """Write the map PNG and return its path (None when there is nothing to draw)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import FancyArrow
    except ImportError:
        return None
    tmap = {t.id: t for t in d.targets}
    ids = [tid for tid in d.chosen if tid in tmap]
    if not ids:
        return None
    ids.sort(key=lambda x: (tmap[x].seqid, d.chosen[x].start))
    seqids = list(dict.fromkeys(tmap[x].seqid for x in ids))
    n_pools = max(d.pools[x] for x in ids)
    lane_h = 0.42
    over_h = 0.6 + 0.45 * n_pools
    fig_h = over_h * len(seqids) + lane_h * len(ids) + 1.2
    fig = plt.figure(figsize=(11, min(fig_h, 60)))
    gs = fig.add_gridspec(len(seqids) + 1, 1, height_ratios=[over_h] * len(seqids) + [lane_h * len(ids)],
                          hspace=0.9 if len(seqids) > 1 else 0.35)

    lengths = {}
    try:
        import pysam
        with pysam.FastaFile(d.genome["fasta"], filepath_index=d.genome.get("fai")) as fa:
            lengths = dict(zip(fa.references, fa.lengths))
    except Exception:
        pass

    # ---- overview per sequence
    for i, sid in enumerate(seqids):
        ax = fig.add_subplot(gs[i, 0])
        length = lengths.get(sid)
        here = [x for x in ids if tmap[x].seqid == sid]
        for x in here:
            pr, p = d.chosen[x], d.pools[x]
            mid = (pr.start + pr.end) / 2
            ax.plot(mid, p, "o", ms=9, color=POOL_COLOURS[(p - 1) % len(POOL_COLOURS)], zorder=3)
            ax.annotate(x, (mid, p), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=6.5,
                        rotation=60, rotation_mode="anchor")
        ax.set_yticks(range(1, n_pools + 1))
        ax.set_yticklabels([f"pool {p}" for p in range(1, n_pools + 1)], fontsize=8)
        ax.set_ylim(0.4, n_pools + 1.1)
        if length:
            ax.set_xlim(0, length)
        ax.set_xlabel(f"{sid} position (Mb)", fontsize=8)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v / 1e6:g}"))
        ax.tick_params(axis="x", labelsize=7)
        ax.grid(axis="x", alpha=0.25)
        ax.set_title(f"Amplicons along {sid}" if len(seqids) > 1 else "Genome overview", fontsize=9, loc="left")

    # ---- one lane per amplicon
    ax = fig.add_subplot(gs[len(seqids), 0])
    ax.set_xlim(-0.02, 1.08)
    ax.set_ylim(-0.6, len(ids) - 0.4)
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks(range(len(ids)))
    labels = []
    for k, x in enumerate(ids):
        pr, p, t = d.chosen[x], d.pools[x], tmap[x]
        a, b = pr.start, pr.end
        span = max(b - a, 1)
        sx = lambda v: (v - a) / span
        col = POOL_COLOURS[(p - 1) % len(POOL_COLOURS)]
        # amplicon bar, insert, deletion
        ax.plot([0, 1], [k, k], color="#bbbbbb", lw=1.2, zorder=1)
        i0, i1 = pr.insert
        ax.plot([sx(i0), sx(i1)], [k, k], color=col, lw=4, alpha=0.55, solid_capstyle="butt", zorder=2)
        if pr.gap:
            ax.plot([sx(pr.gap[0]), sx(pr.gap[1])], [k, k], color="white", lw=4, zorder=3)
            ax.plot([sx(pr.gap[0]), sx(pr.gap[1])], [k, k], color="#f05050", lw=1, ls="--", zorder=4)
        # primers as arrows
        fl = sx(pr.fwd.end) - sx(pr.fwd.start)
        rl = sx(pr.rev.end) - sx(pr.rev.start)
        ax.add_patch(FancyArrow(sx(pr.fwd.start), k, fl, 0, width=0.28, head_width=0.5, head_length=min(fl, 0.02),
                                length_includes_head=True, color="#222222", zorder=5))
        ax.add_patch(FancyArrow(sx(pr.rev.end), k, -rl, 0, width=0.28, head_width=0.5, head_length=min(rl, 0.02),
                                length_includes_head=True, color="#222222", zorder=5))
        # requested positions
        for s0, s1 in _requested(t):
            lo, hi = max(sx(s0), 0), min(sx(s1 + 1), 1)
            ax.plot([lo, hi], [k - 0.3, k - 0.3], color="#f05050", lw=2.2, solid_capstyle="butt", zorder=6)
        kind = f", junction {pr.gap_length} bp del" if pr.gap else ""
        ax.text(1.015, k, f"{len(pr.fwd.seq) + len(pr.rev.seq) + (i1 - i0 + 1)} bp{kind}", va="center", fontsize=6.5)
        labels.append(f"{x}  pool {p}  {t.seqid}:{a}-{b}")
    ax.set_yticklabels(labels, fontsize=6.5, family="monospace")
    for sp in ("top", "right", "bottom"):
        ax.spines[sp].set_visible(False)
    ax.set_title("Each amplicon: primers (arrows), insert the reads cover (coloured), requested positions "
                 "(red ticks), spanned deletion (dashed)", fontsize=8, loc="left")
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return path
