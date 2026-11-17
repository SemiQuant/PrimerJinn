"""A static amplicon map for the report: `<out>_map.png`.

Top: the genome overview, one strip per sequence, each amplicon a marker at its
position coloured by pool and labelled with its id (labels alternate above and below
so neighbours stay readable). Below: one lane per amplicon, all drawn to the same
bp scale from the forward primer's 5' end, showing the primers (arrows), the insert
the reads cover (coloured bar), the requested positions (red marks on the insert) and,
for a junction amplicon, the deletion it spans (hatched gap). A legend at the bottom.
"""
from __future__ import annotations

from typing import Dict, Optional

POOL_COLOURS = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#9085e9", "#008300", "#e66767"]
INK = "#2b2b2b"
MUTED = "#7a7a7a"
RED = "#e0393e"


def _colour(pool: int) -> str:
    return POOL_COLOURS[(pool - 1) % len(POOL_COLOURS)]


def _lengths(d) -> Dict[str, int]:
    try:
        import pysam
        with pysam.FastaFile(d.genome["fasta"], filepath_index=d.genome.get("fai")) as fa:
            return dict(zip(fa.references, fa.lengths))
    except Exception:
        return {}


def _short(text: str, n: int) -> str:
    return text if len(text) <= n else text[: n - 1] + "…"


def draw(d, path: str, dpi: int = 150) -> Optional[str]:
    """Write the map PNG and return its path (None when there is nothing to draw)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        from matplotlib.patches import FancyArrow, Patch, Rectangle
    except ImportError:
        return None
    tmap = {t.id: t for t in d.targets}
    ids = [tid for tid in d.chosen if tid in tmap]
    if not ids:
        return None
    ids.sort(key=lambda x: (tmap[x].seqid, d.chosen[x].start))
    seqids = list(dict.fromkeys(tmap[x].seqid for x in ids))
    pools_used = sorted({d.pools[x] for x in ids})
    lengths = _lengths(d)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": "#cccccc",
                         "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "text.color": INK})
    lane_h = 0.36
    over_h = 1.5
    fig_h = min(over_h * len(seqids) + lane_h * len(ids) + 1.6, 60)
    fig = plt.figure(figsize=(11.5, fig_h), facecolor="white")
    gs = fig.add_gridspec(len(seqids) + 1, 1, height_ratios=[over_h] * len(seqids) + [lane_h * len(ids) + 0.9],
                          hspace=0.3)

    # ------------------------------------------------------------- genome overview
    for i, sid in enumerate(seqids):
        ax = fig.add_subplot(gs[i, 0])
        here = [x for x in ids if tmap[x].seqid == sid]
        glen = lengths.get(sid) or max(d.chosen[x].end for x in here)
        ax.add_patch(Rectangle((0, -0.08), glen, 0.16, color="#e4e7ea", lw=0))
        # amplicons closer than ~1.2 % of the axis share one label (T03-T07) instead of overlapping ones
        groups, cur = [], []
        for x in here:
            mid = (d.chosen[x].start + d.chosen[x].end) / 2
            if cur and mid - cur[-1][1] > glen * 0.012:
                groups.append(cur)
                cur = []
            cur.append((x, mid))
        if cur:
            groups.append(cur)
        for k, grp in enumerate(groups):
            # four label levels (above, below, higher above, lower below) so neighbours stay readable
            lvl = [0.38, -0.38, 0.78, -0.78][k % 4]
            up = lvl > 0
            for x, mid in grp:
                ax.plot([mid, mid], [-0.08, 0.08], color=_colour(d.pools[x]), lw=2.2, solid_capstyle="butt",
                        zorder=3)
            gmid = sum(m for _, m in grp) / len(grp)
            ax.plot([gmid, gmid], [0.08 if up else -0.08, lvl], color="#b0b4b8", lw=0.6, zorder=2)
            pools_here = sorted({d.pools[x] for x, _ in grp})
            for j, pp in enumerate(pools_here):
                ax.plot(gmid, lvl, marker="o", ms=6.5 if len(pools_here) == 1 else 5, color=_colour(pp),
                        mec="white", mew=0.8, zorder=4 + j, alpha=1 if j == 0 else 0.85)
            label = grp[0][0] if len(grp) == 1 else f"{grp[0][0]}\u2013{grp[-1][0]}"
            ax.annotate(label, (gmid, lvl), xytext=(0, 8 if up else -8), textcoords="offset points",
                        ha="center", va="bottom" if up else "top", fontsize=6.8, color=INK)
        ax.set_xlim(-glen * 0.01, glen * 1.01)
        ax.set_ylim(-1.35, 1.35)
        ax.set_yticks([])
        for sp in ("left", "right", "top"):
            ax.spines[sp].set_visible(False)
        unit, div = ("Mb", 1e6) if glen >= 2e6 else ("kb", 1e3)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _, div=div: f"{v / div:g}"))
        ax.set_xlabel(f"{sid} ({unit})", fontsize=7.5, labelpad=2)
        ax.tick_params(axis="x", labelsize=7, length=3)
        ax.set_title(f"{len(here)} amplicon{'s' if len(here) != 1 else ''} along {sid}"
                     + (f", {glen / div:.2f} {unit}" if glen else ""), fontsize=8.5, loc="left", color=INK, pad=4)

    # --------------------------------------------------------- one lane per amplicon
    ax = fig.add_subplot(gs[len(seqids), 0])
    maxlen = max(d.chosen[x].end - d.chosen[x].start + 1 for x in ids)
    ax.set_xlim(-maxlen * 0.01, maxlen * 1.26)
    ax.set_ylim(len(ids) - 0.5, -0.5)
    ax.set_yticks(range(len(ids)))
    labels = []
    for k, x in enumerate(ids):
        pr, p, t = d.chosen[x], d.pools[x], tmap[x]
        a, b = pr.start, pr.end
        rel = lambda v: v - a
        col = _colour(p)
        if k % 2 == 0:
            ax.add_patch(Rectangle((-maxlen * 0.01, k - 0.5), maxlen * 1.27, 1, color="#f5f6f7", lw=0, zorder=0))
        # amplicon backbone and the insert the reads cover
        ax.plot([0, rel(b)], [k, k], color="#c9ccd0", lw=1, zorder=1)
        i0, i1 = pr.insert
        ax.add_patch(Rectangle((rel(i0), k - 0.16), i1 - i0 + 1, 0.32, color=col, alpha=0.45, lw=0, zorder=2))
        if pr.gap:
            g0, g1 = pr.gap
            ax.add_patch(Rectangle((rel(g0), k - 0.16), g1 - g0 + 1, 0.32, facecolor="white", edgecolor=RED,
                                   hatch="////", lw=0.6, zorder=3))
        # primers
        fl, rl = pr.fwd.end - pr.fwd.start + 1, pr.rev.end - pr.rev.start + 1
        hl = min(maxlen * 0.012, fl * 0.6)
        ax.add_patch(FancyArrow(rel(pr.fwd.start), k, fl, 0, width=0.34, head_width=0.62, head_length=hl,
                                length_includes_head=True, color=INK, lw=0, zorder=5))
        ax.add_patch(FancyArrow(rel(pr.rev.end) + 1, k, -rl, 0, width=0.34, head_width=0.62, head_length=hl,
                                length_includes_head=True, color=INK, lw=0, zorder=5))
        # requested positions on the insert
        for s in t.sites:
            s0, s1 = max(s.start, i0), min(s.end, i1)
            if s1 < s0:
                continue
            if s1 - s0 + 1 >= maxlen * 0.004:
                ax.add_patch(Rectangle((rel(s0), k - 0.11), s1 - s0 + 1, 0.22, color=RED, lw=0, zorder=4))
            else:
                ax.plot([rel(s0) + 0.5, rel(s0) + 0.5], [k - 0.24, k + 0.24], color=RED, lw=1.3, zorder=4)
        prod = fl + rl + (i1 - i0 + 1) - (pr.gap_length if pr.gap else 0)
        info = f"{prod} bp"
        if pr.gap:
            info += f", spans {pr.gap_length} bp deletion"
        info += f"  ·  {len(t.sites)} site{'s' if len(t.sites) != 1 else ''}: {_short(t.name, 46)}"
        ax.text(rel(b) + maxlen * 0.015, k, info, va="center", fontsize=6.3, color=MUTED, zorder=6, clip_on=False)
        labels.append(f"{x:<5} pool {p}   {t.seqid}:{a:,}-{b:,}")
    ax.set_yticklabels(labels, fontsize=6.6, family="DejaVu Sans Mono")
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=7, length=3)
    ax.set_xlabel("bp from the forward primer's 5' end (all amplicons to the same scale)", fontsize=7.5)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.set_title("Amplicons", fontsize=8.5, loc="left", color=INK, pad=4)
    handles = [Line2D([0], [0], marker=">", color=INK, ms=7, lw=0, label="forward / reverse primer"),
               Patch(color="#9bb9e0", label="insert the reads cover (colour = pool)"),
               Patch(color=RED, label="requested position"),
               Patch(facecolor="white", edgecolor=RED, hatch="////", label="deletion spanned by a junction amplicon")]
    handles += [Patch(color=_colour(p), label=f"pool {p}") for p in pools_used]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, -0.06 - 1.3 / max(len(ids), 1)), ncol=4,
              fontsize=6.6, frameon=False, handlelength=1.6, columnspacing=1.4)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path
