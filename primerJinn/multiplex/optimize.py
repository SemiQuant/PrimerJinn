"""Pool assignment (graph colouring) and choice of one primer pair per target.

Energy of a choice x (x[t] = candidate index for target t):

    E(x) = sum_t U[t][x[t]] + sum_{t<u, same pool} B[t,u][x[t], x[u]]

U holds the per-pair costs, B the costs between pairs of two targets that share
a tube (primer dimers, unintended products). Minimised by exhaustive search for
small problems, otherwise randomised greedy starts, simulated annealing and a
final 1-opt descent.
"""
from __future__ import annotations

import itertools
import math
import random
from typing import Optional, Dict, List, Sequence, Tuple

import numpy as np


# ----------------------------------------------------------------------------
# pools
# ----------------------------------------------------------------------------
def conflict_graph(targets, cfg) -> Dict[int, set]:
    """Edge when two targets on one sequence are so close that their amplicons
    must overlap, or that the outer product (first forward primer to second
    reverse primer) would be shorter than cfg.outer_min."""
    edges = {i: set() for i in range(len(targets))}
    flank = cfg.margin + cfg.primer_min
    order = sorted(range(len(targets)), key=lambda i: (targets[i].seqid, targets[i].start))
    for a_pos, i in enumerate(order):
        ti = targets[i]
        for j in order[a_pos + 1:]:
            tj = targets[j]
            if tj.seqid != ti.seqid:
                break
            outer_min_len = (tj.end + flank) - (ti.start - flank) + 1
            gap = tj.start - ti.end - 1
            must_overlap = gap < 2 * flank
            # later targets start even further right; stop once they can neither overlap ti (a long
            # target, such as a junction amplicon spanning a deletion, overlaps many that start inside
            # it) nor give a short outer product
            if not must_overlap and outer_min_len >= cfg.outer_min + 2 * cfg.amp_ceiling():
                break
            # tiles of one requested range: only neighbours that must overlap conflict
            # (ARTIC-style alternating pools); their outer products are tolerated
            same_range = bool(ti.origins & tj.origins)
            if must_overlap or (not same_range and outer_min_len < cfg.outer_min):
                edges[i].add(j)
                edges[j].add(i)
    return edges


def dsatur(edges: Dict[int, set], preset: Optional[Dict[int, int]] = None) -> Dict[int, int]:
    """DSATUR greedy colouring. Returns node -> colour (0-based). `preset` seeds
    colours (for example pool numbers of a previous design); a preset colour is
    used unless a neighbour already holds it."""
    colour: Dict[int, int] = {}
    nodes = list(edges)
    for v, c in sorted((preset or {}).items()):
        if v in edges and c >= 0 and all(colour.get(u) != c for u in edges[v]):
            colour[v] = c
    while len(colour) < len(nodes):
        best, key = None, None
        for v in nodes:
            if v in colour:
                continue
            sat = len({colour[u] for u in edges[v] if u in colour})
            k = (sat, len(edges[v]), -v)
            if key is None or k > key:
                best, key = v, k
        used = {colour[u] for u in edges[best] if u in colour}
        c = 0
        while c in used:
            c += 1
        colour[best] = c
    return colour


# ----------------------------------------------------------------------------
# optimisation
# ----------------------------------------------------------------------------
class Problem:
    def __init__(self, unary: Sequence[np.ndarray], pairwise: Dict[Tuple[int, int], np.ndarray]):
        self.U = [np.asarray(u, dtype=float) for u in unary]
        self.n = len(self.U)
        self.nbr: List[List[Tuple[int, np.ndarray]]] = [[] for _ in range(self.n)]
        for (t, u), m in pairwise.items():
            m = np.asarray(m, dtype=float)
            assert m.shape == (len(self.U[t]), len(self.U[u])), (t, u, m.shape)
            self.nbr[t].append((u, m))
            self.nbr[u].append((t, m.T))

    def energy(self, x: Sequence[int]) -> float:
        e = sum(self.U[t][x[t]] for t in range(self.n))
        for t in range(self.n):
            for u, m in self.nbr[t]:
                if u > t:
                    e += m[x[t], x[u]]
        return float(e)

    def local(self, t: int, x: Sequence[int]) -> np.ndarray:
        """Energy contribution of every candidate of t given the others."""
        v = self.U[t].copy()
        for u, m in self.nbr[t]:
            v += m[:, x[u]]
        return v

    # --------------------------------------------------------------
    def greedy(self, rng: random.Random) -> List[int]:
        order = list(range(self.n))
        rng.shuffle(order)
        x = [-1] * self.n
        for t in order:
            v = self.U[t].copy()
            for u, m in self.nbr[t]:
                if x[u] >= 0:
                    v += m[:, x[u]]
            # randomised: pick among the 3 best with weights
            idx = np.argsort(v)[:3]
            w = np.exp(-(v[idx] - v[idx[0]]))
            x[t] = int(rng.choices(list(idx), weights=list(w))[0])
        return x

    def descend(self, x: List[int]) -> List[int]:
        improved = True
        while improved:
            improved = False
            for t in range(self.n):
                v = self.local(t, x)
                b = int(np.argmin(v))
                if v[b] < v[x[t]] - 1e-12:
                    x[t] = b
                    improved = True
        return x

    def anneal(self, x: List[int], rng: random.Random, steps: int, t0: float = 50.0, t1: float = 0.05) -> List[int]:
        e = self.energy(x)
        best, best_e = list(x), e
        movable = [t for t in range(self.n) if len(self.U[t]) > 1]
        if not movable:
            return x
        alpha = (t1 / t0) ** (1.0 / max(steps, 1))
        temp = t0
        for _ in range(steps):
            t = rng.choice(movable)
            c = rng.randrange(len(self.U[t]) - 1)
            if c >= x[t]:
                c += 1
            d = self.U[t][c] - self.U[t][x[t]]
            for u, m in self.nbr[t]:
                d += m[c, x[u]] - m[x[t], x[u]]
            if d <= 0 or rng.random() < math.exp(-d / temp):
                x[t] = c
                e += d
                if e < best_e - 1e-12:
                    best, best_e = list(x), e
            temp *= alpha
        return best

    def solve(self, restarts: int = 12, seed: int = 1, exhaustive_limit: int = 100_000) -> Tuple[List[int], float, str]:
        sizes = [len(u) for u in self.U]
        if any(s == 0 for s in sizes):
            raise ValueError("a target has no candidates")
        space = 1
        for s in sizes:
            space *= s
            if space > exhaustive_limit:
                break
        if space <= exhaustive_limit:
            best, best_e = None, math.inf
            for x in itertools.product(*[range(s) for s in sizes]):
                e = self.energy(x)
                if e < best_e:
                    best, best_e = list(x), e
            return best, best_e, "exhaustive"
        rng = random.Random(seed)
        starts = [self.descend(self.greedy(rng)) for _ in range(30)]
        starts.sort(key=self.energy)
        best, best_e = list(starts[0]), self.energy(starts[0])
        steps = 2000 * self.n
        for r in range(restarts):
            x = list(starts[r % len(starts)])
            x = self.descend(self.anneal(x, rng, steps))
            e = self.energy(x)
            if e < best_e - 1e-12:
                best, best_e = x, e
        return best, best_e, "annealing"
