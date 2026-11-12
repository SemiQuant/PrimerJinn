import itertools
import random

import numpy as np

from primerJinn.multiplex.config import Config
from primerJinn.multiplex.optimize import Problem, conflict_graph, dsatur
from primerJinn.multiplex.targets import Target


def _random_problem(n_t, n_c, rng, density=1.0):
    U = [np.array([rng.uniform(0, 5) for _ in range(n_c)]) for _ in range(n_t)]
    B = {}
    for t, u in itertools.combinations(range(n_t), 2):
        if rng.random() < density:
            m = np.array([[rng.choice([0, 0, 0, rng.uniform(0, 3), 1000]) for _ in range(n_c)] for _ in range(n_c)])
            B[(t, u)] = m
    return Problem(U, B)


def test_exhaustive_matches_bruteforce():
    rng = random.Random(3)
    p = _random_problem(4, 5, rng)
    x, e, how = p.solve()
    assert how == "exhaustive"
    brute = min(p.energy(c) for c in itertools.product(range(5), repeat=4))
    assert abs(e - brute) < 1e-9


def test_annealing_finds_optimum_on_small_hard_case():
    rng = random.Random(7)
    p = _random_problem(7, 8, rng)  # 8^7 = 2.1e6 > exhaustive limit
    x, e, how = p.solve(restarts=8, seed=2)
    assert how == "annealing"
    brute = min(p.energy(c) for c in itertools.product(range(8), repeat=7))
    assert e <= brute + 1e-9
    assert abs(p.energy(x) - e) < 1e-9


def test_conflict_graph_and_colouring():
    cfg = Config(outer_min=1500).resolve(0.65)
    ts = [Target("a", "c1", 1000, 1100), Target("b", "c1", 1101, 1300),   # adjacent tiles: must overlap
          Target("c", "c1", 2000, 2000),                                   # outer product with b < 1500
          Target("d", "c1", 50000, 50000), Target("e", "c2", 1050, 1050)]
    g = conflict_graph(ts, cfg)
    assert 1 in g[0] and 2 in g[1]
    assert not g[3] and not g[4]
    col = dsatur(g)
    for v, nb in g.items():
        for u in nb:
            assert col[v] != col[u]
    # a, b and c are pairwise closer than outer_min (1500 bp): three pools
    assert 2 in g[0]
    assert max(col.values()) + 1 == 3
    ts[2] = Target("c", "c1", 3500, 3500)   # now only b-c... no: 3500 - 1101 > 1500, so a-c and b-c drop
    g = conflict_graph(ts, cfg)
    assert g[2] == set() and max(dsatur(g).values()) + 1 == 2


def test_tiles_of_one_range_alternate_in_two_pools():
    from primerJinn.multiplex.targets import Site, group_sites, span_limits
    cfg = Config().resolve(0.38)
    comfort, hard, gap = span_limits(300, 375, primer_len=22, margin=3, flex=20)
    tiles = group_sites([Site("c", 21563, 25384, "S")], comfort, hard, gap)
    assert len(tiles) >= 10
    g = conflict_graph(tiles, cfg)
    for i, t in enumerate(tiles):   # only neighbours conflict
        assert g[i] <= {i - 1, i + 1}
    assert max(dsatur(g).values()) + 1 == 2
    # an unrelated target 600 bp from the end of the range still conflicts with the last tile
    tiles.append(type(tiles[0])("X", "c", 25984, 25984, [Site("c", 25984, 25984, "x")],
                                frozenset({("c", 25984, 25984, "x")})))
    g = conflict_graph(tiles, cfg)
    assert len(tiles) - 2 in g[len(tiles) - 1]
