"""--one_pool: targets that may not share a tube merge into one amplicon; only the
merged targets that need more room get a larger size class, the rest keep 225-375 bp."""
import pytest

from primerJinn.multiplex.config import Config
from primerJinn.multiplex.optimize import conflict_graph
from primerJinn.multiplex.targets import OnePoolError, Site, group_for_config, one_pool_layout


def _cfg(**kw):
    return Config(one_pool=True, threads=1, **kw).resolve(0.65)


def test_far_targets_keep_the_normal_window():
    targets, windows, note = one_pool_layout([Site("c", 10_000, 10_000, "a"), Site("c", 20_000, 20_000, "b")],
                                             _cfg())
    assert len(targets) == 2 and windows == {} and "no target needed merging" in note


def test_only_merged_targets_get_a_larger_class():
    cfg = _cfg()
    sites = [Site("c", 6735, 6735, "gyrB"), Site("c", 7570, 7570, "gyrA"), Site("c", 20_000, 20_000, "far")]
    assert any(conflict_graph(group_for_config(sites, cfg), cfg).values())   # two pools without --one_pool
    targets, windows, note = one_pool_layout(sites, cfg)
    assert len(targets) == 2
    big = targets[0]
    assert (big.start, big.end) == (6735, 7570) and list(windows) == [big.id]
    lo, opt, hi = windows[big.id]
    assert hi - lo == 150 and hi - opt == 75 and 875 <= hi <= 925
    assert cfg.amp(targets[1].id) == (225, 300, 375)                        # the far target stays short
    cfg.amp_by_target = windows
    assert not any(conflict_graph(targets, cfg).values())
    assert "gyrB" in note and "size classes" in note


def test_classes_follow_the_needs():
    """A 561 bp range and an 835 bp pair: two separate large classes, not one."""
    sites = [Site("c", 50_000, 50_560, "pncA"), Site("c", 6735, 6735, "gyrB"), Site("c", 7570, 7570, "gyrA"),
             Site("c", 90_000, 90_000, "snp")]
    targets, windows, _ = one_pool_layout(sites, _cfg())
    his = sorted({w[2] for w in windows.values()})
    assert len(his) == 2 and 600 <= his[0] <= 650 and his[1] >= 875
    assert sum(1 for t in targets if t.id not in windows) == 1               # the SNP keeps 225-375


def test_cap_too_small_names_the_targets():
    sites = [Site("c", 6735, 6735, "gyrB"), Site("c", 7570, 7570, "gyrA")]
    with pytest.raises(OnePoolError, match=r"gyrB.*gyrA.*about [89]\d\d bp"):
        one_pool_layout(sites, _cfg(one_pool_max=600))


def test_lower_outer_min_avoids_merging():
    sites = [Site("c", 6735, 6735, "gyrB"), Site("c", 7570, 7570, "gyrA")]
    targets, windows, _ = one_pool_layout(sites, _cfg(outer_min=700))
    assert len(targets) == 2 and windows == {}


def test_one_pool_max_below_amp_max_is_rejected():
    with pytest.raises(ValueError, match="one_pool_max"):
        Config(one_pool=True, one_pool_max=300).resolve(0.5)
