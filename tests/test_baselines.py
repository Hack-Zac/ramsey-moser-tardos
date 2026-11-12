"""Checks for the shared local-search helpers and the baseline methods.
Tests for unfinished baselines are skipped until they are implemented."""
import numpy as np
import pytest

from ramsey_mt.core import RamseyInstance, init_state, flip_delta, flip_edge
from ramsey_mt import baselines

KEYS = {"success", "steps", "min_bad_seen", "step_of_min", "time_ms", "final_bad", "colours"}


def test_flip_delta_matches_recount():
    inst = RamseyInstance(10, 4)
    rng = np.random.default_rng(0)
    col = rng.integers(0, 2, inst.num_edges).astype(np.int8)
    red, vlist, vpos, nv = init_state(inst.set_edges, col)
    base = inst.count_monochromatic(col)
    for e in rng.integers(0, inst.num_edges, 30):
        d = flip_delta(e, col, red, inst.inc_ptr, inst.inc_idx, inst.edges_per_set)
        flipped = col.copy(); flipped[e] = 1 - flipped[e]
        assert d == inst.count_monochromatic(flipped) - base


def test_flip_edge_keeps_bad_set_index_consistent():
    inst = RamseyInstance(11, 4)
    rng = np.random.default_rng(1)
    col = rng.integers(0, 2, inst.num_edges).astype(np.int8)
    red, vlist, vpos, nv = init_state(inst.set_edges, col)
    for e in rng.integers(0, inst.num_edges, 200):
        flip_edge(e, col, red, vlist, vpos, nv, inst.inc_ptr, inst.inc_idx, inst.edges_per_set)
    assert nv[0] == inst.count_monochromatic(col)
    reds = col[inst.set_edges].sum(axis=1)
    bad = set(np.flatnonzero((reds == 0) | (reds == inst.edges_per_set)))
    assert set(vlist[: nv[0]]) == bad
    assert all(vpos[s] >= 0 for s in bad)


METHODS = [baselines.random_restart, baselines.greedy_local_search,
           baselines.simulated_annealing]


@pytest.mark.parametrize("method", METHODS, ids=lambda m: m.__name__)
def test_method_interface_and_validity(method):
    n, k = (6, 4) if method is baselines.random_restart else (10, 4)
    inst = RamseyInstance(n, k)
    try:
        r = method(n, k, budget=200000, seed=0, instance=inst)
    except NotImplementedError:
        pytest.skip(f"{method.__name__} not implemented yet")
    assert KEYS <= set(r)
    assert r["success"], f"{method.__name__} should solve K_{n}, k={k} easily"
    assert inst.count_monochromatic(r["colours"]) == 0
    assert r["min_bad_seen"] == 0 and r["final_bad"] == 0


@pytest.mark.parametrize("method", METHODS, ids=lambda m: m.__name__)
def test_method_reports_failure_when_impossible(method):
    # R(3,3) = 6, so every colouring of K_6 has a monochromatic triangle
    inst = RamseyInstance(6, 3)
    try:
        r = method(6, 3, budget=2000, seed=0, instance=inst)
    except NotImplementedError:
        pytest.skip(f"{method.__name__} not implemented yet")
    assert not r["success"] and r["steps"] == 2000
    assert r["min_bad_seen"] >= 2  # Goodman: K_6 always has at least 2 monochromatic triangles
