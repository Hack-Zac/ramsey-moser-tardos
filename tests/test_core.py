"""Correctness checks for the Moser–Tardos Ramsey core. Run: python3 -m pytest -q"""
from itertools import combinations
import numpy as np
from core import RamseyInstance, solve, edge_id, n_lll


def naive_mono(colours, n, k):
    cnt = 0
    for S in combinations(range(n), k):
        cols = {colours[edge_id(u, v, n)] for u, v in combinations(S, 2)}
        cnt += len(cols) == 1
    return cnt


def paley_colouring(q):
    """Red iff u - v is a nonzero square mod q (q prime, q = 1 mod 4)."""
    squares = {(x * x) % q for x in range(1, q)}
    col = np.zeros(q * (q - 1) // 2, dtype=np.int8)
    for u, v in combinations(range(q), 2):
        col[edge_id(u, v, q)] = 1 if (v - u) % q in squares else 0
    return col


def test_edge_ids_are_a_bijection():
    n = 9
    ids = sorted(edge_id(u, v, n) for u, v in combinations(range(n), 2))
    assert ids == list(range(n * (n - 1) // 2))


def test_paley17_has_no_monochromatic_K4():
    inst = RamseyInstance(17, 4)
    col = paley_colouring(17)
    assert inst.count_monochromatic(col) == 0
    assert naive_mono(col, 17, 4) == 0


def test_paley17_has_monochromatic_K3():
    # sanity: R(3,3) = 6 <= 17, so triangles must exist
    assert RamseyInstance(17, 3).count_monochromatic(paley_colouring(17)) > 0


def test_vectorised_count_matches_naive():
    rng = np.random.default_rng(0)
    for n, k in [(8, 3), (9, 4), (10, 4)]:
        inst = RamseyInstance(n, k)
        for _ in range(5):
            col = rng.integers(0, 2, inst.num_edges).astype(np.int8)
            assert inst.count_monochromatic(col) == naive_mono(col, n, k)


def test_incremental_counts_match_recount():
    # stop early so the run ends with violations; the incremental count must agree
    for rule in ["random", "last", "min_index"]:
        inst = RamseyInstance(16, 4)
        for seed in range(5):
            r = solve(16, 4, budget=3, seed=seed, rule=rule, instance=inst)
            assert r["final_bad"] == inst.count_monochromatic(r["colours"])


def test_success_means_valid_colouring():
    inst = RamseyInstance(12, 4)
    for rule in ["random", "last", "min_index"]:
        r = solve(12, 4, budget=10**6, seed=1, rule=rule, instance=inst)
        assert r["success"]
        assert inst.count_monochromatic(r["colours"]) == 0


def test_impossible_instance_never_succeeds():
    # R(3,3) = 6: every 2-colouring of K_6 has a monochromatic triangle
    r = solve(6, 3, budget=20000, seed=0)
    assert not r["success"] and r["steps"] == 20000


def test_n_lll_values():
    assert [n_lll(k) for k in (4, 5, 6)] == [5, 9, 15]
