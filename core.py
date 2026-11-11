"""Moser–Tardos for Ramsey colourings."""

from itertools import combinations
from math import comb
import time

import numpy as np
from numba import njit

RULE_RANDOM = 0      # pick a uniformly random violated event
RULE_LAST = 1        # most recently violated (LIFO)
RULE_MIN_INDEX = 2   # smallest index (deterministic, O(#violated) per step)
RULES = {"random": RULE_RANDOM, "last": RULE_LAST, "min_index": RULE_MIN_INDEX}


def edge_id(u, v, n):
    """Index of edge {u, v} (u < v) in lexicographic order of pairs."""
    if u > v:
        u, v = v, u
    return u * n - u * (u + 1) // 2 + (v - u - 1)


class RamseyInstance:
    """Event structure for "no monochromatic K_k in a 2-colouring of K_n"."""

    def __init__(self, n, k):
        if k < 2 or n < k:
            raise ValueError("need 2 <= k <= n")
        self.n, self.k = n, k
        self.num_edges = comb(n, 2)
        self.edges_per_set = comb(k, 2)
        self.num_sets = comb(n, k)

        ksets = np.array(list(combinations(range(n), k)), dtype=np.int32)
        pairs = list(combinations(range(k), 2))
        set_edges = np.empty((len(ksets), len(pairs)), dtype=np.int32)
        for j, (a, b) in enumerate(pairs):
            u, v = ksets[:, a].astype(np.int64), ksets[:, b].astype(np.int64)
            set_edges[:, j] = u * n - u * (u + 1) // 2 + (v - u - 1)

        # CSR incidence edge -> sets
        counts = np.bincount(set_edges.ravel(), minlength=self.num_edges)
        inc_ptr = np.zeros(self.num_edges + 1, dtype=np.int64)
        np.cumsum(counts, out=inc_ptr[1:])
        order = np.argsort(set_edges.ravel(), kind="stable")
        inc_idx = (order // len(pairs)).astype(np.int32)

        self.ksets, self.set_edges = ksets, set_edges
        self.inc_ptr, self.inc_idx = inc_ptr, inc_idx

    def count_monochromatic(self, colours):
        """Number of monochromatic k-subsets under a 0/1 edge colouring."""
        red = colours[self.set_edges].sum(axis=1)
        return int(np.count_nonzero((red == 0) | (red == self.edges_per_set)))


@njit(cache=True)
def _mt_run(set_edges, inc_ptr, inc_idx, num_edges, rule, budget, seed):
    np.random.seed(seed)
    M, E = set_edges.shape
    colours = np.random.randint(0, 2, num_edges).astype(np.int8)

    red = np.zeros(M, dtype=np.int32)
    for s in range(M):
        c = 0
        for j in range(E):
            c += colours[set_edges[s, j]]
        red[s] = c

    vlist = np.empty(M, dtype=np.int32)
    vpos = np.full(M, -1, dtype=np.int32)
    nv = 0
    for s in range(M):
        if red[s] == 0 or red[s] == E:
            vpos[s] = nv
            vlist[nv] = s
            nv += 1

    best = nv
    step_of_best = 0
    steps = 0
    while nv > 0 and steps < budget:
        # --- choose a violated event
        if rule == 0:
            s = vlist[np.random.randint(0, nv)]
        elif rule == 1:
            s = vlist[nv - 1]
        else:
            s = vlist[0]
            for t in range(1, nv):
                if vlist[t] < s:
                    s = vlist[t]
        # --- resample its edges
        for j in range(E):
            e = set_edges[s, j]
            new = np.int8(np.random.randint(0, 2))
            if new == colours[e]:
                continue
            delta = 1 if new == 1 else -1
            colours[e] = new
            for q in range(inc_ptr[e], inc_ptr[e + 1]):
                t = inc_idx[q]
                was_bad = red[t] == 0 or red[t] == E
                red[t] += delta
                is_bad = red[t] == 0 or red[t] == E
                if was_bad and not is_bad:
                    p = vpos[t]
                    last = vlist[nv - 1]
                    vlist[p] = last
                    vpos[last] = p
                    vpos[t] = -1
                    nv -= 1
                elif is_bad and not was_bad:
                    vpos[t] = nv
                    vlist[nv] = t
                    nv += 1
        steps += 1
        if nv < best:
            best = nv
            step_of_best = steps
    return nv == 0, steps, best, step_of_best, colours, nv


def solve(n, k, budget, seed, rule="random", instance=None):
    """Run Moser–Tardos once.

    Returns a dict matching the group's interface:
      success, steps, min_bad_seen, step_of_min, time_ms, final_bad, colours
    `instance` may be a prebuilt RamseyInstance for (n, k) to avoid rebuilding.
    """
    inst = instance if instance is not None else RamseyInstance(n, k)
    t0 = time.perf_counter()
    ok, steps, best, step_best, colours, final_bad = _mt_run(
        inst.set_edges, inst.inc_ptr, inst.inc_idx, inst.num_edges,
        RULES[rule], int(budget), int(seed))
    return {
        "success": bool(ok),
        "steps": int(steps),
        "min_bad_seen": int(best),
        "step_of_min": int(step_best),
        "time_ms": (time.perf_counter() - t0) * 1000.0,
        "final_bad": int(final_bad),
        "colours": colours,
    }


def n_lll(k):
    """Largest n for which the symmetric LLL certifies R(k,k) > n."""
    from math import e
    p = 2.0 ** (1 - comb(k, 2))
    n = k
    while True:
        m = n + 1
        d = sum(comb(k, i) * comb(m - k, k - i) for i in range(2, k))
        if e * p * (d + 1) <= 1:
            n = m
        else:
            return n


# ---------------------------------------------------------------------------
# Shared helpers for local-search baselines (see baselines.py).
#
# A search "state" is five arrays:
#   colours (int8, one per edge), red (int32, red-edge count per k-set),
#   vlist / vpos (indexed set of bad k-sets), nv (int64 array of length 1:
#   current number of bad k-sets).
# All helpers are Numba-compiled, so call them from @njit code.
# ---------------------------------------------------------------------------

@njit(cache=True)
def init_state(set_edges, colours):
    """Build red counts and the bad-set index for a given colouring."""
    M, E = set_edges.shape
    red = np.zeros(M, dtype=np.int32)
    for s in range(M):
        c = 0
        for j in range(E):
            c += colours[set_edges[s, j]]
        red[s] = c
    vlist = np.empty(M, dtype=np.int32)
    vpos = np.full(M, -1, dtype=np.int32)
    nv = np.zeros(1, dtype=np.int64)
    for s in range(M):
        if red[s] == 0 or red[s] == E:
            vpos[s] = nv[0]
            vlist[nv[0]] = s
            nv[0] += 1
    return red, vlist, vpos, nv


@njit(cache=True)
def flip_delta(e, colours, red, inc_ptr, inc_idx, E):
    """Change in the number of bad k-sets if edge e were flipped (no mutation)."""
    step = -1 if colours[e] == 1 else 1
    delta = 0
    for q in range(inc_ptr[e], inc_ptr[e + 1]):
        t = inc_idx[q]
        was_bad = red[t] == 0 or red[t] == E
        c = red[t] + step
        is_bad = c == 0 or c == E
        delta += (1 if is_bad else 0) - (1 if was_bad else 0)
    return delta


@njit(cache=True)
def flip_edge(e, colours, red, vlist, vpos, nv, inc_ptr, inc_idx, E):
    """Flip edge e, updating red counts and the bad-set index in place."""
    step = -1 if colours[e] == 1 else 1
    colours[e] = 1 - colours[e]
    for q in range(inc_ptr[e], inc_ptr[e + 1]):
        t = inc_idx[q]
        was_bad = red[t] == 0 or red[t] == E
        red[t] += step
        is_bad = red[t] == 0 or red[t] == E
        if was_bad and not is_bad:
            p = vpos[t]
            last = vlist[nv[0] - 1]
            vlist[p] = last
            vpos[last] = p
            vpos[t] = -1
            nv[0] -= 1
        elif is_bad and not was_bad:
            vpos[t] = nv[0]
            vlist[nv[0]] = t
            nv[0] += 1
