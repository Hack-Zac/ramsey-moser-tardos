"""Baseline searches sharing core.solve's interface and return keys."""
import time

import numpy as np
from numba import njit

from .core import RamseyInstance, init_state, flip_delta, flip_edge


def _result(ok, steps, best, step_best, colours, final_bad, t0):
    return {
        "success": bool(ok),
        "steps": int(steps),
        "min_bad_seen": int(best),
        "step_of_min": int(step_best),
        "time_ms": (time.perf_counter() - t0) * 1000.0,
        "final_bad": int(final_bad),
        "colours": colours,
    }


# One step = one fresh colouring.
@njit(cache=True)
def _restart_run(set_edges, num_edges, budget, seed):
    np.random.seed(seed)
    M, E = set_edges.shape
    best = M + 1
    step_of_best = 0
    colours = np.zeros(num_edges, dtype=np.int8)
    bad = M
    steps = 0
    while steps < budget:
        steps += 1
        for e in range(num_edges):
            colours[e] = np.random.randint(0, 2)
        bad = 0
        for s in range(M):
            c = 0
            for j in range(E):
                c += colours[set_edges[s, j]]
            if c == 0 or c == E:
                bad += 1
        if bad < best:
            best = bad
            step_of_best = steps
        if bad == 0:
            break
    return bad == 0, steps, best, step_of_best, colours, bad


def random_restart(n, k, budget, seed, instance=None):
    inst = instance if instance is not None else RamseyInstance(n, k)
    t0 = time.perf_counter()
    ok, steps, best, sb, colours, final = _restart_run(
        inst.set_edges, inst.num_edges, int(budget), int(seed))
    return _result(ok, steps, best, sb, colours, final, t0)


@njit(cache=True)
def _lowest_delta_edge(bad_set, set_edges, colours, red, inc_ptr, inc_idx):
    E = set_edges.shape[1]
    best_edge = set_edges[bad_set, 0]
    lowest_delta = flip_delta(best_edge, colours, red, inc_ptr, inc_idx, E)
    ties = 1
    for j in range(1, E):
        edge = set_edges[bad_set, j]
        delta = flip_delta(edge, colours, red, inc_ptr, inc_idx, E)
        if delta < lowest_delta:
            best_edge = edge
            lowest_delta = delta
            ties = 1
        elif delta == lowest_delta:
            ties += 1
            # Reservoir sampling: uniform tie-break.
            if np.random.randint(0, ties) == 0:
                best_edge = edge
    return best_edge


# One step = one flip.
@njit(cache=True)
def _greedy_run(set_edges, inc_ptr, inc_idx, num_edges, noise, budget, seed):
    np.random.seed(seed)
    E = set_edges.shape[1]
    colours = np.random.randint(0, 2, num_edges).astype(np.int8)
    red, vlist, vpos, nv = init_state(set_edges, colours)
    best, step_of_best, steps = nv[0], 0, 0
    while nv[0] > 0 and steps < budget:
        bad_set = vlist[np.random.randint(0, nv[0])]
        if np.random.random() < noise:
            edge = set_edges[bad_set, np.random.randint(0, E)]
        else:
            edge = _lowest_delta_edge(bad_set, set_edges, colours, red, inc_ptr, inc_idx)
        flip_edge(edge, colours, red, vlist, vpos, nv, inc_ptr, inc_idx, E)
        steps += 1
        if nv[0] < best:
            best, step_of_best = nv[0], steps
    return nv[0] == 0, steps, best, step_of_best, colours, nv[0]


def greedy_local_search(n, k, budget, seed, instance=None, noise=0.2):
    inst = instance if instance is not None else RamseyInstance(n, k)
    t0 = time.perf_counter()
    ok, steps, best, sb, colours, final = _greedy_run(
        inst.set_edges, inst.inc_ptr, inst.inc_idx, inst.num_edges,
        float(noise), int(budget), int(seed))
    return _result(ok, steps, best, sb, colours, final, t0)


# One step = one proposed flip.
@njit(cache=True)
def _anneal_run(set_edges, inc_ptr, inc_idx, num_edges, t0, cooling, t_min, budget, seed):
    np.random.seed(seed)
    E = set_edges.shape[1]
    colours = np.random.randint(0, 2, num_edges).astype(np.int8)
    red, vlist, vpos, nv = init_state(set_edges, colours)
    best, step_of_best, steps = nv[0], 0, 0
    temperature = t0
    while nv[0] > 0 and steps < budget:
        edge = np.random.randint(0, num_edges)
        delta = flip_delta(edge, colours, red, inc_ptr, inc_idx, E)
        if delta <= 0 or np.random.random() < np.exp(-delta / temperature):
            flip_edge(edge, colours, red, vlist, vpos, nv, inc_ptr, inc_idx, E)
        temperature = max(temperature * cooling, t_min)
        steps += 1
        if nv[0] < best:
            best, step_of_best = nv[0], steps
    return nv[0] == 0, steps, best, step_of_best, colours, nv[0]


def simulated_annealing(n, k, budget, seed, instance=None,
                        t0=2.0, cooling=0.99999, t_min=0.05):
    inst = instance if instance is not None else RamseyInstance(n, k)
    start = time.perf_counter()
    ok, steps, best, sb, colours, final = _anneal_run(
        inst.set_edges, inst.inc_ptr, inst.inc_idx, inst.num_edges,
        float(t0), float(cooling), float(t_min), int(budget), int(seed))
    return _result(ok, steps, best, sb, colours, final, start)
