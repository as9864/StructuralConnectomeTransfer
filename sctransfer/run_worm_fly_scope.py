"""Scope control for the worm vs fly comparison (run_worm_fly.py).

C. elegans is a whole nervous system (sensory, inter and motor neurons in one
graph), while each fly graph in run_worm_fly.py comes from one neuropil. Does
the remaining worm-fly gap in clustering and hub structure come from that
difference in scope? Sample the FlyWire whole central brain (optic lobes
excluded) with k-seed breadth-first growth: k seeds in k different home
neuropils each grow n/k neurons, so larger k spans more of the brain. The worm
is sampled the same way for n below its size. The larval fly brain (larva.py)
is a whole nervous system like the worm, so it gives a scope match from the
data itself; it is sampled with one seed and with random node subsets.

Reports, per (n, density, k), the null-normalized statistics and how many
distinct neuropils the sample spans.

Usage:
    python -m sctransfer.run_worm_fly_scope --n-samples 8
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.sparse.csgraph  # noqa: F401  (registers sp.csgraph)

from sctransfer import flywire, larva
from sctransfer import normalized as N
from sctransfer.celegans import build_weight_matrix as celegans_weight_matrix

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
SIZES = [200, 279]
DENSITIES = [0.02, 0.04]
SEEDS_PER_SAMPLE = [1, 3, 9]
FLY_MIN_SYNAPSES = 3


def flywire_brain():
    """Sparse directed W over central-brain neurons, plus each neuron's home neuropil."""
    conn = flywire.load_connections()
    nps = flywire.hemisphere_neuropils(conn, "both")
    conn = conn[conn["neuropil"].isin(nps)]
    counts = flywire.neuron_neuropil_counts(conn)
    home = counts.loc[counts.groupby("body")["w"].idxmax()].set_index("body")["neuropil"]
    pairs = conn.groupby(["pre", "post"], observed=True)["w"].sum().reset_index()
    pairs = pairs[pairs["w"] >= FLY_MIN_SYNAPSES]
    bodies = np.unique(np.concatenate([pairs["pre"].to_numpy(), pairs["post"].to_numpy()]))
    i = np.searchsorted(bodies, pairs["pre"].to_numpy())
    j = np.searchsorted(bodies, pairs["post"].to_numpy())
    W = sp.csr_matrix((pairs["w"].to_numpy().astype(float), (i, j)), shape=(len(bodies), len(bodies)))
    labels = home.reindex(bodies).fillna("none").to_numpy()
    return W, labels


def grow(adj, n, seeds, rng):
    """k-seed breadth-first growth on an undirected sparse adjacency: seeds
    take turns adding one unseen neighbour of their own frontier."""
    k = len(seeds)
    quota = [n // k + (1 if s < n % k else 0) for s in range(k)]
    chosen, seen = [], set()
    for s in seeds:
        seen.add(int(s))
    frontiers = [[int(s)] for s in seeds]
    queues = [[] for _ in seeds]
    taken = [[int(s)] for s in seeds]
    while True:
        progressed = False
        for s in range(k):
            if len(taken[s]) >= quota[s]:
                continue
            while True:  # take the next unseen neighbour, refilling from the frontier as needed
                if not queues[s]:
                    if not frontiers[s]:
                        break
                    v = frontiers[s].pop(0)
                    nbrs = adj.indices[adj.indptr[v]:adj.indptr[v + 1]].copy()
                    rng.shuffle(nbrs)
                    queues[s].extend(int(u) for u in nbrs)
                    continue
                u = queues[s].pop(0)
                if u not in seen:
                    seen.add(u)
                    taken[s].append(u)
                    frontiers[s].append(u)
                    progressed = True
                    break
        if not progressed:
            break
    for t in taken:
        chosen += t
    return np.array(chosen[:n]) if len(chosen) >= n else None


def pick_seeds(labels, k, rng, in_giant):
    """k seeds in k different home neuropils (or anywhere if labels is None)."""
    pool = np.flatnonzero(in_giant)
    if labels is None:
        return rng.choice(pool, size=k, replace=False)
    regions = rng.choice(np.unique(labels[pool]), size=k, replace=False)
    return np.array([rng.choice(pool[labels[pool] == r]) for r in regions])


def giant_component(adj):
    _, comp = sp.csgraph.connected_components(adj, directed=False)
    return comp == np.bincount(comp).argmax()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-samples", type=int, default=8)
    parser.add_argument("--n-nulls", type=int, default=N.N_NULLS_DEFAULT)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    RESULTS_DIR.mkdir(exist_ok=True)
    rng = np.random.default_rng(args.seed)

    Wf, labels = flywire_brain()
    print(f"[scope] FlyWire central brain: {Wf.shape[0]:,} neurons, {Wf.nnz:,} connections (>= {FLY_MIN_SYNAPSES} synapses), "
          f"{len(np.unique(labels))} home neuropils", flush=True)
    Ww = sp.csr_matrix(celegans_weight_matrix())
    Wl, _ = larva.load()
    print(f"[scope] larval fly brain: {Wl.shape[0]:,} neurons, {Wl.nnz:,} connections", flush=True)
    graphs = {"fly_brain": (Wf, labels), "larva": (Wl, None), "worm": (Ww, None)}

    rows = []
    for name, (W, lab) in graphs.items():
        adj = ((W + W.T) > 0).astype(np.int8).tocsr()
        giant = giant_component(adj)
        for n in SIZES:
            for k in SEEDS_PER_SAMPLE + ["random"]:
                whole = name == "worm" and n >= W.shape[0]
                if whole and k != 1:
                    continue
                if name == "larva" and k not in (1, "random"):
                    continue  # larva: one seed (a whole nervous system needs no stitching) or random
                for s in range(1 if whole else args.n_samples):
                    if whole:
                        keep = np.arange(W.shape[0])
                    elif k == "random":
                        keep = rng.choice(W.shape[0], size=n, replace=False)
                    else:
                        keep = grow(adj, n, pick_seeds(lab, k, rng, giant), rng)
                        if keep is None:
                            continue
                    sub = W[keep][:, keep].toarray()
                    native = N.density_of(N.binarize(sub))
                    span = len(np.unique(lab[keep])) if lab is not None else np.nan
                    for d in DENSITIES:
                        if native < 0.95 * d:
                            continue
                        sig = N.normalized_signature(sub, d, args.n_nulls, seed=args.seed + s)
                        rows.append({"graph": name, "n": n, "seeds": "whole" if whole else k, "density": d,
                                     "sample": s, "native_density": native, "neuropils_spanned": span,
                                     **{key: sig[key] for key in N.STATS}})
            print(f"[scope] {name} n={n} done", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS_DIR / "worm_fly_scope.csv", index=False)
    pd.set_option("display.width", 220)
    summary = df.groupby(["n", "density", "graph", "seeds"])[
        ["native_density", "neuropils_spanned"] + N.STATS].agg(["mean", "std"])
    print("\n=== Statistics by sampling scope (mean, std over samples) ===")
    print(summary.round(2).to_string())
    print("\nSaved results/worm_fly_scope.csv")


if __name__ == "__main__":
    main()
