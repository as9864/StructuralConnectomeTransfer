"""Worm vs fly at matched size and density, against within-fly variation.

Both species are electron-microscopy connectomes at single-neuron resolution,
so the measurement-modality problem of the human comparison does not arise.
The question is whether the C. elegans nervous system looks like "one more fly
circuit" once node count and density are matched, judged against how much the
same fly circuit varies between replicates:

  replicates   FlyWire left vs right hemisphere (same animal), and hemibrain
               vs FlyWire (different animals, different reconstruction)
  circuits     mushroom body, lateral horn, antennal lobe, lateral complex,
               central complex

Every graph is sampled at n nodes and thresholded to density d, then scored
with the null-normalized statistics (normalized.py). Two samplers:
  bfs      breadth-first growth from a random seed neuron: keeps local wiring
           intact and works for the very sparse whole-circuit graphs
  random   random node-induced subgraph: unbiased in degree, but too sparse
           for large circuits, which are then excluded at that density

Usage:
    python -m sctransfer.run_worm_fly --n-samples 8
"""
import argparse
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

from sctransfer import flywire, hemibrain
from sctransfer import normalized as N
from sctransfer.celegans import build_weight_matrix as celegans_weight_matrix
from sctransfer.topology import subsample_subgraph

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
SIZES = [100, 200, 279]
DENSITIES = [0.02, 0.04]
FLY_MIN_SYNAPSES = 3  # same rule as the committed hemibrain circuits


def bfs_subgraph(W, n, seed, G=None):
    """n neurons grown breadth-first (neighbours in random order) from a
    random seed in the largest connected component."""
    rng = np.random.default_rng(seed)
    G = G if G is not None else nx.from_numpy_array(N.binarize(W))
    comp = max(nx.connected_components(G), key=len)
    if len(comp) < n:
        return None
    start = rng.choice(sorted(comp))
    order, seen, frontier = [start], {start}, [start]
    while frontier and len(order) < n:
        nxt = []
        for v in frontier:
            nbrs = [u for u in G.neighbors(v) if u not in seen]
            rng.shuffle(nbrs)
            for u in nbrs:
                seen.add(u)
                order.append(u)
                nxt.append(u)
                if len(order) == n:
                    break
            if len(order) == n:
                break
        frontier = nxt
    keep = np.array(order[:n])
    return W[np.ix_(keep, keep)]


def sample(W, n, sampler, seed, G=None):
    if W.shape[0] < n:
        return None
    if W.shape[0] == n:
        return W
    if sampler == "bfs":
        return bfs_subgraph(W, n, seed, G)
    return subsample_subgraph(W, frac_nodes=n / W.shape[0], seed=seed)


def sources():
    """{(species, circuit, replicate): W}"""
    out = {("worm", "whole_animal", "celegans"): celegans_weight_matrix()}
    for name in hemibrain.available_circuits():
        out[("fly", name, "hemibrain_R")] = hemibrain.load_circuit(name)
    conn = flywire.load_connections()
    for name in flywire.CIRCUITS:
        sides = ["C"] if name in flywire.UNPAIRED else ["L", "R"]
        for side in sides:
            out[("fly", name, f"flywire_{side}")] = flywire.circuit_graph(conn, name, side, FLY_MIN_SYNAPSES)
    return out


def summarize(df):
    """Per (sampler, n, d, stat): worm value, fly circuit means, and three gaps."""
    rows = []
    for (sampler, n, d), g in df.groupby(["sampler", "n", "density"]):
        for stat in N.STATS:
            rep = g.groupby(["species", "circuit", "replicate"])[stat].mean()
            if ("worm", "whole_animal", "celegans") not in rep.index:
                continue
            worm = rep[("worm", "whole_animal", "celegans")]
            fly = rep.loc["fly"]
            circ = fly.groupby(level="circuit").mean()
            rep_gaps = []
            for c, reps in fly.groupby(level="circuit"):
                v = reps.to_numpy()
                rep_gaps += [abs(a - b) for i, a in enumerate(v) for b in v[i + 1:]]
            lr = [abs(fly[(c, "flywire_L")] - fly[(c, "flywire_R")]) for c in circ.index
                  if (c, "flywire_L") in fly.index and (c, "flywire_R") in fly.index]
            rows.append({
                "sampler": sampler, "n": n, "density": d, "stat": stat,
                "worm": worm, "fly_min": circ.min(), "fly_max": circ.max(), "n_circuits": len(circ),
                "gap_worm_vs_nearest_circuit": float(np.min(np.abs(circ - worm))),
                "gap_worm_vs_fly_mean": float(abs(circ.mean() - worm)),
                "gap_between_circuits_sd": float(circ.std(ddof=1)) if len(circ) > 1 else np.nan,
                "gap_replicates_median": float(np.median(rep_gaps)) if rep_gaps else np.nan,
                "gap_flywire_LR_median": float(np.median(lr)) if lr else np.nan,
            })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-samples", type=int, default=8)
    parser.add_argument("--n-nulls", type=int, default=N.N_NULLS_DEFAULT)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    RESULTS_DIR.mkdir(exist_ok=True)

    rows = []
    for (species, circuit, replicate), W in sources().items():
        print(f"[worm_fly] {species} {circuit} {replicate}: {W.shape[0]} neurons, "
              f"density {N.density_of(N.binarize(W)):.4f}", flush=True)
        G = nx.from_numpy_array(N.binarize(W))
        for sampler in ["bfs", "random"]:
            for n in SIZES:
                n_draws = 1 if W.shape[0] == n else args.n_samples
                for s in range(n_draws):
                    sub = sample(W, n, sampler, seed=args.seed + 1000 * s + n, G=G)
                    if sub is None:
                        continue
                    native = N.density_of(N.binarize(sub))
                    for d in DENSITIES:
                        if native < 0.95 * d:
                            continue
                        sig = N.normalized_signature(sub, d, args.n_nulls, seed=args.seed + s)
                        rows.append({"species": species, "circuit": circuit, "replicate": replicate,
                                     "sampler": sampler, "n": n, "density": d, "sample": s,
                                     "native_density": native, **{k: sig[k] for k in N.STATS}})

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS_DIR / "worm_fly_samples.csv", index=False)
    summary = summarize(df)
    summary.to_csv(RESULTS_DIR / "worm_fly_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print("\n=== Mean statistic per graph (bfs sampler) ===")
    print(df[df.sampler == "bfs"].groupby(["n", "density", "species", "circuit", "replicate"])[N.STATS]
          .mean().round(2).to_string())
    print("\n=== Worm vs fly gaps (in statistic units) ===")
    print(summary.round(2).to_string(index=False))
    print("\nSaved results/worm_fly_samples.csv, results/worm_fly_summary.csv")


if __name__ == "__main__":
    main()
