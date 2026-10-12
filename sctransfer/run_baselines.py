"""Within-species baselines for the worm vs larva vs adult fly comparison.

run_worm_fly_scope.py found the worm and the larval fly indistinguishable and
the adult fly different, but with one animal per group. Here every group gets
replicates, so a between-group difference can be judged against the
difference between replicates of the same kind:

  worm        eight C. elegans brains (Witvliet et al. 2021: L1 x4, L2, L3,
              adult x2) and the whole-animal hermaphrodite graph
  larva       larval fly brain, left and right hemisphere separately
  adult fly   FlyWire central brain, left and right hemisphere separately
              (neurons assigned by the neuropil holding most of their synapses)

Every graph is sampled the same way (one-seed breadth-first growth to n
neurons, then thresholded to density d) and scored with the null-normalized
statistics. Effect sizes are Cohen's d between the sample distributions of two
graphs.

Usage:
    python -m sctransfer.run_baselines --n-samples 30
"""
import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from sctransfer import flywire, larva, witvliet
from sctransfer import normalized as N
from sctransfer.celegans import build_weight_matrix as celegans_weight_matrix
from sctransfer.run_worm_fly_scope import flywire_brain, giant_component, grow

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
N_NODES = 150   # smallest graph (Witvliet L1 brain) has 187 neurons
DENSITIES = [0.02, 0.04]


def graphs(with_ids=False):
    """{name: (group, sparse W)}, or (group, W, node ids) with `with_ids`
    (worm: neuron names, larva: CATMAID skeleton ids, adult: FlyWire root ids)."""
    out = {}
    for k, stage in witvliet.DATASETS.items():
        Ww, names = witvliet.load(k)
        out[f"worm_w{k}_{stage}"] = ("worm", sp.csr_matrix(Ww), np.array(names))
    out["worm_whole_animal"] = ("worm", sp.csr_matrix(celegans_weight_matrix()), None)

    Wl, meta = larva.load()
    for side in ("L", "R"):
        keep = np.flatnonzero(meta["hemisphere"].to_numpy() == side)
        out[f"larva_{side}"] = ("larva", Wl[keep][:, keep], meta.index.to_numpy()[keep])

    Wf, labels, bodies = flywire_brain(return_ids=True)
    sides = np.array([flywire.neuropil_side(str(x)) for x in labels])
    for side in ("L", "R"):
        keep = np.flatnonzero(sides == side)
        out[f"adult_{side}"] = ("adult_fly", Wf[keep][:, keep].tocsr(), bodies[keep])
    return out if with_ids else {k: v[:2] for k, v in out.items()}


def cohens_d(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    sd = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    return (a.mean() - b.mean()) / sd if sd > 0 else np.nan


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-samples", type=int, default=30)
    parser.add_argument("--n-nulls", type=int, default=N.N_NULLS_DEFAULT)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    RESULTS_DIR.mkdir(exist_ok=True)
    rng = np.random.default_rng(args.seed)

    rows = []
    gs = graphs()
    for name, (group, W) in gs.items():
        adj = ((W + W.T) > 0).astype(np.int8).tocsr()
        giant = np.flatnonzero(giant_component(adj))
        for s in range(args.n_samples):
            keep = grow(adj, N_NODES, [rng.choice(giant)], rng)
            if keep is None:
                continue
            sub = W[keep][:, keep].toarray()
            native = N.density_of(N.binarize(sub))
            for d in DENSITIES:
                if native < 0.95 * d:
                    continue
                sig = N.normalized_signature(sub, d, args.n_nulls, seed=args.seed + s)
                rows.append({"graph": name, "group": group, "density": d, "sample": s,
                             "native_density": native, **{k: sig[k] for k in N.STATS}})
        print(f"[baselines] {name}: {W.shape[0]:,} neurons", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS_DIR / "baselines_samples.csv", index=False)

    pd.set_option("display.width", 220)
    print("\n=== Mean statistic per graph ===")
    print(df.groupby(["density", "group", "graph"])[N.STATS].mean().round(2).to_string())

    # pairwise effect sizes, labelled by comparison kind
    pair_rows = []
    names = list(gs)
    for d, sub in df.groupby("density"):
        for a, b in itertools.combinations(names, 2):
            ga, gb = gs[a][0], gs[b][0]
            kind = "within " + ga if ga == gb else "between " + "-".join(sorted([ga, gb]))
            xa, xb = sub[sub.graph == a], sub[sub.graph == b]
            if len(xa) < 3 or len(xb) < 3:
                continue
            pair_rows.append({"density": d, "a": a, "b": b, "kind": kind,
                              **{k: abs(cohens_d(xa[k], xb[k])) for k in N.STATS}})
    pairs = pd.DataFrame(pair_rows)
    pairs.to_csv(RESULTS_DIR / "baselines_pairs.csv", index=False)
    print("\n=== |Cohen's d| between graphs: median (max) by comparison kind ===")
    agg = pairs.groupby(["density", "kind"])[N.STATS].agg(["median", "max"]).round(2)
    print(agg.to_string())
    print("\nSaved results/baselines_samples.csv, results/baselines_pairs.csv")


if __name__ == "__main__":
    main()
