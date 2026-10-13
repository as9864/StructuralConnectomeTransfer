"""How much does the node definition alone change the normalized statistics?

Invertebrate connectomes have neurons as nodes; human SC has brain regions.
Instead of only listing this as a limitation, measure it inside each species
by coarse-graining (summing a weighted graph over groups of nodes):

  fly    neuron (circuit)  -> cell type (circuit)  -> brain region (whole hemibrain)
  worm   neuron (279)      -> neuron class (left/right partners merged)
  human  Schaefer 400 parcel -> Schaefer 400 sub-network label (parcels sharing
                                e.g. "LH_Default_PFC" merged)

Two ways to build the fly region graph:
  home    each neuron is assigned to the region holding most of its synapses
          and the neuron graph is summed over those groups (pure coarse-graining)
  bridge  regions A and B are linked by the number of neurons with at least
          BRIDGE_MIN_SYNAPSES synapses in both - closer to what tractography
          measures (a fibre running between two regions)

Each graph is scored with normalized.normalized_signature at its native density
and at a grid of matched densities (only where the native density allows).

Needs the public hemibrain export (python -m sctransfer.hemibrain --build, or
just the download it triggers) and the ENIGMA HCP matrices (hcp_data.py).

Usage:
    python -m sctransfer.run_coarse_grain
"""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from sctransfer import hcp_data, hemibrain
from sctransfer import normalized as N
from sctransfer.celegans import DATA_DIR as CELEGANS_DIR, build_weight_matrix as celegans_weight_matrix
from sctransfer.invertebrate_prior import FLY_SUBGRAPH_NODES
from sctransfer.topology import subsample_subgraph

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
DENSITIES = [0.03, 0.1, 0.2, 0.3]
MIN_REGION_SYNAPSES = 10_000   # drop hemibrain regions only fragmentarily inside the volume
BRIDGE_MIN_SYNAPSES = 10
N_NEURON_SAMPLES = 3


def coarse_grain(W, groups):
    """Sum a weighted adjacency over node groups; returns (W_coarse, group labels)."""
    labels, idx = np.unique(np.asarray(groups), return_inverse=True)
    M = np.zeros((len(labels), W.shape[0]))
    M[idx, np.arange(W.shape[0])] = 1
    C = M @ W @ M.T
    np.fill_diagonal(C, 0)
    return C, labels


def celegans_classes():
    names = pd.read_csv(CELEGANS_DIR / "celegans_connectome_nodes.txt", sep=r"\s+")["nodeLabel"].tolist()
    present = set(names)
    out = []
    for name in names:
        base = re.sub(r"\d+$", "", name)  # numbered motor neurons: VA01 -> VA
        if base != name:
            out.append(base)
        elif name[-1] in "LR" and name[:-1] + ("R" if name[-1] == "L" else "L") in present:
            out.append(name[:-1])  # bilateral pair: AVAL/AVAR -> AVA
        else:
            out.append(name)
    return out


def fly_region_graphs(export_dir):
    roi = pd.read_csv(export_dir / "traced-roi-connections.csv")
    roi = roi[roi["roi"] != "NotPrimary"]
    totals = roi.groupby("roi")["weight"].sum()
    regions = sorted(totals[totals >= MIN_REGION_SYNAPSES].index)
    roi = roi[roi["roi"].isin(regions)]

    # synapses per (neuron, region), counting both sides of each synapse
    per_neuron = pd.concat([
        roi.rename(columns={"bodyId_pre": "body"})[["body", "roi", "weight"]],
        roi.rename(columns={"bodyId_post": "body"})[["body", "roi", "weight"]],
    ]).groupby(["body", "roi"])["weight"].sum().reset_index()
    home = per_neuron.loc[per_neuron.groupby("body")["weight"].idxmax()].set_index("body")["roi"]

    r_index = {r: i for i, r in enumerate(regions)}
    edges = roi.groupby(["bodyId_pre", "bodyId_post"])["weight"].sum().reset_index()
    a = edges["bodyId_pre"].map(home).map(r_index).to_numpy()
    b = edges["bodyId_post"].map(home).map(r_index).to_numpy()
    W_home = np.zeros((len(regions), len(regions)))
    np.add.at(W_home, (a, b), edges["weight"].to_numpy())
    np.fill_diagonal(W_home, 0)

    member = per_neuron[per_neuron["weight"] >= BRIDGE_MIN_SYNAPSES]
    X = np.zeros((member["body"].nunique(), len(regions)))
    rows = pd.factorize(member["body"])[0]
    X[rows, member["roi"].map(r_index).to_numpy()] = 1
    W_bridge = X.T @ X
    np.fill_diagonal(W_bridge, 0)
    return {"home": W_home, "bridge": W_bridge}, regions


def schaefer_subnetworks(labels):
    """7Networks_LH_Default_PFC_3 -> LH_Default_PFC"""
    return [re.sub(r"^7Networks_|_\d+$", "", str(lab)) for lab in labels]


def collect_graphs(export_dir, seed):
    """{(species, level, name): [W, ...]}; neuron-level graphs get node subsamples."""
    graphs = {}
    W = celegans_weight_matrix()
    graphs[("worm", "neuron", "celegans")] = [W]
    graphs[("worm", "class", "celegans")] = [coarse_grain(W, celegans_classes())[0]]

    for name in hemibrain.available_circuits():
        d = np.load(hemibrain.CIRCUIT_DIR / f"{name}.npz")
        Wc = hemibrain.load_circuit(name)
        frac = min(1.0, FLY_SUBGRAPH_NODES / Wc.shape[0])
        graphs[("fly", "neuron", name)] = [subsample_subgraph(Wc, frac, seed=seed + i)
                                           for i in range(N_NEURON_SAMPLES)]
        types = [t if t else f"untyped_{b}" for t, b in zip(d["types"], d["body_ids"])]
        graphs[("fly", "type", name)] = [coarse_grain(Wc, types)[0]]

    region, regions = fly_region_graphs(export_dir)
    for how, Wr in region.items():
        graphs[("fly", f"region_{how}", "hemibrain")] = [Wr]
    print(f"[coarse_grain] fly regions: {len(regions)} hemibrain ROIs with >= {MIN_REGION_SYNAPSES} synapses")

    for parc in ["aparc", "schaefer_100", "schaefer_200", "schaefer_400", "glasser_360"]:
        loaded = hcp_data.load_group_connectome(parc)
        if loaded is None:
            continue
        sc, _, labels = loaded
        graphs[("human", "region", parc)] = [sc]
        if parc == "schaefer_400" and len(labels) == sc.shape[0]:
            graphs[("human", "subnetwork", parc)] = [coarse_grain(sc, schaefer_subnetworks(labels))[0]]
    return graphs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-nulls", type=int, default=N.N_NULLS_DEFAULT)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--export-dir", help="extracted hemibrain export (default: download into data/hemibrain/raw)")
    args = parser.parse_args()
    RESULTS_DIR.mkdir(exist_ok=True)
    export_dir = Path(args.export_dir) if args.export_dir else hemibrain._download_export()

    rows = []
    for (species, level, name), Ws in collect_graphs(export_dir, args.seed).items():
        for i, W in enumerate(Ws):
            native = N.density_of(N.binarize(W))
            for d in [None] + DENSITIES:
                if d is not None and native < 0.95 * d:
                    continue
                sig = N.normalized_signature(W, d, args.n_nulls, seed=args.seed + i)
                rows.append({"species": species, "level": level, "graph": name, "sample": i,
                             "target_density": "native" if d is None else d,
                             **{k: sig[k] for k in N.STATS + ["density", "n_nodes"]}})
        print(f"[coarse_grain] {species:5s} {level:14s} {name:16s} n={W.shape[0]:5d} native density {native:.3f}")

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS_DIR / "coarse_grain.csv", index=False)
    summary = (df.groupby(["target_density", "species", "level", "graph"], sort=False)
               [["n_nodes", "density"] + N.STATS].mean())
    pd.set_option("display.width", 200)
    print("\n=== Null-normalized statistics by node definition (1 = random) ===")
    print(summary.round(2).to_string())
    print("\nSaved results/coarse_grain.csv")


if __name__ == "__main__":
    main()
