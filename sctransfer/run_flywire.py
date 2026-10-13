"""Region-level fly vs human comparison with uncertainty (FlyWire whole brain).

Repeats the hemibrain region pilot (run_coarse_grain.py) on FlyWire, which
covers both hemispheres: each hemisphere's central-brain neuropils form one
region graph, and neurons are bootstrapped to give confidence intervals. The
human side is split by hemisphere the same way (cortex only, ENIGMA HCP group
SC), so both species contribute left/right near-replicates at a similar node
count. All statistics are null-normalized at matched densities.

Usage:
    python -m sctransfer.run_flywire --n-boot 20
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from sctransfer import flywire, hcp_data
from sctransfer import normalized as N
from sctransfer.run_coarse_grain import coarse_grain, schaefer_subnetworks

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
DENSITIES = [0.1, 0.2, 0.3]
HUMAN_PARCELLATIONS = ["aparc", "schaefer_100", "schaefer_200"]


def human_hemisphere_graphs():
    """{(level, parcellation, side): W} from the ENIGMA group matrices."""
    out = {}
    for parc in HUMAN_PARCELLATIONS + ["schaefer_400"]:
        loaded = hcp_data.load_group_connectome(parc)
        if loaded is None:
            continue
        sc, _, labels = loaded
        labels = np.array([str(x) for x in labels])
        for side, tag in (("L", ("L_", "LH_", "7Networks_LH")), ("R", ("R_", "RH_", "7Networks_RH"))):
            keep = np.array([lab.startswith(tag) for lab in labels])
            if not keep.any():
                continue
            W, lab = sc[np.ix_(keep, keep)], labels[keep]
            if parc == "schaefer_400":
                out[("subnetwork", parc, side)] = coarse_grain(W, schaefer_subnetworks(lab))[0]
            else:
                out[("region", parc, side)] = W
        if parc == "schaefer_400":
            out[("subnetwork", parc, "both")] = coarse_grain(sc, schaefer_subnetworks(labels))[0]
        else:
            out[("region", parc, "both")] = sc
    for parc in ["aparc", "schaefer_100"]:  # the fly whole brain includes non-cortical hubs too
        loaded = hcp_data.load_group_connectome(parc, with_subcortex=True)
        if loaded is not None:
            out[("region+subcortex", parc, "both")] = loaded[0]
    return out


def score(W, densities, n_nulls, seed):
    native = N.density_of(N.binarize(W))
    for d in densities:
        if native < 0.95 * d:
            continue
        sig = N.normalized_signature(W, d, n_nulls, seed=seed)
        yield d, sig


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-boot", type=int, default=20)
    parser.add_argument("--n-nulls", type=int, default=N.N_NULLS_DEFAULT)
    parser.add_argument("--include-optic", action="store_true", help="keep optic-lobe neuropils")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    RESULTS_DIR.mkdir(exist_ok=True)
    rng = np.random.default_rng(args.seed)

    rows = []
    conn = flywire.load_connections()
    print(f"[flywire] {len(conn):,} neuron-pair-neuropil rows")
    for side in ("L", "R", "both"):
        nps = flywire.hemisphere_neuropils(conn, side, include_optic=args.include_optic)
        builder = flywire.RegionGraphBuilder(conn, nps)
        print(f"[flywire] hemisphere {side}: {len(nps)} neuropils, {len(builder.bodies):,} neurons")
        for b in range(args.n_boot + 1):
            graphs = builder.build(None if b == 0 else rng)
            for how, W in graphs.items():
                for d, sig in score(W, DENSITIES, args.n_nulls, seed=args.seed + b):
                    rows.append({"species": "fly", "level": f"region_{how}", "graph": "flywire", "side": side,
                                 "boot": b, "density": d, "n_nodes": W.shape[0], **{k: sig[k] for k in N.STATS}})
            if b % 5 == 0:
                print(f"[flywire] hemisphere {side}: bootstrap {b}/{args.n_boot}")

    for (level, parc, side), W in human_hemisphere_graphs().items():
        for d, sig in score(W, DENSITIES, args.n_nulls, seed=args.seed):
            rows.append({"species": "human", "level": level, "graph": parc, "side": side, "boot": 0,
                         "density": d, "n_nodes": W.shape[0], **{k: sig[k] for k in N.STATS}})
        print(f"[human] {parc} {level} {side}: {W.shape[0]} nodes")

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS_DIR / "flywire_regions.csv", index=False)

    def lo(x):
        return np.nanpercentile(x, 2.5)

    def hi(x):
        return np.nanpercentile(x, 97.5)

    pd.set_option("display.width", 220)
    df["scope"] = np.where(df["side"] == "both", "whole brain", "hemisphere")
    keys = ["scope", "density", "species", "level", "graph"]
    print("\n=== Point estimates (unresampled; hemisphere rows average L and R) ===")
    point = df[df.boot == 0].groupby(keys)[["n_nodes"] + N.STATS].mean()
    print(point.round(2).to_string())
    print("\n=== Fly bootstrap 95% intervals (neurons resampled; hemisphere rows pool L and R) ===")
    boot = df[(df.species == "fly") & (df.boot > 0)]
    print(boot.groupby(keys)[N.STATS].agg([lo, hi]).round(2).to_string())
    print("\nSaved results/flywire_regions.csv")


if __name__ == "__main__":
    main()
