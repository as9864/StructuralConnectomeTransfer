"""FlyWire whole adult brain (both hemispheres), coarse-grained to neuropils.

Source: FlyWire connectome v783 on Zenodo (record 10676866, CC BY 4.0),
file proofread_connections_783.feather - one row per (pre neuron, post neuron,
neuropil) with the synapse count. Download it into the gitignored
data/flywire/raw/ (about 850 MB):

    https://zenodo.org/api/records/10676866/files/proofread_connections_783.feather/content

Cite: Dorkenwald, S. et al. (2024). Neuronal wiring diagram of an adult brain.
Nature 634, 124-138; Schlegel, P. et al. (2024). Whole-brain annotation and
multi-connectome cell typing of Drosophila. Nature 634, 139-152.

Region graphs are built as in run_coarse_grain.py ("home": each neuron is
assigned to the neuropil holding most of its synapses and the neuron graph is
summed over those groups; "bridge": two neuropils are linked by the number of
neurons with >= BRIDGE_MIN_SYNAPSES synapses in both). Each hemisphere gives
one graph, which makes left and right two near-replicates of the same brain.
"""
from pathlib import Path

import numpy as np
import pandas as pd

RAW = Path(__file__).resolve().parent.parent / "data" / "flywire" / "raw" / "proofread_connections_783.feather"
MIN_SYNAPSES = 5          # FlyWire's own recommended threshold for a connection
BRIDGE_MIN_SYNAPSES = 10
OPTIC = ("ME", "AME", "LO", "LOP", "LA")


def load_connections(path=RAW):
    df = pd.read_feather(path, columns=["pre_pt_root_id", "post_pt_root_id", "neuropil", "syn_count"])
    df = df.rename(columns={"pre_pt_root_id": "pre", "post_pt_root_id": "post", "syn_count": "w"})
    return df[(df["neuropil"] != "UNASGD") & (df["pre"] != df["post"])]


def neuropil_side(name):
    for side in ("_L", "_R"):
        if name.endswith(side):
            return side[1]
    return "C"  # central, unpaired (e.g. PB, EB, FB, GNG)


def neuron_neuropil_counts(conn):
    """Synapses per (neuron, neuropil), both sides of each synapse counted."""
    a = conn[["pre", "neuropil", "w"]].rename(columns={"pre": "body"})
    b = conn[["post", "neuropil", "w"]].rename(columns={"post": "body"})
    return pd.concat([a, b]).groupby(["body", "neuropil"], observed=True)["w"].sum().reset_index()


class RegionGraphBuilder:
    """Precomputes neuron-level tables once so bootstrap resamples are cheap."""

    def __init__(self, conn, neuropils):
        self.neuropils = list(neuropils)
        idx = {p: i for i, p in enumerate(self.neuropils)}
        conn = conn[conn["neuropil"].isin(self.neuropils)]
        counts = neuron_neuropil_counts(conn)
        home = counts.loc[counts.groupby("body")["w"].idxmax()].set_index("body")["neuropil"]
        self.bodies = home.index.to_numpy()
        body_idx = pd.Series(np.arange(len(self.bodies)), index=self.bodies)
        self.home = home.map(idx).to_numpy()

        pairs = conn.groupby(["pre", "post"], observed=True)["w"].sum().reset_index()
        pairs = pairs[pairs["w"] >= MIN_SYNAPSES]
        self.pre = body_idx[pairs["pre"]].to_numpy()
        self.post = body_idx[pairs["post"]].to_numpy()
        self.w = pairs["w"].to_numpy().astype(float)

        member = counts[counts["w"] >= BRIDGE_MIN_SYNAPSES]
        self.m_body = body_idx[member["body"]].to_numpy()
        self.m_np = member["neuropil"].map(idx).to_numpy()

    def build(self, rng=None):
        """(home, bridge) region graphs; with `rng`, neurons are resampled with
        replacement (a neuron drawn k times contributes k-fold)."""
        n_b, n_r = len(self.bodies), len(self.neuropils)
        mult = (np.bincount(rng.integers(0, n_b, n_b), minlength=n_b).astype(float)
                if rng is not None else np.ones(n_b))
        H = np.zeros((n_r, n_r))
        np.add.at(H, (self.home[self.pre], self.home[self.post]), self.w * mult[self.pre] * mult[self.post])
        np.fill_diagonal(H, 0)
        X = np.zeros((n_b, n_r))
        X[self.m_body, self.m_np] = 1
        B = (X * mult[:, None]).T @ X
        np.fill_diagonal(B, 0)
        return {"home": H, "bridge": B}


def hemisphere_neuropils(conn, side, include_optic=False, min_synapses=1_000):
    """Neuropils of one hemisphere ("L"/"R"), or of the whole brain ("both":
    both hemispheres plus the unpaired central neuropils)."""
    totals = conn.groupby("neuropil", observed=True)["w"].sum()
    out = []
    for p, total in totals.items():
        if (side != "both" and neuropil_side(p) != side) or total < min_synapses or p == "OCG":
            continue
        if not include_optic and p.rsplit("_", 1)[0] in OPTIC:
            continue
        out.append(p)
    return sorted(out)
