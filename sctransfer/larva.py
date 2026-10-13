"""Larval Drosophila whole-brain connectome (Winding, Pedigo et al. 2023).

Like C. elegans, this is a whole nervous system of a few thousand neurons
(brain plus sensory inputs), so it matches the worm in scope without any
sampling trick. Files come from the analysis repository of Pedigo et al.
(2023, eLife), which redistributes the Winding et al. (2023) connectome:

  https://raw.githubusercontent.com/neurodata/bilateral-connectome/main/data/elife/G_edgelist.txt
  https://raw.githubusercontent.com/neurodata/bilateral-connectome/main/data/elife/meta_data.csv

Download both into the gitignored data/larva/raw/ (about 5 MB). See
data/larva/README.md for citations.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

RAW = Path(__file__).resolve().parent.parent / "data" / "larva" / "raw"


def load(selected_lcc=True):
    """Sparse directed W[pre, post] = synapse count, and the node metadata.
    `selected_lcc` keeps the neurons Pedigo et al. analysed (largest
    connected component of the curated set)."""
    edges = pd.read_csv(RAW / "G_edgelist.txt", sep=" ", header=None, names=["pre", "post", "w"])
    meta = pd.read_csv(RAW / "meta_data.csv", index_col=0)
    nodes = np.unique(np.concatenate([edges["pre"].to_numpy(), edges["post"].to_numpy()]))
    if selected_lcc:
        nodes = np.intersect1d(nodes, meta.index[meta["selected_lcc"].astype(bool)])
    edges = edges[edges["pre"].isin(nodes) & edges["post"].isin(nodes) & (edges["pre"] != edges["post"])]
    i = np.searchsorted(nodes, edges["pre"].to_numpy())
    j = np.searchsorted(nodes, edges["post"].to_numpy())
    W = sp.csr_matrix((edges["w"].to_numpy().astype(float), (i, j)), shape=(len(nodes), len(nodes)))
    return W, meta.reindex(nodes)
