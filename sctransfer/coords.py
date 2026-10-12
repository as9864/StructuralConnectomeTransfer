"""Neuron positions for the spatial term of the wiring-rule model (E4).

worm        3D neuron atlas of Skuhersky et al. (2022, BMC Bioinformatics,
            data CC0): https://github.com/bluevex/elegans-atlas
            (LowResAtlasWithHighResHeadsAndTails.csv) -> data/coords/raw/
adult fly   FlyWire soma positions (Schlegel et al. 2024 annotations):
            https://github.com/flyconnectome/flywire_annotations
            (supplemental_files/Supplemental_file1_neuron_annotations.tsv)
            -> data/coords/raw/
larva       soma node of each skeleton on the public L1EM CATMAID server
            hosted by Virtual Fly Brain, root node where no soma is labelled
            (sensory neurons enter through nerves). `python -m sctransfer.coords
            --fetch-larva` caches them in data/larva/raw/larva_positions.csv;
            `--fetch-larva-synapses` caches each neuron's synapse centroid.

Fly positions are in nm (FlyWire voxels are scaled by 4x4x40), the worm atlas
in its own units;
the generative model only uses distances relative to the graph's own scale.
"""
import argparse
import json
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent / "data"
WORM_ATLAS = ROOT / "coords" / "raw" / "LowResAtlasWithHighResHeadsAndTails.csv"
FLYWIRE_ANN = ROOT / "coords" / "raw" / "Supplemental_file1_neuron_annotations.tsv"
LARVA_STATS = ROOT / "larva" / "raw" / "catmaid_label_stats.json"
LARVA_POS = ROOT / "larva" / "raw" / "larva_positions.csv"
CATMAID = "https://l1em.catmaid.virtualflybrain.org/1"


def worm_positions():
    """{neuron name: (x, y, z)}"""
    df = pd.read_csv(WORM_ATLAS, header=None, names=["name", "x", "y", "z"])
    return {str(r.name).strip(): (r.x, r.y, r.z) for r in df.itertuples(index=False)}


FLYWIRE_VOXEL_NM = (4.0, 4.0, 40.0)  # annotation coordinates are in 4x4x40 nm voxels


def flywire_positions(kind="soma"):
    """{root_id: (x, y, z)} in nm. kind="soma": cell body, falling back to the
    anchor; kind="anchor": the annotated anchor point, "typically on the
    backbone of the neuron" (inside the neuropil/tracts rather than the cell-body
    rind), our proxy for where the neuron's synapses are."""
    df = pd.read_csv(FLYWIRE_ANN, sep="\t", usecols=["root_id", "pos_x", "pos_y", "pos_z",
                                                    "soma_x", "soma_y", "soma_z"], low_memory=False)
    for c, scale in zip(("x", "y", "z"), FLYWIRE_VOXEL_NM):
        if kind == "soma":
            df[c] = df[f"soma_{c}"].where(df[f"soma_{c}"].notna() & (df[f"soma_{c}"] != 0), df[f"pos_{c}"])
        else:
            df[c] = df[f"pos_{c}"]
        df[c] = df[c].astype(float) * scale
    return {int(rid): (x, y, z) for rid, x, y, z in df[["root_id", "x", "y", "z"]].itertuples(index=False)}


FLYWIRE_SYN_POS = ROOT / "coords" / "raw" / "flywire_synapse_centroids_approx.csv"
MIN_ANCHOR_SOMA_NM = 10_000  # anchors this far from the soma are on the backbone, not the cell body


def flywire_synapse_centroids():
    """Approximate synapse centroid per FlyWire neuron, without the 9.5 GB
    synapse table: (1) estimate each neuropil's centre as the synapse-weighted
    mean of backbone anchors (anchors >= 10 um from the soma) of neurons with
    synapses there; (2) place each neuron at the synapse-weighted mean of the
    centres of the neuropils it has synapses in. Neurons confined to one
    neuropil share that neuropil's centre, so within-neuropil distances are
    lost; cached in data/coords/raw/."""
    if FLYWIRE_SYN_POS.exists():
        df = pd.read_csv(FLYWIRE_SYN_POS)
        return {int(r.root_id): (r.x, r.y, r.z) for r in df.itertuples(index=False)}
    from sctransfer import flywire
    counts = flywire.neuron_neuropil_counts(flywire.load_connections())
    soma, anchor = flywire_positions("soma"), flywire_positions("anchor")
    backbone = {k: v for k, v in anchor.items()
                if k in soma and np.linalg.norm(np.subtract(v, soma[k])) >= MIN_ANCHOR_SOMA_NM}
    c = counts[counts["body"].isin(backbone)].copy()
    xyz = np.array([backbone[b] for b in c["body"]])
    for k, axis in enumerate("xyz"):
        c[axis] = xyz[:, k] * c["w"]
    g = c.groupby("neuropil", observed=True)[["x", "y", "z", "w"]].sum()
    centre = {p: (r.x / r.w, r.y / r.w, r.z / r.w) for p, r in g.iterrows()}
    counts = counts[counts["neuropil"].isin(centre)].copy()
    cxyz = np.array([centre[p] for p in counts["neuropil"]])
    for k, axis in enumerate("xyz"):
        counts[axis] = cxyz[:, k] * counts["w"]
    pos = counts.groupby("body")[["x", "y", "z", "w"]].sum()
    out = pd.DataFrame({"root_id": pos.index, "x": pos.x / pos.w, "y": pos.y / pos.w, "z": pos.z / pos.w,
                        "n_synapses": pos.w})
    out.to_csv(FLYWIRE_SYN_POS, index=False)
    return {int(r.root_id): (r.x, r.y, r.z) for r in out.itertuples(index=False)}


LARVA_SYN_POS = ROOT / "larva" / "raw" / "larva_synapse_centroids.csv"


def larva_positions(kind="soma"):
    """{skeleton_id: (x, y, z)} in nm. kind="soma": soma node (root where no soma
    is labelled); kind="synapse": centroid of the neuron's synapse (connector)
    positions, see fetch_larva_synapses."""
    df = pd.read_csv(LARVA_POS if kind == "soma" else LARVA_SYN_POS)
    return {int(r.skeleton_id): (r.x, r.y, r.z) for r in df.itertuples(index=False)}


def fetch_larva_synapses(skeleton_ids, batch=50, pause=0.5):
    """Centroid of every connector linked to each skeleton (pre or post), from
    the public L1EM CATMAID: one GET per batch of skeletons."""
    sums = {}
    skeleton_ids = list(skeleton_ids)
    for i in range(0, len(skeleton_ids), batch):
        part = skeleton_ids[i:i + batch]
        query = "&".join(f"skeleton_ids%5B{k}%5D={s}" for k, s in enumerate(part))
        try:
            links = _get(f"{CATMAID}/connectors/links/?{query}")["links"]
        except Exception as exc:  # noqa: BLE001 - report and continue; missing ones fall back to soma
            print(f"[coords] batch {i}: {exc}")
            continue
        seen = set()
        for skid, connector, x, y, z, *_ in links:
            if (skid, connector) in seen:  # one point per connector per neuron
                continue
            seen.add((skid, connector))
            s = sums.setdefault(int(skid), [0.0, 0.0, 0.0, 0])
            s[0] += x
            s[1] += y
            s[2] += z
            s[3] += 1
        print(f"[coords] larval synapses {i + len(part)}/{len(skeleton_ids)}", flush=True)
        time.sleep(pause)
    rows = [dict(skeleton_id=k, x=v[0] / v[3], y=v[1] / v[3], z=v[2] / v[3], n_synapses=v[3])
            for k, v in sums.items()]
    pd.DataFrame(rows).to_csv(LARVA_SYN_POS, index=False)
    print(f"[coords] saved {len(rows)} larval synapse centroids to {LARVA_SYN_POS}")


def _get(url):
    with urllib.request.urlopen(url, timeout=60) as r:
        return json.loads(r.read())


def fetch_larva(skeleton_ids, pause=0.15):
    """Soma (or root) position of every skeleton, one polite GET per neuron."""
    stats = json.load(open(LARVA_STATS))
    soma = {}
    for _, label, skid, node in stats:
        if label == "soma":
            soma.setdefault(int(skid), int(node))
    done = {}
    if LARVA_POS.exists():
        done = {int(r.skeleton_id): r for r in pd.read_csv(LARVA_POS).itertuples(index=False)}
    rows = [dict(skeleton_id=k, x=v.x, y=v.y, z=v.z, source=v.source) for k, v in done.items()]
    todo = [s for s in skeleton_ids if s not in done]
    for i, skid in enumerate(todo):
        try:
            if skid in soma:
                d = _get(f"{CATMAID}/treenodes/{soma[skid]}/compact-detail")
                rows.append(dict(skeleton_id=skid, x=d[2], y=d[3], z=d[4], source="soma"))
            else:
                d = _get(f"{CATMAID}/skeletons/{skid}/root")
                rows.append(dict(skeleton_id=skid, x=d["x"], y=d["y"], z=d["z"], source="root"))
        except Exception as exc:  # noqa: BLE001 - skip and report; rerun resumes
            print(f"[coords] skeleton {skid}: {exc}")
        if i % 200 == 0:
            pd.DataFrame(rows).to_csv(LARVA_POS, index=False)
            print(f"[coords] larva {i}/{len(todo)}", flush=True)
        time.sleep(pause)
    pd.DataFrame(rows).to_csv(LARVA_POS, index=False)
    print(f"[coords] saved {len(rows)} larval positions to {LARVA_POS}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fetch-larva", action="store_true", help="soma/root positions")
    parser.add_argument("--fetch-larva-synapses", action="store_true", help="synapse centroids")
    args = parser.parse_args()
    from sctransfer import larva
    if args.fetch_larva:
        _, meta = larva.load(selected_lcc=False)
        fetch_larva([int(s) for s in meta.index])
    if args.fetch_larva_synapses:
        _, meta = larva.load()
        fetch_larva_synapses([int(s) for s in meta.index])


if __name__ == "__main__":
    main()
