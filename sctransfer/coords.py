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
            --fetch-larva` caches them in data/larva/raw/larva_positions.csv.

All positions are in the dataset's own units (worm: atlas units, fly: nm);
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


def flywire_positions():
    """{root_id: (x, y, z)} soma position, falling back to the annotated point on the neuron."""
    df = pd.read_csv(FLYWIRE_ANN, sep="\t", usecols=["root_id", "pos_x", "pos_y", "pos_z",
                                                    "soma_x", "soma_y", "soma_z"], low_memory=False)
    out = {}
    for c in ("x", "y", "z"):
        df[c] = df[f"soma_{c}"].where(df[f"soma_{c}"].notna() & (df[f"soma_{c}"] != 0), df[f"pos_{c}"])
    for rid, x, y, z in df[["root_id", "x", "y", "z"]].itertuples(index=False):
        out[int(rid)] = (float(x), float(y), float(z))
    return out


def larva_positions():
    df = pd.read_csv(LARVA_POS)
    return {int(r.skeleton_id): (r.x, r.y, r.z) for r in df.itertuples(index=False)}


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
    parser.add_argument("--fetch-larva", action="store_true")
    args = parser.parse_args()
    if args.fetch_larva:
        from sctransfer import larva
        _, meta = larva.load(selected_lcc=False)
        fetch_larva([int(s) for s in meta.index])


if __name__ == "__main__":
    main()
