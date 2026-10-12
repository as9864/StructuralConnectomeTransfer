"""C. elegans brain (nerve ring) connectomes of eight isogenic animals across
development (Witvliet et al. 2021, Nature 596, 257-261): datasets 1-4 L1,
5 L2, 6 L3, 7-8 adult. Gives the worm side a within-species baseline.

The spreadsheets (edge lists: pre, post, type, synapses) come from WormWiring
via the OpenWorm ConnectomeToolbox repository (MIT):
  https://github.com/openworm/ConnectomeToolbox/tree/main/cect/data
  ("witvliet_2020_<k> <stage>.xlsx"). Download into the gitignored
data/witvliet/raw/ (about 300 KB in total). Read with the standard library
(xlsx is zipped XML), so no spreadsheet package is needed.
"""
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

import numpy as np

RAW = Path(__file__).resolve().parent.parent / "data" / "witvliet" / "raw"
DATASETS = {1: "L1", 2: "L1", 3: "L1", 4: "L1", 5: "L2", 6: "L3", 7: "adult", 8: "adult"}
_NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def _read_rows(path):
    with zipfile.ZipFile(path) as z:
        strings = [''.join(t.itertext()) for t in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", _NS)]
        sheet = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
    rows = []
    for r in sheet.findall(".//m:row", _NS):
        vals = []
        for c in r.findall("m:c", _NS):
            v = c.find("m:v", _NS)
            val = None if v is None else v.text
            if c.get("t") == "s" and val is not None:
                val = strings[int(val)]
            vals.append(val)
        rows.append(vals)
    return rows


def load(dataset, synapse_type="chemical"):
    """Directed W[pre, post] = synapse count for one animal, and node names."""
    path = RAW / f"witvliet_2020_{dataset} {DATASETS[dataset]}.xlsx"
    rows = _read_rows(path)
    header, body = rows[0], rows[1:]
    i_pre, i_post, i_type, i_syn = (header.index(k) for k in ("pre", "post", "type", "synapses"))
    edges = [(r[i_pre], r[i_post], float(r[i_syn])) for r in body
             if len(r) > i_syn and r[i_type] == synapse_type and r[i_pre] != r[i_post]]
    names = sorted({e[0] for e in edges} | {e[1] for e in edges})
    idx = {n: k for k, n in enumerate(names)}
    W = np.zeros((len(names), len(names)))
    for pre, post, w in edges:
        W[idx[pre], idx[post]] += w
    return W, names
