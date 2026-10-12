# FlyWire whole adult *Drosophila* brain (v783)

Nothing in this folder is committed except this README. `sctransfer/flywire.py` reads
`raw/proofread_connections_783.feather` (gitignored, about 850 MB), downloaded from:

<https://zenodo.org/api/records/10676866/files/proofread_connections_783.feather/content>

- **Source**: FlyWire Whole-brain Connectome Connectivity Data, Zenodo record 10676866
- **License**: CC BY 4.0
- **Cite**: Dorkenwald, S. et al. (2024). Neuronal wiring diagram of an adult brain. *Nature* 634,
  124–138. Schlegel, P. et al. (2024). Whole-brain annotation and multi-connectome cell typing of
  *Drosophila*. *Nature* 634, 139–152.

One row per (pre neuron, post neuron, neuropil) with the synapse count. `flywire.py` drops
unassigned synapses, keeps neuron pairs with at least 5 synapses (FlyWire's recommended
threshold), and coarse-grains the central brain (optic lobes excluded) to neuropil graphs.
