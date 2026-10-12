# Larval *Drosophila* whole-brain connectome

Nothing here is committed except this README. `sctransfer/larva.py` reads `raw/G_edgelist.txt`
(pre, post, synapse count) and `raw/meta_data.csv` (gitignored, about 5 MB), downloaded from the
analysis repository of Pedigo et al. (2023):
<https://github.com/neurodata/bilateral-connectome/tree/main/data/elife>

- **Cite**: Winding, M., Pedigo, B.D. et al. (2023). The connectome of an insect brain. *Science*
  379, eadd9330. doi:10.1126/science.add9330 — and for the processed files, Pedigo, B.D. et al.
  (2023). Generative network modeling reveals quantitative definitions of bilateral symmetry
  exhibited by a whole insect brain connectome. *eLife* 12, e83739.
- 3,554 neurons in the edge list (brain and sensory inputs of a first-instar larva); 2,952 in the
  curated largest connected component used by default.
