# Neuron positions (spatial term of the wiring-rule model)

Nothing here is committed except this README; `sctransfer/coords.py` reads the files in `raw/`.

- `raw/LowResAtlasWithHighResHeadsAndTails.csv` — 3D *C. elegans* neuron atlas from
  <https://github.com/bluevex/elegans-atlas>. Cite Skuhersky, M. et al. (2022). Toward a more accurate
  3D atlas of *C. elegans* neurons. *BMC Bioinformatics* 23, 195 (data CC0 per the article).
- `raw/Supplemental_file1_neuron_annotations.tsv` — FlyWire v783 neuron annotations incl. soma
  positions, from <https://github.com/flyconnectome/flywire_annotations>. Cite Schlegel, P. et al. (2024),
  *Nature* 634, 139–152.
- Larval positions are cached in `data/larva/raw/larva_positions.csv` by
  `python -m sctransfer.coords --fetch-larva` (one request per neuron to the public L1EM CATMAID
  server hosted by Virtual Fly Brain).
