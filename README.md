# omni-tenxbrain-ladder

A cell-count ladder over **TENxBrainData** — 1.3M E18 mouse brain cells, 10x
HDF5 — built as an omnibenchmark of its own. The dataset never went to GEO/SRA;
10x published it directly and Bioconductor wraps it, so its one real accession
is the ExperimentHub id **EH1039** (`TENxBrainData`).

```
RAW      omni-data  --hapiq--> EH1039, checksum-verified, once
LADDER   run.py     --scsampler draw + streaming gather--> one rung per module
```

Each rung is a publishable DATA directory:

```
tenx-0020k.h5ad                 counts in layers/counts (CSR, f32), no X,
tenx-0020k.clusters_truth.tsv   flat string obs columns
tenx-0020k.clusters_truth_num.txt
tenx-0020k_properties.yaml
```

## Why a ladder

Rungs double, so the whole ladder costs ~2x its top rung and the points sit
evenly in log n — what a power-law fit of time/memory vs n wants. Draws are
**nested** (rung k ⊂ rung k+1) and composition-preserving (proportional within
`Library`, the 133 sequencing libraries), so a change in a downstream number is
"more cells", not "different cells". Both properties come from
[omni-scsampler](https://github.com/btraven00/omni-scsampler).

## Run

```sh
ob run benchmark -b benchmark.yaml            # every rung
./run.py --selftest                           # blocks() + nesting, no data
```

The top rung is 1.31M cells / 2.6e9 nonzeros. Only the *selection* is delegated
to scsampler — the row gather streams, because the top rungs do not fit in
memory as scipy CSR and anndata's backed mode cannot materialize them either.
Peak RAM is bounded by `SPAN` (~200MB of source data+indices), not by rung size.

## Consuming a rung

Rungs land in `<out>/tenxbrain/RAW/tenxbrain/<hash>/LADDER/<rung>/<hash>/tenxbrain/`.
Point a downstream DATA stage at that directory with omni-data, using the rung
id as the module id (the file stems must match it):

```yaml
- id: tenx-0020k
  software_environment: omni_data
  repository:
    url: https://github.com/btraven00/omni-data
    commit: eb4d368b5e83b26389e6047dd8bad6badede78f8
  parameters:
    - uri: file:///…/LADDER/tenx-0020k/<hash>/tenxbrain
```

## Caveats

- **No cell-type labels exist for this dataset.** The truth files are emitted
  header-only on purpose — enough for the stage contract, useless as a label
  source. Any metric stage that consumes them is invalid downstream of here.
- 640k is the last **R**-readable rung: at 1.31M cells nnz passes 2^31, and
  anndataR silently drops `layers/counts` rather than failing.
- The rung modules point at a local checkout in `benchmark.yaml`; swap to the
  GitHub URL and pin a commit before publishing. Same for the scsampler pin in
  `envs/tenx-ladder.yml`.
