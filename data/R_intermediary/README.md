# R_intermediary — frozen phenograph clusterings

Phenograph clustering in this project runs through `Rphenoannoy`, which uses
**Annoy** (approximate nearest-neighbour search, C++) followed by
`igraph::cluster_louvain`. Neither step is bit-reproducible across `igraph`
versions, compilers or CPUs, so the **same seed and same input data can produce
different clusters on a different machine or container** — the same caveat the
`uwot`/UMAP documentation gives for `calculateUMAP`.

To make the published results reproducible everywhere, the clusterings computed
in the reference (docker) environment are shipped as `.rds` files.
`analysis(_AD)/02_QC_1.Rmd` and `analysis(_AD)/04_phenotyping.Rmd` **load** these
files when present (keyed by unique cell id, `sample_id_ObjectNumber`) instead of
recomputing, and only fall back to computing (and saving) when they are absent.

## Files

Each file is a **named character vector** (`names` = unique cell id → cluster):

| file | used by | contents |
|------|---------|----------|
| `pg_clusters_lostcells_AD.rds`   | `analysis_AD/02_QC_1.Rmd` | lost-cells cluster (k=45), retained cells; ghosts rebuilt at load |
| `pg_clusters_AD.rds`             | `analysis_AD/04_phenotyping.Rmd` | phenograph cluster (k=30) → `pg_clusters` |
| `pg_clusters_lostcells_PTBM.rds` | `analysis/02_QC_1.Rmd` | lost-cells cluster (k=45), retained cells; ghosts rebuilt at load |
| `pg_clusters_PTBM.rds`           | `analysis/04_phenotyping.Rmd` | phenograph cluster (k=30, before GMM subclustering) |

With these loaded, the hard-coded lost-cells cluster is **17** for AD and **26**
for PT/BM, and the manual celltype annotation stays valid because the phenograph
cluster numbers are frozen.

## How to generate

The labels already exist in the docker SPE objects (`pg_clusters_lostcells`,
`pg_clusters_k30`), so nothing is re-clustered. Run once in the reference docker
image (see `code/generate_phenograph_intermediary.R`):

```sh
Rscript code/generate_phenograph_intermediary.R      # reads spe_*.rds, writes here
```

It reads:
* `pg_clusters_k30` from `spe_04_phenotyping*.rds` → `pg_clusters_*.rds`
* `pg_clusters_lostcells` (retained cells only) from `spe_02_QC_1*.rds` →
  `pg_clusters_lostcells_*.rds`. The excluded ghost cells are rebuilt at load
  time in `02_QC_1.Rmd` (any post-QC1 cell missing from the file is a ghost and
  is labelled with the lost-cells cluster number, 17 for AD / 26 for PT/BM).

The load lookup order is `params$input/R_intermediary/` first, then
`params$output/R_intermediary/`. Ship the files next to the input data
(`params$input/R_intermediary/`) or download them from Zenodo alongside the
other `R_intermediary` objects.

> These `.rds` files are **not** committed here (they depend on the dataset).
> This folder and README document the scheme and give the loader a stable path.
