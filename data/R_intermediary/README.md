# R_intermediary — frozen phenograph clusterings

Phenograph clustering in this project runs through `Rphenoannoy`, which uses
**Annoy** (approximate nearest-neighbour search, C++) followed by
`igraph::cluster_louvain`. Neither step is bit-reproducible across `igraph`
versions, compilers or CPUs, so the **same seed and same input data can produce
different clusters on a different machine or container** — the same caveat the
`uwot`/UMAP documentation gives for `calculateUMAP`.

To make the published results reproducible everywhere, the clusterings are
computed once in the reference (docker) environment and shipped as `.rds` files.
`analysis(_AD)/02_QC_1.Rmd` and `analysis(_AD)/04_phenotyping.Rmd` **load** these
files when present (keyed by unique cell id, `sample_id_ObjectNumber`) instead of
recomputing, and only fall back to computing (and saving) when they are absent.

## Files

| file | produced by | contents |
|------|-------------|----------|
| `pg_clusters_lostcells_AD.rds`     | `analysis_AD/02_QC_1.Rmd` | named character vector: cell id → lost-cells cluster (k=45) |
| `pg_clusters_phenotyping_AD.rds`   | `analysis_AD/04_phenotyping.Rmd` | data.frame (rownames = cell id) with `pg_clusters_k15/30/45/60` |
| `pg_clusters_lostcells_PTBM.rds`   | `analysis/02_QC_1.Rmd` | named character vector: cell id → lost-cells cluster (k=45) |
| `pg_clusters_phenotyping_PTBM.rds` | `analysis/04_phenotyping.Rmd` | data.frame (rownames = cell id) with `pg_clusters_k15/30/45/60` |

With these loaded, the hard-coded lost-cells cluster is **17** for AD and **26**
for PT/BM, and the manual celltype annotation stays valid because the cluster
numbers are frozen.

## How to (re)generate

Run once in the reference docker image (see the header of
`code/generate_phenograph_intermediary.R` for details):

```sh
Rscript code/generate_phenograph_intermediary.R
```

The lookup order at load time is `params$input/R_intermediary/` first, then
`params$output/R_intermediary/`. Place the shipped files next to the input data
(`params$input/R_intermediary/`), or download them from Zenodo alongside the
other `R_intermediary` objects.

> These `.rds` files are **not** committed here (they depend on the dataset).
> This folder and README document the scheme and give the loader a stable path.
