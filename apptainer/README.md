# MapMetSC — Apptainer image

A single [Apptainer](https://apptainer.org/) (formerly Singularity) image that
runs the **entire** R analysis for Lazic et al.:

- `analysis/` — `01_read_data` … `06.*_correlation`
- `analysis_AD/` — the adrenal-gland reference pipeline
- `analysis/07_spatial_analysis.Rmd` — the spatial analysis

It is built **from scratch** on the same pinned base the project's `Dockerfile`
uses (`rocker/rstudio:4.4.0`), so it is reproducible and self-contained rather
than a wrapper around the pre-built Docker Hub image.

## Why one image instead of two

The project historically needed a second image
([`lazdaria/mapmetsc_spatial:v1.0`](https://hub.docker.com/repository/docker/lazdaria/mapmetsc_spatial/general))
purely because of a bug in `imcRtools::testInteractions` in versions **< 1.5.5**.

That fix first ships in the released **imcRtools 1.6.0** (Bioconductor 3.17).
Because this image is built on `rocker/rstudio:4.4.0` → **Bioconductor 3.19**, it
gets **imcRtools 1.10.0**, which already contains the fix. The only `imcRtools`
function the non-spatial scripts use is `read_steinbock()`, whose API is
unchanged across these versions — so moving everything onto the fixed
`imcRtools` does not alter the rest of the pipeline. The build **asserts**
`imcRtools >= 1.6.0` (and that `testInteractions` is exported), so it can never
silently produce the buggy version. One image therefore covers both pipelines.

`ggdendro` (used by `05.1_PT_DE_cellularcomm.Rmd` and `07_spatial_analysis.Rmd`),
along with the other packages the original `Dockerfile` was missing —
`EBImage`, `batchelor`, `harmony`, `ggrastr`, `hexbin`, `lme4`, `writexl`,
`workflowr` — is included.

## Build

```bash
# from the repository root
./apptainer/build.sh            # -> ./mapmetsc.sif
# or explicitly
apptainer build mapmetsc.sif apptainer/mapmetsc.def
```

Building pulls `rocker/rstudio:4.4.0` from Docker Hub and downloads the
CRAN/Bioconductor packages, so run it somewhere with outbound network access.
`apptainer build` needs root or `--fakeroot` (the helper script adds
`--fakeroot` automatically when you are not root — this is available on most HPC
systems). Expect a multi-GB `.sif` and a build time dominated by compiling the
Bioconductor stack.

Bioconductor is pinned to release **3.19** (the one matching the rocker image's
frozen CRAN snapshot, which avoids CRAN/Bioc version conflicts). The 3.19
repositories are added to `options(repos)` and installed with plain
`install.packages()` — deliberately *not* via `BiocManager::install(version=)`,
whose version-switch handshake needs a `bioconductor.org` download and package
downgrades that fail behind a strict HPC egress proxy. Packages are installed
with a retry-only-missing loop (many rounds, short sleeps), so a cluster where
`bioconductor.org` returns intermittent `504`s still builds. If a package
genuinely can't be fetched, the build fails loudly at the completeness check
rather than shipping a half-installed image.

No local Apptainer? A manual GitHub Actions workflow
([`.github/workflows/build-apptainer.yml`](../.github/workflows/build-apptainer.yml))
builds the image and uploads it as an artifact — run it from the repo's
**Actions** tab.

## Run

The scripts read/write everything under `/mnt/data` (see the main
[README](../README.md) for the expected `input_all` / `input_AD` /
`public_datasets` / `output` layout). Bind your extracted single-cell data
there:

```bash
# render the whole pipeline in order (headless)
apptainer exec --bind /path/to/singlecelldata:/mnt/data \
    mapmetsc.sif Rscript run_analysis.R
```

```bash
# render a single script (e.g. the spatial analysis)
apptainer exec --bind /path/to/singlecelldata:/mnt/data mapmetsc.sif \
    Rscript -e 'rmarkdown::render("analysis/07_spatial_analysis.Rmd",
                  output_dir="docs",
                  params=list(input="/mnt/data/input_all",
                              output="/mnt/data/output"))'
```

```bash
# drop into an interactive R session
apptainer run --bind /path/to/singlecelldata:/mnt/data mapmetsc.sif
```

### RStudio Server (optional)

The base image ships RStudio Server. Apptainer containers are read-only, so
redirect the server's writable directories to bind mounts:

```bash
mkdir -p rstudio-run/{tmp,var,run}
PASSWORD=mapmetsc apptainer exec \
    --bind rstudio-run/var:/var/lib/rstudio-server \
    --bind rstudio-run/run:/var/run/rstudio-server \
    --bind rstudio-run/tmp:/tmp \
    --bind /path/to/MapMetSC:/home/rstudio/MapMetSC \
    --bind /path/to/singlecelldata:/mnt/data \
    mapmetsc.sif \
    /usr/lib/rstudio-server/bin/rserver --www-port 8787 \
        --auth-none 0 --auth-pam-helper-path=pam-helper --server-user=$(whoami)
# then open http://localhost:8787  (user: your shell user, password: mapmetsc)
```

## Notes

- `07_spatial_analysis.Rmd` also expects a **tumour-only cellular-community**
  column (see the note at the top of that script and in `run_analysis.R`). That
  is a data/analysis prerequisite, independent of the container — this image
  provides the software (fixed `imcRtools` + all packages) needed to run it.
- Built `.sif` files are git-ignored (`*.sif`).
- The image records provenance at `/opt/mapmetsc-sessionInfo.txt` and
  `/opt/mapmetsc-packages.csv`.
