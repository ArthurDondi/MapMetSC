#!/bin/bash
#SBATCH --output logs/MapMetSC_%j.log
#SBATCH --error  logs/MapMetSC_%j.err
#SBATCH --job-name=MapMetSC
#SBATCH --partition=longq      # 30d limit: comfortably covers the whole workflow
#SBATCH --qos=longq            # qos must match the partition
#SBATCH --time=2-00:00:00      # --time=days-hours:minutes:seconds
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16     # SingleR steps fork MulticoreParam(workers=10..16)
#SBATCH --mem=192000           # parallel workers multiply memory (see NOTE below)
#SBATCH --mail-type=end
#SBATCH --mail-user=arthur.dondi@cemm.at
#
# Run the whole MapMetSC R analysis (Lazic et al.) on the SLURM cluster.
#
# Unlike the Snakemake template this is NOT a lightweight controller: the job
# runs the pipeline itself by executing `Rscript run_analysis.R` inside the
# Apptainer image. run_analysis.R renders, in order, the adrenal-gland (AD)
# reference -> main PT/BM pipeline -> the four correlation steps. It is one
# heavy compute job, hence the cpus/mem above (tune them to your data size) and
# the long walltime.
#
# Submit from the MapMetSC/ repo root (where run_analysis.R lives):
#   mkdir -p logs                       # once: SLURM needs the log dir to exist
#   sbatch run_mapmetsc_slurm.sh
#
# Paths can be overridden at submit time, e.g.:
#   MAPMET_DATA=/nobackup/.../mapmet MAPMET_SIF=/path/to/mapmetsc.sif \
#       sbatch run_mapmetsc_slurm.sh
#
# NOTE on parallelism / memory: the analysis hard-codes its worker counts -
# 04_phenotyping uses SnowParam(workers=30) for the UMAP, and the SingleR
# correlation steps use MulticoreParam(workers=10..16). With cpus-per-task=16
# the 30-worker step just oversubscribes (slower, not fatal), but every worker
# holds a copy of the data, so memory is the real risk - hence the generous
# --mem. If a step is OOM-killed, raise --mem (or --cpus-per-task toward 30);
# if it queues too long, lower both and accept the phenotyping step running
# slower. RNGseed is fixed in those calls, so results do not depend on the
# worker count.

set -euo pipefail

# Some clusters expose apptainer/singularity through modules; uncomment if so:
module load apptainer/1.1.9

# ----------------------------- CONFIG (edit me) -----------------------------
# The built image (see apptainer/build.sh). Defaults to the repo root.
SIF="${MAPMET_SIF:-$SLURM_SUBMIT_DIR/mapmetsc.sif}"

# The directory that contains Publication/ (your extracted data root).
DATA="${MAPMET_DATA:-/nobackup/lab_taschner-mandl/arthurdondi/projects/mapmet}"

# Derived data paths (match the Zenodo layout; edit if yours differs).
INPUT="${MAPMET_INPUT:-$DATA/Publication/20240811_Zenodo-Upload/MapMetIP_ProcessedDataset}"
PUBLIC="${MAPMET_PUBLIC:-$DATA/Publication/public_datasets}"
OUTPUT="${MAPMET_OUTPUT:-$DATA/Publication/output}"
# Frozen phenograph intermediaries (kept_cells_*.rds, pg_clusters_*.rds). In the
# Zenodo layout these live in an R_intermediary/ folder that sits NEXT TO
# MapMetIP_ProcessedDataset, i.e. one level above $INPUT. Override if yours
# differs; if you point it outside $DATA, add a matching --bind below.
INTERMEDIARY="${MAPMET_INTERMEDIARY:-$(dirname "$INPUT")/R_intermediary}"
# ----------------------------------------------------------------------------

echo "======================"
echo "submit dir : $SLURM_SUBMIT_DIR"
echo "job name   : $SLURM_JOB_NAME"
echo "partition  : $SLURM_JOB_PARTITION"
echo "job id     : $SLURM_JOB_ID"
echo "cpus       : ${SLURM_CPUS_PER_TASK:-1}   mem(MB): ${SLURM_MEM_PER_NODE:-?}"
echo "sif        : $SIF"
echo "input      : $INPUT"
echo "public     : $PUBLIC"
echo "output     : $OUTPUT"
echo "intermed.  : $INTERMEDIARY"
echo "started    : $(date)"
echo "======================"

# Sanity checks (fail early with a clear message rather than deep in R).
[[ -f "$SIF" ]]                 || { echo "ERROR: image not found: $SIF (build it with apptainer/build.sh)"; exit 1; }
[[ -d "$INPUT/regionprops" ]]   || { echo "ERROR: '$INPUT' is not a MapMetIP_ProcessedDataset (no regionprops/)"; exit 1; }
[[ -d "$PUBLIC" ]]              || { echo "ERROR: public datasets dir not found: $PUBLIC"; exit 1; }
[[ -d "$INTERMEDIARY" ]]       || echo "WARNING: intermediary dir not found: $INTERMEDIARY (02/04 will stop if the frozen kept_cells_*/pg_clusters_* lists are missing; generate them with code/generate_phenograph_intermediary.R)"
mkdir -p "$OUTPUT"

# Keep BLAS/OpenMP threads inside the SLURM allocation (APPTAINERENV_* is how
# apptainer forwards env vars into the container).
export APPTAINERENV_OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export APPTAINERENV_OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

# --cleanenv: do NOT inherit the host shell env (an active conda `base` sets
# LD_LIBRARY_PATH to conda libs, which can shadow the image's own libraries and
# break e.g. git2r's libgit2). --env / APPTAINERENV_* still pass through.
apptainer exec \
    --cleanenv \
    --bind "$DATA" \
    --bind "$SLURM_SUBMIT_DIR" \
    --pwd "$SLURM_SUBMIT_DIR" \
    --env MAPMET_INPUT="$INPUT" \
    --env MAPMET_PUBLIC="$PUBLIC" \
    --env MAPMET_OUTPUT="$OUTPUT" \
    --env MAPMET_INTERMEDIARY="$INTERMEDIARY" \
    "$SIF" \
    Rscript run_analysis.R

echo "======================"
echo "finished   : $(date)"
echo "objects & figures in: $OUTPUT"
echo "rendered HTML in    : $SLURM_SUBMIT_DIR/docs"
echo "======================"
