#!/usr/bin/env bash
# =============================================================================
# Build the MapMetSC Apptainer image.
#
#   ./apptainer/build.sh                 # -> mapmetsc.sif in repo root
#   ./apptainer/build.sh my.sif          # custom output path
#
# Requires Apptainer (>= 1.1) or SingularityCE. Building needs to pull the
# rocker/rstudio:4.4.0 base from Docker Hub and download CRAN/Bioconductor
# packages, so run it somewhere with outbound network access. `apptainer build`
# needs root or `--fakeroot` (configured on most HPC systems).
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEF="${HERE}/mapmetsc.def"
OUT="${1:-${HERE}/../mapmetsc.sif}"

# Pick apptainer or singularity, whichever is on PATH.
if command -v apptainer >/dev/null 2>&1; then
    RUNTIME=apptainer
elif command -v singularity >/dev/null 2>&1; then
    RUNTIME=singularity
else
    echo "ERROR: neither 'apptainer' nor 'singularity' found on PATH." >&2
    echo "Install Apptainer: https://apptainer.org/docs/admin/main/installation.html" >&2
    exit 1
fi

# Use --fakeroot when not root (the usual case on HPC login nodes).
FAKEROOT=""
if [ "$(id -u)" -ne 0 ]; then
    FAKEROOT="--fakeroot"
fi

echo ">> ${RUNTIME} build ${FAKEROOT} ${OUT} ${DEF}"
exec "${RUNTIME}" build ${FAKEROOT} "${OUT}" "${DEF}"
