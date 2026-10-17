#!/usr/bin/env bash
# Runs cellbender.sh inside the CellBender 0.4.0 container (Apptainer by default, or Docker)
#
# Input, steps, output and settings: see cellbender.sh
#
# Requires one of:
#   Apptainer (default): container image at containers/cellbender-0.4.0.sif, created with
#     apptainer pull containers/cellbender-0.4.0.sif docker://us.gcr.io/broad-dsde-methods/cellbender:v0.4.0
#     (another location can be given with CELLBENDER_SIF=/path/to/image.sif)
#   --docker: Docker (Docker Desktop on macOS / Windows); the image is pulled on first use
#     (the image is linux/amd64 only, so it runs emulated on Apple Silicon)
#
# Usage:  bash 2.Scripts/cellbender_container.sh                          # all samples, Apptainer
#         bash 2.Scripts/cellbender_container.sh 24h-blastoids            # selected samples, Apptainer
#         bash 2.Scripts/cellbender_container.sh --docker                 # all samples, Docker
#         bash 2.Scripts/cellbender_container.sh --docker 24h-blastoids   # selected samples, Docker

set -euo pipefail

Root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
Script="${Root}/2.Scripts/cellbender.sh"
Image="us.gcr.io/broad-dsde-methods/cellbender:v0.4.0"
Container="${CELLBENDER_SIF:-${Root}/containers/cellbender-0.4.0.sif}"

Engine=apptainer
if [[ "${1:-}" == "--docker" ]]; then
    Engine=docker
    shift
fi

if [[ "$Engine" == "apptainer" ]]; then
    if ! command -v apptainer > /dev/null; then
        echo "ERROR: apptainer not found; install Apptainer or rerun with --docker" >&2
        exit 1
    fi
    if [[ ! -f "$Container" ]]; then
        echo "ERROR: container not found at ${Container}; create it with" >&2
        echo "  apptainer pull ${Container} docker://${Image}" >&2
        exit 1
    fi
    Gpu_flag=""
    command -v nvidia-smi > /dev/null && Gpu_flag="--nv"
    echo "Engine: Apptainer (${Container})"
    apptainer exec $Gpu_flag --bind "$Root" "$Container" bash "$Script" "$@"
else
    if ! command -v docker > /dev/null; then
        echo "ERROR: docker not found; install Docker or rerun without --docker" >&2
        exit 1
    fi
    if ! docker info > /dev/null 2>&1; then
        echo "ERROR: the Docker daemon is not running (start Docker Desktop)" >&2
        exit 1
    fi
    if ! docker image inspect "$Image" > /dev/null 2>&1; then
        docker pull --platform linux/amd64 "$Image"
    fi
    Gpu_flag=""
    command -v nvidia-smi > /dev/null && Gpu_flag="--gpus all"
    echo "Engine: Docker (${Image})"
    # Runs as the host user so outputs are not owned by root; HOME must be writable for caches
    docker run --rm --platform linux/amd64 $Gpu_flag \
        -u "$(id -u):$(id -g)" -e HOME=/tmp -e SC_RESULTS \
        -v "${Root}:${Root}" -w "$Root" \
        "$Image" bash "$Script" "$@"
fi
