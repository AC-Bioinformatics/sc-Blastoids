#!/usr/bin/env bash
# CellBender remove-background on each sample's raw 10X matrix
#
# Input:  1.Data/<Sample>-raw_feature_bc_matrix/ (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)
# Steps:  Run CellBender 0.4.0 with the same settings on every sample found in 1.Data/
#         (GPU is used when PyTorch can see one, otherwise CPU)
# Output: 3.Results/2.CellBender/
#           <Sample>_cellbender.h5, <Sample>_cellbender_filtered.h5, <Sample>_cellbender_posterior.h5
#           <Sample>_cellbender.pdf, <Sample>_cellbender_report.html
#           <Sample>_cellbender_metrics.csv, <Sample>_cellbender_cell_barcodes.csv
#           <Sample>_cellbender.log
#           ckpt/<Sample>/ckpt.tar.gz (training checkpoint)
#
# Requires: CellBender 0.4.0 on PATH (pip install cellbender==0.4.0, e.g. in a conda environment)
#           To run in the Apptainer or Docker container instead, use cellbender_container.sh
# Usage:    bash 2.Scripts/cellbender.sh                  # all samples
#           bash 2.Scripts/cellbender.sh 24h-blastoids    # selected samples

set -euo pipefail

Root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
Data_dir="${Root}/1.Data"
# Results folder: 3.Results, or SC_RESULTS (relative to the repository root unless absolute)
Results="${SC_RESULTS:-3.Results}"
[[ "$Results" == /* ]] || Results="${Root}/${Results}"
Out_dir="${Results}/2.CellBender"

# CellBender settings (identical for all samples)
Total_droplets=65000
Fpr=0.01
Epochs=150
Checkpoint_mins=30

if ! command -v cellbender > /dev/null; then
    echo "ERROR: cellbender not found on PATH (pip install cellbender==0.4.0)" >&2
    exit 1
fi

if python -c "import torch, sys; sys.exit(0 if torch.cuda.is_available() else 1)" 2> /dev/null; then
    Device_flag="--cuda"
    echo "Device: GPU ($(python -c "import torch; print(torch.cuda.get_device_name(0))"))"
else
    Device_flag=""
    echo "Device: CPU"
fi

# Samples from the command line, otherwise every raw matrix in 1.Data/
if [[ $# -gt 0 ]]; then
    Samples=("$@")
else
    Samples=()
    for dir in "$Data_dir"/*-raw_feature_bc_matrix; do
        [[ -d "$dir" ]] || continue
        name="$(basename "$dir")"
        Samples+=("${name%-raw_feature_bc_matrix}")
    done
fi

if [[ ${#Samples[@]} -eq 0 ]]; then
    echo "ERROR: no <Sample>-raw_feature_bc_matrix directories in ${Data_dir}" >&2
    exit 1
fi

mkdir -p "$Out_dir"
echo "Samples: ${Samples[*]}"
echo "Params:  total_droplets=${Total_droplets}, fpr=${Fpr}, epochs=${Epochs}"

for sample in "${Samples[@]}"; do
    input="${Data_dir}/${sample}-raw_feature_bc_matrix"
    output="${Out_dir}/${sample}_cellbender.h5"

    if [[ ! -d "$input" ]]; then
        echo "ERROR: ${input} not found" >&2
        exit 1
    fi

    # Skip samples that already have a filtered output
    if [[ -f "${Out_dir}/${sample}_cellbender_filtered.h5" ]]; then
        echo "Skipping ${sample}: output already exists"
        continue
    fi

    # CellBender writes its checkpoint to the working directory; one per sample
    # lets an interrupted run resume from it
    ckpt_dir="${Out_dir}/ckpt/${sample}"
    mkdir -p "$ckpt_dir"

    echo "$(date '+%H:%M:%S')  ${sample}: started"
    (
        cd "$ckpt_dir"
        cellbender remove-background \
            --input "$input" \
            --output "$output" \
            --total-droplets-included "$Total_droplets" \
            --fpr "$Fpr" \
            --epochs "$Epochs" \
            --checkpoint-mins "$Checkpoint_mins" \
            $Device_flag
    )

    if [[ ! -f "${Out_dir}/${sample}_cellbender_filtered.h5" ]]; then
        echo "ERROR: ${sample}_cellbender_filtered.h5 was not produced; see ${sample}_cellbender.log" >&2
        exit 1
    fi
    echo "$(date '+%H:%M:%S')  ${sample}: done"
done

echo "All samples processed"
