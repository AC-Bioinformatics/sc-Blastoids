#!/usr/bin/env bash
# Runs the sc-Blastoids analysis step by step, from the raw 10X matrices to the trajectories
#
# Steps (each reads the previous step's output; results in 3.Results/):
#   cellbender     cellbender.sh                  -> 2.CellBender/
#   h5_to_seurat   h5_to_seurat.R                 -> 2.CellBender/merged_filtered_seurat.rds
#   qc             QC.R                           -> 3.QC_filtering/
#   normalisation  Normalisation.R                -> 5.Normalisation/
#   clustering     Clustering.R                   -> 6.Clustering/Louvain/
#   annotation     Annotation_SingleR.R           -> 7.Annotation/singler_zhao2024/
#   clustassess    ClustAssess.R, ClustAssess_postprocess.R -> 6.Clustering/ClustAssess/
#   flufftail      Flufftail.R                    -> 8.Pseudotime/flufftail/res0.3/
#   slingshot      Slingshot.R                    -> 8.Pseudotime/slingshot/
#
# R packages: Rscript 2.Scripts/install_packages.R; CellBender: see cellbender.sh / cellbender_container.sh
#
# Input data (see 1.Data/README.md):
#   1.Data/<Sample>-raw_feature_bc_matrix/
#   1.Data/references/zhao2024_human_embryo/psd.R3.6.em.seurat.ob.rds, zhao2024_celltype_markers_top25.csv
#
# Usage:  bash 2.Scripts/run_pipeline.sh                   # all steps
#         bash 2.Scripts/run_pipeline.sh --from qc         # from one step to the end
#         N_CORES=8 bash 2.Scripts/run_pipeline.sh         # workers for SingleR and ClustAssess (default 4)
#         CELLBENDER=apptainer bash 2.Scripts/run_pipeline.sh   # CellBender in the container (apptainer or docker)
#         SC_RESULTS=test bash 2.Scripts/run_pipeline.sh --from qc   # write a separate run to test/ instead of 3.Results
#         SC_RECLUSTER=true bash 2.Scripts/run_pipeline.sh       # recompute the Louvain clusters (default: reuse
#                                                                # 3.Results/6.Clustering/Louvain/cluster_metadata.rds)
# Logs:   <results folder>/logs/<step>.log

set -euo pipefail

Root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$Root"
# Results folder: 3.Results, or SC_RESULTS (relative to the repository root unless absolute)
Results="${SC_RESULTS:-3.Results}"
[[ "$Results" == /* ]] || Results="${Root}/${Results}"
export SC_RESULTS="$Results"
Log_dir="${Results}/logs"
mkdir -p "$Log_dir"

export N_CORES="${N_CORES:-4}"
CELLBENDER="${CELLBENDER:-local}"

Steps=(cellbender h5_to_seurat qc normalisation clustering annotation clustassess flufftail slingshot)

From="${Steps[0]}"
if [[ "${1:-}" == "--from" ]]; then
    From="${2:?--from needs a step name}"
fi
if [[ ! " ${Steps[*]} " == *" ${From} "* ]]; then
    echo "ERROR: unknown step '${From}'; steps: ${Steps[*]}" >&2
    exit 1
fi

run_step() {
    case "$1" in
        cellbender)
            case "$CELLBENDER" in
                local)     bash 2.Scripts/cellbender.sh ;;
                apptainer) bash 2.Scripts/cellbender_container.sh ;;
                docker)    bash 2.Scripts/cellbender_container.sh --docker ;;
                *) echo "ERROR: CELLBENDER must be local, apptainer or docker" >&2; return 1 ;;
            esac ;;
        h5_to_seurat)  Rscript 2.Scripts/h5_to_seurat.R ;;
        qc)            Rscript 2.Scripts/QC.R ;;
        normalisation) Rscript 2.Scripts/Normalisation.R ;;
        clustering)    Rscript 2.Scripts/Clustering.R ;;
        annotation)    Rscript 2.Scripts/Annotation_SingleR.R ;;
        clustassess)   Rscript 2.Scripts/ClustAssess.R && Rscript 2.Scripts/ClustAssess_postprocess.R ;;
        flufftail)     Rscript 2.Scripts/Flufftail.R ;;
        slingshot)     Rscript 2.Scripts/Slingshot.R ;;
    esac
}

started=false
for step in "${Steps[@]}"; do
    [[ "$step" == "$From" ]] && started=true
    $started || continue
    echo "$(date '+%F %H:%M:%S')  ${step}: started (log: ${Log_dir}/${step}.log)"
    if ! run_step "$step" > "${Log_dir}/${step}.log" 2>&1; then
        echo "ERROR: ${step} failed; see ${Log_dir}/${step}.log" >&2
        tail -20 "${Log_dir}/${step}.log" >&2
        exit 1
    fi
    echo "$(date '+%F %H:%M:%S')  ${step}: done"
done

echo "All steps finished"
