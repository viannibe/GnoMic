#!/bin/bash
#SBATCH --job-name=gco_01_parquet
#SBATCH --account=CNHPC_2116672
#SBATCH --partition=boost_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=192G
#SBATCH --time=08:00:00
#SBATCH --output=../logs/01_build_parquet_%j.out
#SBATCH --error=../logs/01_build_parquet_%j.err

# ============================================================
# Step 01 — Build gnomAD parquet from merged_annotated.csv.gz
# DEVE completare prima di eseguire gli step 02/03/04
# ============================================================

set -euo pipefail

PROJECT_DIR="/leonardo/home/userexternal/viannibe/GnoMic"
mkdir -p "${PROJECT_DIR}/logs"

MINI="$HOME/miniconda3"
source "${MINI}/etc/profile.d/conda.sh"
conda activate gnomad_venv

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MPLBACKEND=Agg

echo "========== STEP 01: Build gnomAD parquet =========="
echo "Job ID  : ${SLURM_JOB_ID}"
echo "Start   : $(date)"
echo "Node    : $(hostname)"
echo "==================================================="

ANNOTATED_CSV="/leonardo_work/CNHPC_2116672/gnomAD/processed/annotation_2/merged_annotated.csv.gz"
OUTPUT_PARQ="${PROJECT_DIR}/results/cache/gnomad_corrected.parquet"

echo "[INFO] Input  : ${ANNOTATED_CSV}"
echo "[INFO] Output : ${OUTPUT_PARQ}"
echo ""

mkdir -p "${PROJECT_DIR}/results/cache"

/usr/bin/time -v python3 "${PROJECT_DIR}/scripts/01_build_parquet.py" \
    --annotated-csv "${ANNOTATED_CSV}" \
    --output        "${OUTPUT_PARQ}"

EXITCODE=$?

echo ""
echo "==================================================="
if [ ${EXITCODE} -eq 0 ]; then
    echo "SUCCESS | End: $(date)"
    ls -lh "${OUTPUT_PARQ}"
else
    echo "FAILED (exit code ${EXITCODE}) | End: $(date)"
fi
echo "==================================================="
exit ${EXITCODE}
