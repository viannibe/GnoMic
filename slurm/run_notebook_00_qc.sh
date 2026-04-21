#!/bin/bash
#SBATCH --job-name=gco_00_nbqc
#SBATCH --account=CNHPC_2116672
#SBATCH --partition=boost_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/leonardo/home/userexternal/viannibe/GnoMic/logs/00_dataset_qc_%j.out
#SBATCH --error=/leonardo/home/userexternal/viannibe/GnoMic/logs/00_dataset_qc_%j.err

set -euo pipefail

PROJECT_DIR="/leonardo/home/userexternal/viannibe/GnoMic"
NOTEBOOK="${PROJECT_DIR}/notebooks/00_dataset_qc.ipynb"
OUTDIR="${PROJECT_DIR}/results/notebook_runs"
OUTNOTEBOOK="${OUTDIR}/00_dataset_qc_executed_${SLURM_JOB_ID}.ipynb"

mkdir -p "${PROJECT_DIR}/logs"
mkdir -p "${OUTDIR}"

MINI="$HOME/miniconda3"
source "${MINI}/etc/profile.d/conda.sh"
conda activate gnomad_venv

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MPLBACKEND=Agg

echo "========== RUN NOTEBOOK: 00_dataset_qc =========="
echo "Job ID   : ${SLURM_JOB_ID}"
echo "Start    : $(date)"
echo "Notebook : ${NOTEBOOK}"
echo "Output   : ${OUTNOTEBOOK}"
echo "================================================"

if [ ! -f "${NOTEBOOK}" ]; then
    echo "ERROR: notebook non trovato: ${NOTEBOOK}"
    exit 1
fi

cd "${PROJECT_DIR}"

# Esegue il notebook e salva un notebook eseguito con output
/usr/bin/time -v jupyter nbconvert \
    --to notebook \
    --execute "${NOTEBOOK}" \
    --output "${OUTNOTEBOOK}" \
    --ExecutePreprocessor.timeout=-1 \
    --ExecutePreprocessor.kernel_name=python3

EXITCODE=$?

echo ""
echo "================================================"
if [ ${EXITCODE} -eq 0 ]; then
    echo "SUCCESS | End: $(date)"
    ls -lh "${OUTDIR}"
else
    echo "FAILED (exit code ${EXITCODE}) | End: $(date)"
fi
echo "================================================"

exit ${EXITCODE}