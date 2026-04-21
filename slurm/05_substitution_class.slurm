#!/bin/bash
#SBATCH --job-name=gco_05_subst
#SBATCH --account=CNHPC_2116672
#SBATCH --partition=boost_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=../logs/05_substitution_class_%j.out
#SBATCH --error=../logs/05_substitution_class_%j.err

# ============================================================
# Step 05 — stratificazione per classe di sostituzione (SBS6)
# Dipende da Step 01 rigenerato con --force (parquet deve
# contenere colonne REF e ALT)
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

echo "========== STEP 05: substitution class analysis =========="
echo "Job ID  : ${SLURM_JOB_ID}"
echo "Start   : $(date)"
echo "==================================================="

PARQUET="${PROJECT_DIR}/results/cache/gnomad_corrected.parquet"
COSMIC_TSV="/leonardo_work/CNHPC_2116672/COSMIC_v103/processed/cosmic_final.tsv"
OUTDIR="${PROJECT_DIR}/results/substitution_class"

if [ ! -f "${PARQUET}" ]; then
    echo "ERROR: parquet non trovato: ${PARQUET}"
    echo "Esegui prima lo step 01."
    exit 1
fi

echo "[INFO] Parquet  : ${PARQUET}"
echo "[INFO] COSMIC   : ${COSMIC_TSV}"
echo "[INFO] Output   : ${OUTDIR}"
echo ""

mkdir -p "${OUTDIR}"

/usr/bin/time -v python3 "${PROJECT_DIR}/scripts/05_substitution_class.py" \
    --gnomad-parquet "${PARQUET}" \
    --cosmic-tsv     "${COSMIC_TSV}" \
    -o               "${OUTDIR}"

EXITCODE=$?

echo ""
echo "==================================================="
if [ ${EXITCODE} -eq 0 ]; then
    echo "SUCCESS | End: $(date)"
    ls -lh "${OUTDIR}/"
else
    echo "FAILED (exit code ${EXITCODE}) | End: $(date)"
fi
echo "==================================================="
exit ${EXITCODE}
