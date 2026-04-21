#!/bin/bash
# =============================================================
# submit_all.sh — Sottomette l'intera pipeline con dipendenze
#
# Uso:
#   cd /leonardo/home/userexternal/viannibe/GnoMic/slurm
#   bash submit_all.sh
#
# Flusso:
#   01_build_parquet  ──┬── 02_kde_region
#                       ├── 03_cpg_split
#                       └── 04_driver_passenger
#
# Step 02/03/04 partono solo dopo che 01 termina con successo.
# =============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"

echo "Project dir : ${PROJECT_DIR}"
echo "Script dir  : ${SCRIPT_DIR}"
echo ""

# --- Step 01: build parquet (lungo, 192 GB, 8h) ---
JOB01=$(sbatch --parsable "${SCRIPT_DIR}/01_build_parquet.slurm")
echo "Submitted Step 01 (build parquet) → Job ID: ${JOB01}"

# --- Step 02/03/04: partono solo dopo 01 (--dependency=afterok) ---
JOB02=$(sbatch --parsable --dependency=afterok:${JOB01} "${SCRIPT_DIR}/02_kde_region.slurm")
echo "Submitted Step 02 (kde_region)     → Job ID: ${JOB02}  [depends on ${JOB01}]"

JOB03=$(sbatch --parsable --dependency=afterok:${JOB01} "${SCRIPT_DIR}/03_cpg_split.slurm")
echo "Submitted Step 03 (cpg_split)      → Job ID: ${JOB03}  [depends on ${JOB01}]"

JOB04=$(sbatch --parsable --dependency=afterok:${JOB01} "${SCRIPT_DIR}/04_driver_passenger.slurm")
echo "Submitted Step 04 (driver_pass)    → Job ID: ${JOB04}  [depends on ${JOB01}]"

echo ""
echo "=== Pipeline queued ==="
echo "Monitor con:  squeue -u \$USER"
echo "Logs in     : ${PROJECT_DIR}/logs/"
echo ""
echo "Output attesi:"
echo "  results/cache/gnomad_corrected.parquet"
echo "  results/kde_region/A_kde_region_all.png"
echo "  results/kde_region/A_kde_region_all.csv"
echo "  results/cpg_split/A_exon_cpg_split.png"
echo "  results/cpg_split/A_exon_cpg_split_summary.csv"
echo "  results/cpg_split/B_exon_noncpg_paired.png"
echo "  results/cpg_split/B_exon_noncpg_paired.csv"
echo "  results/cpg_split/C_chrom_consistency.png"
echo "  results/driver_passenger/A_driver_passenger_kde.png"
echo "  results/driver_passenger/A_driver_passenger_summary.csv"
