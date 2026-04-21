# gnomAD-COSMIC Clean Pipeline

Pipeline pulita e riproducibile per l'analisi comparativa dell'energia
mutazionale (ΔE) tra varianti somatiche COSMIC e varianti germinali gnomAD.

---

## Struttura

```
GnoMic/
├── config.yaml                        # Tutti i percorsi e parametri
├── scripts/
│   ├── 01_build_parquet.py            # gnomAD CSV → parquet cache
│   ├── 02_kde_region.py               # → A_kde_region_all.png
│   ├── 03_cpg_split.py                # → A_exon_cpg_split + B_paired + C_consistency
│   └── 04_driver_passenger.py         # → A_driver_passenger_kde.png
├── notebooks/
│   └── 00_dataset_qc.ipynb            # Esplorazione: VQSLOD, CpG, stats
├── results/
│   ├── cache/                         # gnomad_corrected.parquet (generato da step 01)
│   ├── kde_region/                    # Output step 02
│   ├── cpg_split/                     # Output step 03
│   └── driver_passenger/              # Output step 04
├── slurm/
│   ├── 01_build_parquet.slurm
│   ├── 02_kde_region.slurm
│   ├── 03_cpg_split.slurm
│   ├── 04_driver_passenger.slurm
│   └── submit_all.sh                  # Sottomette tutta la pipeline con dipendenze
└── logs/                              # Log SLURM
```

---

## Input datasets

| File | Path | Usato da |
|------|------|----------|
| gnomAD annotated CSV | `/leonardo_work/CNHPC_2116672/gnomAD/processed/annotation_2/merged_annotated.csv.gz` | step 01 |
| COSMIC full (con ANN_TOP_ANNOTATION) | `/leonardo_work/CNHPC_2116672/COSMIC_v103/processed_ws2/cosmic_full_annotato.csv` | step 02 |
| COSMIC final (con E_wt7, in_cgc) | `/leonardo_work/CNHPC_2116672/COSMIC_v103/processed/cosmic_final.tsv` | step 03, 04 |

---

## Come eseguire

### Opzione A — SLURM (raccomandato su Leonardo)

```bash
cd /leonardo/home/userexternal/viannibe/GnoMic/slurm
bash submit_all.sh
```

Questo sottomette 4 job con la dipendenza corretta:
- Step 01 parte subito (8h, 192 GB)
- Step 02/03/04 partono in parallelo solo dopo che Step 01 termina con successo

Monitoraggio:
```bash
squeue -u $USER
# Logs in: logs/
```

### Opzione B — Esecuzione manuale sequenziale

```bash
conda activate gnomad_venv
cd /leonardo/home/userexternal/viannibe/GnoMic

# Step 01: build parquet (lento, ~2–4h)
python scripts/01_build_parquet.py

# Step 02/03/04: veloci (~30–60 min ciascuno), si possono eseguire in parallelo
python scripts/02_kde_region.py
python scripts/03_cpg_split.py
python scripts/04_driver_passenger.py
```

### Rigenerare solo un subset

```bash
# Rigenera solo il parquet (es. se l'annotated CSV è cambiato)
python scripts/01_build_parquet.py --force

# Cambia output dir
python scripts/02_kde_region.py -o /path/custom/output/
```

---

## Output prodotti

| Figura | Script | Descrizione |
|--------|--------|-------------|
| `A_kde_region_all.png` | 02 | KDE ΔE per regione (exon/intron/extragenic), COSMIC vs gnomAD. µ(exon): COSMIC +0.246, gnomAD +0.042 |
| `A_kde_region_all.csv` | 02 | Statistiche: mean, delta_mu, p-value MWU per regione |
| `A_exon_cpg_split.png` | 03 | KDE esoni, CpG vs non-CpG separati |
| `A_exon_cpg_split_summary.csv` | 03 | mean_dE, sem, n per (dataset × CpG_class) |
| `B_exon_noncpg_paired.png` | 03 | Scatter + strip plot per-chrom, non-CpG esoni. Δµ=+0.178, p=3e-28, d=18.3 |
| `B_exon_noncpg_paired.csv` | 03 | Tabella per cromosoma: gnomad_mean_dE, cosmic_mean_dE, diff |
| `C_chrom_consistency.png` | 03 | Strip plot standalone: tutti 22 autosomi con diff positiva (+0.163–+0.200) |
| `A_driver_passenger_kde.png` | 04 | KDE driver/passenger/noncoding vs gnomAD. d(driver)=8.4, d(passenger)=18.4 |
| `A_driver_passenger_summary.csv` | 04 | Statistiche per gruppo: t, p, Cohen's d |

---

## Metodo

**ΔE = E_ref − E_mut = −ln(WT_7mer_freq) + ln(MUT_7mer_freq)**

- Positivo: il contesto WT è più raro del contesto mutato → mutazione "destabilizzante"
- Aggregazione per cromosoma (n=22 medie autosomiche), poi KDE su queste medie
- Test statistici: Mann-Whitney U (confronto distribuzioni), paired t-test (per-chromosome), Cohen's d (effect size)
- CpG identificato dal 7-mer centrale (posizione 4=C, posizione 5=G, o complementare)

**Coerenza con gnomAD_vi_clean:**
I calcoli sono identici agli script `09_gnomad_cosmic_comparison.py`, `13c_combined_AB.py` e `14_driver_passenger_control.py`. Il parquet cache (`gnomad_corrected.parquet`) è interscambiabile con quello in `gnomAD_vi_clean/results/09_cosmic_vs_gnomad/`.

---

## Notebook

`notebooks/00_dataset_qc.ipynb` — apri con JupyterLab su Leonardo:

```bash
conda activate gnomad_venv
jupyter lab notebooks/00_dataset_qc.ipynb
```

Contenuto:
1. Overview dataset gnomAD (n varianti, distribuzione per cromosoma, maf_bin)
2. AS_VQSLOD: distribuzione e effetto del filtro su ΔE
3. CpG fraction: confronto gnomAD vs COSMIC per regione
4. ΔE overview: distribuzione per regione e CpG class
5. COSMIC overview: n varianti, distribuzione coding/noncoding, driver stats
6. Tabella riassuntiva esportata in `results/dataset_qc_summary.csv`

---

## Ambiente

```bash
conda activate gnomad_venv   # già configurato su Leonardo
# Dipendenze: numpy, pandas, matplotlib, seaborn, scipy, pyarrow
```
