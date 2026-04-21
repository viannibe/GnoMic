# Struttura completa di GnoMic

## Directory tree

```
GnoMic/
├── config.yaml                    ← Tutti i percorsi e parametri centralizzati
├── DATA_PROVENANCE.md             ← Dove vengono i file di input (alto livello)
├── PREPROCESSING_GUIDE.md         ← COMANDI ESATTI per rigenerare tutto da zero
├── STRUCTURE.md                   ← Questo file (overview struttura)
├── README.md                      ← Come eseguire la pipeline (entrypoint)
│
├── scripts/                       ← Pipeline principale (01–04)
│   ├── 01_build_parquet.py        → gnomad_corrected.parquet
│   ├── 02_kde_region.py           → A_kde_region_all.png + CSV
│   ├── 03_cpg_split.py            → A_exon_cpg_split + B_paired + C_consistency PNGs + CSVs
│   └── 04_driver_passenger.py     → A_driver_passenger_kde.png + CSV
│
├── notebooks/
│   └── 00_dataset_qc.ipynb        ← Notebook Jupyter: VQSLOD, CpG, dataset stats
│                                    (apri con: jupyter lab notebooks/00_dataset_qc.ipynb)
│
├── slurm/                         ← Job submission scripts
│   ├── 01_build_parquet.slurm     (8h, 192 GB)
│   ├── 02_kde_region.slurm        (2h, 64 GB)
│   ├── 03_cpg_split.slurm         (1h, 64 GB)
│   ├── 04_driver_passenger.slurm  (1h, 64 GB)
│   └── submit_all.sh              ← Sottometti tutto con dipendenze SLURM
│
├── results/                       ← Output (generati dai scripts)
│   ├── cache/
│   │   └── gnomad_corrected.parquet  (generato da step 01, ~20-40 GB)
│   ├── kde_region/
│   │   ├── A_kde_region_all.png      (μ_COSMIC=+0.246, μ_gnomAD=+0.042)
│   │   └── A_kde_region_all.csv      (statistiche MWU)
│   ├── cpg_split/
│   │   ├── A_exon_cpg_split.png
│   │   ├── A_exon_cpg_split_summary.csv
│   │   ├── B_exon_noncpg_paired.png  (Δμ=+0.178, p=3e-28, d=18.3)
│   │   ├── B_exon_noncpg_paired.csv
│   │   ├── C_chrom_consistency.png   (tutti 22 chroms positivi)
│   └── driver_passenger/
│       ├── A_driver_passenger_kde.png (d_driver=8.4, d_passenger=18.4)
│       └── A_driver_passenger_summary.csv
│
├── logs/                          ← Log SLURM (creati da submit_all.sh)
│   ├── 01_build_parquet_[jobid].out
│   ├── 02_kde_region_[jobid].out
│   ├── 03_cpg_split_[jobid].out
│   └── 04_driver_passenger_[jobid].out
│
└── [questi file di documentazione]
    ├── config.yaml                # Centralizza tutti i path
    ├── DATA_PROVENANCE.md         # Donde provengo i dati (overview)
    ├── PREPROCESSING_GUIDE.md     # Comandi esatti per rigenerare (step-by-step)
    ├── STRUCTURE.md               # Questo file
    └── README.md                  # Istruzioni uso
```

---

## Flusso di esecuzione

### Opzione A: Batch SLURM (raccomandato)

```bash
cd slurm/
bash submit_all.sh

# Monitoraggio:
squeue -u $USER
```

Sequenza di esecuzione (con dipendenze SLURM):
1. **Step 01** (8h): build_parquet → `gnomad_corrected.parquet`
2. **Step 02/03/04** (paralleli, dopo step 01):
   - Step 02 (2h): kde_region
   - Step 03 (1h): cpg_split
   - Step 04 (1h): driver_passenger

### Opzione B: Manuale

```bash
conda activate gnomad_venv
cd /leonardo/home/userexternal/viannibe/GnoMic

# Step 01 (lento, ~2-8h)
python scripts/01_build_parquet.py

# Step 02/03/04 (veloci, ~1h ciascuno, possono partire in parallelo)
python scripts/02_kde_region.py &
python scripts/03_cpg_split.py &
python scripts/04_driver_passenger.py &
wait
```

---

## Input file (read-only, su scratch)

| File | Path | Generato da | Usa in |
|------|------|-------------|--------|
| `merged_annotated.csv.gz` | `/leonardo_work/CNHPC_2116672/gnomAD/processed/annotation_2/` | gnomAD preprocessing (vedi PREPROCESSING_GUIDE) | step 01 |
| `cosmic_full_annotato.csv` | `/leonardo_work/CNHPC_2116672/COSMIC_v103/processed_ws2/` | COSMIC preprocessing (vedi PREPROCESSING_GUIDE) | step 02 |
| `cosmic_final.tsv` | `/leonardo_work/CNHPC_2116672/COSMIC_v103/processed/` | COSMIC preprocessing (vedi PREPROCESSING_GUIDE) | step 03, 04 |

---

## Output file (CSV + PNG)

| File | Step | Descrizione | Righe/Comandi |
|------|------|-------------|---------------|
| `A_kde_region_all.png` | 02 | KDE ΔE per regione, COSMIC vs gnomAD | Region, mean_cosmic, mean_gnomad, delta_mu, p_mwu |
| `A_exon_cpg_split.png` | 03A | KDE esoni, CpG vs non-CpG | Dataset, cpg_class, mean_dE, n_variants, sem_dE |
| `B_exon_noncpg_paired.png` | 03B | Scatter + strip per chrom, non-CpG | Chrom, gnomad_mean_dE, cosmic_mean_dE, diff, n_* |
| `C_chrom_consistency.png` | 03C | Strip plot standalone, 22 chroms positivi | [stessi dati di B] |
| `A_driver_passenger_kde.png` | 04 | KDE driver/passenger/noncoding | Group, dataset, mean_dE, p, cohens_d, t |

---

## Metodo (matematica)

### Energia mutazionale

```
E_ref = −ln(freq_7mer_reference_in_hg38)
E_mut = −ln(freq_7mer_mutant_in_hg38)
ΔE = E_ref − E_mut

Interpretazione:
  ΔE > 0  →  contesto WT è più raro  →  mutazione destabilizzante
  ΔE < 0  →  contesto MUT è più raro  →  mutazione stabilizzante
```

### Aggregazione per cromosoma

Non si usa il raw ΔE di ogni variante.
Invece:
1. Per ogni cromosoma, calcola mean(ΔE) di tutte le varianti
2. Ottieni 22 valori (uno per cromosoma)
3. Fai KDE su questi 22 valori

Vantaggi:
- Riduce rumore
- Highlights sistematica differenza tra COSMIC e gnomAD
- Permite paired statistics (paired t-test)

### Test statistici

- **Mann-Whitney U**: confronto distribuzioni KDE (COSMIC vs gnomAD)
- **Paired t-test**: differenza per-chrom (COSMIC − gnomAD)
- **Cohen's d**: effect size (d = mean_diff / std_diff)

---

## Configurazione

Vedi `config.yaml` per:
- Path assoluti di tutti i file
- Parametri di analisi (chunksizes, DPI, bandwidth)
- Mappature funzionali (VEP category → region)

---

## Notebook interattivo

```bash
conda activate gnomad_venv
jupyter lab notebooks/00_dataset_qc.ipynb
```

Celle:
1. Dataset overview (n righe, n chroms, dist maf_bin)
2. AS_VQSLOD (distribuzione, effetto filtro)
3. CpG statistics (frazione per regione)
4. ΔE overview (KDE per regione)
5. COSMIC overview (n varianti, driver stats)
6. Tabella riepilogativa (salvata in CSV)

---

## Riproducibilità

### Scenario 1: File input già presenti

Se i tre file esisono già:
```
✓ /leonardo_work/CNHPC_2116672/gnomAD/processed/annotation_2/merged_annotated.csv.gz
✓ /leonardo_work/CNHPC_2116672/COSMIC_v103/processed_ws2/cosmic_full_annotato.csv
✓ /leonardo_work/CNHPC_2116672/COSMIC_v103/processed/cosmic_final.tsv
```

Vai diretto a `cd slurm && bash submit_all.sh`

Tempo totale: ~12h (step 01 è il collo di bottiglia)

### Scenario 2: Rigenerare da zero

Se i file **non existono**, segui `PREPROCESSING_GUIDE.md`:

1. Scarica COSMIC v103 + gnomAD v4.1 raw VCF
2. Esegui step di preprocessing (VCF → clean → 7mer → energies → labels)
3. Quindi esegui questa pipeline

Tempo totale: ~1-2 giorni (dipende dal networking e risorse compute)

---

## Contact / Issues

Se i file di input sono corrotti o mancanti:
1. Verifica i path in `config.yaml`
2. Leggi `DATA_PROVENANCE.md` e `PREPROCESSING_GUIDE.md`
3. Contatta: vedi README.md
