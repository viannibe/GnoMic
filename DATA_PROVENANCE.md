# Data Provenance — Tracciabilità completa dei dataset

Questo file descrive come sono stati generati tutti i file di input usati dalla pipeline, 
partendo dai dati pubblici scaricati.

**Per i comandi esatti per rigenerare tutto da zero**, vedi: [`PREPROCESSING_GUIDE.md`](PREPROCESSING_GUIDE.md)

---

## Panoramica end-to-end

```
DOWNLOAD PUBBLICO
    gnomAD v4.1 VCF (GRCh38)           COSMIC v103 VCF (GRCh38)
           │                                     │
           ▼  PREPROCESSING                      ▼  PREPROCESSING
    merged_annotated.csv.gz             cosmic_full_annotato.csv
                                        cosmic_final.tsv
           │                                     │
           ▼  QUESTA PIPELINE (scripts/01–04)
    gnomad_corrected.parquet  ──────────────────►  Figure + CSV statistiche
```

---

## 1. gnomAD v4.1

### 1a. Sorgente raw

```
/leonardo_work/CNHPC_2116672/gnomAD/gnomad_genomes_4.1/
    gnomad.genomes.v4.1.sites.chrN.vcf.bgz   (N = 1..22)
    gnomad.genomes.v4.1.sites.chrN.vcf.bgz.tbi
    all_chromosomes.vcf.bgz    (merge di tutti i cromosomi)
```

Scaricati con:
```bash
# /leonardo_work/CNHPC_2116672/gnomAD/download_gnomAD.sh
# (account gnomAD, accesso HTTPS / AWS S3)
```

### 1b. Estrazione campi → gnomad_exonic.csv

Comando bcftools (vedere `download_gnomAD.sh`):
```bash
bcftools query -H \
  -f "%CHROM,%POS,%REF,%ALT,%ID,%FILTER,%AC,%AF,%AN,%allele_type,%AS_VQSLOD,\
$(for p in afr nfe eas mid sas grpmax; do echo -n "%AC_$p,%AF_$p,%AN_$p,"; done)%vep\n" \
  gnomad.exonic.vcf.gz > gnomad_exonic.csv
```

File prodotto:
```
/leonardo_work/CNHPC_2116672/gnomAD/gnomad_exonic.csv
```

Colonne chiave: `CHROM, POS, REF, ALT, FILTER, AF, AN, AS_VQSLOD, allele_type,
AF_afr/nfe/eas/mid/sas/grpmax, vep` (stringa VEP raw da parsare)

### 1c. Aggiunta energie 7-mer → merged_with_energies.csv.gz

Script: `kmer7.py` (nell'ambiente di lavoro originale)

Operazioni:
- Estrae il contesto 7-mer attorno a ogni variante usando `hg38.fa`
- Calcola `E_ref = −ln(WT_7mer_freq)` e `E_mut = −ln(MUT_7mer_freq)`
- Frequenze di background: `hg38_k7_freq_withPos.csv`
  (`/leonardo_work/CNHPC_2116672/gnomAD/hg38_k7_freq_withPos.csv`)
  (frequenza osservata di ogni k-mero 7 nel genoma hg38, senza cromosoma M)

File prodotto:
```
/leonardo_work/CNHPC_2116672/gnomAD/processed/merged_with_energies.csv.gz
```

Colonne aggiunte: `WT_7mer, MUT_7mer, E_ref, E_mut`

### 1d. Annotazione VEP + features → merged_annotated.csv.gz

Script: `gnomAD_vi_clean/scripts/06_gnomad_annotated.py`

SLURM: `gnomAD_vi_clean/run_06_annotated.slurm`

Operazioni:
- Parsifica la stringa VEP grezza (campo `vep`)
- Seleziona la conseguenza più grave per gene canonico (severity order)
- Aggiunge `vep_category` (forma corta: missense/synonymous/intron/ecc.)
- Calcola `is_CpG` dal 7-mer (posizione centrale: C→G o G←C)
- Calcola `maf_bin` (ultra_rare/very_rare/rare/low_freq/common_5/ecc.)
- Aggiunge `mutation_class` (transition/transversion), `is_coding`

File prodotto (INPUT di `01_build_parquet.py`):
```
/leonardo_work/CNHPC_2116672/gnomAD/processed/annotation_2/merged_annotated.csv.gz
```

Colonne chiave aggiunte:
`vep_category, vep_consequence, vep_impact, vep_symbol,
is_CpG, maf_bin, mutation_class, mutation_type, is_coding, functional_impact`

---

## 2. COSMIC v103

### 2a. Sorgente raw

```
/leonardo_work/CNHPC_2116672/COSMIC_v103/mutation/
    Cosmic_GenomeScreensMutant_v103_GRCh38.vcf
    README_Cosmic_GenomeScreensMutant_v103_GRCh38.txt
```

Scaricato dal portale COSMIC (account richiesto):
`https://cancer.sanger.ac.uk/cosmic/download`
→ "Genome Screens Mutant" → GRCh38 → v103

### 2b. Filtro SNV → cosmic_snv_clean.tsv

Filtra:
- Solo SNV (REF e ALT singolo nucleotide)
- Aggiunge `is_coding` dal campo `consequence` nel VCF INFO

File:
```
/leonardo_work/CNHPC_2116672/COSMIC_v103/processed/cosmic_snv_clean.tsv
```

Colonne: `chrom, pos, ref, alt, cosmic_id, gene, so_term, aa, hgvsp,
          sample_count, transcript, is_canonical, gene_tier, consequence, is_coding`

### 2c. Aggiunta contesto 7-mer → cosmic_snv_with_7mer.tsv

Estrae il 7-mer WT e MUT usando:
- Genoma di riferimento: `/leonardo_work/CNHPC_2116672/COSMIC_v103/ref/hg38.fa`

File:
```
/leonardo_work/CNHPC_2116672/COSMIC_v103/processed/cosmic_snv_with_7mer.tsv
```

Colonne aggiunte: `WT_7mer, MUT_7mer`

### 2d. Aggiunta frequenze background → cosmic_snv_with_background.tsv

Associa ogni 7-mer alla sua frequenza nel genoma hg38 usando:
```
/leonardo_work/CNHPC_2116672/COSMIC_v103/processed/hg38_k7_freq_noM.csv
```
(stesso file di background di gnomAD, senza cromosoma M)

File:
```
/leonardo_work/CNHPC_2116672/COSMIC_v103/processed/cosmic_snv_with_background.tsv
```

Colonne aggiunte: `p_wt7, p_mut7` (frequenze osservate)

### 2e. Aggiunta energie + CpG + CGC → cosmic_final.tsv

Calcola:
- `E_wt7 = −ln(p_wt7)`
- `E_mut7 = −ln(p_mut7)`
- `deltaE7 = E_wt7 − E_mut7`
- `is_cpg` (dal 7-mer centrale)
- `in_cgc` (1 se il gene è nel Cancer Gene Census)
  Fonte CGC: `census/` directory (lista ufficiale COSMIC dei driver genes)

File prodotto (INPUT di `03_cpg_split.py` e `04_driver_passenger.py`):
```
/leonardo_work/CNHPC_2116672/COSMIC_v103/processed/cosmic_final.tsv
```

### 2f. Aggiunta annotazione funzionale completa → cosmic_full_annotato.csv

Aggiunge `ANN_TOP_ANNOTATION` con la conseguenza funzionale completa
(formato SnpEff/VEP lungo: `missense_variant`, `intron_variant`, ecc.)
necessaria per mappare exon/intron/extragenic in `02_kde_region.py`.

File prodotto (INPUT di `02_kde_region.py`):
```
/leonardo_work/CNHPC_2116672/COSMIC_v103/processed_ws2/cosmic_full_annotato.csv
```

---

## 3. Tabelle di riferimento condivise

| File | Path | Descrizione |
|------|------|-------------|
| `hg38_k7_freq_withPos.csv` | `/leonardo_work/CNHPC_2116672/gnomAD/hg38_k7_freq_withPos.csv` | Frequenza 7-mer in hg38 per gnomAD |
| `hg38_k7_freq_noM.csv` | `/leonardo_work/CNHPC_2116672/COSMIC_v103/processed/hg38_k7_freq_noM.csv` | Stessa tabella (senza chrM) per COSMIC |
| `hg38.fa` | `/leonardo_work/CNHPC_2116672/gnomAD/hg38.fa` | Genoma di riferimento hg38 |
| `gencode.v46.exons.bed` | `/leonardo_work/CNHPC_2116672/gnomAD/gencode.v46.exons.bed` | Coordinate esoni (GENCODE v46) |
| CGC (Cancer Gene Census) | `/leonardo_work/CNHPC_2116672/COSMIC_v103/census/` | Lista driver genes COSMIC |

La tabella di frequenza 7-mer è il cuore del metodo:
```
E = −ln(frequenza del 7-mer nel genoma hg38)
ΔE = E_ref − E_mut  (positivo = mutazione verso contesto più raro)
```

---

## 4. Riepilogo percorsi file

```
# INPUT (read-only)
GNOMAD_ANNOTATED = /leonardo_work/CNHPC_2116672/gnomAD/processed/annotation_2/merged_annotated.csv.gz
COSMIC_FULL      = /leonardo_work/CNHPC_2116672/COSMIC_v103/processed_ws2/cosmic_full_annotato.csv
COSMIC_FINAL     = /leonardo_work/CNHPC_2116672/COSMIC_v103/processed/cosmic_final.tsv

# CACHE (generato da scripts/01_build_parquet.py)
GNOMAD_PARQUET   = results/cache/gnomad_corrected.parquet

# OUTPUT (generati da scripts/02–04)
results/kde_region/A_kde_region_all.png
results/cpg_split/A_exon_cpg_split.png
results/cpg_split/B_exon_noncpg_paired.png
results/cpg_split/C_chrom_consistency.png
results/driver_passenger/A_driver_passenger_kde.png
```

---

## 5. Nota su riproducibilità

Il parquet cache `gnomad_corrected.parquet` è identico a quello già prodotto in:
```
/leonardo/home/userexternal/viannibe/gnomAD_vi_clean/results/09_cosmic_vs_gnomad/gnomad_corrected.parquet
```
Se esiste già puoi evitare lo step 01 (8h) con un symlink:
```bash
ln -s /leonardo/home/userexternal/viannibe/gnomAD_vi_clean/results/09_cosmic_vs_gnomad/gnomad_corrected.parquet \
      /leonardo/home/userexternal/viannibe/GnoMic/results/cache/gnomad_corrected.parquet
```
