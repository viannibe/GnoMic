# Preprocessing Guide — Riproducibilità da zero

Questa guida descrive come rigenerare i file di input **da zero**, partendo dai VCF pubblici.

Se i file preprocessati **esistono già** su `/leonardo_work/CNHPC_2116672/`, puoi saltare questa fase e andare diretto a `scripts/01_build_parquet.py`.

---

## Overview: cosa viene generato

```
Cosmic_GenomeScreensMutant_v103_GRCh38.vcf  (raw)
  ↓
cosmic_snv_clean.tsv                     [step 01]
  ↓
cosmic_snv_with_7mer.tsv                 [step 02]
  ↓
cosmic_snv_with_energies.tsv             [step 03]
  ↓
cosmic_final.tsv                         [step 04] ← INPUT per scripts/03, 04
cosmic_full_annotato.csv                 [step 05] ← INPUT per scripts/02

gnomad.genomes.v4.1.sites.chr*.vcf.bgz (raw)
  ↓
gnomad_exonic.csv                        [bcftools query]
  ↓
merged_with_energies.csv.gz              [kmer7.py]
  ↓
merged_annotated.csv.gz                  [06_gnomad_annotated.py] ← INPUT per scripts/01
```

---

## COSMIC Preprocessing

### Prerequisiti

```bash
# Accertati di avere:
# - Ambiente conda: gnomad_venv
# - Python 3.9+
# - pysam
# - pandas
# - numpy
```

### Step 0: Scarica COSMIC v103 raw VCF

Accedi a: `https://cancer.sanger.ac.uk/cosmic/download`
- Seleziona "Genome Screens Mutant"
- Versione: v103
- Genoma: GRCh38
- Scarica: `Cosmic_GenomeScreensMutant_v103_GRCh38.vcf.gz`

Salva in:
```
/leonardo_work/CNHPC_2116672/COSMIC_v103/mutation/
```

Decomprimi se necessario:
```bash
gunzip Cosmic_GenomeScreensMutant_v103_GRCh38.vcf.gz
```

### Step 1: VCF → Clean TSV

Script: `gnomAD_vi_clean/scripts/cosmic_pipeline/01_vcf_to_clean.py`

```bash
cd /leonardo/home/userexternal/viannibe/GnoMic

python ../../gnomAD_vi_clean/scripts/cosmic_pipeline/01_vcf_to_clean.py
```

Output: `/leonardo_work/CNHPC_2116672/COSMIC_v103/processed/cosmic_snv_clean.tsv`

**Cosa fa:**
- Filtra solo SNV (no INDEL)
- Esclude chrM
- Valida basi (solo A/T/C/G)
- Esclude automutazioni (REF ≠ ALT)
- Conserva annotazioni COSMIC (gene, consequence, sample_count, ecc.)

### Step 2: Aggiungi contesto 7-mer

Script: `gnomAD_vi_clean/scripts/cosmic_pipeline/02_add_7mer.py`

```bash
python ../../gnomAD_vi_clean/scripts/cosmic_pipeline/02_add_7mer.py \
  --in  /leonardo_work/CNHPC_2116672/COSMIC_v103/processed/cosmic_snv_clean.tsv \
  --ref /leonardo_work/CNHPC_2116672/COSMIC_v103/ref/hg38.fa \
  --out /leonardo_work/CNHPC_2116672/COSMIC_v103/processed/cosmic_snv_with_7mer.tsv \
  --k 7 \
  --chunksize 200000
```

Output: `/leonardo_work/CNHPC_2116672/COSMIC_v103/processed/cosmic_snv_with_7mer.tsv`

**Cosa fa:**
- Estrae il contesto 7-mer attorno a ogni posizione (dal genoma hg38.fa)
- Calcola caratteristiche sequenza: GC%, Shannon entropy, homopolymer run
- Identifica CpG nel contesto centrale
- Canonicalizza il 7-mer (revcomp se occorre)

**Prerequisiti:**
- Genoma: `/leonardo_work/CNHPC_2116672/COSMIC_v103/ref/hg38.fa` (indicizzato con samtools)

### Step 3: Calcola energie

Script: `gnomAD_vi_clean/scripts/cosmic_pipeline/03_calc_energies.py`

```bash
python ../../gnomAD_vi_clean/scripts/cosmic_pipeline/03_calc_energies.py
```

Output: `/leonardo_work/CNHPC_2116672/COSMIC_v103/processed/cosmic_snv_with_energies.tsv`

**Cosa fa:**
- Legge la tabella di frequenza 7-mer genome-wide: `hg38_k7_freq_noM.csv`
- Calcola `E_wt7 = −ln(freq_WT_7mer)`
- Calcola `E_mut7 = −ln(freq_MUT_7mer)`
- Calcola `deltaE = E_wt − E_mut`
- Aggiunge pseudocount (1e-12) per k-mer mai osservati

**Prerequisiti:**
- Tabella di frequenza: `/leonardo_work/CNHPC_2116672/COSMIC_v103/ref/processed/hg38_k7_freq_noM.csv`
  (scaricata da COSMIC o precomputata dal genoma hg38)

### Step 4: Aggiungi etichette (CGC + driver/passenger)

Script: `gnomAD_vi_clean/scripts/cosmic_pipeline/04_label.py`

```bash
python ../../gnomAD_vi_clean/scripts/cosmic_pipeline/04_label.py \
  --in /leonardo_work/CNHPC_2116672/COSMIC_v103/processed/cosmic_snv_with_energies.tsv \
  --cgc /leonardo_work/CNHPC_2116672/COSMIC_v103/census/cancer_gene_census.tsv \
  --out /leonardo_work/CNHPC_2116672/COSMIC_v103/processed/cosmic_final.tsv \
  --min_recurrence 3 \
  --canonical_only
```

Output: `/leonardo_work/CNHPC_2116672/COSMIC_v103/processed/cosmic_final.tsv`

**Cosa fa:**
- Legge Cancer Gene Census (lista driver genes)
- Aggiunge flag `in_cgc` (1 se il gene è un driver noto)
- Classifica `is_coding` (coding vs non-coding)
- Etichetta `is_driver` se: in_cgc=1 AND sample_count≥3 AND is_coding="coding"

**Prerequisiti:**
- Cancer Gene Census: `/leonardo_work/CNHPC_2116672/COSMIC_v103/census/cancer_gene_census.tsv`
  (scaricato da: `https://cancer.sanger.ac.uk/census`)

---

## gnomAD Preprocessing

### Step 0: Scarica gnomAD v4.1

Scarica i VCF per ogni cromosoma da:
```
https://gnomad.broadinstitute.org/downloads
→ v4.1
→ Genomes (GRCh38)
→ gnomad.genomes.v4.1.sites.chrN.vcf.bgz (N=1..22)
```

Salva in:
```
/leonardo_work/CNHPC_2116672/gnomAD/gnomad_genomes_4.1/
```

### Step 1: Estrai campi + VEP con bcftools

```bash
# Merge tutti i cromosomi (una tantum)
bcftools concat \
  gnomad.genomes.v4.1.sites.chr{1..22}.vcf.bgz \
  -o gnomad.genomes.all_chroms.vcf.bgz -O z

# Indizza
bcftools index gnomad.genomes.all_chroms.vcf.bgz

# Estrai campi
bcftools query -H \
  -f "%CHROM,%POS,%REF,%ALT,%ID,%FILTER,%AC,%AF,%AN,%allele_type,%AS_VQSLOD,\
$(for p in afr nfe eas mid sas grpmax; do echo -n "%AC_$p,%AF_$p,%AN_$p,"; done)%vep\n" \
  gnomad.genomes.all_chroms.vcf.bgz \
  > /leonardo_work/CNHPC_2116672/gnomAD/gnomad_exonic.csv
```

Output: `gnomad_exonic.csv`

**Colonne estratte:**
- Base: CHROM, POS, REF, ALT, FILTER, AC, AF, AN, allele_type, AS_VQSLOD
- PopAF: AC_afr, AF_afr, AN_afr, ... (per ogni popolazione)
- VEP: vep (stringa raw VEP SnpEff format)

### Step 2: Aggiungi energie 7-mer

Script: `kmer7.py` (in `gnomAD_vi_clean/` o in home utente di lavoro)

Questo script **non è nel repo**, ma l'output è salvato in:
```
/leonardo_work/CNHPC_2116672/gnomAD/processed/merged_with_energies.csv.gz
```

**Operazioni (dal codice del progetto):**
- Legge gnomad_exonic.csv
- Estrae 7-mer WT e MUT dal genoma hg38.fa
- Cerca la frequenza in hg38_k7_freq_withPos.csv
- Calcola E_ref = −ln(WT_freq), E_mut = −ln(MUT_freq)
- Salva in CSV.gz

**Se non esiste**, puoi rigenerare con:
```bash
# Pseudo-codice (si trova in gnomAD_vi_clean/data/ se disponibile)
python kmer7.py \
  --input gnomad_exonic.csv \
  --ref hg38.fa \
  --kmer_freq hg38_k7_freq_withPos.csv \
  --output merged_with_energies.csv.gz
```

### Step 3: VEP parsing + CpG + maf_bin

Script: `gnomAD_vi_clean/scripts/06_gnomad_annotated.py`

```bash
python ../../gnomAD_vi_clean/scripts/06_gnomad_annotated.py
```

(Path hardcodati nel file; modifica se necessario)

Output: `/leonardo_work/CNHPC_2116672/gnomAD/processed/annotation_2/merged_annotated.csv.gz`

**Cosa fa:**
- Parsa la stringa VEP raw (SnpEff format)
- Seleziona la conseguenza più grave (severity order)
- Mappa a `vep_category` (forma corta: missense/synonymous/intron/ecc.)
- Calcola `is_CpG` dal 7-mer (central C-G o G-C)
- Calcola `maf_bin` dall'AF (ultra_rare/very_rare/rare/ecc.)
- Aggiunge `mutation_class` (transition/transversion)
- Aggiunge `is_coding` (infeito dall'HGVS protein o dalle conseguenze)

**Colonne output:**
```
CHROM, POS, REF, ALT, AF, AN, AS_VQSLOD,
E_ref, E_mut,
vep_category, vep_consequence, vep_impact, vep_symbol,
is_CpG, maf_bin, mutation_class, mutation_type, is_coding, 
functional_impact
```

---

## Tabelle di riferimento (pre-scaricate)

Questi file devono già esistere (scaricati separatamente):

| File | Path | Fonte | Uso |
|------|------|-------|-----|
| `hg38.fa` | `/leonardo_work/CNHPC_2116672/COSMIC_v103/ref/hg38.fa` | NCBI/Ensembl | Estrazione 7-mer |
| `hg38_k7_freq_noM.csv` | `/leonardo_work/CNHPC_2116672/COSMIC_v103/ref/processed/hg38_k7_freq_noM.csv` | precomputed | Frequenze 7-mer COSMIC |
| `hg38_k7_freq_withPos.csv` | `/leonardo_work/CNHPC_2116672/gnomAD/hg38_k7_freq_withPos.csv` | precomputed | Frequenze 7-mer gnomAD |
| `cancer_gene_census.tsv` | `/leonardo_work/CNHPC_2116672/COSMIC_v103/census/` | COSMIC | Driver genes |

---

## Shortcut: usa file preesistenti

Se i file di input **esistono già** su scratch:

```bash
# Non serve rigenerare; usa direttamente:
merged_annotated.csv.gz      → scripts/01_build_parquet.py
cosmic_full_annotato.csv     → scripts/02_kde_region.py
cosmic_final.tsv             → scripts/03_cpg_split.py, scripts/04_driver_passenger.py
```

Verifica che i file siano presenti:
```bash
ls /leonardo_work/CNHPC_2116672/gnomAD/processed/annotation_2/merged_annotated.csv.gz
ls /leonardo_work/CNHPC_2116672/COSMIC_v103/processed_ws2/cosmic_full_annotato.csv
ls /leonardo_work/CNHPC_2116672/COSMIC_v103/processed/cosmic_final.tsv
```

Se tutti e tre **esistono**, puoi iniziare da `scripts/01_build_parquet.py`.

---

## Note sulla riproducibilità

1. **Script originali**: Tutti gli script di preprocessing sono in
   ```
   /leonardo/home/userexternal/viannibe/gnomAD_vi_clean/scripts/cosmic_pipeline/
   /leonardo/home/userexternal/viannibe/gnomAD_vi_clean/scripts/06_gnomad_annotated.py
   ```
   Sono script **non modificati** da quelli originali del progetto.

2. **Hardcoding di percorsi**: Gli script hanno percorsi hardcodati. Se vuoi modificarli,
   edita i file `.py` alle prime righe (sezione `PATHS` o `INPUT_FILE`/`OUTPUT_FILE`).

3. **Versioni software**: Verifica che sia disponibile
   - Python 3.9+
   - pandas, numpy, pysam
   - bcftools (per gnomAD extraction)

4. **Dipendenze di file**: Ogni step dipende dal precedente. Non puoi saltare step!
