# Preprocessing

The pipeline requires three input files that are not included in this repository.
They must be generated from publicly available raw data before running any script.

---

## Full data flow

```
PUBLIC DATA                     PREPROCESSING                         PIPELINE INPUT
-----------                     ------------                          --------------

gnomAD v4.1 VCF (chr1–22)
  → bcftools query              extract CHROM, POS, REF, ALT,
                                AF, AN, AS_VQSLOD, VEP string
  → kmer7.py                    extract 7-mer context from hg38.fa
                                E_ref = −ln(WT_7mer_freq)
                                E_mut = −ln(MUT_7mer_freq)
  → 06_gnomad_annotated.py      parse VEP, add vep_category,
                                is_CpG, maf_bin, mutation_class
                                                              →  merged_annotated.csv.gz
                                                                   (step 01)

COSMIC v103 VCF
  → filter SNVs                 remove indels, validate bases
  → kmer7.py                    extract 7-mer context from hg38.fa
  → calc_energies.py            E_wt7 = −ln(WT_7mer_freq)
                                E_mut7 = −ln(MUT_7mer_freq)
  → label.py                    add is_cpg, in_cgc (Cancer Gene
                                Census), is_driver, is_coding
                                                              →  cosmic_final.tsv
                                                                   (steps 03–05)
  → annotate_regions.py         add region label
                                (exon / intron / extragenic)
                                                              →  cosmic_full_annotato.csv
                                                                   (step 02)
```

---

## Reference files required

| File | Source |
|------|--------|
| `hg38.fa` | NCBI / Ensembl GRCh38 reference genome |
| `hg38_k7_freq_noM.csv` | 7-mer frequencies computed from hg38 (no chrM) |
| `cancer_gene_census.tsv` | COSMIC Cancer Gene Census |

The 7-mer frequency table is central to the method:
`E = −ln(observed 7-mer frequency in hg38)` → `ΔE = E_ref − E_mut`

---

## Data sources

| Dataset | Version | URL |
|---------|---------|-----|
| gnomAD genomes | v4.1 (GRCh38) | https://gnomad.broadinstitute.org/downloads |
| COSMIC Genome Screens Mutant | v103 (GRCh38) | https://cancer.sanger.ac.uk/cosmic/download |
| Cancer Gene Census | v103 | https://cancer.sanger.ac.uk/census |
