#!/usr/bin/env python3
"""
01_build_parquet.py
-------------------
Legge merged_annotated.csv.gz e genera il parquet gnomAD corretto.

Operazioni:
  - Filtra autosomi 1-22
  - Calcola delta_E = E_ref - E_mut  (positivo → contesto mutato più raro)
  - Mappa vep_category → func_class → region (exon/intron/extragenic)
  - Standardizza is_cpg (0/1 int)
  - Conserva maf_bin, AS_VQSLOD, AF per analisi downstream

Output:
  results/cache/gnomad_corrected.parquet

Uso:
  python 01_build_parquet.py \\
      --annotated-csv  /path/to/merged_annotated.csv.gz \\
      --output         results/cache/gnomad_corrected.parquet

  Se il parquet esiste già, lo script lo salta (usa --force per rigenerarlo).
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
PROJECT_DIR  = Path(__file__).resolve().parents[1]
DEFAULT_CSV  = "/leonardo_work/CNHPC_2116672/gnomAD/processed/annotation_2/merged_annotated.csv.gz"
DEFAULT_PARQ = str(PROJECT_DIR / "results/cache/gnomad_corrected.parquet")
CHUNKSIZE    = 500_000
AUTOSOMES    = {str(i) for i in range(1, 23)}

# ---------------------------------------------------------------------------
# Functional class mappings  (identici a gnomAD_vi_clean 09_*)
# ---------------------------------------------------------------------------
FUNC_MAP_GN = {
    "missense":         "missense",
    "synonymous":       "synonymous",
    "nonsense":         "stop",
    "splice_region":    "splice",
    "splice_essential": "splice",
    "intron":           "intron",
    "utr_3":            "utr",
    "utr_5":            "utr",
    "upstream":         "extragenic",
    "downstream":       "extragenic",
    "noncoding_exon":   "extragenic",
    "non_coding":       "extragenic",
    "start_lost":       "stop",
    "stop_lost":        "stop",
    "frameshift":       "frameshift",
    "inframe_indel":    "inframe_indel",
    "other":            "extragenic",
}

REGION_MAP = {
    "missense":      "exon",
    "synonymous":    "exon",
    "stop":          "exon",
    "frameshift":    "exon",
    "inframe_indel": "exon",
    "utr":           "exon",
    "splice":        "intron",
    "intron":        "intron",
    "extragenic":    "extragenic",
}

WANT_COLS = {
    "CHROM", "REF", "ALT", "E_ref", "E_mut",
    "vep_category", "is_CpG", "is_cpg",
    "mutation_type_normalized",
    "is_coding", "maf_bin", "AS_VQSLOD", "AF",
}


# ---------------------------------------------------------------------------

def build_parquet(annotated_csv: str, out_parquet: str, force: bool = False) -> None:
    out_path = Path(out_parquet)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.exists() and not force:
        print(f"[skip] Parquet già esistente: {out_path}")
        print("       Usa --force per rigenerarlo.")
        return

    print(f"[01] Building parquet from:\n     {annotated_csv}\n")

    chunks_out = []
    total_in = total_out = 0

    for i, chunk in enumerate(pd.read_csv(
        annotated_csv,
        chunksize=CHUNKSIZE,
        low_memory=False,
        usecols=lambda c: c in WANT_COLS,
    )):
        total_in += len(chunk)

        # Autosomi 1-22
        chunk["CHROM"] = (
            chunk["CHROM"].astype(str)
            .str.replace("chr", "", regex=False).str.strip()
        )
        chunk = chunk[chunk["CHROM"].isin(AUTOSOMES)].copy()
        if chunk.empty:
            continue

        # Numerici
        for col in ("E_ref", "E_mut", "AS_VQSLOD", "AF"):
            if col in chunk.columns:
                chunk[col] = pd.to_numeric(chunk[col], errors="coerce")

        # delta_E = E_ref - E_mut  (positivo = WT più raro = mutazione destabilizzante)
        chunk["delta_E"] = chunk["E_ref"] - chunk["E_mut"]

        # is_cpg → int 0/1
        cpg_src = next((c for c in ("is_CpG", "is_cpg") if c in chunk.columns), None)
        if cpg_src:
            chunk["is_cpg"] = (
                chunk[cpg_src].astype(str)
                .map({"True": 1, "False": 0, "1": 1, "0": 0,
                      "true": 1, "false": 0})
                .fillna(0).astype(int)
            )
        else:
            chunk["is_cpg"] = 0

        # mut_class: classe di sostituzione normalizzata alla pirimidina (C>A, C>G, C>T, T>A, T>C, T>G)
        VALID_CLASSES = {"C>A", "C>G", "C>T", "T>A", "T>C", "T>G"}
        if "mutation_type_normalized" in chunk.columns:
            chunk["mut_class"] = chunk["mutation_type_normalized"].where(
                chunk["mutation_type_normalized"].isin(VALID_CLASSES), other=pd.NA
            )
        else:
            chunk["mut_class"] = pd.NA

        # func_class da vep_category
        if "vep_category" in chunk.columns:
            chunk["func_class"] = (
                chunk["vep_category"].str.strip().str.lower().map(FUNC_MAP_GN)
            )
            chunk["func_class"] = chunk["func_class"].fillna(
                chunk["vep_category"].str.strip().str.lower()
            )
        else:
            chunk["func_class"] = "other"

        chunk["region"] = chunk["func_class"].map(REGION_MAP).fillna("extragenic")

        chunk = chunk.dropna(subset=["delta_E", "E_ref", "E_mut"])
        if chunk.empty:
            continue

        keep = [c for c in (
            "CHROM", "REF", "ALT", "E_ref", "E_mut", "delta_E",
            "is_cpg", "func_class", "region",
            "mut_class",
            "maf_bin", "AS_VQSLOD", "AF",
        ) if c in chunk.columns]

        chunks_out.append(chunk[keep])
        total_out += len(chunk)

        if (i + 1) % 10 == 0:
            print(f"  chunk {i+1:3d}: letto={total_in:>12,}  salvato={total_out:>12,}")

    if not chunks_out:
        sys.exit("ERRORE: nessuna riga sopravvissuta al parsing.")

    print(f"\n  Concatenazione {len(chunks_out)} chunks...")
    df = pd.concat(chunks_out, ignore_index=True)

    print(f"\n  === SANITY CHECK ===")
    print(f"  righe totali   : {len(df):,}")
    print(f"  mean(delta_E)  : {df['delta_E'].mean():+.4f}  (atteso > 0)")
    print(f"  frac_CpG       : {df['is_cpg'].mean():.4f}")
    print(f"  region:\n{df['region'].value_counts().to_string()}")
    if "maf_bin" in df.columns:
        print(f"  maf_bin:\n{df['maf_bin'].value_counts().to_string()}")
    if "AS_VQSLOD" in df.columns:
        print(f"  AS_VQSLOD:\n{df['AS_VQSLOD'].describe().to_string()}")

    df.to_parquet(out_path, index=False)
    print(f"\n  ✓ Parquet salvato: {out_path}  ({len(df):,} righe)")


# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Costruisce il parquet gnomAD da merged_annotated.csv.gz"
    )
    parser.add_argument("--annotated-csv", default=DEFAULT_CSV)
    parser.add_argument("--output",        default=DEFAULT_PARQ)
    parser.add_argument("--force", action="store_true",
                        help="Rigenera il parquet anche se esiste già")
    args = parser.parse_args()
    build_parquet(args.annotated_csv, args.output, args.force)


if __name__ == "__main__":
    main()
