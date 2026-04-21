#!/usr/bin/env python3
"""
02_kde_region.py
----------------
Genera A_kde_region_all.png: confronto COSMIC vs gnomAD per compartimento
genomico (esone / introne / extragenic).

Metodo:
  - Aggrega per cromosoma (n=22 medie) sia COSMIC che gnomAD
  - KDE sulle 22 medie cromosomiche
  - Mann-Whitney U test (two-sided) per ogni regione
  - delta_mu = mean(COSMIC) - mean(gnomAD)

Input:
  --gnomad-parquet   results/cache/gnomad_corrected.parquet  (da script 01)
  --cosmic-csv       cosmic_full_annotato.csv  (ha ANN_TOP_ANNOTATION → region)

Output (in --outdir, default results/kde_region/):
  A_kde_region_all.png
  A_kde_region_all.csv     ← statistiche per regione
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
PROJECT_DIR   = Path(__file__).resolve().parents[1]
DEFAULT_PARQ  = str(PROJECT_DIR / "results/cache/gnomad_corrected.parquet")
DEFAULT_CSV   = "/leonardo_work/CNHPC_2116672/COSMIC_v103/processed_ws2/cosmic_full_annotato.csv"
DEFAULT_OUT   = str(PROJECT_DIR / "results/kde_region")
CHUNKSIZE     = 500_000
AUTOSOMES     = {str(i) for i in range(1, 23)}
REGIONS       = ["exon", "intron", "extragenic"]

# ---------------------------------------------------------------------------
# Functional class mappings — COSMIC ANN_TOP_ANNOTATION (identici a 09_*)
# ---------------------------------------------------------------------------
FUNC_MAP_CO = {
    "missense_variant":                               "missense",
    "missense_variant&splice_region_variant":         "missense",
    "synonymous_variant":                             "synonymous",
    "stop_retained_variant":                          "synonymous",
    "splice_region_variant&synonymous_variant":       "synonymous",
    "stop_gained":                                    "stop",
    "stop_gained&splice_region_variant":              "stop",
    "stop_lost":                                      "stop",
    "start_lost":                                     "stop",
    "start_lost&splice_region_variant":               "stop",
    "splice_acceptor_variant":                        "splice",
    "splice_donor_variant":                           "splice",
    "splice_acceptor_variant&intron_variant":         "splice",
    "splice_donor_variant&intron_variant":            "splice",
    "splice_region_variant":                          "splice",
    "splice_region_variant&intron_variant":           "splice",
    "intron_variant":                                 "intron",
    "3_prime_utr_variant":                            "utr",
    "5_prime_utr_variant":                            "utr",
    "5_prime_utr_premature_start_codon_gain_variant": "utr",
    "upstream_gene_variant":                          "extragenic",
    "downstream_gene_variant":                        "extragenic",
    "non_coding_transcript_exon_variant":             "extragenic",
    "intergenic_region":                              "extragenic",
    "nmd_transcript_variant":                         "extragenic",
    "inframe_deletion":                               "inframe_indel",
    "inframe_insertion":                              "inframe_indel",
    "frameshift_variant":                             "frameshift",
    "disruptive_inframe_deletion":                    "inframe_indel",
    "disruptive_inframe_insertion":                   "inframe_indel",
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

COLORS = {"COSMIC": "#e63946", "gnomAD": "#1d3557"}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _cpg_7mer(seq) -> int:
    """1 se il 7-mer ha CpG in posizione centrale o complementare."""
    if pd.isna(seq) or len(str(seq)) != 7:
        return 0
    s = str(seq)
    return 1 if (s[3] == "C" and s[4] == "G") or (s[3] == "G" and s[2] == "C") else 0


def chrom_mean(df: pd.DataFrame, val_col: str, chrom_col: str = "CHROM") -> np.ndarray:
    """Media per cromosoma; ritorna array di 22 valori (o meno se qualche chrom manca)."""
    df2 = df.copy()
    df2["_c"] = df2[chrom_col].astype(str).str.replace("chr", "", regex=False)
    df2 = df2[df2["_c"].isin(AUTOSOMES)]
    agg = df2.groupby("_c")[val_col].mean()
    return agg.values


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_gnomad(parquet_path: str) -> pd.DataFrame:
    print(f"[gnomAD] reading {parquet_path}", flush=True)
    df = pd.read_parquet(parquet_path, columns=["CHROM", "delta_E", "region", "is_cpg"])
    df["CHROM"] = df["CHROM"].astype(str).str.replace("chr", "", regex=False)
    df = df[df["CHROM"].isin(AUTOSOMES)].copy()
    df["is_cpg"] = pd.to_numeric(df["is_cpg"], errors="coerce").fillna(0).astype(int)
    print(f"  {len(df):,} righe  |  mean(delta_E)={df['delta_E'].mean():+.4f}", flush=True)
    return df


def load_cosmic(csv_path: str) -> pd.DataFrame:
    print(f"[COSMIC] reading {csv_path}", flush=True)
    want = {"CHROM", "WT_7mer", "WT_freq", "MUT_freq", "ANN_TOP_ANNOTATION"}
    chunks = []
    eps = 1e-12
    for chunk in pd.read_csv(
        csv_path, chunksize=CHUNKSIZE, low_memory=False,
        usecols=lambda c: c in want,
    ):
        chunk["WT_freq"]  = pd.to_numeric(chunk["WT_freq"],  errors="coerce")
        chunk["MUT_freq"] = pd.to_numeric(chunk["MUT_freq"], errors="coerce")
        chunk["E_wt"]     = -np.log(chunk["WT_freq"]  + eps)
        chunk["E_mut"]    = -np.log(chunk["MUT_freq"] + eps)
        chunk["delta"]    = chunk["E_wt"] - chunk["E_mut"]
        chunk["is_cpg"]   = chunk["WT_7mer"].apply(_cpg_7mer)

        ann = chunk["ANN_TOP_ANNOTATION"].str.strip().str.lower()
        chunk["func_class"] = ann.map(FUNC_MAP_CO).fillna(ann)
        chunk["region"]     = chunk["func_class"].map(REGION_MAP).fillna("extragenic")

        chunk["CHROM"] = (
            chunk["CHROM"].astype(str).str.replace("chr", "", regex=False)
        )
        chunks.append(
            chunk[["CHROM", "delta", "is_cpg", "func_class", "region"]]
            .dropna(subset=["delta"])
        )
    df = pd.concat(chunks, ignore_index=True)
    df = df[df["CHROM"].isin(AUTOSOMES)]
    print(f"  {len(df):,} righe  |  mean(delta)={df['delta'].mean():+.4f}", flush=True)
    return df


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot_kde_region(cosmic: pd.DataFrame, gnomad: pd.DataFrame, out_dir: Path) -> None:
    sns.set_theme(style="whitegrid")
    plt.rcParams["font.size"] = 10

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    rows = []

    for ax, reg in zip(axes, REGIONS):
        co_sub = cosmic[cosmic["region"] == reg]
        gn_sub = gnomad[gnomad["region"] == reg]

        co_vals = chrom_mean(co_sub, "delta",   chrom_col="CHROM")
        gn_vals = chrom_mean(gn_sub, "delta_E", chrom_col="CHROM")

        n_co = len(co_sub)
        n_gn = len(gn_sub)

        # KDE
        for vals, label, color in [
            (co_vals, "COSMIC",  COLORS["COSMIC"]),
            (gn_vals, "gnomAD",  COLORS["gnomAD"]),
        ]:
            if len(vals) < 2:
                continue
            mu = vals.mean()
            sns.kdeplot(
                vals, ax=ax, fill=True, alpha=0.25, lw=2.5,
                color=color, bw_adjust=0.8,
                label=f"{label}  μ = {mu:+.3f}  (n={len(vals)})",
            )
            ax.axvline(mu, color=color, ls=":", lw=2)

        # Mann-Whitney
        pval = np.nan
        if len(co_vals) >= 3 and len(gn_vals) >= 3:
            _, pval = mannwhitneyu(co_vals, gn_vals, alternative="two-sided")

        delta_mu = (co_vals.mean() - gn_vals.mean()
                    if len(co_vals) and len(gn_vals) else np.nan)

        ax.axvline(0, color="#aaa", ls="--", lw=1, alpha=0.5)
        title = f"{reg.upper()}\nco={n_co:,}  gn={n_gn:,}"
        if not np.isnan(pval):
            title += f"  p={pval:.2e}"
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_xlabel("per-chromosome mean ΔE  (E_ref − E_mut)", fontsize=10)
        ax.set_ylabel("Density", fontsize=10)
        ax.legend(fontsize=9, frameon=True)
        ax.grid(axis="y", alpha=0.2)

        rows.append({
            "region":       reg,
            "n_cosmic_var": n_co,
            "n_gnomad_var": n_gn,
            "n_chrom_cosmic": len(co_vals),
            "n_chrom_gnomad": len(gn_vals),
            "mean_cosmic":  round(co_vals.mean(), 5) if len(co_vals) else np.nan,
            "mean_gnomad":  round(gn_vals.mean(), 5) if len(gn_vals) else np.nan,
            "delta_mu":     round(delta_mu, 5),
            "pval_mwu":     pval,
        })

    fig.suptitle(
        "COSMIC vs gnomAD  ΔE distributions by genomic region\n"
        "(per-chromosome means, n=22 autosomi)",
        fontsize=13, fontweight="bold",
    )
    fig.tight_layout()

    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / "A_kde_region_all.png"
    out_csv = out_dir / "A_kde_region_all.csv"

    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {out_png}", flush=True)

    stats = pd.DataFrame(rows)
    stats.to_csv(out_csv, index=False, float_format="%.6g")
    print(f"  ✓ {out_csv}", flush=True)

    print("\n  === STATISTICS ===")
    print(stats[["region", "mean_cosmic", "mean_gnomad", "delta_mu", "pval_mwu"]].to_string(index=False))


# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="KDE ΔE per regione genomica: COSMIC vs gnomAD"
    )
    parser.add_argument("--gnomad-parquet", default=DEFAULT_PARQ)
    parser.add_argument("--cosmic-csv",     default=DEFAULT_CSV)
    parser.add_argument("-o", "--outdir",   default=DEFAULT_OUT)
    args = parser.parse_args()

    gnomad = load_gnomad(args.gnomad_parquet)
    cosmic = load_cosmic(args.cosmic_csv)
    plot_kde_region(cosmic, gnomad, Path(args.outdir))
    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
