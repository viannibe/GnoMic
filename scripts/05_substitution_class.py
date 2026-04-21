#!/usr/bin/env python3
"""
05_substitution_class.py
------------------------
Stratificazione per classe di sostituzione nucleotidica (6 classi SBS).

Filtra esoni non-CpG autosomici in entrambi i dataset.
Normalizza al pirimidinico: se REF ∈ {A, G} si prende il complemento di
REF e ALT, ottenendo le 6 classi canoniche:
  C>A  C>G  C>T  T>A  T>C  T>G

Per ogni classe: media ΔE per cromosoma (22 valori), paired t-test
COSMIC vs gnomAD, Cohen's d.

ATTENZIONE: richiede che il parquet gnomad_corrected.parquet sia stato
rigenerato con 01_build_parquet.py --force (deve contenere REF e ALT).

Input:
  --gnomad-parquet   results/cache/gnomad_corrected.parquet
  --cosmic-tsv       cosmic_final.tsv

Output:
  results/substitution_class/A_subst_kde.png
  results/substitution_class/subst_class_summary.tsv
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde, ttest_rel

# ---------------------------------------------------------------------------
PROJECT_DIR  = Path(__file__).resolve().parents[1]
DEFAULT_PARQ = str(PROJECT_DIR / "results/cache/gnomad_corrected.parquet")
DEFAULT_TSV  = "/leonardo_work/CNHPC_2116672/COSMIC_v103/processed/cosmic_final.tsv"
DEFAULT_OUT  = str(PROJECT_DIR / "results/substitution_class")
CHUNKSIZE    = 500_000
AUTOSOMES    = {str(i) for i in range(1, 23)}

GN_COLOR = "#1d3557"
CO_COLOR = "#e63946"
BW       = 0.8

SBS6 = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]

_COMPL = str.maketrans("ACGT", "TGCA")


def _to_pyrimidine(ref: str, alt: str) -> str | None:
    """
    Normalizza una sostituzione alla base pirimidinica.
    Restituisce la classe 'X>Y' oppure None se non è un SNV valido.
    """
    if not (isinstance(ref, str) and isinstance(alt, str)):
        return None
    ref = ref.strip().upper()
    alt = alt.strip().upper()
    if len(ref) != 1 or len(alt) != 1 or ref == alt:
        return None
    if ref in ("C", "T"):
        return f"{ref}>{alt}"
    # purinica → complemento
    return f"{ref.translate(_COMPL)}>{alt.translate(_COMPL)}"


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_gnomad(parquet_path: str) -> pd.DataFrame:
    """
    Carica gnomAD parquet, filtra esoni non-CpG autosomici.
    Richiede colonne REF e ALT (parquet rigenerato con --force).
    Ritorna aggregato per (sub_class, CHROM): mean_dE, n.
    """
    print(f"[gnomAD] reading {parquet_path}", flush=True)

    needed = {"CHROM", "REF", "ALT", "E_ref", "E_mut", "region", "is_cpg"}
    available = set(pd.read_parquet(parquet_path, columns=[]).columns
                    if False else  # workaround: leggi schema
                    pd.read_parquet(parquet_path).columns)
    missing = needed - available
    if missing:
        raise SystemExit(
            f"[ERROR] Il parquet manca delle colonne: {missing}\n"
            "  → Rigenera con:  python 01_build_parquet.py --force"
        )

    df = pd.read_parquet(
        parquet_path,
        columns=["CHROM", "REF", "ALT", "E_ref", "E_mut", "region", "is_cpg"],
    )
    df["dE"]     = df["E_ref"] - df["E_mut"]
    df["CHROM"]  = df["CHROM"].astype(str).str.replace("chr", "", regex=False)
    df["is_cpg"] = pd.to_numeric(df["is_cpg"], errors="coerce").fillna(0).astype(int)
    df = df[
        (df["region"] == "exon") &
        (df["is_cpg"] == 0) &
        df["CHROM"].isin(AUTOSOMES)
    ].copy()

    df["sub_class"] = [
        _to_pyrimidine(r, a) for r, a in zip(df["REF"], df["ALT"])
    ]
    df = df[df["sub_class"].isin(SBS6)]

    agg = (
        df.groupby(["sub_class", "CHROM"])["dE"]
        .agg(mean_dE="mean", n="count")
        .reset_index()
    )
    total = int(agg["n"].sum())
    cpg_note = "(non-CpG exon SNVs)"
    print(f"  {total:,} varianti {cpg_note}", flush=True)
    return agg   # columns: sub_class, CHROM, mean_dE, n


def load_cosmic(tsv_path: str) -> pd.DataFrame:
    """
    Carica cosmic_final.tsv, filtra coding non-CpG autosomici.
    Normalizza sub_class da ref/alt (colonne lowercase nel TSV).
    Ritorna aggregato per (sub_class, chrom): mean_dE, n.
    """
    print(f"[COSMIC] reading {tsv_path}", flush=True)
    records = []
    for chunk in pd.read_csv(
        tsv_path, sep="\t", chunksize=CHUNKSIZE,
        usecols=["chrom", "ref", "alt", "E_wt7", "E_mut7", "is_coding", "is_cpg"],
        dtype={"chrom": str, "is_coding": str, "ref": str, "alt": str,
               "is_cpg": "Int64"},
    ):
        chunk = chunk[chunk["is_coding"] == "coding"].copy()
        chunk["chrom"]  = chunk["chrom"].str.replace("chr", "", regex=False)
        chunk["is_cpg"] = pd.to_numeric(chunk["is_cpg"], errors="coerce").fillna(0).astype(int)
        chunk = chunk[
            (chunk["is_cpg"] == 0) &
            chunk["chrom"].isin(AUTOSOMES)
        ].copy()
        if chunk.empty:
            continue

        chunk["dE"] = (
            pd.to_numeric(chunk["E_wt7"],  errors="coerce")
            - pd.to_numeric(chunk["E_mut7"], errors="coerce")
        )
        chunk = chunk.dropna(subset=["dE"])

        chunk["sub_class"] = [
            _to_pyrimidine(r, a) for r, a in zip(chunk["ref"], chunk["alt"])
        ]
        chunk = chunk[chunk["sub_class"].isin(SBS6)]
        if chunk.empty:
            continue

        records.append(
            chunk.groupby(["sub_class", "chrom"])["dE"]
            .agg(sum_dE="sum", n="count")
            .reset_index()
        )

    raw = pd.concat(records, ignore_index=True)
    agg = (
        raw.groupby(["sub_class", "chrom"])
        .agg(sum_dE=("sum_dE", "sum"), n=("n", "sum"))
        .reset_index()
    )
    agg["mean_dE"] = agg["sum_dE"] / agg["n"]
    total = int(agg["n"].sum())
    print(f"  {total:,} coding non-CpG SNV rows", flush=True)
    return agg[["sub_class", "chrom", "mean_dE", "n"]]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _kde_x_y(vals: np.ndarray):
    kde = gaussian_kde(vals, bw_method="silverman")
    kde.set_bandwidth(kde.factor * BW)
    x = np.linspace(vals.min() - 0.15, vals.max() + 0.15, 500)
    return x, kde(x)


def paired_stats(gn_ser: pd.Series, co_ser: pd.Series):
    """t-test paired e Cohen's d su cromosomi condivisi."""
    shared = gn_ser.index.intersection(co_ser.index)
    if len(shared) < 3:
        return np.nan, np.nan, np.nan, np.array([])
    a = gn_ser[shared].values
    b = co_ser[shared].values
    diff = b - a
    _, p = ttest_rel(a, b)
    d = diff.mean() / diff.std(ddof=1) if diff.std(ddof=1) > 0 else np.nan
    return p, d, diff.mean(), diff


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

def build_summary(gn_agg: pd.DataFrame, co_agg: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cls in SBS6:
        gn_sub = gn_agg[gn_agg["sub_class"] == cls]
        co_sub = co_agg[co_agg["sub_class"] == cls]

        gn_ser = gn_sub.set_index("CHROM")["mean_dE"]
        co_ser = co_sub.set_index("chrom")["mean_dE"]
        p, d, delta_mu, _ = paired_stats(gn_ser, co_ser)

        rows.append({
            "sub_class":      cls,
            "n_gnomad":       int(gn_sub["n"].sum()),
            "n_cosmic":       int(co_sub["n"].sum()),
            "mean_dE_gnomad": round(float(gn_ser.mean()), 5) if len(gn_ser) else np.nan,
            "mean_dE_cosmic": round(float(co_ser.mean()), 5) if len(co_ser) else np.nan,
            "delta_mu":       round(float(delta_mu), 5) if not np.isnan(delta_mu) else np.nan,
            "p":              float(p)   if not np.isnan(p) else np.nan,
            "cohens_d":       round(float(d), 2) if not np.isnan(d) else np.nan,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Figure A: KDE 2×3
# ---------------------------------------------------------------------------

def plot_A(gn_agg: pd.DataFrame, co_agg: pd.DataFrame,
           summary: pd.DataFrame, out_dir: Path) -> None:

    fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharey=False)

    for ax, cls in zip(axes.flat, SBS6):
        row = summary[summary["sub_class"] == cls].iloc[0]

        gn_ser = gn_agg[gn_agg["sub_class"] == cls].set_index("CHROM")["mean_dE"]
        co_ser = co_agg[co_agg["sub_class"] == cls].set_index("chrom")["mean_dE"]
        shared = gn_ser.index.intersection(co_ser.index)

        pval_str = f"p = {row['p']:.2e}" if not np.isnan(row["p"]) else "p = n/a"
        d_str    = f"d = {row['cohens_d']:.1f}" if not np.isnan(row["cohens_d"]) else ""

        for vals, color, ds_label, n_key in [
            (gn_ser[shared].values, GN_COLOR, "gnomAD germline", "n_gnomad"),
            (co_ser[shared].values, CO_COLOR, "COSMIC somatic",  "n_cosmic"),
        ]:
            if len(vals) < 3:
                continue
            mu  = vals.mean()
            n_v = int(row[n_key])
            x, y = _kde_x_y(vals)
            ax.plot(x, y, color=color, lw=2.5,
                    label=f"{ds_label}  μ = {mu:+.3f}  (n = {n_v:,})")
            ax.fill_between(x, y, alpha=0.25, color=color)
            ax.axvline(mu, color=color, lw=1.2, ls="--", alpha=0.8)

        ax.axvline(0, color="gray", lw=1.0, ls="--", alpha=0.4, zorder=0)
        ax.set_title(f"{cls}   {pval_str}   {d_str}", fontsize=11)
        ax.set_xlabel("Mean ΔE per chromosome  (E_ref − E_mut)", fontsize=9)
        ax.set_ylabel("Density", fontsize=9)
        ax.legend(frameon=False, fontsize=8.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle(
        "Exonic non-CpG variants — substitution class KDE\n"
        "COSMIC somatic vs gnomAD germline  (per-chromosome means, n=22 autosomi)",
        fontsize=13, fontweight="bold",
    )
    fig.tight_layout()
    out_png = out_dir / "A_subst_kde.png"
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {out_png}", flush=True)


# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Stratificazione per classe di sostituzione (SBS6) — esoni non-CpG"
    )
    parser.add_argument("--gnomad-parquet", default=DEFAULT_PARQ)
    parser.add_argument("--cosmic-tsv",     default=DEFAULT_TSV)
    parser.add_argument("-o", "--outdir",   default=DEFAULT_OUT)
    args = parser.parse_args()

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gn_agg = load_gnomad(args.gnomad_parquet)
    co_agg = load_cosmic(args.cosmic_tsv)

    summary = build_summary(gn_agg, co_agg)

    print("\n=== SUMMARY ===")
    print(summary.to_string(index=False))

    out_tsv = out_dir / "subst_class_summary.tsv"
    summary.to_csv(out_tsv, sep="\t", index=False, float_format="%.6g")
    print(f"  ✓ {out_tsv}", flush=True)

    # Segnale diagnostico rapido
    sig_pos = summary[(summary["delta_mu"] > 0) & (summary["p"] < 0.05)]
    print(f"\n  Classi con delta_mu > 0 e p < 0.05 : {len(sig_pos)}/6")
    print(f"  {sig_pos['sub_class'].tolist()}")
    if len(sig_pos) == 6:
        print("  → segnale positivo e significativo in tutte e 6 le classi ✓")

    plot_A(gn_agg, co_agg, summary, out_dir)

    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
