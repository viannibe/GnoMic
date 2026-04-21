#!/usr/bin/env python3
"""
03_cpg_split.py
---------------
Produce tre figure per le varianti esoniche, confronto COSMIC vs gnomAD:

  A_exon_cpg_split.png
      KDE a 2 pannelli: non-CpG | CpG
      Dimostra che il segnale ΔE COSMIC > gnomAD NON è spiegato dalla
      deaminazione CpG: persiste anche nei soli non-CpG.
      → A_exon_cpg_split_summary.csv

  B_exon_noncpg_paired.png
      Scatter (gnomAD vs COSMIC, per-chrom) + strip plot (differenza paired)
      per le sole varianti non-CpG esoniche.
      → B_exon_noncpg_paired.csv

  C_chrom_consistency.png
      Strip plot standalone: differenza COSMIC−gnomAD per cromosoma
      (tutti 22 autosomi positivi).  Identica al pannello destro di B
      ma come figura indipendente ad alta risoluzione.
      → stessi dati di B

Metodo (identico a gnomAD_vi_clean 13c_combined_AB.py):
  - Per ogni (is_cpg, chrom) calcola mean ΔE aggregando varianti
  - Paired t-test sulle 22 medie cromosomiche condivise
  - Cohen's d = mean(diff) / std(diff)

Input:
  --gnomad-parquet   results/cache/gnomad_corrected.parquet
  --cosmic-tsv       cosmic_final.tsv   (ha E_wt7, E_mut7, is_coding, is_cpg)

Output:
  results/cpg_split/A_exon_cpg_split.png
  results/cpg_split/A_exon_cpg_split_summary.csv
  results/cpg_split/B_exon_noncpg_paired.png
  results/cpg_split/B_exon_noncpg_paired.csv
  results/cpg_split/C_chrom_consistency.png
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
DEFAULT_OUT  = str(PROJECT_DIR / "results/cpg_split")
CHUNKSIZE    = 500_000
AUTOSOMES    = {str(i) for i in range(1, 23)}

GN_COLOR = "#1d3557"
CO_COLOR = "#e63946"
BW       = 0.8   # bandwidth multiplier (Silverman × BW)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_gnomad(parquet_path: str) -> pd.DataFrame:
    """
    Carica gnomAD parquet, filtra esoni autosomici.
    Ritorna aggregato per (is_cpg, CHROM): mean_dE, n.
    """
    print(f"[gnomAD] reading {parquet_path}", flush=True)
    df = pd.read_parquet(
        parquet_path,
        columns=["CHROM", "E_ref", "E_mut", "region", "is_cpg"],
    )
    df["dE"]     = df["E_ref"] - df["E_mut"]
    df["CHROM"]  = df["CHROM"].astype(str).str.replace("chr", "", regex=False)
    df["is_cpg"] = pd.to_numeric(df["is_cpg"], errors="coerce").fillna(0).astype(int)
    df = df[(df["region"] == "exon") & df["CHROM"].isin(AUTOSOMES)].copy()

    agg = (
        df.groupby(["is_cpg", "CHROM"])["dE"]
        .agg(mean_dE="mean", n="count")
        .reset_index()
    )
    print(f"  {len(df):,} exon rows  |  CpG frac={df['is_cpg'].mean():.3f}", flush=True)
    return agg   # columns: is_cpg, CHROM, mean_dE, n


def load_cosmic(tsv_path: str) -> pd.DataFrame:
    """
    Carica cosmic_final.tsv, filtra coding autosomici.
    Ritorna aggregato per (is_cpg, chrom): mean_dE, n.
    """
    print(f"[COSMIC] reading {tsv_path}", flush=True)
    records = []
    for chunk in pd.read_csv(
        tsv_path, sep="\t", chunksize=CHUNKSIZE,
        usecols=["chrom", "E_wt7", "E_mut7", "is_coding", "is_cpg"],
        dtype={"chrom": str, "is_coding": str, "is_cpg": "Int64"},
    ):
        chunk = chunk[chunk["is_coding"] == "coding"].copy()
        chunk["chrom"]   = chunk["chrom"].str.replace("chr", "", regex=False)
        chunk["is_cpg"]  = pd.to_numeric(chunk["is_cpg"], errors="coerce").fillna(0).astype(int)
        chunk = chunk[chunk["chrom"].isin(AUTOSOMES)].copy()
        if chunk.empty:
            continue
        chunk["dE"] = (
            pd.to_numeric(chunk["E_wt7"],  errors="coerce")
            - pd.to_numeric(chunk["E_mut7"], errors="coerce")
        )
        chunk = chunk.dropna(subset=["dE"])
        records.append(
            chunk.groupby(["is_cpg", "chrom"])["dE"]
            .agg(sum_dE="sum", n="count")
            .reset_index()
        )

    raw = pd.concat(records, ignore_index=True)
    raw["sum_dE_w"] = raw["sum_dE"]
    agg = (
        raw.groupby(["is_cpg", "chrom"])
        .agg(sum_dE=("sum_dE_w", "sum"), n=("n", "sum"))
        .reset_index()
    )
    agg["mean_dE"] = agg["sum_dE"] / agg["n"]
    total = int(agg["n"].sum())
    print(f"  {total:,} coding rows totali", flush=True)
    return agg[["is_cpg", "chrom", "mean_dE", "n"]]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _kde_x_y(vals: np.ndarray):
    kde = gaussian_kde(vals, bw_method="silverman")
    kde.set_bandwidth(kde.factor * BW)
    x = np.linspace(vals.min() - 0.15, vals.max() + 0.15, 500)
    return x, kde(x)


def paired_stats(a: np.ndarray, b: np.ndarray):
    """t-test rel e Cohen's d su b - a."""
    diff = b - a
    t, p = ttest_rel(a, b)
    d    = diff.mean() / diff.std(ddof=1) if diff.std(ddof=1) > 0 else np.nan
    return t, p, d, diff.mean(), diff.std(ddof=1)


# ---------------------------------------------------------------------------
# Figure A: KDE CpG split
# ---------------------------------------------------------------------------

def build_summary(gn_agg: pd.DataFrame, co_agg: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cpg_val, cpg_label in [(0, "non-CpG"), (1, "CpG")]:
        for ds_label, agg, chrom_col in [
            ("gnomAD", gn_agg, "CHROM"),
            ("COSMIC", co_agg, "chrom"),
        ]:
            sub  = agg[agg["is_cpg"] == cpg_val]
            vals = sub["mean_dE"].dropna().values
            rows.append({
                "dataset":    ds_label,
                "cpg_class":  cpg_label,
                "n_variants": int(sub["n"].sum()),
                "n_chrom":    len(vals),
                "mean_dE":    round(float(vals.mean()), 5) if len(vals) else np.nan,
                "std_dE":     round(float(vals.std(ddof=1)), 5) if len(vals) > 1 else np.nan,
                "sem_dE":     round(float(vals.std(ddof=1) / np.sqrt(len(vals))), 6)
                              if len(vals) > 1 else np.nan,
            })
    return pd.DataFrame(rows)


def plot_A(gn_agg: pd.DataFrame, co_agg: pd.DataFrame,
           summary: pd.DataFrame, out_dir: Path) -> None:

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax, cpg_val, cpg_label in [
        (axes[0], 0, "non-CpG variants"),
        (axes[1], 1, "CpG variants"),
    ]:
        gn_vals = gn_agg[gn_agg["is_cpg"] == cpg_val]["mean_dE"].dropna().values
        co_vals = co_agg[co_agg["is_cpg"] == cpg_val]["mean_dE"].dropna().values

        # paired t-test su cromosomi condivisi
        gn_ser  = gn_agg[gn_agg["is_cpg"] == cpg_val].set_index("CHROM")["mean_dE"]
        co_ser  = co_agg[co_agg["is_cpg"] == cpg_val].set_index("chrom")["mean_dE"]
        shared  = gn_ser.index.intersection(co_ser.index)
        _, pval = ttest_rel(gn_ser[shared].values, co_ser[shared].values) \
                  if len(shared) >= 3 else (np.nan, np.nan)

        for vals, color, ds_label in [
            (gn_vals, GN_COLOR, "gnomAD germline"),
            (co_vals, CO_COLOR, "COSMIC somatic"),
        ]:
            if len(vals) < 3:
                continue
            mu    = vals.mean()
            n_var = int(
                summary.loc[
                    (summary["dataset"] == ("gnomAD" if color == GN_COLOR else "COSMIC"))
                    & (summary["cpg_class"] == ("non-CpG" if cpg_val == 0 else "CpG")),
                    "n_variants",
                ].values[0]
            )
            x, y = _kde_x_y(vals)
            ax.plot(x, y, color=color, lw=2.5,
                    label=f"{ds_label}  μ = {mu:+.3f}  (n = {n_var:,})")
            ax.fill_between(x, y, alpha=0.25, color=color)
            ax.axvline(mu, color=color, lw=1.2, ls="--", alpha=0.8)

        ax.axvline(0, color="gray", lw=1.0, ls="--", alpha=0.4, zorder=0)
        pval_str = f"p = {pval:.2e}" if not np.isnan(pval) else "p = n/a"
        ax.set_title(f"{cpg_label}  ({pval_str})", fontsize=11)
        ax.set_xlabel("Mean ΔE per chromosome  (E_ref − E_mut)", fontsize=10)
        ax.set_ylabel("Density", fontsize=10)
        ax.legend(frameon=False, fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle(
        "Exonic variants — CpG split  (COSMIC somatic vs gnomAD germline)\n"
        "per-chromosome means, n=22 autosomi",
        fontsize=12, fontweight="bold",
    )
    fig.tight_layout()
    out_png = out_dir / "A_exon_cpg_split.png"
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {out_png}", flush=True)


# ---------------------------------------------------------------------------
# Paired table (non-CpG only)
# ---------------------------------------------------------------------------

def build_paired_table(gn_agg: pd.DataFrame, co_agg: pd.DataFrame) -> pd.DataFrame:
    gn = (gn_agg[gn_agg["is_cpg"] == 0]
          .rename(columns={"CHROM": "chrom", "mean_dE": "gnomad_mean_dE", "n": "n_gnomad"}))
    co = (co_agg[co_agg["is_cpg"] == 0]
          .rename(columns={"mean_dE": "cosmic_mean_dE", "n": "n_cosmic"}))
    m = pd.merge(
        gn[["chrom", "gnomad_mean_dE", "n_gnomad"]],
        co[["chrom", "cosmic_mean_dE", "n_cosmic"]],
        on="chrom", how="inner",
    )
    m["diff"]      = m["cosmic_mean_dE"] - m["gnomad_mean_dE"]
    m["chrom_int"] = m["chrom"].astype(int)
    m = m.sort_values("chrom_int").drop(columns="chrom_int")
    return m[["chrom", "gnomad_mean_dE", "cosmic_mean_dE", "diff", "n_gnomad", "n_cosmic"]]


# ---------------------------------------------------------------------------
# Figure B: scatter + strip
# ---------------------------------------------------------------------------

def plot_B(tbl: pd.DataFrame, out_dir: Path) -> None:
    gn_v   = tbl["gnomad_mean_dE"].values
    co_v   = tbl["cosmic_mean_dE"].values
    diff_v = tbl["diff"].values
    chroms = tbl["chrom"].values

    t_stat, pval = ttest_rel(gn_v, co_v)
    mean_diff    = diff_v.mean()
    sd_diff      = diff_v.std(ddof=1)
    cohens_d     = mean_diff / sd_diff if sd_diff > 0 else np.nan

    cmap   = plt.cm.tab20
    colors = [cmap(i % 20) for i in range(len(chroms))]

    fig, (ax_sc, ax_st) = plt.subplots(1, 2, figsize=(12, 5))

    # Scatter
    ax_sc.scatter(gn_v, co_v, color=GN_COLOR, s=70, zorder=3)
    lo = min(gn_v.min(), co_v.min()) - 0.01
    hi = max(gn_v.max(), co_v.max()) + 0.01
    ax_sc.plot([lo, hi], [lo, hi], color="gray", ls="--", lw=1.2, alpha=0.6)
    for x, y, ch in zip(gn_v, co_v, chroms):
        ax_sc.annotate(ch, (x, y), textcoords="offset points",
                       xytext=(5, 3), fontsize=7.5, color=GN_COLOR)
    ax_sc.set_xlabel("gnomAD germline  mean ΔE", fontsize=10)
    ax_sc.set_ylabel("COSMIC somatic  mean ΔE",  fontsize=10)
    ax_sc.set_title(
        f"Non-CpG exonic  (per-chromosome means)\n"
        f"p = {pval:.2e}  |  d = {cohens_d:.1f}",
        fontsize=10,
    )
    ax_sc.spines["top"].set_visible(False)
    ax_sc.spines["right"].set_visible(False)

    # Strip plot
    y_pos = np.arange(len(diff_v), 0, -1)
    ax_st.axvline(mean_diff,            color=CO_COLOR, ls="--", lw=1.8,
                  zorder=1, label=f"μ = +{mean_diff:.3f}")
    ax_st.axvline(mean_diff - sd_diff,  color="gray",   ls="--", lw=1.0, alpha=0.6)
    ax_st.axvline(mean_diff + sd_diff,  color="gray",   ls="--", lw=1.0, alpha=0.6)
    for y, d, col in zip(y_pos, diff_v, colors):
        ax_st.scatter(d, y, color=col, s=70, zorder=3, clip_on=False)
    ax_st.set_yticks(y_pos)
    ax_st.set_yticklabels(chroms, fontsize=8)
    ax_st.set_ylim(0.0, len(diff_v) + 0.8)
    ax_st.set_xlabel("COSMIC − gnomAD  mean ΔE  (non-CpG exons)", fontsize=10)
    ax_st.set_ylabel("Chromosome", fontsize=10)
    ax_st.set_title("Per-chromosome paired difference", fontsize=10)
    ax_st.legend(fontsize=9, frameon=False, loc="lower right")
    ax_st.text(0.97, 0.02, "all 22 chromosomes above zero",
               transform=ax_st.transAxes, ha="right", va="bottom",
               fontsize=8.5, style="italic", color="#444")
    ax_st.spines["top"].set_visible(False)
    ax_st.spines["right"].set_visible(False)

    fig.tight_layout()
    out_png = out_dir / "B_exon_noncpg_paired.png"
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {out_png}", flush=True)


# ---------------------------------------------------------------------------
# Figure C: standalone chromosome consistency strip plot
# ---------------------------------------------------------------------------

def plot_C(tbl: pd.DataFrame, out_dir: Path) -> None:
    """
    Strip plot standalone (pannello destro di B come figura dedicata).
    Mostra tutti 22 autosomi con differenza COSMIC−gnomAD positiva.
    Range valori tipici: +0.163 – +0.200
    """
    diff_v = tbl["diff"].values
    chroms = tbl["chrom"].values
    mean_diff = diff_v.mean()
    sd_diff   = diff_v.std(ddof=1)

    cmap   = plt.cm.tab20
    colors = [cmap(i % 20) for i in range(len(chroms))]

    fig, ax = plt.subplots(figsize=(7, 7))

    y_pos = np.arange(len(diff_v), 0, -1)

    ax.axvline(mean_diff,            color=CO_COLOR, ls="--", lw=2.0,
               zorder=1, label=f"μ = +{mean_diff:.3f}")
    ax.axvline(mean_diff - sd_diff,  color="gray",   ls="--", lw=1.2, alpha=0.5)
    ax.axvline(mean_diff + sd_diff,  color="gray",   ls="--", lw=1.2, alpha=0.5)
    ax.axvline(0, color="#aaa", ls="-", lw=0.8, alpha=0.3)

    for y, d, col, ch in zip(y_pos, diff_v, colors, chroms):
        ax.scatter(d, y, color=col, s=80, zorder=3, clip_on=False)
        ax.text(d + 0.001, y, f"{d:+.3f}", va="center", fontsize=7.5, color="#333")

    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"chr{c}" for c in chroms], fontsize=9)
    ax.set_ylim(0.0, len(diff_v) + 1)
    ax.set_xlabel("COSMIC − gnomAD  mean ΔE  (non-CpG exonic variants)", fontsize=11)
    ax.set_title(
        "Per-chromosome consistency: COSMIC somatic > gnomAD germline\n"
        "(non-CpG exonic variants, n=22 autosomi)",
        fontsize=11, fontweight="bold",
    )
    ax.legend(fontsize=10, frameon=False, loc="lower right")
    ax.text(0.97, 0.02, "all 22 chromosomes above zero",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=9, style="italic", color="#444")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    out_png = out_dir / "C_chrom_consistency.png"
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {out_png}", flush=True)


# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="CpG split + paired analysis (non-CpG exonic variants)"
    )
    parser.add_argument("--gnomad-parquet", default=DEFAULT_PARQ)
    parser.add_argument("--cosmic-tsv",     default=DEFAULT_TSV)
    parser.add_argument("-o", "--outdir",   default=DEFAULT_OUT)
    args = parser.parse_args()

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load
    gn_agg = load_gnomad(args.gnomad_parquet)
    co_agg = load_cosmic(args.cosmic_tsv)

    # Figure A — CpG split
    summary = build_summary(gn_agg, co_agg)
    print("\n=== SUMMARY ===")
    print(summary[["dataset", "cpg_class", "mean_dE", "n_variants"]].to_string(index=False))
    summary.to_csv(out_dir / "A_exon_cpg_split_summary.csv", index=False, float_format="%.6g")
    print(f"  ✓ {out_dir}/A_exon_cpg_split_summary.csv", flush=True)
    plot_A(gn_agg, co_agg, summary, out_dir)

    # Paired table (non-CpG only) — usata da B e C
    tbl = build_paired_table(gn_agg, co_agg)
    gn_v, co_v, diff_v = (tbl["gnomad_mean_dE"].values,
                           tbl["cosmic_mean_dE"].values,
                           tbl["diff"].values)
    t_stat, pval = ttest_rel(gn_v, co_v)
    cohens_d = diff_v.mean() / diff_v.std(ddof=1)

    # Aggiungi statistiche globali al CSV
    tbl["mean_diff_global"] = diff_v.mean()
    tbl["pval_paired_t"]    = pval
    tbl["cohens_d"]         = cohens_d
    tbl.to_csv(out_dir / "B_exon_noncpg_paired.csv", index=False, float_format="%.6f")
    print(f"  ✓ {out_dir}/B_exon_noncpg_paired.csv", flush=True)

    print(f"\n  Δμ (COSMIC−gnomAD) non-CpG exonic : {diff_v.mean():+.3f}")
    print(f"  p (paired t-test)                  : {pval:.2e}")
    print(f"  Cohen's d                          : {cohens_d:.1f}")
    print(f"  all 22 chroms positive?            : {(diff_v > 0).all()}")
    print(f"\n{tbl[['chrom','gnomad_mean_dE','cosmic_mean_dE','diff']].to_string(index=False)}")

    # Figure B — scatter + strip
    plot_B(tbl, out_dir)

    # Figure C — standalone strip (chromosome consistency)
    plot_C(tbl, out_dir)

    print("\n=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
