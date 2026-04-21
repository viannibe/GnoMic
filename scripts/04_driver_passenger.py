#!/usr/bin/env python3
"""
04_driver_passenger.py
----------------------
Stratification control: il segnale ΔE COSMIC > gnomAD è genome-wide o
è dovuto all'ascertainment bias sui cancer driver genes?

Tre pannelli (solo varianti non-CpG):
  A  COSMIC driver    (in_cgc == 1, is_coding == 'coding')  vs gnomAD esonico
  B  COSMIC passenger (in_cgc == 0, is_coding == 'coding')  vs gnomAD esonico
  C  COSMIC non-coding (is_coding == 'non_coding')           vs gnomAD extragenico
     [controllo negativo]

Metodo (identico a gnomAD_vi_clean 14_driver_passenger_control.py):
  - Aggregazione per cromosoma (n=22)
  - Paired t-test + Cohen's d su diff = COSMIC − gnomAD

Input:
  --gnomad-parquet   results/cache/gnomad_corrected.parquet
  --cosmic-tsv       cosmic_final.tsv   (ha in_cgc per driver/passenger)

Output:
  results/driver_passenger/A_driver_passenger_kde.png
  results/driver_passenger/A_driver_passenger_summary.csv
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
DEFAULT_OUT  = str(PROJECT_DIR / "results/driver_passenger")
CHUNKSIZE    = 500_000
AUTOSOMES    = {str(i) for i in range(1, 23)}

GN_COLOR = "#1d3557"
CO_COLOR = "#e63946"
BW       = 0.8


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_gnomad(parquet_path: str) -> dict:
    """
    Ritorna dict con due DataFrame (non-CpG, autosomici, aggregati per chrom):
      'exon'       → usato per confronto con driver/passenger
      'extragenic' → usato per confronto con non-coding (controllo)
    """
    print(f"[gnomAD] reading {parquet_path}", flush=True)
    df = pd.read_parquet(
        parquet_path,
        columns=["CHROM", "E_ref", "E_mut", "region", "is_cpg"],
    )
    df["dE"]     = df["E_ref"] - df["E_mut"]
    df["CHROM"]  = df["CHROM"].astype(str).str.replace("chr", "", regex=False)
    df["is_cpg"] = pd.to_numeric(df["is_cpg"], errors="coerce").fillna(0).astype(int)
    df = df[(df["is_cpg"] == 0) & df["CHROM"].isin(AUTOSOMES)]

    result = {}
    for reg in ("exon", "extragenic"):
        sub = df[df["region"] == reg]
        agg = (
            sub.groupby("CHROM")["dE"]
            .agg(mean_dE="mean", n="count")
            .reset_index()
            .rename(columns={"CHROM": "chrom"})
        )
        result[reg] = agg
        print(f"  gnomAD {reg}: {len(sub):,} righe, {len(agg)} chroms", flush=True)
    return result


def load_cosmic(tsv_path: str) -> dict:
    """
    Ritorna dict con tre DataFrame (non-CpG, autosomici, aggregati per chrom):
      'driver'    → in_cgc==1, coding
      'passenger' → in_cgc==0, coding
      'noncoding' → is_coding=='non_coding'
    """
    print(f"[COSMIC] reading {tsv_path}", flush=True)
    want = ["chrom", "E_wt7", "E_mut7", "is_coding", "is_cpg", "in_cgc"]
    acc  = {"driver": [], "passenger": [], "noncoding": []}

    for chunk in pd.read_csv(
        tsv_path, sep="\t", chunksize=CHUNKSIZE,
        usecols=want,
        dtype={"chrom": str, "is_coding": str, "in_cgc": "Int64", "is_cpg": "Int64"},
    ):
        chunk["chrom"]  = chunk["chrom"].str.replace("chr", "", regex=False)
        chunk["is_cpg"] = pd.to_numeric(chunk["is_cpg"], errors="coerce").fillna(0).astype(int)
        chunk["in_cgc"] = pd.to_numeric(chunk["in_cgc"], errors="coerce").fillna(0).astype(int)
        chunk = chunk[(chunk["is_cpg"] == 0) & chunk["chrom"].isin(AUTOSOMES)].copy()
        if chunk.empty:
            continue

        chunk["dE"] = (
            pd.to_numeric(chunk["E_wt7"],  errors="coerce")
            - pd.to_numeric(chunk["E_mut7"], errors="coerce")
        )
        chunk = chunk.dropna(subset=["dE"])

        coding     = chunk[chunk["is_coding"] == "coding"]
        non_coding = chunk[chunk["is_coding"] == "non_coding"]

        for key, sub in [
            ("driver",    coding[coding["in_cgc"] == 1]),
            ("passenger", coding[coding["in_cgc"] == 0]),
            ("noncoding", non_coding),
        ]:
            if sub.empty:
                continue
            acc[key].append(
                sub.groupby("chrom")["dE"]
                .agg(sum_dE="sum", n="count")
                .reset_index()
            )

    result = {}
    for key, records in acc.items():
        if not records:
            result[key] = pd.DataFrame(columns=["chrom", "mean_dE", "n"])
            continue
        raw = pd.concat(records, ignore_index=True)
        agg = (
            raw.groupby("chrom")
            .agg(sum_dE=("sum_dE", "sum"), n=("n", "sum"))
            .reset_index()
        )
        agg["mean_dE"] = agg["sum_dE"] / agg["n"]
        result[key] = agg[["chrom", "mean_dE", "n"]]
        print(f"  COSMIC {key}: {int(agg['n'].sum()):,} varianti, {len(agg)} chroms", flush=True)
    return result


# ---------------------------------------------------------------------------
# Stats
# ---------------------------------------------------------------------------

def paired_stats(gn_ser: pd.Series, co_ser: pd.Series) -> dict:
    shared = sorted(gn_ser.index.intersection(co_ser.index), key=int)
    if len(shared) < 3:
        return {}
    a = gn_ser[shared].values
    b = co_ser[shared].values
    diff = b - a
    t, p = ttest_rel(a, b)
    d    = diff.mean() / diff.std(ddof=1) if diff.std(ddof=1) > 0 else np.nan
    return {
        "n_chrom":   len(shared),
        "t":         round(t, 3),
        "p":         p,
        "cohens_d":  round(d, 3),
        "mean_diff": round(diff.mean(), 5),
    }


# ---------------------------------------------------------------------------
# KDE panel helper
# ---------------------------------------------------------------------------

def _kde_panel(ax, gn_df, co_df, title):
    gn_ser = gn_df.set_index("chrom")["mean_dE"]
    co_ser = co_df.set_index("chrom")["mean_dE"] if not co_df.empty else pd.Series(dtype=float)

    gn_vals = gn_ser.values
    co_vals = co_ser.values if len(co_ser) else np.array([])

    stats = paired_stats(gn_ser, co_ser) if len(co_vals) else {}

    for vals, color, lab, n_var in [
        (gn_vals, GN_COLOR, "gnomAD germline", int(gn_df["n"].sum())),
        (co_vals, CO_COLOR, "COSMIC somatic",  int(co_df["n"].sum()) if not co_df.empty else 0),
    ]:
        if len(vals) < 3:
            continue
        mu = vals.mean()
        kde = gaussian_kde(vals, bw_method="silverman")
        kde.set_bandwidth(kde.factor * BW)
        x = np.linspace(vals.min() - 0.12, vals.max() + 0.12, 500)
        ax.plot(x, kde(x), color=color, lw=2.5,
                label=f"{lab}  μ = {mu:+.3f}  (n = {n_var:,})")
        ax.fill_between(x, kde(x), alpha=0.25, color=color)
        ax.axvline(mu, color=color, lw=1.2, ls="--", alpha=0.8)

    ax.axvline(0, color="gray", lw=1.0, ls="--", alpha=0.4, zorder=0)

    subtitle = ""
    if stats:
        subtitle = f"p = {stats['p']:.2e}  |  d = {stats['cohens_d']:.1f}"
    ax.set_title(f"{title}\n{subtitle}", fontsize=10)
    ax.set_xlabel("Mean ΔE per chromosome  (E_ref − E_mut)", fontsize=9)
    ax.set_ylabel("Density", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return stats, gn_vals, co_vals


# ---------------------------------------------------------------------------
# Main plot
# ---------------------------------------------------------------------------

def plot(gn: dict, co: dict, out_dir: Path) -> pd.DataFrame:
    panels = [
        # (co_key, gn_key, title, group_label)
        ("driver",    "exon",       "Driver genes  (CGC,  in_cgc = 1)",   "driver"),
        ("passenger", "exon",       "Passenger genes  (in_cgc = 0)",       "passenger"),
        ("noncoding", "extragenic", "Non-coding / extragenic  [control]",  "noncoding"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    summary_rows = []

    for ax, (lbl, (co_key, gn_key, title, grp)) in zip(
        axes,
        zip(["A", "B", "C"], panels),
    ):
        stats, gn_vals, co_vals = _kde_panel(ax, gn[gn_key], co[co_key], title)
        ax.text(-0.08, 1.05, lbl, transform=ax.transAxes,
                fontsize=13, fontweight="bold", va="top")

        for ds, vals, df_ref in [
            ("gnomAD", gn_vals, gn[gn_key]),
            ("COSMIC", co_vals, co[co_key]),
        ]:
            summary_rows.append({
                "group":      grp,
                "dataset":    ds,
                "n_variants": int(df_ref["n"].sum()) if not df_ref.empty else 0,
                "n_chrom":    len(vals),
                "mean_dE":    round(float(vals.mean()), 5) if len(vals) else np.nan,
                "sem_dE":     round(float(vals.std(ddof=1) / np.sqrt(len(vals))), 6)
                              if len(vals) > 1 else np.nan,
                "t":          stats.get("t", np.nan),
                "p":          stats.get("p", np.nan),
                "cohens_d":   stats.get("cohens_d", np.nan),
                "mean_diff_cosmic_minus_gnomad": stats.get("mean_diff", np.nan),
            })

    fig.suptitle(
        "Driver vs Passenger vs Non-coding: COSMIC somatic vs gnomAD germline\n"
        "(non-CpG only, per-chromosome means)",
        fontsize=12, fontweight="bold",
    )
    fig.tight_layout(w_pad=3.0)

    out_png = out_dir / "A_driver_passenger_kde.png"
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {out_png}", flush=True)

    summary = pd.DataFrame(summary_rows)
    out_csv = out_dir / "A_driver_passenger_summary.csv"
    summary.to_csv(out_csv, index=False, float_format="%.6g")
    print(f"  ✓ {out_csv}", flush=True)
    return summary


# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Driver/passenger/noncoding stratification (non-CpG)"
    )
    parser.add_argument("--gnomad-parquet", default=DEFAULT_PARQ)
    parser.add_argument("--cosmic-tsv",     default=DEFAULT_TSV)
    parser.add_argument("-o", "--outdir",   default=DEFAULT_OUT)
    args = parser.parse_args()

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gn = load_gnomad(args.gnomad_parquet)
    co = load_cosmic(args.cosmic_tsv)

    summary = plot(gn, co, out_dir)

    print("\n=== SUMMARY ===")
    print(summary[["group", "dataset", "mean_dE", "p", "cohens_d"]].to_string(index=False))
    print("=== DONE ===", flush=True)


if __name__ == "__main__":
    main()
