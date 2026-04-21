#!/usr/bin/env python3
"""
STEP 3: Calculate k-mer energies - GENOME-WIDE aggregation

Aggregate per-chromosome frequencies to genome-wide!

Formula:
  E = -log(freq_genome_wide)
  delta_E = E_WT - E_MUT

Output: cosmic_snv_with_energies.tsv
Adds: E_WT, E_MUT, delta_E, abs_delta_E
"""

from pathlib import Path
import argparse
import pandas as pd
import numpy as np
import time

print("="*70)
print("STEP 3: Calculate energies (genome-wide)")
print("="*70)

def main():
    ap = argparse.ArgumentParser(description="STEP 3: Calculate k-mer energies")
    ap.add_argument("--in", dest="inp", required=True, help="Input TSV from step2 (with 7-mers)")
    ap.add_argument("--kmer-freq", required=True, help="K-mer frequency table (genome-wide aggregated)")
    ap.add_argument("--out", required=True, help="Output TSV with energies")
    ap.add_argument("--chunksize", type=int, default=200_000, help="Chunk size (default 200k)")
    args = ap.parse_args()

    # ============================================================================
    # PARAMETERS
    # ============================================================================

    PSEUDOCOUNT = 1e-12  # For k-mers never observed
    CHUNKSIZE = args.chunksize

    inp_path = Path(args.inp)
    freq_path = Path(args.kmer_freq)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"\n Input:     {inp_path}")
    print(f" Freq:      {freq_path}")
    print(f" Output:    {out_path}")
    print(f"\n Parameters:")
    print(f"  Pseudocount: {PSEUDOCOUNT:.0e} (for unseen k-mers)")
    print(f"  Chunk size: {CHUNKSIZE:,}")

    # ============================================================================
    # LOAD AND AGGREGATE K-MER FREQUENCIES
    # ============================================================================

    print(f"\n Loading per-chromosome k-mer frequencies: {freq_path}")
    freq_df = pd.read_csv(freq_path)
    print(f"  ✓ Loaded {len(freq_df):,} entries")

    # CRITICAL: Aggregate to GENOME-WIDE frequencies!
    print(f"\n Aggregating to genome-wide frequencies...")

    # Step 1: Sum counts across all chromosomes for each k-mer
    print(f"  → Grouping by k-mer...")
    kmer_grouped = freq_df.groupby('kmer')['count'].sum()

    # Step 2: Calculate genome-wide frequency
    print(f"  → Calculating genome-wide frequencies...")
    total_kmers_genome = kmer_grouped.sum()
    kmer_freq_dict = (kmer_grouped / total_kmers_genome).to_dict()

    print(f"  ✓ Genome-wide frequency dict: {len(kmer_freq_dict):,} unique k-mers")
    print(f"  ✓ Total k-mers genome-wide: {total_kmers_genome:,}")

    # Check coverage
    expected = 4**7
    observed = len(kmer_freq_dict)
    print(f"  ✓ Expected k-mers (4^7): {expected:,}")
    print(f"  ✓ Observed k-mers: {observed:,}")

    if observed < expected:
        missing = expected - observed
        print(f"  {missing} k-mers never observed (will use pseudocount)")

    # ============================================================================
    # CALCULATE ENERGIES (CHUNKED)
    # ============================================================================

    start_time = time.time()

    stats = {
        'total': 0,
        'wt_unseen': 0,
        'mut_unseen': 0
    }

    first_chunk = True

    print(f"\n Calculating energies in chunks...")

    for chunk_num, chunk in enumerate(pd.read_csv(inp_path, sep="\t", chunksize=CHUNKSIZE, low_memory=False), 1):

        print(f" Chunk {chunk_num}: {len(chunk):,} mutations...", end='', flush=True)

        wt_kmers = chunk["WT_7mer"].astype(str).values
        mut_kmers = chunk["MUT_7mer"].astype(str).values

        # Get frequencies (with pseudocount for unseen)
        wt_freq = np.array([kmer_freq_dict.get(k, PSEUDOCOUNT) for k in wt_kmers], dtype=float)
        mut_freq = np.array([kmer_freq_dict.get(k, PSEUDOCOUNT) for k in mut_kmers], dtype=float)

        # Track unseen k-mers
        stats['total'] += len(chunk)
        stats['wt_unseen'] += np.sum(wt_freq == PSEUDOCOUNT)
        stats['mut_unseen'] += np.sum(mut_freq == PSEUDOCOUNT)

        # Calculate energies: E = -log(frequency)
        E_WT = -np.log(wt_freq)
        E_MUT = -np.log(mut_freq)

        # Delta E = E_WT - E_MUT
        delta_E = E_WT - E_MUT

        # Add to chunk
        chunk["E_WT"] = E_WT
        chunk["E_MUT"] = E_MUT
        chunk["delta_E"] = delta_E
        chunk["abs_delta_E"] = np.abs(delta_E)

        # Write to output
        chunk.to_csv(
            out_path,
            sep="\t",
            index=False,
            mode="w" if first_chunk else "a",
            header=first_chunk
        )

        first_chunk = False

        print(f" Done")

    elapsed = time.time() - start_time

    # ============================================================================
    # LOAD FULL DATASET FOR SUMMARY (Sample if too large)
    # ============================================================================

    print(f"\n Loading output for summary statistics...")

    try:
        # Try loading full dataset
        cosmic_final = pd.read_csv(out_path, sep="\t")
    except:
        # If too large, sample
        print(f"  (Dataset too large, sampling 1M rows)")
        cosmic_final = pd.read_csv(out_path, sep="\t", nrows=1_000_000)

    # ============================================================================
    # SUMMARY STATISTICS
    # ============================================================================

    print(f"\n{'='*70}")
    print("ENERGY CALCULATION SUMMARY")
    print("="*70)

    print(f"\n Processed:")
    print(f"  Total mutations: {stats['total']:,}")
    print(f"  WT unseen k-mers:  {stats['wt_unseen']:8,} ({100*stats['wt_unseen']/stats['total']:.3f}%)")
    print(f"  MUT unseen k-mers: {stats['mut_unseen']:8,} ({100*stats['mut_unseen']/stats['total']:.3f}%)")

    print(f"\n Energy Statistics (from sample):")

    print(f"\n  E_WT (Wild-Type Energy):")
    print(f"    Mean:   {cosmic_final['E_WT'].mean():.4f}")
    print(f"    Median: {cosmic_final['E_WT'].median():.4f}")
    print(f"    Std:    {cosmic_final['E_WT'].std():.4f}")
    print(f"    Range:  [{cosmic_final['E_WT'].min():.4f}, {cosmic_final['E_WT'].max():.4f}]")

    print(f"\n  E_MUT (Mutant Energy):")
    print(f"    Mean:   {cosmic_final['E_MUT'].mean():.4f}")
    print(f"    Median: {cosmic_final['E_MUT'].median():.4f}")
    print(f"    Std:    {cosmic_final['E_MUT'].std():.4f}")
    print(f"    Range:  [{cosmic_final['E_MUT'].min():.4f}, {cosmic_final['E_MUT'].max():.4f}]")

    print(f"\n  ΔE (Delta Energy = E_WT - E_MUT):")
    print(f"    Mean:   {cosmic_final['delta_E'].mean():+.4f}")
    print(f"    Median: {cosmic_final['delta_E'].median():+.4f}")
    print(f"    Std:    {cosmic_final['delta_E'].std():.4f}")
    print(f"    Range:  [{cosmic_final['delta_E'].min():+.4f}, {cosmic_final['delta_E'].max():+.4f}]")

    # Interpretation
    mean_delta = cosmic_final['delta_E'].mean()
    if mean_delta > 0.01:
        print(f"\n  ✓ Positive mean ΔE: mutations tend to move to MORE common contexts")
    elif mean_delta < -0.01:
        print(f"\n    Negative mean ΔE: mutations tend to move to RARER contexts")
    else:
        print(f"\n    Mean ΔE ≈ 0: no clear directional bias")

    print(f"\n Time: {elapsed/60:.1f} minutes")

    print(f"\n{'='*70}")
    print("✓ STEP 3 COMPLETE")
    print("="*70)
    print(f"Output: {out_path}")


if __name__ == "__main__":
    main()
