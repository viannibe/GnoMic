#!/usr/bin/env python3
"""
STEP 0: Build k-mer frequencies from reference genome

 Exclude chrM (mitochondrial chromosome)
 Canonical chromosomes only (chr1-22, chrX, chrY)
All 16,384 k-mers output for easy lookup

Output: hg38_k7_freq_withPos.csv
Columns: chrom, kmer, count, total_windows, freq

"""

from pathlib import Path
import csv
import pysam
import argparse
import time

print("="*70)
print("STEP 0: Build k-mer frequencies from reference genome")
print("="*70)

# ============================================================================
# PARAMETERS
# ============================================================================

K = 7
ALPH = "ACGT"
BASE2INT = {b: i for i, b in enumerate(ALPH)}

# DESIGN DECISION 1 & 5: Canonical chromosomes only (exclude chrM)
CANONICAL_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def int_to_kmer(idx: int, k: int) -> str:
    """
    Convert integer index to k-mer string (base-4 decoding)
    Example: 0 → AAAAAAA, 16383 → TTTTTTT
    """
    out = []
    for _ in range(k):
        out.append(ALPH[idx & 3])
        idx >>= 2
    return "".join(reversed(out))


def count_kmers_in_seq(seq: str, k: int):
    """
    Count k-mers using rolling hash (efficient!)

    DESIGN DECISION 4: Skips k-mers containing 'N'

    Returns: (counts_array, total_windows)
    - counts_array: length = 4^k with counts for each k-mer
    - total_windows: total valid k-mers counted
    """
    size = 4 ** k
    counts = [0] * size
    mask = (1 << (2 * k)) - 1

    rolling = 0
    valid_len = 0
    total = 0

    for ch in seq:
        ch = ch.upper()
        v = BASE2INT.get(ch, None)

        if v is None:  # DECISION 4: N or ambiguous base → reset
            rolling = 0
            valid_len = 0
            continue

        # Rolling hash: shift left 2 bits, add new base, mask to k bases
        rolling = ((rolling << 2) | v) & mask
        valid_len += 1

        if valid_len >= k:
            counts[rolling] += 1
            total += 1

    return counts, total


# ============================================================================
# MAIN
# ============================================================================

def main():
    ap = argparse.ArgumentParser(
        description="Build k-mer frequencies from reference genome"
    )
    ap.add_argument("--fasta", required=True,
                    help="Path to hg38 fasta (indexed with .fai)")
    ap.add_argument("--out", required=True,
                    help="Output CSV path")
    ap.add_argument("--include-contigs", action="store_true",
                    help="Include non-canonical contigs (NOT RECOMMENDED)")
    args = ap.parse_args()

    # Open FASTA
    fa = pysam.FastaFile(args.fasta)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Select chromosomes
    chroms = list(fa.references)
    if not args.include_contigs:
        chroms = [c for c in chroms if c in CANONICAL_CHROMS]
        print(f"\n✓ Using canonical chromosomes only: {len(chroms)}")
        print(f"  Excluded: chrM (mitochondrial)")
    else:
        print(f"\n Including all contigs: {len(chroms)}")

    print(f"\n Parameters:")
    print(f"  K-mer size: {K}")
    print(f"  Total possible k-mers: {4**K:,} (4^{K})")
    print(f"  FASTA: {args.fasta}")
    print(f"  Output: {out_path}")

    # Create output file with header
    print(f"\n Creating output file...")
    with out_path.open("w", newline="") as fout:
        w = csv.writer(fout)
        w.writerow(["chrom", "kmer", "count", "total_windows", "freq"])

    # Process each chromosome
    start_time = time.time()
    grand_total_kmers = 0

    print(f"\n Processing chromosomes...")
    print("─" * 70)

    for chrom in chroms:
        chrom_start = time.time()

        print(f"\n {chrom}")

        # Fetch entire chromosome sequence (forward strand)
        seq = fa.fetch(chrom)
        seq_len = len(seq)
        print(f"  Length: {seq_len:,} bp")

        # Count k-mers
        print(f"  Counting k-mers...", end='', flush=True)
        counts, total = count_kmers_in_seq(seq, K)

        chrom_elapsed = time.time() - chrom_start
        print(f" Done ({chrom_elapsed:.1f}s)")

        if total == 0:
            print(f"  No valid k-mers (all ambiguous bases)")
            continue

        # Statistics
        unique_present = sum(1 for x in counts if x > 0)
        print(f"  Total k-mers: {total:,}")
        print(f"  Unique observed: {unique_present:,} / {4**K:,}")

        grand_total_kmers += total

        # Write ALL 4^7 k-mers for this chromosome
        # (Even count=0, for easy lookup later)
        print(f"  Writing to CSV...", end='', flush=True)

        with out_path.open("a", newline="") as fout:
            w = csv.writer(fout)
            for idx, cnt in enumerate(counts):
                kmer = int_to_kmer(idx, K)
                freq = cnt / total if total > 0 else 0.0
                w.writerow([chrom, kmer, cnt, total, freq])

        print(f" Done")

    elapsed = time.time() - start_time

    # Summary
    print(f"\n{'='*70}")
    print("K-MER FREQUENCY CALCULATION COMPLETE")
    print("="*70)

    print(f"\n Summary:")
    print(f"  Chromosomes processed: {len(chroms)}")
    print(f"  K-mer size: {K}")
    print(f"  Total possible k-mers: {4**K:,}")
    print(f"  Total k-mers counted: {grand_total_kmers:,}")
    print(f"  Output rows: {len(chroms) * (4**K):,}")
    print(f"  Output file: {out_path}")
    print(f"  File size: {out_path.stat().st_size / 1e6:.1f} MB")
    print(f"  Total time: {elapsed/60:.1f} minutes")


    fa.close()


if __name__ == "__main__":
    main()
