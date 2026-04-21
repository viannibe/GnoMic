#!/usr/bin/env python3
"""
STEP 0: Build reference tables from genome FASTA

Outputs (genome-wide):
  1) k7_table.tsv
     columns: kmer7, count, freq
  2) k7_table_canonical.tsv
     columns: kmer7_canon, count, freq
  3) flank6_table.tsv
     columns: flank6, count, freq
  4) flank6_center_table.tsv
     columns: flank6, center, count7, total_flank, p_center_given_flank

- Counts are computed on forward strand of reference (FASTA). Canonical table
  collapses kmer with its reverse complement: canon = min(kmer, revcomp(kmer)).
- Skips windows containing non-ACGT (e.g., N) by resetting rolling hash.
- Excludes chrM by default, and uses canonical chromosomes only unless --include-contigs.
"""

from __future__ import annotations
from pathlib import Path
import argparse
import csv
import pysam
import time
from collections import defaultdict

ALPH = "ACGT"
BASE2INT = {b: i for i, b in enumerate(ALPH)}

CANONICAL_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]

def canonical_kmer(seq: str) -> str:
    rc = revcomp(seq)
    return seq if seq <= rc else rc

def int_to_kmer(idx: int, k: int) -> str:
    out = []
    for _ in range(k):
        out.append(ALPH[idx & 3])
        idx >>= 2
    return "".join(reversed(out))

def count_kmers_in_seq(seq: str, k: int) -> tuple[list[int], int]:
    size = 4 ** k
    counts = [0] * size
    mask = (1 << (2 * k)) - 1

    rolling = 0
    valid_len = 0
    total = 0

    for ch in seq:
        ch = ch.upper()
        v = BASE2INT.get(ch, None)
        if v is None:
            rolling = 0
            valid_len = 0
            continue

        rolling = ((rolling << 2) | v) & mask
        valid_len += 1
        if valid_len >= k:
            counts[rolling] += 1
            total += 1

    return counts, total

def main():
    ap = argparse.ArgumentParser(description="Build reference k-mer tables")
    ap.add_argument("--fasta", required=True, help="Path to reference FASTA (indexed with .fai)")
    ap.add_argument("--outdir", required=True, help="Output directory for tables")
    ap.add_argument("--k", type=int, default=7, help="k-mer length (default 7)")
    ap.add_argument("--include-contigs", action="store_true", help="Include non-canonical contigs")
    args = ap.parse_args()

    k = args.k
    if k % 2 == 0:
        print("WARNING: k is even; canonical collapsing still works but symmetry differs.")
    if k < 3:
        raise ValueError("k must be >= 3")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    out_k7 = outdir / "k7_table.tsv"
    out_k7_canon = outdir / "k7_table_canonical.tsv"
    out_flank = outdir / "flank6_table.tsv"
    out_flank_center = outdir / "flank6_center_table.tsv"

    print("=" * 70)
    print("STEP 0: Build reference tables")
    print("=" * 70)
    print(f"FASTA:  {args.fasta}")
    print(f"OUTDIR: {outdir}")
    print(f"k:      {k}")

    fa = pysam.FastaFile(args.fasta)
    chroms = list(fa.references)

    if not args.include_contigs:
        chroms = [c for c in chroms if c in CANONICAL_CHROMS]
        chroms = [c for c in chroms if c not in ("chrM", "chrMT", "MT")]
        print(f"Using canonical chromosomes only: {len(chroms)}")
    else:
        print(f"Including all contigs: {len(chroms)}")

    size = 4 ** k
    global_counts = [0] * size
    global_total = 0

    start = time.time()
    for chrom in chroms:
        t0 = time.time()
        seq = fa.fetch(chrom)
        counts, total = count_kmers_in_seq(seq, k)

        # accumulate
        for i, c in enumerate(counts):
            global_counts[i] += c
        global_total += total

        print(f"{chrom:>5}  len={len(seq):,}  windows={total:,}  time={time.time()-t0:.1f}s")

    print(f"Genome-wide total k-mers counted: {global_total:,}")
    print(f"Elapsed: {time.time()-start:.1f}s")

    # Build k7 freq dict
    k7_count_dict: dict[str, int] = {}
    for i, c in enumerate(global_counts):
        kmer = int_to_kmer(i, k)
        k7_count_dict[kmer] = c

    # Canonical collapsing
    canon_counts: dict[str, int] = defaultdict(int)
    for kmer, c in k7_count_dict.items():
        canon_counts[canonical_kmer(kmer)] += c

    # flank6 counts and conditional center counts derived from 7-mer counts
    flank_counts: dict[str, int] = defaultdict(int)
    flank_center_counts: dict[tuple[str, str], int] = defaultdict(int)

    for kmer, c in k7_count_dict.items():
        if c == 0:
            continue
        flank6 = kmer[:3] + kmer[4:]
        center = kmer[3]
        flank_counts[flank6] += c
        flank_center_counts[(flank6, center)] += c

    # Write tables
    print("Writing outputs...")

    with out_k7.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["kmer7", "count", "freq"])
        for kmer in sorted(k7_count_dict.keys()):
            c = k7_count_dict[kmer]
            freq = c / global_total if global_total else 0.0
            w.writerow([kmer, c, f"{freq:.12e}"])

    canon_total = sum(canon_counts.values())
    with out_k7_canon.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["kmer7_canon", "count", "freq"])
        for kmer in sorted(canon_counts.keys()):
            c = canon_counts[kmer]
            freq = c / canon_total if canon_total else 0.0
            w.writerow([kmer, c, f"{freq:.12e}"])

    flank_total = sum(flank_counts.values())
    with out_flank.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["flank6", "count", "freq"])
        for flank6 in sorted(flank_counts.keys()):
            c = flank_counts[flank6]
            freq = c / flank_total if flank_total else 0.0
            w.writerow([flank6, c, f"{freq:.12e}"])

    with out_flank_center.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["flank6", "center", "count7", "total_flank", "p_center_given_flank"])
        for flank6 in sorted(flank_counts.keys()):
            total_f = flank_counts[flank6]
            for center in "ACGT":
                c = flank_center_counts.get((flank6, center), 0)
                p = c / total_f if total_f else 0.0
                w.writerow([flank6, center, c, total_f, f"{p:.12e}"])

    print(f"✓ Wrote: {out_k7}")
    print(f"✓ Wrote: {out_k7_canon}")
    print(f"✓ Wrote: {out_flank}")
    print(f"✓ Wrote: {out_flank_center}")

    fa.close()

if __name__ == "__main__":
    main()
