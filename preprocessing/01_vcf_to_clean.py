#!/usr/bin/env python3
"""
STEP 1: VCF → Clean TSV (Robust Filtering)

Exclude chrM (mitochondrial)
SNVs only (no INDELs)
Only A/T/C/G (unambiguous bases)
Canonical chromosomes only
No self-mutations (REF ≠ ALT)

Output: cosmic_snv_clean.tsv
"""

from pathlib import Path
import argparse
import csv
import time

# ============================================================================
# DESIGN DECISIONS (CONSTANTS)
# ============================================================================

# Decision 1 & 5: Canonical chromosomes (exclude chrM)
CANONICAL_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

# Decision 3: Valid bases only
VALID_BASES = {'A', 'T', 'C', 'G'}

# ============================================================================
# HELPER FUNCTION
# ============================================================================

def parse_info(info_str: str) -> dict:
    """Parse VCF INFO field into dictionary"""
    d = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
        else:
            d[item] = True
    return d

# ============================================================================
# MAIN
# ============================================================================

def main():
    ap = argparse.ArgumentParser(description="STEP 1: VCF → Clean TSV (Robust Filtering)")
    ap.add_argument("--vcf", required=True, help="Input COSMIC VCF file")
    ap.add_argument("--out", required=True, help="Output clean TSV")
    args = ap.parse_args()

    vcf_file = Path(args.vcf)
    out_tsv = Path(args.out)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    print("="*70)
    print("STEP 1: VCF → TSV (Robust Filtering)")
    print("="*70)

    print(f"\n Filtering Rules:")
    print(f"  1. Canonical chromosomes: {len(CANONICAL_CHROMS)} (chr1-22, chrX, chrY)")
    print(f"  2. Exclude chrM: YES")
    print(f"  3. SNVs only: YES")
    print(f"  4. Valid bases: {VALID_BASES}")
    print(f"  5. No self-mutations: YES")
    print(f"  6. Require gene annotation: YES")

    print(f"\n Input:  {vcf_file}")
    print(f" Output: {out_tsv}")

    start_time = time.time()

    # Statistics
    stats = {
        'total': 0,
        'excluded_chrM': 0,
        'excluded_non_canonical': 0,
        'excluded_indel': 0,
        'excluded_invalid_base': 0,
        'excluded_self_mutation': 0,
        'excluded_no_gene': 0,
        'excluded_multiallelic': 0,
        'kept': 0
    }

    # Output columns
    output_cols = [
        "chrom", "pos", "ref", "alt", "cosmic_id",
        "gene", "so_term", "aa", "hgvsp",
        "sample_count", "transcript", "is_canonical",
        "gene_tier", "legacy_id"
    ]

    print(f"\n⏳ Parsing VCF...")

    with vcf_file.open("r") as fin, out_tsv.open("w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(output_cols)

        for line in fin:
            # Skip header
            if line.startswith("#"):
                continue

            stats['total'] += 1

            # Progress indicator
            if stats['total'] % 1_000_000 == 0:
                print(f"  Processed {stats['total']:,} variants...", end='\r')

            fields = line.rstrip("\n").split("\t")

            # Parse VCF columns
            chrom = fields[0]
            pos = int(fields[1])
            cosmic_id = fields[2]
            ref = fields[3].upper()
            alt = fields[4].upper()
            info_str = fields[7]

            # FILTER 1: Exclude chrM
            if chrom in ("MT", "chrMT", "chrM"):
                stats['excluded_chrM'] += 1
                continue

            # FILTER 2: Canonical chromosomes only
            if chrom not in CANONICAL_CHROMS:
                stats['excluded_non_canonical'] += 1
                continue

            # FILTER 3: No multi-allelic sites
            if "," in alt:
                stats['excluded_multiallelic'] += 1
                continue

            # FILTER 4: SNVs only - no INDELs
            if len(ref) != 1 or len(alt) != 1:
                stats['excluded_indel'] += 1
                continue

            # FILTER 5: Only A/T/C/G
            if ref not in VALID_BASES or alt not in VALID_BASES:
                stats['excluded_invalid_base'] += 1
                continue

            # FILTER 6: No self-mutations
            if ref == alt:
                stats['excluded_self_mutation'] += 1
                continue

            # Parse INFO field
            info = parse_info(info_str)

            # Extract fields
            gene = info.get("GENE", "")

            # FILTER 7: Require gene annotation
            if not gene or gene == ".":
                stats['excluded_no_gene'] += 1
                continue

            # Extract other INFO fields
            so_term = info.get("SO_TERM", "")
            aa = info.get("AA", "")
            hgvsp = info.get("HGVSP", "")
            sample_count = info.get("GENOME_SCREEN_SAMPLE_COUNT", "")
            transcript = info.get("TRANSCRIPT", "")
            is_canonical = info.get("IS_CANONICAL", "")
            gene_tier = info.get("TIER", "")
            legacy_id = info.get("LEGACY_ID", "")

            # Write row
            writer.writerow([
                chrom, pos, ref, alt, cosmic_id,
                gene, so_term, aa, hgvsp,
                sample_count, transcript, is_canonical,
                gene_tier, legacy_id
            ])

            stats['kept'] += 1

    elapsed = time.time() - start_time

    print(f"\n")  # Clear progress line

    # ============================================================================
    # SUMMARY STATISTICS
    # ============================================================================

    print(f"\n{'='*70}")
    print("FILTERING SUMMARY")
    print("="*70)

    print(f"\n Input:")
    print(f"  Total variants: {stats['total']:,}")

    print(f"\n Excluded:")
    print(f"  chrM (mitochondrial):     {stats['excluded_chrM']:8,} ({100*stats['excluded_chrM']/stats['total']:.2f}%)")
    print(f"  Non-canonical chr:        {stats['excluded_non_canonical']:8,} ({100*stats['excluded_non_canonical']/stats['total']:.2f}%)")
    print(f"  Multi-allelic sites:      {stats['excluded_multiallelic']:8,} ({100*stats['excluded_multiallelic']/stats['total']:.2f}%)")
    print(f"  INDELs:                   {stats['excluded_indel']:8,} ({100*stats['excluded_indel']/stats['total']:.2f}%)")
    print(f"  Invalid bases:            {stats['excluded_invalid_base']:8,} ({100*stats['excluded_invalid_base']/stats['total']:.2f}%)")
    print(f"  Self-mutations (REF=ALT): {stats['excluded_self_mutation']:8,} ({100*stats['excluded_self_mutation']/stats['total']:.2f}%)")
    print(f"  No gene annotation:       {stats['excluded_no_gene']:8,} ({100*stats['excluded_no_gene']/stats['total']:.2f}%)")

    total_excluded = stats['total'] - stats['kept']
    print(f"  {'─'*50}")
    print(f"  Total excluded:           {total_excluded:8,} ({100*total_excluded/stats['total']:.2f}%)")

    print(f"\n Output:")
    print(f"  Clean SNVs:               {stats['kept']:8,} ({100*stats['kept']/stats['total']:.2f}%)")

    print(f"\n Time: {elapsed/60:.1f} minutes")
    print(f"✓ Output: {out_tsv}")


if __name__ == "__main__":
    main()
