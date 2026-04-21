#!/usr/bin/env python3
"""
STEP 2: Add 7-mer contexts (WT and MUT) + flank6 + explainable sequence features

Input : cosmic_snv_clean.tsv
Output: cosmic_snv_with_7mer.tsv

Adds:
- WT_7mer, MUT_7mer
- flank6 = left3 + right3
- WT_7mer_canon, MUT_7mer_canon (canonical by revcomp)
- mut_type = "REF>ALT"
Explainable features (from WT_7mer):
- wt_gc (GC fraction)
- wt_entropy (Shannon entropy, base-level)
- wt_max_homopolymer (max run length)
- is_CpG (central base is C and next base is G, or central is G and prev is C)
"""

from __future__ import annotations
from pathlib import Path
import argparse
import pandas as pd
import pysam
import time
import math

def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]

def canonical_kmer(seq: str) -> str:
    rc = revcomp(seq)
    return seq if seq <= rc else rc

def gc_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    seq = seq.upper()
    gc = sum(1 for c in seq if c in ("G", "C"))
    return gc / len(seq)

def shannon_entropy(seq: str) -> float:
    # entropy in bits (0..2 for DNA)
    if not seq:
        return 0.0
    seq = seq.upper()
    n = len(seq)
    counts = {}
    for c in seq:
        counts[c] = counts.get(c, 0) + 1
    ent = 0.0
    for c, k in counts.items():
        p = k / n
        ent -= p * math.log2(p)
    return ent

def max_homopolymer_run(seq: str) -> int:
    seq = seq.upper()
    best = 1
    cur = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            cur += 1
            best = max(best, cur)
        else:
            cur = 1
    return best if seq else 0

def main():
    ap = argparse.ArgumentParser(description="Add 7-mer context + explainable sequence features")
    ap.add_argument("--in", dest="inp", required=True, help="Input TSV from step1")
    ap.add_argument("--ref", required=True, help="Reference FASTA (indexed)")
    ap.add_argument("--out", required=True, help="Output TSV")
    ap.add_argument("--k", type=int, default=7, help="k-mer length (default 7)")
    ap.add_argument("--chunksize", type=int, default=200_000, help="Chunk size (default 200k)")
    args = ap.parse_args()

    K = args.k
    HALF = K // 2
    if K != 7:
        print("WARNING: scripts are designed for K=7 context, but will work for odd K>=3.")
        if K % 2 == 0:
            raise ValueError("K must be odd for centered SNV context")

    in_path = Path(args.inp)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("STEP 2: Add context + explainable sequence features")
    print("=" * 70)
    print(f"Input:  {in_path}")
    print(f"Ref:    {args.ref}")
    print(f"Output: {out_path}")
    print(f"K={K} HALF={HALF} chunksize={args.chunksize:,}")

    fa = pysam.FastaFile(args.ref)

    stats = {
        "total": 0,
        "excluded_wrong_chrom": 0,
        "excluded_boundary": 0,
        "excluded_ref_mismatch": 0,
        "excluded_N": 0,
        "kept": 0
    }

    first = True
    t0 = time.time()

    for chunk_i, chunk in enumerate(pd.read_csv(in_path, sep="\t", chunksize=args.chunksize, low_memory=False), 1):
        print(f"Chunk {chunk_i}: {len(chunk):,}")

        wt_list = []
        mut_list = []
        flank_list = []
        wt_canon_list = []
        mut_canon_list = []
        mut_type_list = []

        wt_gc_list = []
        wt_ent_list = []
        wt_hpoly_list = []
        is_cpg_list = []

        for chrom, pos, ref, alt in zip(
            chunk["chrom"].astype(str),
            chunk["pos"].astype(int),
            chunk["ref"].astype(str),
            chunk["alt"].astype(str)
        ):
            stats["total"] += 1
            ref = ref.upper()
            alt = alt.upper()

            if chrom not in fa.references:
                stats["excluded_wrong_chrom"] += 1
                wt_list.append(None); mut_list.append(None); flank_list.append(None)
                wt_canon_list.append(None); mut_canon_list.append(None); mut_type_list.append(None)
                wt_gc_list.append(None); wt_ent_list.append(None); wt_hpoly_list.append(None); is_cpg_list.append(None)
                continue

            pos0 = pos - 1
            start = pos0 - HALF
            end = pos0 + HALF + 1

            if start < 0 or end > fa.get_reference_length(chrom):
                stats["excluded_boundary"] += 1
                wt_list.append(None); mut_list.append(None); flank_list.append(None)
                wt_canon_list.append(None); mut_canon_list.append(None); mut_type_list.append(None)
                wt_gc_list.append(None); wt_ent_list.append(None); wt_hpoly_list.append(None); is_cpg_list.append(None)
                continue

            # verify REF
            ref_fa = fa.fetch(chrom, pos0, pos0 + 1).upper()
            if ref_fa != ref:
                stats["excluded_ref_mismatch"] += 1
                wt_list.append(None); mut_list.append(None); flank_list.append(None)
                wt_canon_list.append(None); mut_canon_list.append(None); mut_type_list.append(None)
                wt_gc_list.append(None); wt_ent_list.append(None); wt_hpoly_list.append(None); is_cpg_list.append(None)
                continue

            wt = fa.fetch(chrom, start, end).upper()
            if len(wt) != K:
                stats["excluded_boundary"] += 1
                wt_list.append(None); mut_list.append(None); flank_list.append(None)
                wt_canon_list.append(None); mut_canon_list.append(None); mut_type_list.append(None)
                wt_gc_list.append(None); wt_ent_list.append(None); wt_hpoly_list.append(None); is_cpg_list.append(None)
                continue

            if "N" in wt:
                stats["excluded_N"] += 1
                wt_list.append(None); mut_list.append(None); flank_list.append(None)
                wt_canon_list.append(None); mut_canon_list.append(None); mut_type_list.append(None)
                wt_gc_list.append(None); wt_ent_list.append(None); wt_hpoly_list.append(None); is_cpg_list.append(None)
                continue

            if wt[HALF] != ref:
                stats["excluded_ref_mismatch"] += 1
                wt_list.append(None); mut_list.append(None); flank_list.append(None)
                wt_canon_list.append(None); mut_canon_list.append(None); mut_type_list.append(None)
                wt_gc_list.append(None); wt_ent_list.append(None); wt_hpoly_list.append(None); is_cpg_list.append(None)
                continue

            mut = wt[:HALF] + alt + wt[HALF + 1:]
            flank6 = wt[:HALF] + wt[HALF + 1:]

            wt_list.append(wt)
            mut_list.append(mut)
            flank_list.append(flank6)
            wt_canon_list.append(canonical_kmer(wt))
            mut_canon_list.append(canonical_kmer(mut))
            mut_type_list.append(f"{ref}>{alt}")

            # explainable features
            wt_gc_list.append(gc_fraction(wt))
            wt_ent_list.append(shannon_entropy(wt))
            wt_hpoly_list.append(max_homopolymer_run(wt))
            # CpG around center (dinucleotide)
            is_cpg = 1 if ((wt[HALF] == "C" and wt[HALF+1] == "G") or (wt[HALF] == "G" and wt[HALF-1] == "C")) else 0
            is_cpg_list.append(is_cpg)

            stats["kept"] += 1

        chunk["WT_7mer"] = wt_list
        chunk["MUT_7mer"] = mut_list
        chunk["flank6"] = flank_list
        chunk["WT_7mer_canon"] = wt_canon_list
        chunk["MUT_7mer_canon"] = mut_canon_list
        chunk["mut_type"] = mut_type_list

        chunk["wt_gc"] = wt_gc_list
        chunk["wt_entropy"] = wt_ent_list
        chunk["wt_max_homopolymer"] = wt_hpoly_list
        chunk["is_cpg"] = is_cpg_list

        chunk = chunk[chunk["WT_7mer"].notna()]

        chunk.to_csv(out_path, sep="\t", index=False, mode="w" if first else "a", header=first)
        first = False

    print("=" * 70)
    print("STEP 2 SUMMARY")
    print("=" * 70)
    for k, v in stats.items():
        print(f"{k:>22}: {v:,}")
    print(f"Time: {(time.time()-t0)/60:.1f} min")

    fa.close()
    print(f"✓ Output: {out_path}")

if __name__ == "__main__":
    main()
