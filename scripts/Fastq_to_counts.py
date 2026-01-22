#!/usr/bin/env python
"""
01_fastq_to_counts.py
Build a variant-level count table from amplicon FASTQ files.

Example:
  python scripts/01_fastq_to_counts.py \
    --exon 23 \
    --barcode_csv data/ALK_E23_aa.csv \
    --sample D10_R1=/path/E23R1D10.extendedFrags.fastq \
    --sample D10_R2=/path/E23R2D10.extendedFrags.fastq \
    --sample REF=/path/H3122_E23.extendedFrags.fastq \
    --sample A=/path/23R2A.extendedFrags.fastq \
    --sample L=/path/23R2L.extendedFrags.fastq \
    --sample T=/path/23R2T.extendedFrags.fastq \
    --sample U=/path/23R2U.extendedFrags.fastq \
    --out results/ALK_E23_R2_D21_counts.csv
"""

from __future__ import annotations

import argparse
from typing import Dict, Tuple

import pandas as pd

from src import endo_ref
from src import endo_read


def parse_sample_kv(kv: str) -> Tuple[str, str]:
    if "=" not in kv:
        raise argparse.ArgumentTypeError("--sample must be in LABEL=PATH format")
    label, path = kv.split("=", 1)
    label = label.strip()
    path = path.strip()
    if not label:
        raise argparse.ArgumentTypeError("Empty LABEL in --sample")
    if not path:
        raise argparse.ArgumentTypeError("Empty PATH in --sample")
    return label, path


def main() -> None:
    ap = argparse.ArgumentParser(description="FASTQ -> variant-level count table")
    ap.add_argument("--exon", type=int, required=True, help="ALK exon number (20–28)")
    ap.add_argument("--barcode_csv", required=True, help="Barcode CSV (e.g., ALK_E23_aa.csv)")
    ap.add_argument("--out", required=True, help="Output CSV path")
    ap.add_argument("--fp_len", type=int, default=10, help="Forward primer seed length (default: 10)")
    ap.add_argument("--trim_left", type=int, default=55, help="Trim barcode seq left (default: 55)")
    ap.add_argument("--trim_right", type=int, default=50, help="Trim barcode seq right (default: 50)")
    ap.add_argument(
        "--sample",
        action="append",
        default=[],
        help="Sample in LABEL=FASTQ_PATH format (repeatable)",
    )

    args = ap.parse_args()

    if len(args.sample) == 0:
        raise SystemExit("ERROR: provide at least one --sample LABEL=FASTQ_PATH")

    ref_seq = endo_ref.exon(args.exon)

    barcode_df = endo_read.load_barcode_table(
        args.barcode_csv,
        trim_spec=endo_read.BarcodeTrimSpec(left=args.trim_left, right=args.trim_right),
        seq_col="seq",
    )

    # Build counts per sample
    sample_counts: Dict[str, Dict[str, int]] = {}
    for kv in args.sample:
        label, path = parse_sample_kv(kv)
        dd = endo_read.read_fastq_trimmed_counts(path, ref_seq, fp_len=args.fp_len)
        sample_counts[label] = dd

    # Create matrix
    mat = endo_read.build_count_matrix(sample_counts, barcode_df, ref_seq, fp_len=args.fp_len)

    # Save
    mat.to_csv(args.out, index=False)
    print(f"[OK] wrote: {args.out}  (rows={mat.shape[0]}, cols={mat.shape[1]})")


if __name__ == "__main__":
    main()
