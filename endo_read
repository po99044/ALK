"""
endo_read.py
FASTQ -> trimmed read counting -> variant-level count table builder.

Cleaned & CLI-friendly core functions derived from original Endo_read.py (custom code).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import pandas as pd


@dataclass
class BarcodeTrimSpec:
    """How to trim barcode/oligo sequence to match the amplicon target region."""
    left: int = 55
    right: int = 50


def read_fastq_trimmed_counts(
    fastq_path: str,
    ref_seq: str,
    fp_len: int = 10,
) -> Dict[str, int]:
    """
    Parse a FASTQ file and count trimmed reads.

    Logic (matches original):
      - Use FP = first fp_len bases of reference
      - For each read line (FASTQ line index % 4 == 1), find FP
      - Extract segment after FP of length len(ref_seq[fp_len:])
      - Count exact matches
      - Filter out indel-length reads (keep only exact expected length)
    """
    fp = ref_seq[:fp_len]
    expected_len = len(ref_seq[fp_len:])

    counts: Dict[str, int] = {}
    with open(fastq_path, "r") as f:
        for i, line in enumerate(f):
            if i % 4 != 1:
                continue
            seq = line.strip()
            pos = seq.find(fp)
            if pos == -1:
                continue
            trimmed = seq[pos + fp_len : pos + fp_len + expected_len]
            counts[trimmed] = counts.get(trimmed, 0) + 1

    # Remove indel reads by enforcing length
    filtered = {k: v for k, v in counts.items() if len(k) == expected_len}
    return filtered


def load_barcode_table(
    barcode_csv: str,
    trim_spec: BarcodeTrimSpec = BarcodeTrimSpec(),
    seq_col: str = "seq",
) -> pd.DataFrame:
    """
    Load barcode table (e.g., ALK_E23_aa.csv) and trim the seq field as in original code:
      seq = seq[55:-50]
    """
    df = pd.read_csv(barcode_csv)
    if seq_col not in df.columns:
        raise ValueError(f"Barcode CSV must contain column '{seq_col}'")

    df = df.copy()
    df[seq_col] = df[seq_col].astype(str).apply(
        lambda s: s[trim_spec.left : len(s) - trim_spec.right]
    )
    return df


def attach_readcounts_individual(
    counts: Dict[str, int],
    barcode_df: pd.DataFrame,
    ref_seq: str,
    fp_len: int = 10,
    seq_col: str = "seq",
    id_col: str = "ID",
    aa_col: str = "AA_Change",
) -> pd.DataFrame:
    """
    Map per-sequence counts onto barcode table and add WT row.
    This corresponds to conv_indi() in the original Endo_read.py.

    Returns a dataframe including:
      - original barcode columns
      - ReadCounts
      - WT row appended
    """
    df = barcode_df.copy()

    # Map counts
    df["ReadCounts"] = df[seq_col].str.upper().map(counts).fillna(0).astype(int)

    # Remove duplicates by ID+seq
    df["IDseq"] = df[id_col].astype(str) + df[seq_col].astype(str)
    df = df.drop_duplicates(subset=["IDseq"])

    # Append WT
    wt_seq = ref_seq[fp_len:]
    wt_count = int(counts.get(wt_seq, 0))
    wt_row = {id_col: "WT", seq_col: wt_seq, "ReadCounts": wt_count, "IDseq": "WT" + wt_seq}
    if aa_col in df.columns:
        wt_row[aa_col] = "WT"
    df = pd.concat([df, pd.DataFrame([wt_row])], ignore_index=True)

    return df


def build_count_matrix(
    samples: Dict[str, Dict[str, int]],
    barcode_df: pd.DataFrame,
    ref_seq: str,
    fp_len: int = 10,
) -> pd.DataFrame:
    """
    Create a wide count matrix with one row per variant ID (including WT) and one column per sample label.

    samples: dict of {label: counts_dict}, where counts_dict is output of read_fastq_trimmed_counts.
    """
    merged: Optional[pd.DataFrame] = None

    for label, dd in samples.items():
        one = attach_readcounts_individual(dd, barcode_df, ref_seq, fp_len=fp_len)
        keep_cols = ["ID", "AA_Change", "ReadCounts"] if "AA_Change" in one.columns else ["ID", "ReadCounts"]
        one = one[keep_cols].copy()
        one = one.rename(columns={"ReadCounts": label})

        if merged is None:
            merged = one
        else:
            # merge on ID (and AA_Change if available)
            on_cols = ["ID", "AA_Change"] if "AA_Change" in merged.columns and "AA_Change" in one.columns else ["ID"]
            merged = pd.merge(merged, one, on=on_cols, how="outer")

    if merged is None:
        raise ValueError("No samples provided.")
    merged = merged.fillna(0)
    # ensure integer counts for sample columns
    for c in merged.columns:
        if c not in ("ID", "AA_Change"):
            merged[c] = merged[c].astype(int)
    return merged
