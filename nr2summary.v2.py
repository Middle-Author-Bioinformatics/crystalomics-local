#!/usr/bin/env python3
import argparse
import csv
import sys
from collections import OrderedDict
from typing import List

import pandas as pd


def detect_key_column(df: pd.DataFrame) -> str:
    """
    Try to infer which column in the summary CSV corresponds to the BLAST qseqid.
    Preference order: 'qseqid', 'query', 'peptide', first column.
    """
    candidates = ["qseqid", "query", "peptide"]
    for c in candidates:
        if c in df.columns:
            return c
    # fallback: first column
    return df.columns[0]


def read_blast_nr_table(path: str) -> pd.DataFrame:
    """
    Read DIAMOND/BLAST tabular (-outfmt 6) with the exact columns used in the user's command:
      qseqid sseqid pident length evalue bitscore stitle qseq sseq
    Extra columns will be retained if present; missing ones will raise a helpful error.
    """
    # Read as TSV without header, then assign if needed.
    # We'll try to detect if a header exists by checking for 'qseqid' in first line.
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        first = f.readline().strip().split("\t")
    has_header = first and first[0].lower() in {"qseqid", "#qseqid"}
    if has_header:
        df = pd.read_csv(path, sep="\t", dtype=str, quoting=csv.QUOTE_NONE, engine="python")
    else:
        df = pd.read_csv(path, sep="\t", header=None, dtype=str, quoting=csv.QUOTE_NONE, engine="python")
        # Assign the expected headers when column count matches 9; otherwise create generic names.
        expected = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "stitle", "qseq", "sseq"]
        if df.shape[1] >= len(expected):
            df.columns = expected + [f"extra_{i}" for i in range(df.shape[1] - len(expected))]
        else:
            # If fewer columns, pad with generic names to avoid crashes but mark missing expected cols.
            df.columns = [f"col{i+1}" for i in range(df.shape[1])]

    # Ensure key columns exist
    required = ["qseqid", "stitle"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"BLAST file '{path}' is missing required columns: {missing}. "
            "Expected at least 'qseqid' and 'stitle'. If your file includes headers, make sure they are present."
        )

    # Coerce bitscore to float for sorting if present
    if "bitscore" in df.columns:
        with pd.option_context('mode.use_inf_as_na', True):
            df["bitscore_num"] = pd.to_numeric(df["bitscore"], errors="coerce")
    else:
        df["bitscore_num"] = pd.NA

    return df


def collapse_nr_hits(df_nr: pd.DataFrame, n: int) -> pd.DataFrame:
    """
    For each qseqid, take the top-n hits by bitscore (descending; NaNs last) and join stitle(s).
    Returns a DataFrame with columns ['qseqid', 'nr_stitle'].
    """
    # Sort by bitscore desc (NaNs at the end by using na_position='last')
    if "bitscore_num" in df_nr.columns:
        df_nr_sorted = df_nr.sort_values(["qseqid", "bitscore_num"], ascending=[True, False], na_position="last")
    else:
        df_nr_sorted = df_nr.copy()

    # For stability, drop duplicate stitles per qseqid while preserving order
    def unique_in_order(series: pd.Series) -> List[str]:
        seen = OrderedDict()
        for s in series:
            if pd.isna(s):
                continue
            if s not in seen:
                seen[s] = True
        return list(seen.keys())

    grouped = df_nr_sorted.groupby("qseqid", sort=False)["stitle"].apply(unique_in_order).reset_index(name="titles")
    grouped["nr_stitle"] = grouped["titles"].apply(lambda lst: " || ".join(lst[:n]) if lst else "")
    result = grouped[["qseqid", "nr_stitle"]]
    return result


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Append NR 'stitle' information to an existing summary CSV.\n"
            "Reads a DIAMOND/BLAST outfmt 6 file against NR (with at least qseqid and stitle), "
            "collapses the top-N titles per query, and adds a new 'nr_stitle' column to the summary."
        )
    )
    ap.add_argument("-b1", "--blast-nr", required=True, help="Path to NR BLAST/DIAMOND tabular file (outfmt 6).")
    ap.add_argument("-b2", "--summary-csv", required=True, help="Path to the existing summary CSV to augment.")
    ap.add_argument("-o", "--out", required=True, help="Output CSV path for the augmented summary.")
    ap.add_argument("--nr-max", type=int, default=1, help="Maximum number of NR titles to include per query (default: 1).")
    ap.add_argument("--csv-delim", default=",", help="CSV delimiter for reading and writing the summary (default: ',').")
    ap.add_argument("--encoding", default="utf-8", help="Text encoding for I/O (default: utf-8).")

    args = ap.parse_args()

    # Load inputs
    df_summary = pd.read_csv(args.summary_csv, sep=args.csv_delim, dtype=str, encoding=args.encoding)
    key_col = detect_key_column(df_summary)

    df_nr = read_blast_nr_table(args.blast_nr)
    df_collapsed = collapse_nr_hits(df_nr, args.nr_max)

    # Merge
    df_out = df_summary.merge(df_collapsed, how="left", left_on=key_col, right_on="qseqid")
    if "qseqid" in df_out.columns and "qseqid" != key_col:
        df_out = df_out.drop(columns=["qseqid"])

    # Ensure the new column exists even if no matches
    if "nr_stitle" not in df_out.columns:
        df_out["nr_stitle"] = pd.NA

    # Write output
    df_out.to_csv(args.out, sep=args.csv_delim, index=False, encoding=args.encoding)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(f"[nr2summary.py] ERROR: {e}\n")
        sys.exit(1)