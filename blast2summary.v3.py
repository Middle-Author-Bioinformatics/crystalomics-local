#!/usr/bin/env python3
import argparse
import csv
import os
from typing import Dict, Iterable
from Bio import SeqIO  # Requires Biopython

def _is_better_by_bitscore(new_row: dict, cur_row: dict) -> bool:
    """Return True if new_row is better (higher bitscore; tie-breaker = lower evalue)."""
    try:
        bn = float(new_row["bitscore"])
    except Exception:
        bn = float("-inf")
    try:
        bc = float(cur_row["bitscore"])
    except Exception:
        bc = float("-inf")

    if bn > bc:
        return True
    if bn < bc:
        return False

    # tie-breaker: smaller evalue wins
    try:
        en = float(new_row["evalue"])
        ec = float(cur_row["evalue"])
        return en < ec
    except Exception:
        return False

def parse_args():
    p = argparse.ArgumentParser(
        description="Summarize one or more BLAST (outfmt 6) files; include peptide sequence and BLAST file prefix. "
                    "Selects the best hit per query across all files, then optionally filters by minimum bitscore."
    )

    p.add_argument(
        "-b", "--blast",
        nargs="+",
        required=True,
        help="One or more tabular BLAST files (outfmt 6)."
    )

    p.add_argument(
        "-f", "--fasta", required=True,
        help="FASTA file of peptide query sequences (the BLAST query set)."
    )
    p.add_argument(
        "-o", "--output", required=True,
        help="Output CSV file."
    )

    # New: configurable filters
    p.add_argument(
        "--min-bitscore", type=float, default=None,
        help="If set, only write best-per-query hits with bitscore >= this value."
    )
    p.add_argument(
        "--max-evalue", type=float, default=1.0,
        help="Discard hits with e-value greater than this (default: 1.0)."
    )
    return p.parse_args()

def parse_float_safe(x: str) -> float:
    try:
        return float(x)
    except Exception:
        return float("inf")

def load_fasta_sequences(fasta_file: str) -> Dict[str, str]:
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_file, "fasta")}

def infer_blast_prefix(path: str) -> str:
    """
    Derive a readable prefix from a BLAST filename by stripping common suffixes.
    e.g., 'sample.ref.blast' -> 'sample.ref'
          'foo.tblastn'     -> 'foo'
          'bar.blastp'      -> 'bar'
    """
    base = os.path.basename(path)
    for suf in (".blast", ".blastp", ".tblastn", ".out", ".txt"):
        if base.endswith(suf):
            return base[: -len(suf)]
    return os.path.splitext(base)[0]

def make_subject_coords(sid: str, sstart_raw: str, send_raw: str) -> str:
    """Return 'sid_start_end' with start <= end; handles non-integer gracefully."""
    try:
        sstart = int(sstart_raw)
        send = int(send_raw)
        start, end = (sstart, send) if sstart <= send else (send, sstart)
        return f"{sid}_{start}_{end}"
    except Exception:
        return f"{sid}_{sstart_raw}_{send_raw}"

def parse_outfmt6_lines(lines: Iterable[str], seqs: Dict[str, str], blast_prefix: str, max_evalue: float):
    """
    Yield rows matching the output header:
      query_id, subject_id, subject_locus, subject_coords, ID, aln_len, evalue, peptide_seq, blast_prefix, bitscore
    """
    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            continue

        try:
            qid      = fields[0]
            sid      = fields[1]
            pident   = fields[2]   # kept as 'ID' for backward compatibility
            aln_len  = fields[3]
            sstart   = fields[8]
            send     = fields[9]
            evalue   = fields[10]
            bitscore = fields[11]
        except (IndexError, ValueError):
            continue

        evalue_val = parse_float_safe(evalue)
        if evalue_val > max_evalue:
            continue

        bitscore_val = parse_float_safe(bitscore)
        if bitscore_val == float("inf"):
            continue

        peptide_seq = seqs.get(qid, "NA")
        subject_coords = make_subject_coords(sid, sstart, send)

        yield qid, [
            qid,          # 0
            sid,          # 1
            sid,          # 2 subject_locus (kept as original sid unless you map differently)
            subject_coords,# 3
            pident,       # 4 ID
            aln_len,      # 5
            evalue,       # 6
            peptide_seq,  # 7
            blast_prefix, # 8
            bitscore_val  # 9 (used for ranking + filtering)
        ]

def better_hit(row_a: list, row_b: list) -> list:
    # last element is bitscore
    bits_a, bits_b = row_a[-1], row_b[-1]
    if bits_a != bits_b:
        return row_a if bits_a > bits_b else row_b

    # tie-breakers
    def _pf(x):
        try: return float(x)
        except Exception: return float("inf")

    e_a, e_b = _pf(row_a[6]), _pf(row_b[6])
    if e_a != e_b:
        return row_a if e_a < e_b else row_b

    try:
        pid_a, pid_b = float(row_a[4]), float(row_b[4])
    except Exception:
        pid_a = pid_b = 0.0
    if pid_a != pid_b:
        return row_a if pid_a > pid_b else row_b

    try:
        len_a, len_b = int(row_a[5]), int(row_b[5])
    except Exception:
        len_a = len_b = 0
    if len_a != len_b:
        return row_a if len_a > len_b else row_b

    return row_a

def main():
    args = parse_args()
    query_seqs = load_fasta_sequences(args.fasta)

    # Collect best hit per query across all files
    best_by_query = {}

    for bpath in args.blast:
        prefix = infer_blast_prefix(bpath)
        with open(bpath, "r") as bf:
            for qid, row in parse_outfmt6_lines(bf, query_seqs, prefix, args.max_evalue):
                if qid not in best_by_query:
                    best_by_query[qid] = row
                else:
                    best_by_query[qid] = better_hit(best_by_query[qid], row)

    # Optional global bitscore filtering on the best-per-query winners
    if args.min_bitscore is not None:
        best_by_query = {
            q: r for q, r in best_by_query.items()
            if r[-1] >= args.min_bitscore
        }

    # Write output (omit internal bitscore column; add it if you want)
    with open(args.output, "w", newline="") as outfh:
        writer = csv.writer(outfh)
        writer.writerow([
            "query_id",
            "subject_id",
            "subject_locus",
            "subject_coords",
            "ID",
            "aln_len",
            "evalue",
            "peptide_seq",
            "blast_prefix"
        ])
        for qid in sorted(best_by_query.keys()):
            row = best_by_query[qid]
            writer.writerow(row[:-1])  # drop bitscore

if __name__ == "__main__":
    main()






