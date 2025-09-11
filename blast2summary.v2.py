#!/usr/bin/env python3
import argparse
import csv
import os
from typing import Dict, Iterable
from Bio import SeqIO  # Requires Biopython

def parse_args():
    p = argparse.ArgumentParser(
        description="Summarize one or more BLAST (outfmt 6) files; include peptide sequence and BLAST file prefix."
    )

    p.add_argument(
        "-b", "--blast",
        nargs="+",  # <— accept one or more BLAST files
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
    # Strip known trailing extensions in order of specificity
    for suf in (".blast", ".blastp", ".tblastn", ".out", ".txt"):
        if base.endswith(suf):
            return base[: -len(suf)]
    # Fallback: strip only the final extension
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

def parse_outfmt6_lines(lines: Iterable[str], seqs: Dict[str, str], blast_prefix: str):
    """
    Yield rows matching the output header:
      query_id, subject_id, ID, aln_len, evalue, peptide_seq, blast_prefix
    'ID' is preserved from your original script (percent identity).
    """
    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            # Not a standard outfmt 6 line; skip defensively
            continue
        try:
            qid      = fields[0]
            sid      = fields[1]
            pident   = fields[2]   # keep column name "ID" for backward compat
            aln_len  = fields[3]
            evalue   = fields[10]
        except (IndexError, ValueError):
            continue
        evalue_val = parse_float_safe(evalue)
        if evalue_val > 1.0:
            continue  # skip hits with evalue > 1

        bitscore_val = parse_float_safe(fields[11])
        if bitscore_val == float("inf"):
            continue

        peptide_seq = seqs.get(qid, "NA")
        subject_locus = sid
        sstart = fields[8]
        send = fields[9]
        subject_coords = make_subject_coords(sid, sstart, send)

        yield qid, [
            qid,
            sid,
            subject_locus,
            subject_coords,
            pident,
            aln_len,
            evalue,
            peptide_seq,
            blast_prefix,
            bitscore_val  # keep bitscore at end for ranking
        ]

def better_hit(row_a: list, row_b: list) -> list:
    # last element is bitscore
    bits_a, bits_b = row_a[-1], row_b[-1]
    if bits_a != bits_b:
        return row_a if bits_a > bits_b else row_b

    e_a, e_b = parse_float_safe(row_a[6]), parse_float_safe(row_b[6])
    if e_a != e_b:
        return row_a if e_a < e_b else row_b

    try:
        pid_a, pid_b = float(row_a[4]), float(row_b[4])
    except Exception:
        pid_a = pid_b = 0
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

        best_by_query = {}

        for bpath in args.blast:
            prefix = infer_blast_prefix(bpath)
            with open(bpath, "r") as bf:
                for qid, row in parse_outfmt6_lines(bf, query_seqs, prefix):
                    if qid not in best_by_query:
                        best_by_query[qid] = row
                    else:
                        best_by_query[qid] = better_hit(best_by_query[qid], row)

        for qid in sorted(best_by_query.keys()):
            writer.writerow(best_by_query[qid][:-1])  # drop bitscore before writing

if __name__ == "__main__":
    main()





