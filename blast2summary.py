#!/usr/bin/env python3

import argparse
import csv
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Summarize BLAST output with peptide and subject info.")
    parser.add_argument("-b", "--blast", required=True, help="BLAST output file (outfmt 6)")
    parser.add_argument("-f", "--fasta", required=True, help="FASTA file of peptide query sequences")
    parser.add_argument("-db", "--database", required=True, help="FASTA database used in BLAST (contains subject sequences)")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file")
    return parser.parse_args()

def load_fasta_sequences(fasta_file):
    return {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

def load_fasta_headers(fasta_file):
    return {record.id: record.description for record in SeqIO.parse(fasta_file, "fasta")}

def main():
    args = parse_args()

    query_seqs = load_fasta_sequences(args.fasta)
    subject_headers = load_fasta_headers(args.database)

    with open(args.blast) as blast_file, open(args.output, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["query_id", "subject_id", "ID", "aln_len", "evalue", "peptide_seq"])

        for line in blast_file:
            fields = line.strip().split("\t")
            if len(fields) < 2:
                continue
            query_id, subject_id = fields[0], fields[1]
            peptide_seq = query_seqs.get(query_id, "NA")
            evalue = fields[10]
            ID = fields[2]
            aln = fields[3]

            writer.writerow([query_id, subject_id, ID, aln, evalue, peptide_seq])

if __name__ == "__main__":
    main()
