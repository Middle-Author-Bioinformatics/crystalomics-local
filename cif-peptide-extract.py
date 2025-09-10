#!/usr/bin/env python3
import gemmi
import textwrap
import argparse
import sys


parser = argparse.ArgumentParser(
    prog="cif-peptide-extract.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    *******************************************************
    Developed by Arkadiy Garber;
    *******************************************************
    '''))

parser.add_argument('-cif', type=str, help='input .cif file', default="bowtie2")

parser.add_argument('-faa', type=str, help='output peptide predictions in FASTA format', default="")

parser.add_argument('-txt', type=str, help='output text-file with peptide predictions for each chain', default="NA")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]

# Define a dictionary for three-letter to one-letter amino acid conversion
three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

# Load the CIF file
model = gemmi.read_structure(args.cif)
if args.txt != "NA":
    out = open(args.txt, "w")
out2 = open(args.faa, "w")

# Extract sequences
for chain in model[0]:  # Access the first model in the structure
    sequence = ''.join(three_to_one.get(res.name, 'X') for res in chain if res.name in three_to_one)
    if len(sequence) > 5:
        # print(f"Chain {chain.name}: {sequence}")
        out2.write(">Chain_" + chain.name + "\n")
        out2.write(sequence + "\n")
    if args.txt != "NA":
        out.write(f"Chain {chain.name}: {sequence}\n")

if args.txt != "NA":
    out.close()
out2.close()