#!/usr/bin/env python3
import argparse
import sys
from collections import defaultdict
from statistics import mean

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Combine local BLAST best hit and NR top-N hits into a summary.\n"
            "NR/DIAMOND must be outfmt 6: qseqid sseqid pident length stitle (in this order).\n"
            "REF must be outfmt 6: qseqid sseqid."
        )
    )
    p.add_argument("-b1", "--nr", required=True,
                   help="NR/DIAMOND file (outfmt 6: qseqid sseqid pident length stitle).")
    p.add_argument("-b2", "--ref", required=True,
                   help="Local BLAST file (outfmt 6: qseqid sseqid). Expected 1 hit/query.")
    p.add_argument("-o", "--out", required=True, help="Output TSV path.")
    p.add_argument("--nr-max", type=int, default=10,
                   help="Number of NR hits to report per query (default: 10).")
    p.add_argument("--sort-queries", action="store_true",
                   help="Sort query IDs in the output (default: keep input union order).")
    return p.parse_args()

def read_ref_best(ref_path):
    """Read best ref hit per query (first occurrence). Expect: qseqid sseqid."""
    best = {}
    with open(ref_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            qid, sid = parts[0], parts[1]
            if qid not in best:
                best[qid] = sid
    return best

def read_nr_hits(nr_path, nr_max):
    """
    Read NR hits, keeping:
      - up to nr_max stitles per query (file order)
      - all pident and length values per query (for averaging)

    Expect columns (exact order): qseqid sseqid pident length stitle
    """
    stitles = defaultdict(list)   # qid -> [stitle...]
    pidents = defaultdict(list)   # qid -> [float pident...]
    lengths = defaultdict(list)   # qid -> [int length...]

    with open(nr_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            # Split only the first 4 tabs; leave the rest as stitle (which may contain spaces)
            parts = line.split("\t")
            qid, _sid, pident_str, length_str, stitle = parts[0], parts[1], parts[2], parts[3], parts[6]
            # if len(parts) < 5:
            #     # try strict split as fallback
            #     parts = line.split("\t")
            #     if len(parts) < 5:
            #         continue
            #     qid, _sid, pident_str, length_str, stitle = parts[0], parts[1], parts[2], parts[3], parts[4]
            # else:
            #     qid, _sid, pident_str, length_str, stitle = parts
            #
            # Collect averages (even if we exceed nr_max for stitles)
            try:
                pidents[qid].append(float(pident_str))
            except ValueError:
                pass
            try:
                lengths[qid].append(int(float(length_str)))
            except ValueError:
                pass

            # Keep up to nr_max stitles
            if len(stitles[qid]) < nr_max:
                stitles[qid].append(stitle)

    return stitles, pidents, lengths

def fmt_avg(vals):
    if not vals:
        return ""
    # 4 decimal places for pident; integer for lengths average also 1 decimal?
    # We'll do 4 decimals for pident, and 1 decimal for length to keep it readable.
    return f"{mean(vals):.4f}"

def fmt_avg_len(vals):
    if not vals:
        return ""
    return f"{mean(vals):.1f}"

def main():
    args = parse_args()

    ref_best = read_ref_best(args.ref)
    stitles, pidents, lengths = read_nr_hits(args.nr, args.nr_max)

    # Union of queries present anywhere
    all_queries = list({**{q: None for q in ref_best.keys()},
                        **{q: None for q in stitles.keys()}}.keys())
    if args.sort_queries:
        all_queries = sorted(all_queries)

    with open(args.out, "w", encoding="utf-8") as out:
        # Header:
        # 1) query_id
        # 2) subject_id_ref
        # 3) avg_seqid (mean pident of NR hits)
        # 4) avg_align_len (mean length of NR hits)
        # 5-... stitle_1..stitle_N
        header = (["query_id", "subject_id_ref", "avg_seqid", "avg_align_len"] +
                  [f"stitle_{i}" for i in range(1, args.nr_max + 1)])
        out.write("\t".join(header) + "\n")

        for q in all_queries:
            ref_sid = ref_best.get(q, "")
            avg_pid = fmt_avg(pidents.get(q, []))
            avg_len = fmt_avg_len(lengths.get(q, []))
            titles = stitles.get(q, [])
            if len(titles) < args.nr_max:
                titles = titles + [""] * (args.nr_max - len(titles))
            else:
                titles = titles[:args.nr_max]

            row = [q, ref_sid, avg_pid, avg_len] + titles
            out.write("\t".join(row) + "\n")

if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        sys.exit(0)
