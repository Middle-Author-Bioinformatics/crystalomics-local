#!/usr/bin/env python3
import argparse
import sys
from collections import defaultdict
from statistics import mean

# Column indices (0-based) as requested:
NR_PIDENT_IDX = 2
NR_LENGTH_IDX = 3
NR_EVALUE_IDX = 4      # NR evalue is the 5th column
REF_PIDENT_IDX = 2
REF_LENGTH_IDX = 3
REF_EVALUE_IDX = 10    # local ref evalue in default 12-col outfmt6

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Combine local BLAST best hit and NR top-N hits into a summary.\n"
            "NR/DIAMOND must be -outfmt 6 with: qseqid sseqid pident length evalue stitle qseq sseq\n"
            "  (i.e., evalue at idx 4; last three columns are stitle, qseq, sseq → stitle is idx -3).\n"
            "Local REF is BLAST -outfmt 6; supports 12-col default (evalue@10), 4-col, or 2-col."
        )
    )
    p.add_argument("-b1", "--nr", required=True,
                   help="NR/DIAMOND file (expects pident@2, length@3, evalue@4, and stitle at -3).")
    p.add_argument("-b2", "--ref", required=True,
                   help="Local BLAST file. Supports 12-col default (evalue@10), 4-col, or 2-col.")
    p.add_argument("-o", "--out", required=True, help="Output TSV path.")
    p.add_argument("--nr-max", type=int, default=10,
                   help="Number of NR titles to report per query (default: 10).")
    p.add_argument("--sort-queries", action="store_true",
                   help="Sort query IDs in the output (default: keep input union order).")
    return p.parse_args()

def read_ref_best(ref_path):
    """
    Read best ref hit per query (first occurrence).
    Accepts:
      - 12-col default BLAST outfmt6 (evalue at idx 10)
      - 4-col: qseqid sseqid pident length
      - 2-col: qseqid sseqid
    Returns: qid -> (sid, pident_str, length_str, evalue_str)
    """
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
            if qid in best:
                continue
            pid = parts[REF_PIDENT_IDX] if len(parts) > REF_PIDENT_IDX else ""
            alen = parts[REF_LENGTH_IDX] if len(parts) > REF_LENGTH_IDX else ""
            eval_str = parts[REF_EVALUE_IDX] if len(parts) > REF_EVALUE_IDX else ""
            best[qid] = (sid, pid, alen, eval_str)
    return best

def read_nr_hits(nr_path, nr_max):
    """
    Read NR hits, keeping:
      - up to nr_max stitles per query (file order)
      - all pident, length, evalue values per query (for averaging)

    Expected indices:
      pident@2, length@3, evalue@4, stitle@-3 (with trailing qseq, sseq)
    """
    stitles = defaultdict(list)   # qid -> [stitle...]
    pidents = defaultdict(list)   # qid -> [float ...]
    lengths = defaultdict(list)   # qid -> [int   ...]
    evalues = defaultdict(list)   # qid -> [float ...]

    with open(nr_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            # Need at least the first 5 fields for pident/length/evalue
            if len(parts) <= max(NR_PIDENT_IDX, NR_LENGTH_IDX, NR_EVALUE_IDX):
                continue

            qid = parts[0]
            pident_str = parts[NR_PIDENT_IDX]
            length_str = parts[NR_LENGTH_IDX]
            evalue_str = parts[NR_EVALUE_IDX]

            # stitle handling:
            # Preferred: ... evalue, stitle, qseq, sseq  => stitle is -3
            if len(parts) >= 8:
                stitle = parts[-3]
            # Fallback: if file lacks qseq/sseq but has stitle at fixed idx 5
            elif len(parts) >= 6:
                stitle = parts[5]
            else:
                # Last-resort fallback
                stitle = parts[-1]

            # Collect for averages
            try:
                pidents[qid].append(float(pident_str))
            except ValueError:
                pass
            try:
                lengths[qid].append(int(float(length_str)))
            except ValueError:
                pass
            try:
                evalues[qid].append(float(evalue_str))
            except ValueError:
                pass

            # Keep up to nr_max stitles
            if len(stitles[qid]) < nr_max:
                stitles[qid].append(stitle)

    return stitles, pidents, lengths, evalues

def fmt_avg(vals, places=4):
    return "" if not vals else f"{mean(vals):.{places}f}"

def fmt_avg_len(vals):
    return "" if not vals else f"{mean(vals):.1f}"

def main():
    args = parse_args()

    ref_best = read_ref_best(args.ref)
    stitles, pidents, lengths, evalues = read_nr_hits(args.nr, args.nr_max)

    # Union of queries present anywhere
    all_queries = list({**{q: None for q in ref_best.keys()},
                        **{q: None for q in stitles.keys()}}.keys())
    if args.sort_queries:
        all_queries = sorted(all_queries)

    with open(args.out, "w", encoding="utf-8") as out:
        header = ([
            "query_id", "subject_id_ref",
            "ref_perc_id", "ref_aln_len", "ref_evalue",
            "avg_seqid", "avg_align_len", "avg_evalue"
        ] + [f"stitle_{i}" for i in range(1, args.nr_max + 1)])
        out.write("\t".join(header) + "\n")

        for q in all_queries:
            ref_sid, ref_pid, ref_alen, ref_eval = ref_best.get(q, ("", "", "", ""))
            avg_pid  = fmt_avg(pidents.get(q, []), places=4)
            avg_len  = fmt_avg_len(lengths.get(q, []))
            avg_eval = fmt_avg(evalues.get(q, []), places=6)
            titles = stitles.get(q, [])
            if len(titles) < args.nr_max:
                titles += [""] * (args.nr_max - len(titles))
            else:
                titles = titles[:args.nr_max]

            row = [q, ref_sid, ref_pid, ref_alen, ref_eval, avg_pid, avg_len, avg_eval] + titles
            out.write("\t".join(row) + "\n")

if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        sys.exit(0)
