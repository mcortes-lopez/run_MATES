import pandas as pd
import numpy as np
import pickle
import os
from os.path import join
import scipy
from scipy import sparse
import shutil
import sys
from pathlib import Path

import pysam
import pybedtools

from intervaltree import IntervalTree

from collections import defaultdict

import argparse

def create_TE_vec_df(te_bed):
    """
    pd.DataFrame
        ['chromosome','start','end','TE_Name','te_index','strand','TE_fam','length']
        with each feature extended by 50 bp at its start.
    """
    records = []
    
    for feat in te_bed:
        chrom  = feat.chrom
        start  = int(feat.start)
        end    = int(feat.end)
        name   = feat.name
        idx    = int(feat.fields[4])
        strand = feat.strand
        fam    = feat.fields[6]
        length = int(feat.fields[7])

        records.append({
            'chromosome': chrom,
            'start':      start,
            'end':        end,
            'TE_Name':    name,
            'te_index':      idx,
            'strand':     strand,
            'TE_fam':     fam,
            'length':     length
        })

    df = pd.DataFrame(records, columns=['chromosome','start','end','TE_Name','te_index','strand','TE_fam','length'])
    
    return df

def build_TE_tss_trees(TE_vec_df, window=50):
    """
    Build one IntervalTree per chromosome that covers only the ±window
    around each TE's TSS: start for '+' strands, end for '-' strands.
    Stores (te_index, strand) in each node.
    """
    trees = {}
    for chrom, grp in TE_vec_df.groupby("chromosome"):
        tree = IntervalTree()
        for r in grp.itertuples(index=False):
            ti = r.te_index
            if r.strand == "+":
                ws = max(0, r.start - window)
                we = r.start + window
            else:
                ws = max(0, r.end   - window)
                we = r.end   + window
            tree.addi(ws, we, (ti, r.strand))
        trees[chrom] = tree
    return trees

def generate_unique_matric_fast(sample, barcode, path_to_bam, TE_vec_df):
    from collections import defaultdict
    import pandas as pd
    import os
    from pathlib import Path
    import pysam
    from intervaltree import IntervalTree

    cur_path = os.getcwd()
    base = Path(cur_path) / 'count_long_reads' / sample / barcode
    base.mkdir(parents=True, exist_ok=True)
    
    f1 = base / "read_weights.csv"
    f2 = base / "TE_unique_Info.csv"
    f3 = base / "AS_reads.csv"
    if f1.exists() and f2.exists() and f3.exists():
        print(f"[{barcode}] all three files exist, skipping")
        return

    # --- Pass 1: alignment counts per read ---
    bam1 = pysam.AlignmentFile(path_to_bam, "rb")
    align_counts = defaultdict(int)
    for read in bam1.fetch(until_eof=True):
        if read.is_unmapped or read.is_supplementary or read.is_duplicate:
            continue
        align_counts[read.query_name] += 1
    bam1.close()

    counts_df = pd.DataFrame({
        'read_name': list(align_counts.keys()),
        'count': list(align_counts.values())
    })
    wpath = base / "read_weights.csv"
    counts_df.to_csv(wpath, index=False)
    print(f"[{barcode}] wrote raw-mapping-counts → {wpath} ({counts_df.shape[0]} reads)")

    # --- Build per-chromosome interval trees ---
    tree_by_chrom = build_TE_tss_trees(TE_vec_df, window=50)

    # --- Pass 2: assign alignments to TEs ---
    bam2 = pysam.AlignmentFile(path_to_bam, "rb")
    sense_counts     = defaultdict(float)
    anti_counts      = defaultdict(float)
    seen_sense       = defaultdict(set)
    seen_antisense   = defaultdict(set)
    
    for read in bam2.fetch(until_eof=True):
        # 1) basic filters
        if read.is_unmapped or read.is_supplementary or read.is_duplicate:
            continue
        # drop pure‐skip alignments
        if read.cigartuples and all(op == 3 for op,_ in read.cigartuples):
            continue
    
        name   = read.query_name
        weight = 1.0 / align_counts.get(name, 1)
        chrom  = read.reference_name
        if chrom not in tree_by_chrom:
            continue
    
        # 2) pick anchor position
        pos = read.reference_start if not read.is_reverse else read.reference_end
    
        # 3) lookup in ±50bp TSS windows
        for iv in tree_by_chrom[chrom].at(pos):
            te_idx, te_strand = iv.data
    
            # 4) classify orientation
            if not read.is_reverse and te_strand == "+":
                target, seen = sense_counts, seen_sense
            elif read.is_reverse and te_strand == "-":
                target, seen = sense_counts, seen_sense
            elif read.is_reverse and te_strand == "+":
                target, seen = anti_counts, seen_antisense
            elif not read.is_reverse and te_strand == "-":
                target, seen = anti_counts, seen_antisense
            else:
                continue
    
            # 5) dedupe & tally
            if name in seen[te_idx]:
                continue
            seen[te_idx].add(name)
            target[te_idx] += weight
    
    bam2.close()

    # --- Write output files ---
    df_sense = pd.DataFrame({
        "te_index": list(sense_counts),
        "TE_sense_read_num": [sense_counts[i] for i in sense_counts],
    }).sort_values("te_index", ignore_index=True)
    p1 = base / "TE_unique_Info.csv"
    df_sense.to_csv(p1, index=False)
    print(f"[{barcode}] wrote sense counts → {p1} ({len(sense_counts)} TEs)")

    df_anti = pd.DataFrame({
        "te_index": list(anti_counts),
        "TE_antisense_read_num": [anti_counts[i] for i in anti_counts],
    }).sort_values("te_index", ignore_index=True)
    p2 = base / "AS_reads.csv"
    df_anti.to_csv(p2, index=False)
    print(f"[{barcode}] wrote antisense counts → {p2} ({len(anti_counts)} TEs)")


def main():
    p = argparse.ArgumentParser(
        description="Count sense/antisense TE reads for one barcode."
    )
    p.add_argument("--sample",   required=True,
                   help="Sample name (will write into count_long_reads/<sample>/<barcode>)")
    p.add_argument("--barcode",  required=True,
                   help="Barcode ID (basename of the BAM file, no .bam)")
    p.add_argument("--bam",      required=True,
                   help="Path to the barcode-specific BAM")
    p.add_argument("--te-ref",   required=True,
                   help="Path to the TE reference BED (uncompressed .bed)")
    args = p.parse_args()

    raw_te   = pybedtools.BedTool(args.te_ref)
    TE_vec_df = create_TE_vec_df(raw_te)

    # run the fast per-barcode count
    generate_unique_matric_fast(
        args.sample,
        args.barcode, 
        args.bam, 
        TE_vec_df
    )

if __name__ == "__main__":
    main()
