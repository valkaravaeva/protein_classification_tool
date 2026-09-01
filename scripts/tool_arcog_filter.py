#!/usr/bin/env python3
"""
Filter raw arCOG search hits (BLAST/DIAMOND-style tabular output) down to the
single best hit per query, using bit score, then percent identity, then e-value as
tie-breakers, subject to minimum identity/e-value cutoffs.

Standalone tool (not called from pipeline.nf) — used for the arCOG side analysis.
"""

import os
import sys
import argparse


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments (positional, matching the original script order)."""
    parser = argparse.ArgumentParser(
        description='Filter raw arCOG search hits (BLAST/DIAMOND-style tabular output) down to the single best hit per query, using bit score, then percent identity, then e-value as tie-breakers, subject to minimum identity/e-value cutoffs.'
    )
    parser.add_argument("taxfile", help="taxonomy.tsv")
    parser.add_argument("outfile", help="savefile")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    taxfile = args.taxfile
    outfile = args.outfile

    besthits = {}
    currdir = os.path.dirname(taxfile)
    with open(taxfile) as tf:
        for line in tf:
            currgen = line.strip().split("\t")[0] + "_arcogs.tsv"
            infile = os.path.join(currdir,currgen)
            with open(infile) as ff:
                for line in ff:
                    sp = line.strip().split("\t")
                    q = sp[0]
                    s = sp[1]
                    if "gi|" in s:
                        continue
                    id = float(sp[2])
                    ev = float(sp[10]) #10
                    if id < 25:
                        continue
                    if ev > 0.0000000001:
                        continue
                    bitscore = float(sp[11]) #11
                    if q not in besthits:
                        besthits[q] = [s,id,ev,bitscore]
                    else:
                        curr_bh = besthits[q]
                        curr_id = curr_bh[1]
                        curr_ev = curr_bh[2]
                        curr_bt = curr_bh[3]
                        if bitscore > curr_bt:
                            besthits[q] = [s, id, ev, bitscore]
                        elif bitscore == curr_bt:
                            if id > curr_id:
                                besthits[q] = [s,id,ev,bitscore]
                            elif id == curr_id:
                                if ev <= curr_ev:
                                    besthits[q] = [s, id, ev, bitscore]

    with open(outfile, "w") as svf:
        for k,v in besthits.items():
            strv = [str(vv) for vv in v]
            tw = k + "\t" + "\t".join(strv) + "\n"
            svf.write(tw)


if __name__ == "__main__":
    main()
