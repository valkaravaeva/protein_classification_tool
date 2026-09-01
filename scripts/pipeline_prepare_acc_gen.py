#!/usr/bin/env python3
"""
Build a simple accession-to-genome mapping table by reading the headers of each
genome's protein FASTA (.faa) file.

Nextflow process: createAccgen
"""

import sys
import argparse
import os
from Bio import SeqIO


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments (positional, matching the original script order)."""
    parser = argparse.ArgumentParser(
        description="Build a simple accession-to-genome mapping table by reading the headers of each genome's protein FASTA (.faa) file."
    )
    parser.add_argument("genomelist", help="path to genome list file")
    parser.add_argument("savefile", help="path to output file")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    genomelist = args.genomelist
    savefile = args.savefile

    print("Start creating accession - genome map")
    currdir = os.path.dirname(genomelist)
    with open(savefile, "w") as svf:
        with open(genomelist) as gf:
            for l in gf:
                currgenome = l.strip().split("\t")[0]
                faaname = currgenome + ".faa"
                genomefaa = os.path.join(currdir,faaname)
                with open(genomefaa) as pf:
                    for rec in SeqIO.parse(pf,"fasta"):
                        acc = rec.id
                        tw = currgenome + "\t" + acc + "\n"
                        svf.write(tw)

    print("Finished creating accession - genome map")


if __name__ == "__main__":
    main()
