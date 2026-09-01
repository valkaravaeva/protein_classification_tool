#!/usr/bin/env python3
"""
Reshape either an arCOG or a KEGG-pipeline "uncharacterized percent" table into a
long-format table (genome, percent, taxon) suitable for plotting boxplots per
taxon.

Standalone tool (not called from pipeline.nf).
"""

import sys
import argparse


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments (positional, matching the original script order)."""
    parser = argparse.ArgumentParser(
        description='Reshape either an arCOG or a KEGG-pipeline "uncharacterized percent" table into a long-format table (genome, percent, taxon) suitable for plotting boxplots per taxon.'
    )
    parser.add_argument("typedata", help="options: arcog, pipeline")
    parser.add_argument("table", help="path to arcog/pipeline table counts per genome")
    parser.add_argument("taxonomy", help="path to taxonomy.tsv file")
    parser.add_argument("savefile", help="path to save the file, boxplot_table.tsv or your custom name")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    typedata = args.typedata
    table = args.table
    taxonomy = args.taxonomy
    savefile = args.savefile

    if typedata == "arcog":
        ##arcogs
        with open(savefile, "w") as svf:
            header = "Genome\tUncharacterized_percent\tTaxon\n"
            svf.write(header)
            with open(table) as tf:
                next(tf)
                for line in tf:
                    tsp = line.strip().split("\t")
                    gen = tsp[0]
                    unch = tsp[6]
                    phyl = tsp[8].replace(" ", "_")
                    tw = gen + "\t" + unch + "\t" + phyl + "\n"
                    svf.write(tw)
    else:
        ##pipeline
        taxdict = {}
        with open(taxonomy) as tf:
            for line in tf:
                tsp = line.strip().split("\t")
                genome = tsp[0]
                phyl = tsp[2]
                taxdict[genome] = phyl

        with open(savefile, "w") as svf:
            header = "Genome\tUncharacterized_percent\tTaxon\n"
            svf.write(header)
            with open(table) as tf:
                next(tf)
                for line in tf:
                    tsp = line.strip().split("\t")
                    gen = tsp[0]
                    unch = tsp[1]
                    phyl = taxdict[gen].replace(" ", "_")
                    tw = gen + "\t" + unch + "\t" + phyl + "\n"
                    svf.write(tw)


if __name__ == "__main__":
    main()
