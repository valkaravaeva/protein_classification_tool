#!/usr/bin/env python3
"""
Filter raw HMMER tabular (--tblout/--domtblout style) output against per-KO score
and e-value cutoffs (either full-sequence or domain cutoffs), keeping only the
best-scoring KO assignment per protein accession.

Nextflow process: filterHmmer
"""

import sys
import argparse
import os


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments (positional, matching the original script order)."""
    parser = argparse.ArgumentParser(
        description='Filter raw HMMER tabular (--tblout/--domtblout style) output against per-KO score and e-value cutoffs (either full-sequence or domain cutoffs), keeping only the best-scoring KO assignment per protein accession.'
    )
    parser.add_argument("genomefile", help="path to the genome list file")
    parser.add_argument("savefile", help="path to save the filtered file")
    parser.add_argument("ko_list", help="path to cutoff file, must be in tabular format if custom: ko\tcutoff")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    genomefile = args.genomefile
    savefile = args.savefile
    ko_list = args.ko_list

    currdir = os.path.dirname(genomefile)

    print("Start filtering hmmer")

    evaluecutoff=0.00005 #change if necessary

    full_score_dict = {}  # save cutoffs as model:cutoff pair
    domain_score_dict = {}
    no_score_dict = {} #dictionary for models with no score
    with open(ko_list) as kl:
        next(kl)
        for line in kl:
            ksp = line.strip().split("\t")
            ko = ksp[0]
            score = ksp[1]
            if score == '-':
                score = 50.00
            score_type = ksp[2]
            if score_type == "full":
                full_score_dict[ko] = float(score)
            elif score_type == "domain":
                domain_score_dict[ko] = float(score)

    # parse hmmfile, check if each hit passes the corresponding cutoff, save certain columns to the filtered file
    with open(savefile, "w") as svf:
        tosavedict = {}
        besthit = {}
        with open(genomefile) as gf:
            for l in gf:
                currgen = l.strip().split("\t")[0] + "_hmmer.tsv"
                hmmfile = os.path.join(currdir,currgen)
                with open(hmmfile) as hmmf:
                    for line in hmmf:
                        if line.startswith("#"):
                            continue
                        sp = line.strip().split(" ")
                        spfilt = [x for x in sp if x != ""]
                        acc = spfilt[0]
                        if acc not in besthit:
                            besthit[acc] = 0
                        ko = spfilt[2]
                        if ko in full_score_dict:
                            score = float(spfilt[5])
                            evalue = float(spfilt[4])
                            if evalue > evaluecutoff:
                                continue
                            thresh = full_score_dict[ko]
                        elif ko in domain_score_dict:
                            score = float(spfilt[8])
                            evalue = float(spfilt[7])
                            if evalue > evaluecutoff:
                                continue
                            thresh = domain_score_dict[ko]
                        if score >= thresh:
                            if score > besthit[acc]:
                                besthit[acc] = score
                                tosavedict[acc] = ko
        for k,vs in tosavedict.items():
            tw = k + "\t" + vs + "\n"
            svf.write(tw)

    print("Finished filtering hmmer")


if __name__ == "__main__":
    main()
