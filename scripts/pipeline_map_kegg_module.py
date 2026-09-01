#!/usr/bin/env python3
"""
Map characterized protein hits to their KEGG module(s) using the KEGG module
architecture file, and additionally group hits by module for downstream steps.

Nextflow process: mapModule
"""

import sys
import argparse


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments (positional, matching the original script order)."""
    parser = argparse.ArgumentParser(
        description='Map characterized protein hits to their KEGG module(s) using the KEGG module architecture file, and additionally group hits by module for downstream steps.'
    )
    parser.add_argument("chfile", help="file with characterized hits")
    parser.add_argument("modfile", help="module to ko map file, built into the pipeline")
    parser.add_argument("savefile", help="path to save the uncharacterized hits with mapped KEGG module information")
    parser.add_argument("savefile1", help="path to save the uncharacterized hits GROUPED by KEGG module information")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    chfile = args.chfile
    modfile = args.modfile
    savefile = args.savefile
    savefile1 = args.savefile1

    ko2mod = {}
    with open(modfile) as mf:
        # FIXED: see the identical note (and fix) in pipeline_analyze_characterized.py - the
        # original wrapped this in an extra "for line in mf:" together with "mf.read()" inside
        # the loop, which silently dropped the first line of modfile. Now parses the whole file
        # directly.
        mr = mf.read().split("----------------------------------------------------------------------------------\n")
        for el in mr:
            if el == "":
                continue
            curr_mod_and_kos = el.strip().split("\n")
            mod = curr_mod_and_kos[0].split("|")[0].strip()
            modname = curr_mod_and_kos[0].split("|")[1].strip()
            kos = curr_mod_and_kos[1:]
            kolist = set()
            for ko in kos:
                ks = ko.split(",")
                for k in ks:
                    if "+" not in k:
                        if "-" not in k:
                            kolist.add(k)
                        else:
                            ksplit = k.split("-")
                            for kt in ksplit:
                                kolist.add(kt)
                    else:
                        ksplit = k.split("+")
                        for kt in ksplit:
                            if "-" not in kt:
                                kolist.add(kt)
                            else:
                                ksplitt = kt.split("-")
                                for ktt in ksplitt:
                                    kolist.add(ktt)
            for ko in kolist:
                if ko not in ko2mod:
                    ko2mod[ko] = [mod]
                else:
                    ko2mod[ko].append(mod)

    groupsave = {} #module:entry
    with open(savefile, "w") as svf:
        with open(chfile) as chf:
            for line in chf:
                chsp = line.strip().split("\t")
                acc = chsp[0]
                tax = "\t".join(chsp[19:])
                ko = chsp[1].split("|")[0]
                if ko == "No_KO":
                    mod = "-"
                else:
                    if ko in ko2mod:
                        mods = sorted(ko2mod[ko])
                        mod = "|".join(mods)
                        ts = "\t".join([acc,chsp[1],tax]) + "\n"
                        for m in mods:
                            if m not in groupsave:
                                groupsave[m] = [ts]
                            else:
                                groupsave[m].append(ts)
                    else:
                        mod = "-"
                tw = "\t".join(chsp) + "\t" + mod + "\n"
                svf.write(tw)

    ##module,acc,ko,taxonomy
    with open(savefile1, "w") as sf:
        for k,vs in groupsave.items():
            for v in vs:
                tw = k + "\t" + v
                sf.write(tw)


if __name__ == "__main__":
    main()
