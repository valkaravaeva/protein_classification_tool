#!/usr/bin/env python3
"""
Build a per-KEGG-module table listing the score cutoff used for every KO in that
module (including KOs inside complexes and alternative-KO groups), for reporting
alongside module-completeness results.

Nextflow process: prepCutoffs
"""

import sys
import argparse


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments (positional, matching the original script order)."""
    parser = argparse.ArgumentParser(
        description='Build a per-KEGG-module table listing the score cutoff used for every KO in that module (including KOs inside complexes and alternative-KO groups), for reporting alongside module-completeness results.'
    )
    parser.add_argument("modulesfile", help="path to the module kegg architecture")
    parser.add_argument("thresholdfile", help="path to the ko_list")
    parser.add_argument("save_cutoff", help="where to save")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    modulesfile = args.modulesfile
    thresholdfile = args.thresholdfile
    save_cutoff = args.save_cutoff

    cutoff_dict = {}
    with open(thresholdfile) as tf:
        for line in tf:
            tsp = line.strip().split("\t")
            ko = tsp[0]
            thresh = tsp[1]
            cutoff_dict[ko] = thresh

    mod2thresh = {}
    with open(modulesfile) as mf:
        mr = mf.read().split("----------------------------------------------------------------------------------\n")
        for el in mr:
            if el == "":
                continue
            curr_mod_and_kos = el.strip().split("\n")
            mod = curr_mod_and_kos[0].split("|")[0].strip()
            modname = curr_mod_and_kos[0].split("|")[1].strip()
            kos = curr_mod_and_kos[1:]
            modcutoffs = {}
            for ko in kos:
                if "," in ko:
                    alts = ko.split(",")
                    altcutoffs = []
                    for alt in alts:
                        if "+" in alt:
                            sub_cutoffs = []
                            subs = alt.split("+")
                            for sub in subs:
                                if sub in cutoff_dict:
                                    curr_thr = cutoff_dict[sub]
                                else:
                                    curr_thr = "-"
                                sub_cutoffs.append(curr_thr)
                            altcutoffs.append("+".join(sub_cutoffs))
                        else:
                            if alt in cutoff_dict:
                                altcut = cutoff_dict[alt]
                            else:
                                altcut = "-"
                            altcutoffs.append(altcut)
                    modcutoffs[ko] = ",".join(altcutoffs)
                else:
                    if "+" in ko:
                        sub_cutoffs = []
                        subs = ko.split("+")
                        for sub in subs:
                            if sub in cutoff_dict:
                                curr_thr = cutoff_dict[sub]
                            else:
                                curr_thr = "-"
                            sub_cutoffs.append(curr_thr)
                        modcutoffs[ko] = "+".join(sub_cutoffs)
                    else:
                        if ko in cutoff_dict:
                            modcutoffs[ko] = cutoff_dict[ko]
                        else:
                            modcutoffs[ko] = "-"
            if mod not in mod2thresh:
                mod2thresh[mod] = modcutoffs

    with open(save_cutoff, "w") as svf:
        for m,cuts in mod2thresh.items():
            cutsave = []
            for k,v in cuts.items():
                kv = k + ":" + v
                cutsave.append(kv)
            tw = m + "\t" + "|".join(cutsave) + "\n"
            svf.write(tw)


if __name__ == "__main__":
    main()
