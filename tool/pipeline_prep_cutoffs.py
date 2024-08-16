import sys

##make table with ko cutoffs per module
modulesfile = sys.argv[1] #path to the module kegg architecture
thresholdfile = sys.argv[2] #path to the ko_list
save_cutoff = sys.argv[3] #where to save

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

with open(save_cutoff,"w+") as svf:
    for m,cuts in mod2thresh.items():
        cutsave = []
        for k,v in cuts.items():
            kv = k + ":" + v
            cutsave.append(kv)
        tw = m + "\t" + "|".join(cutsave) + "\n"
        svf.write(tw)
