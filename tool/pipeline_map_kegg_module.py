###map KEGG module information to characterized hits and group hits per module
import sys
chfile = sys.argv[1] #file with characterized hits
modfile = sys.argv[2] #module to ko map file, built into the pipeline
savefile = sys.argv[3] #path to save the uncharacterized hits with mapped KEGG module information
savefile1 = sys.argv[4] #path to save the uncharacterized hits GROUPED by KEGG module information

ko2mod = {}
with open(modfile) as mf:
    for line in mf:
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
with open(savefile,"w+") as svf:
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
with open(savefile1,"w+") as sf:
    for k,vs in groupsave.items():
        for v in vs:
            tw = k + "\t" + v
            sf.write(tw)
