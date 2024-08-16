##calculate presence absence per genome of characterized per module -> fullness of the module / presence-absence
import sys
import os
import shutil

charhits = sys.argv[1] #path to characterized grouped file with module info
mod2kofile = sys.argv[2] #kegg module architecture file
dbfile = sys.argv[3] #path to taxonomy file
modulesfile = sys.argv[4] #path to a unique list of KEGG modules (1 per line, e.g. "M00001	Glycolysis_(Embden-Meyerhof_pathway)_glucose_=>_pyruvate"; tab separated)
konamefile = sys.argv[5] #path to file with names for KOs, inbuilt
savefile_presabs = sys.argv[6] #where to save total output file
savenohits = sys.argv[7] #path to save the list of modules that have no hits
save_with_hits = sys.argv[8] #path to save the list of modules that have hits

savedir_permodule = "hits_per_module"
kickdir = os.path.join(savedir_permodule,"no_hits")
checkdir = os.path.join(savedir_permodule,"yes_hits")

if not os.path.isdir(kickdir):
   os.makedirs(kickdir)

if not os.path.isdir(checkdir):
   os.makedirs(checkdir)

taxdict = {}
with open(dbfile) as dbf:
    for line in dbf:
        dbsp = line.strip().split("\t")
        genome = dbsp[0]
        tax = [x.replace(" ","_") for x in dbsp[1:]]
        taxdict[genome] = "\t".join(tax)

print("taxonomy done")

mod2kodict = {}
modnamedict = {}
with open(mod2kofile) as mf:
    for line in mf:
        mr = mf.read().split("----------------------------------------------------------------------------------\n")
        for el in mr:
            if el == "":
                continue
            curr_mod_and_kos = el.strip().split("\n")
            mod = curr_mod_and_kos[0].split("|")[0].strip()
            modname = curr_mod_and_kos[0].split("|")[1].strip()
            modnamedict[mod] = modname
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
            mod2kodict[mod] = sorted(list(kolist))

knamedict = {}
with open(konamefile) as knf:
    for line in knf:
        knsp = line.strip().split("\t")
        knamedict[knsp[0]] = knsp[1]

print("modules done")

presencedict = {}
knamedict_specific = {}
gen2ko2accdict = {}
with open(charhits) as cf:
    for line in cf:
        csp = line.strip().split("\t")
        acc = csp[1]
        modid = csp[0].split("|")[0]
        kos = csp[2].split("|")
        if kos == "No_KO":
            continue
        if len(kos) > 2:
            actualkos = []
            konames = []
            for idx, element in enumerate(kos):
                if idx % 2 == 0:
                    actualkos.append(element)
                else:
                    konames.append(element)
        else:
            actualkos = [kos[0]]
            konames = [kos[1]]
        for idx, ko in enumerate(actualkos):
            kname = konames[idx]
            knamedict_specific[ko] = kname
        curr_gendict = {}
        gens = csp[3].split("|")
        for gen in gens:
            for kk in actualkos:
                if gen not in presencedict:
                    presencedict[gen] = {modid:{kk:1}}
                else:
                    if modid not in presencedict[gen]:
                        presencedict[gen][modid] = {kk:1}
                    else:
                        if kk not in presencedict[gen][modid]:
                            presencedict[gen][modid][kk] = 1
                        else:
                            presencedict[gen][modid][kk] += 1
            if gen not in gen2ko2accdict:
                gen2ko2accdict[gen] = {ko:[acc]}
            else:
                if ko not in gen2ko2accdict[gen]:
                    gen2ko2accdict[gen][ko] = [acc]
                else:
                    if acc not in gen2ko2accdict[gen][ko]:
                        gen2ko2accdict[gen][ko].append(acc)

modulecompleteness_percent = {}
with open(savefile_presabs,"w+") as svf:
    header = "Genome\tKEGG_Module\tKO\tPresence_(counts)\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tStrain\n"
    svf.write(header)
    for gen,genpres in presencedict.items():
        tax = taxdict[gen]
        genaccs = gen2ko2accdict[gen]
        for module,kos in mod2kodict.items():
            if module in genpres:
                modpres = genpres[module]
                modname = modnamedict[module]
                for ko in kos:
                    curr_presence = 0
                    curr_accs = "-"
                    if ko in genaccs:
                        curr_accs = "|".join(sorted(genaccs[ko]))
                    if ko in knamedict_specific:
                        koname = knamedict_specific[ko]
                    elif ko in knamedict:
                        koname = knamedict[ko]
                    else:
                        koname = ""
                    if ko in modpres:
                        curr_presence = modpres[ko]
                    tw = gen + "\t" + module + "|" + modname + "\t" + ko + "|" + koname + "\t" + str(curr_presence) + "\t" + tax + "\t" + curr_accs + "\n"
                    svf.write(tw)
            else:
                modname = modnamedict[module]
                for ko in kos:
                    curr_presence = 0
                    curr_accs = "-"
                    if ko in knamedict_specific:
                        koname = knamedict_specific[ko]
                    elif ko in knamedict:
                        koname = knamedict[ko]
                    else:
                        koname = ""
                    tw = gen + "\t" + module + "|" + modname + "\t" + ko + "|" + koname + "\t" + str(curr_presence) + "\t" + tax + "\t" + curr_accs + "\n"
                    svf.write(tw)

print("module completeness loaded")

##separate per module
modulesdict = {}
with open(modulesfile) as mf:
    for line in mf:
        msp = line.strip().split("\t")
        modulesdict[msp[0]] = msp[1]

print("separating files, might take a while")

savelines_dict = {}
linecount=0
with open(savefile_presabs) as ff:
    next(ff)
    for line in ff:
        linecount += 1
        fsp = line.strip().split("\t")
        mod = fsp[1].split("|")[0]
        if mod not in savelines_dict:
            savelines_dict[mod] = [line]
        else:
            savelines_dict[mod].append(line)
print("lines loaded")

for mod, lines in savelines_dict.items():
    print("Processing: ", mod)
    curr_savefilename = "characterized_hits_module_fullness_presence_absence_" + mod + ".tsv"
    curr_savefile = os.path.join(savedir_permodule,curr_savefilename)
    with open(curr_savefile,"w+") as currsvf:
        for line in lines:
            currsvf.write(line)

print("files separated") #very time limiting step

##sort separated files by empty modules and non-empty modules, create a list of modules with no hits
nohits = {}
for m,mn in modulesdict.items():
    curr_mod = m
    curr_filename = "characterized_hits_module_fullness_presence_absence_" + curr_mod + ".tsv"
    curr_file = os.path.join(savedir_permodule, curr_filename)
    curr_presences = {"0":0,"1":0}
    with open(curr_file) as ff:
        for line in ff:
            fsp = line.strip().split("\t")
            presence = fsp[3]
            if presence not in curr_presences:
                curr_presences[presence] = 1
            else:
                curr_presences[presence] += 1
    if curr_presences["1"] == 0:
        if len(curr_presences.keys()) <= 2:
            nohits[curr_mod] = 1
            curr_move = os.path.join(kickdir,curr_filename)
            shutil.move(curr_file, curr_move)
        else:
            curr_move = os.path.join(checkdir,curr_filename)
            shutil.move(curr_file, curr_move)
    else:
        curr_move = os.path.join(checkdir,curr_filename)
        shutil.move(curr_file, curr_move)

with open(savenohits,"w+") as svf:
    for k in nohits.keys():
        tw = k + "\n"
        svf.write(tw)

print("files sorted")

# create a list of modules that have hits
from pathlib import Path
with open(save_with_hits,"w+") as svf:
    for m,mn in modulesdict.items():
        curr_mod = m
        curr_filename = "characterized_hits_module_fullness_presence_absence_" + curr_mod + ".tsv"
        curr_file = os.path.join(checkdir, curr_filename)
        my_file = Path(curr_file)
        if my_file.is_file():
            tw = m + "\t" + mn + "\t?\n"
            svf.write(tw)

print("files saved")
