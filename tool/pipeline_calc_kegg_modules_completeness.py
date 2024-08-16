# ##calculate the presence/absence but in module completeness: 0%, <50%, 50%-<75%, >=75%, 100% of KOs from a module present ---> accounting for ALTERNATIVE KOs
# ###1. calculate completeness of complexes per genome
# ##2. calculate general completeness with alternatives for non-complex KOs
# ##3. merge/sum the results
import sys
import os

modulesfile = sys.argv[1] #kegg module architecture file
allgenfile = sys.argv[2] #accgen file
absentmodfile = sys.argv[3] #no hit modules txt file
cutoffs = sys.argv[4] #ko_list file
taxonomy = sys.argv[5]
modulenamefile = sys.argv[6] #module fullness file from previous step
complexpresence = sys.argv[7] #complex_presence_per_genome.tsv
altpresence = sys.argv[8] #alternative_presence_per_genome.tsv
save_pergenome = sys.argv[9] #modules_presence_per_genome.tsv
savefile_perc = sys.argv[10] #module_completeness_percent.tsv
moddir = sys.argv[11] #per module directory
checkdir = os.path.join(moddir,"yes_hits")
kickdir = os.path.join(moddir,"no_hits")

##1.
###calculate presence of complexes per genome --> if subunit occurs more than once -> count it only once

absent_modules = set()
with open(absentmodfile) as absf:
    for line in absf:
        absp = line.strip()
        absent_modules.add(absp)

complexes = {}
modnamedict = {}
with open(modulesfile) as mf:
    mr = mf.read().split("----------------------------------------------------------------------------------\n")
    for el in mr:
        if el == "":
            continue
        curr_mod_and_kos = el.strip().split("\n")
        mod = curr_mod_and_kos[0].split("|")[0].strip()
        modname = curr_mod_and_kos[0].split("|")[1].strip()
        modnamedict[mod] = modname
        kos = curr_mod_and_kos[1:]
        complexes[mod] = set()
        for ko in kos:
            if "," in ko:
                altkos = ko.split(",")
                for altko in altkos:
                    if "+" not in altko:
                        continue
                    else:
                        complexes[mod].add(altko)
            elif "+" in ko:
                complexes[mod].add(ko)

genomes2ko = {}
for m,compls in complexes.items():
    if not compls:
        continue
    mn = modnamedict[m]
    if m in absent_modules:
        continue
    curr_mod = m
    curr_file = os.path.join(checkdir,"characterized_hits_module_fullness_presence_absence_" + curr_mod + ".tsv")
    with open(curr_file) as mf:
        for line in mf:
            msp = line.strip().split("\t")
            genome = msp[0]
            ko = msp[2].split("|")[0]
            presence = int(msp[3])
            if genome not in genomes2ko:
                genomes2ko[genome] = {ko:presence}
            else:
                if ko not in genomes2ko[genome]:
                    genomes2ko[genome][ko] = presence
                else:
                    continue ##the same genome entry for ko can exist in different modules

used = {}
with open(complexpresence,"w+") as svf:
    header = "Genome\tComplex\tPresence_each_subunit\tTotal_presence\n"
    svf.write(header)
    for m,compls in complexes.items():
        if not compls:
            continue
        for curr_compl in compls:
            curr_kos = curr_compl.split("+")
            curr_numsub = len(curr_kos)
            for genome, kod in genomes2ko.items():
                curr_gen_comp = []
                for ko,pres in kod.items():
                    if pres > 1:
                        pres = 1
                    if ko in curr_kos:
                        curr_gen_comp.append(pres)
                if not curr_gen_comp:
                    curr_gen_comp = [0]*curr_numsub
                sumgencomp = sum(curr_gen_comp)
                fracgencomp = sumgencomp / curr_numsub
                if "|".join([genome,curr_compl]) not in used:
                    used["|".join([genome,curr_compl])] = 1
                    twl = [genome,curr_compl,"+".join([str(x) for x in curr_gen_comp]),"{:.2f}".format(fracgencomp)]
                    tw = "\t".join(twl) + "\n"
                    svf.write(tw)

print("Complexes presence per genome calculated")

##2. calculate presence of alternative kos per genome --> per group of alternative KOs -> presence should not be > 1 because alternatives catalyze the same reactions/step in the module

absent_modules = set()
with open(absentmodfile) as absf:
    for line in absf:
        absp = line.strip()
        absent_modules.add(absp)

genome2complex = {}
ko2complex = {}
with open(complexpresence) as cf:
    next(cf)
    for line in cf:
        csp = line.strip().split("\t")
        gen = csp[0]
        complex = csp[1]
        complexsplit = complex.split("+")
        for k in complexsplit:
            ko2complex[k] = complex
        pres = float(csp[-1])
        if gen not in genome2complex:
            genome2complex[gen] = {complex:pres}
        else:
            genome2complex[gen][complex] = pres

altkos_dict = {}
ko2altkos = {}
modnamedict = {}
modulecomposition = {}
real_alts = set()
with open(modulesfile) as mf:
    mr = mf.read().split("----------------------------------------------------------------------------------\n")
    for el in mr:
        if el == "":
            continue
        curr_mod_and_kos = el.strip().split("\n")
        mod = curr_mod_and_kos[0].split("|")[0].strip()
        modname = curr_mod_and_kos[0].split("|")[1].strip()
        modnamedict[mod] = modname
        kos = curr_mod_and_kos[1:]
        altkos_dict[mod] = set()
        modulecomposition[mod] = kos
        for ko in kos:
            if "," in ko:
                altkos_dict[mod].add(ko)
                real_alts.add(ko)
                altkos = ko.split(",")
                for altko in altkos:
                    if altko not in ko2altkos:
                        ko2altkos[altko] = [ko]
                    else:
                        if ko not in ko2altkos[altko]:
                            ko2altkos[altko].append(ko)

genomes2ko = {}
for m,compls in altkos_dict.items():
    if compls == set():
        continue
    mn = modnamedict[m]
    if m in absent_modules:
        continue
    curr_mod = m
    curr_file = os.path.join(checkdir,"characterized_hits_module_fullness_presence_absence_" + curr_mod + ".tsv")
    with open(curr_file) as mf:
        for line in mf:
            msp = line.strip().split("\t")
            genome = msp[0]
            ko = msp[2].split("|")[0]
            presence = float(msp[3])
            if genome not in genomes2ko:
                genomes2ko[genome] = {ko:presence}
            else:
                if ko not in genomes2ko[genome]:
                    genomes2ko[genome][ko] = presence
                else:
                    continue ##the same genome entry for ko can exist in different modules

used = {}
with open(altpresence,"w+") as svf:
    header = "Genome\tAlternatives\tPresence_each_alternative\tTotal_presence\n"
    svf.write(header)
    for m,alts in altkos_dict.items():
        if m in absent_modules:
            continue
        if not alts:
            continue
        for curr_alt in alts:
            curr_kos = curr_alt.split(",")
            curr_comps = {}
            for x in curr_kos:
                if "+" in x:
                    xsplit = x.split("+")
                    for xx in xsplit:
                        if xx not in curr_comps:
                            curr_comps[xx] = [x]
                        else:
                            if x not in curr_comps[xx]:
                                curr_comps[xx].append(x)
            for genome, kod in genomes2ko.items():
                curr_altpresence = {}
                curr_gen_comp = 0
                for ko,pres in kod.items():
                    if ko in curr_kos:
                        curr_altpresence[ko] = pres
                        if pres >= 1:
                            curr_gen_comp = 1
                    elif ko in curr_comps:
                        kocoms = curr_comps[ko]
                        for kocom in kocoms:
                            curr_complexes = genome2complex[genome]
                            curr_comp_pres = curr_complexes[kocom]
                            if kocom not in curr_altpresence:
                                curr_altpresence[kocom] = [curr_comp_pres]
                            else:
                                curr_altpresence[kocom].append(curr_comp_pres)
                            if curr_comp_pres >= 1:
                                curr_gen_comp = 1
                            elif curr_comp_pres > 0:
                                if curr_gen_comp == 0:
                                    curr_gen_comp = curr_comp_pres
                                else:
                                    if curr_comp_pres > curr_gen_comp:
                                        curr_gen_comp = curr_comp_pres
                if "|".join([genome,curr_alt]) not in used:
                    used["|".join([genome,curr_alt])] = 1
                    if not curr_altpresence:
                        curr_altpresencetw = [0]*len(curr_kos)
                    else:
                        curr_altpresencetw = []
                        for a in curr_alt.split(","):
                            if "+" not in a:
                                apres = genomes2ko[genome][a]
                                curr_altpresencetw.append(apres)
                            else:
                                apres = genome2complex[genome][a]
                                curr_altpresencetw.append(apres)
                    twl = [genome,curr_alt,",".join([str(x) for x in curr_altpresencetw]),"{:.2f}".format(curr_gen_comp)]
                    tw = "\t".join(twl) + "\n"
                    svf.write(tw)

print("Alternatives per genome calculated")

###3. get all kos, alternatives and complexes presence per module per genome
genomes = set()
with open(allgenfile) as agf:
    for line in agf:
        agsp = line.strip().split("\t")
        gen = agsp[0]
        genomes.add(gen)

absent_modules = set()
with open(absentmodfile) as absf:
    for line in absf:
        absp = line.strip()
        absent_modules.add(absp)

genome2complex = {}
ko2complex = {}
with open(complexpresence) as cf:
    next(cf)
    for line in cf:
        csp = line.strip().split("\t")
        gen = csp[0]
        complex = csp[1]
        complexsplit = complex.split("+")
        for k in complexsplit:
            if k not in ko2complex:
                ko2complex[k] = [complex]
            else:
                if complex not in ko2complex[k]:
                    ko2complex[k].append(complex)
        pres = float(csp[-1])
        if gen not in genome2complex:
            genome2complex[gen] = {complex:pres}
        else:
            genome2complex[gen][complex] = pres

genome2alts = {}
ko2alt = {}
with open(altpresence) as alf:
    next(alf)
    for line in alf:
        alsp = line.strip().split("\t")
        gen = alsp[0]
        alt = alsp[1]
        pres = float(alsp[-1])
        altsplit = alt.split(",")
        if gen not in genome2alts:
            genome2alts[gen] = {alt:pres}
        else:
            genome2alts[gen][alt] = pres
        for ko in altsplit:
            koplus = ko.split("+")
            for k in koplus:
                if k not in ko2alt:
                    ko2alt[k] = [alt]
                else:
                    if alt not in ko2alt[k]:
                        ko2alt[k].append(alt)

modulecomposition = {}
modnamedict = {}
with open(modulesfile) as mf:
    mr = mf.read().split("----------------------------------------------------------------------------------\n")
    for el in mr:
        if el == "":
            continue
        curr_mod_and_kos = el.strip().split("\n")
        mod = curr_mod_and_kos[0].split("|")[0].strip()
        modname = curr_mod_and_kos[0].split("|")[1].strip()
        modnamedict[mod] = modname
        kos = curr_mod_and_kos[1:]
        modulecomposition[mod] = kos

modulepres = {}
modulepres_prelim = {}
for m,compos in modulecomposition.items():
    mn = modnamedict[m]
    curr_mod = m
    modulepres[m] = {}
    modulepres_prelim[m] = {}
    modulesize = len(compos)
    if m in absent_modules:
        modulepres[m] = {}
        for gen in genomes:
            for ko in compos:
                modulepres[m][gen] = 0
        continue
    for gen in genomes:
        modulepres_prelim[m][gen] = {k:0 for k in compos}
    curr_file = os.path.join(checkdir,"characterized_hits_module_fullness_presence_absence_" + curr_mod + ".tsv")
    with open(curr_file) as mf:
        for line in mf:
            msp = line.strip().split("\t")
            genome = msp[0]
            ko = msp[2].split("|")[0]
            presence = float(msp[3])
            if ko in ko2complex: #if individual ko is in complexes
                ko_comps = ko2complex[ko] #list of all complexes including this ko
                for ko_comp in ko_comps: #iterate over each complex
                    if ko not in ko2alt: #check if individual ko has alternatives
                        if ko_comp in compos: #check if complex is in current module composition
                            ko_pres = genome2complex[genome][ko_comp] #get presence of current complex
                            checkpres = modulepres_prelim[m][genome][ko_comp] #check if presence was already assigned
                            if checkpres > 0:
                                if ko_pres > checkpres:
                                    modulepres_prelim[m][genome][ko_comp] = ko_pres
                            else:
                                modulepres_prelim[m][genome][ko_comp] = ko_pres
                        elif ko in compos: #if complex is not in module composition but individual ko is
                            if presence >= 1:
                                ko_pres = 1
                            elif presence == 0:
                                ko_pres = 0
                            modulepres_prelim[m][genome][ko] = ko_pres
                    else: #if individual ko does not have alternatives
                        ko_alts = ko2alt[ko] #list of all alternatives
                        for ko_alt in ko_alts:
                            if ko_alt in compos: #if alt in module composition
                                ko_pres = genome2alts[genome][ko_alt]
                                checkpres = modulepres_prelim[m][genome][ko_alt]
                                if checkpres > 0:
                                    if ko_pres > checkpres:
                                        modulepres_prelim[m][genome][ko_alt] = ko_pres
                                else:
                                    modulepres_prelim[m][genome][ko_alt] = ko_pres
                            elif ko in compos: #if individual ko in module composition
                                if presence >= 1:
                                    ko_pres = 1
                                elif presence == 0:
                                    ko_pres = 0
                                modulepres_prelim[m][genome][ko] = ko_pres
            if ko in ko2alt: #if ko not in complex but is in alts
                ko_alts = ko2alt[ko]
                for ko_alt in ko_alts:
                    if ko_alt in compos:
                        ko_pres = genome2alts[genome][ko_alt]
                        checkpres = modulepres_prelim[m][genome][ko_alt]
                        if checkpres > 0:
                            if ko_pres > checkpres:
                                modulepres_prelim[m][genome][ko_alt] = ko_pres
                        else:
                            modulepres_prelim[m][genome][ko_alt] = ko_pres
                    elif ko in compos:
                        if presence >= 1:
                            ko_pres = 1
                        elif presence == 0:
                            ko_pres = 0
                        modulepres_prelim[m][genome][ko] = ko_pres
            else:
                if presence >= 1:
                    ko_pres = 1
                elif presence == 0:
                    ko_pres = 0
                else:
                    ko_pres = presence
                if ko in compos:
                    modulepres_prelim[m][genome][ko] = "{:.2f}".format(ko_pres)

for m,gencomp in modulepres_prelim.items():
    for gen,comp in gencomp.items():
        compv = [float(v) for v in comp.values()]
        sumcompv = sum(compv)
        fraccompv = sumcompv / len(compv)
        modulepres[m][gen] = fraccompv


##methanogenesis module M00617 is a combination module of other methanogenesis modules; it is present if any of those are present: M00567,M00357,M00356,M00563 (sum up per genome)
methanfix = modulepres["M00617"]
for gen,frac in methanfix.items():
    allfracs = []
    for m in ["M00567","M00357","M00356","M00563"]:
        currfrac = modulepres[m][gen]
        allfracs.append(currfrac)
    sumfrac = sum(allfracs)
    if sumfrac > 1:
        sumfrac = 1
    modulepres["M00617"][gen] = sumfrac

with open(save_pergenome,"w+") as svf:
    header = "Module\tGenome\tModule_completeness_(%)\n"
    svf.write(header)
    for mod,gencomp in modulepres.items():
        for gen,comp in gencomp.items():
            tw = mod + "\t" + gen + "\t" + "{:.2f}".format(comp) + "\n"
            svf.write(tw)

##methanogenesis module M00617 is a combination module of other methanogenesis modules; it is present if any of those are present: M00567,M00357,M00356,M00563 (sum up per genome)
methanfix = modulepres["M00617"]
for gen,frac in methanfix.items():
    allfracs = []
    for m in ["M00567","M00357","M00356","M00563"]:
        currfrac = modulepres[m][gen]
        allfracs.append(currfrac)
    sumfrac = sum(allfracs)
    modulepres["M00617"][gen] = sumfrac

with open(save_pergenome,"w+") as svf:
    header = "Module\tGenome\tModule_completeness_(%)\n"
    svf.write(header)
    for mod,gencomp in modulepres.items():
        for gen,comp in gencomp.items():
            tw = mod + "\t" + gen + "\t" + "{:.2f}".format(comp) + "\n"
            svf.write(tw)

print("Module completeness per genome calculated")

###4. calculate module completeness in # genomes: 0%, <50%, >=50%-<75%, >=75%-<100%,100% completeness + add taxonomy->class (for each category) + add ko cutoffs to the table
genome2class = {}
with open(taxonomy) as tf:
    for line in tf:
        tsp = line.strip().split("\t")
        gen = tsp[0]
        class_ = tsp[3]
        genome2class[gen] = class_

modnamedict = {}
with open(modulenamefile) as mnf:
    next(mnf)
    for line in mnf:
        mnsp = line.strip().split("\t")
        modfull = mnsp[1].split("|")
        modnamedict[modfull[0]] = modfull[1]

cutoff_dict = {}
stepsnumdict = {}
konumdict = {}
with open(cutoffs) as cf:
    for line in cf:
        csp = line.strip().split("\t")
        mod = csp[0]
        cuts = csp[1]
        cutoff_dict[mod] = cuts
        cs = cuts.split("|")
        steps = []
        kos = set()
        for c in cs:
            csplit = c.split(":")[0]
            steps.append(csplit)
            if "," in csplit:
                ccomma = csplit.split(",")
                for cspp in ccomma:
                    if "+" in cspp:
                        cplus = cspp.split("+")
                        for cp in cplus:
                            kos.add(cp)
                    else:
                        kos.add(cspp)
            else:
                if "+" in csplit:
                    cplus = csplit.split("+")
                    for cp in cplus:
                        kos.add(cp)
                else:
                    kos.add(csplit)
        stepsnumdict[mod] = len(steps)
        konumdict[mod] = len(kos)

modpres_calc = {}
modpres_tax = {}
with open(save_pergenome) as mf:
    next(mf)
    for line in mf:
        msp = line.strip().split("\t")
        mod = msp[0]
        gen = msp[1]
        pres = float(msp[2])
        tax = genome2class[gen]
        if mod not in modpres_calc:
            modpres_tax[mod] = {"0%": [], "<50%": [], "<75%": [], ">=75%": [], "100%": []}
            modpres_calc[mod] = {"0%": 0, "<50%": 0, "<75%": 0, ">=75%": 0, "100%": 0}
            if pres == 0:
                modpres_calc[mod]["0%"] += 1
                modpres_tax[mod]["0%"].append(tax)
            elif pres == 1:
                modpres_calc[mod]["100%"] += 1
                modpres_tax[mod]["100%"].append(tax)
            elif pres < 0.5:
                modpres_calc[mod]["<50%"] += 1
                modpres_tax[mod]["<50%"].append(tax)
            elif pres < 0.75:
                modpres_calc[mod]["<75%"] += 1
                modpres_tax[mod]["<75%"].append(tax)
            elif pres >= 0.75:
                modpres_calc[mod][">=75%"] += 1
                modpres_tax[mod][">=75%"].append(tax)
            else:
                print("PROBLEM: ", msp)
        else:
            if pres == 0:
                modpres_calc[mod]["0%"] += 1
                modpres_tax[mod]["0%"].append(tax)
            elif pres == 1:
                modpres_calc[mod]["100%"] += 1
                modpres_tax[mod]["100%"].append(tax)
            elif pres < 0.5:
                modpres_calc[mod]["<50%"] += 1
                modpres_tax[mod]["<50%"].append(tax)
            elif pres < 0.75:
                modpres_calc[mod]["<75%"] += 1
                modpres_tax[mod]["<75%"].append(tax)
            elif pres >= 0.75:
                modpres_calc[mod][">=75%"] += 1
                modpres_tax[mod][">=75%"].append(tax)
            else:
                print("PROBLEM: ", msp)

mod_tax_count = {}
from collections import Counter
for m,taxpres in modpres_tax.items():
    mod_tax_count[m] = {}
    for perc,taxs in taxpres.items():
        mod_tax_count[m][perc] = "-"
        counttaxs = Counter(taxs)
        sorted_counttaxs = {k: v for k, v in sorted(counttaxs.items(), key=lambda pair: pair[1], reverse=True)}
        mergedcount = []
        for k,v in sorted_counttaxs.items():
            kv = k + ":" + str(v)
            mergedcount.append(kv)
        if mergedcount:
            mod_tax_count[m][perc] = "|".join(mergedcount)

with open(savefile_perc, "w+") as sf:
    header = "Module\t#_KOs_in_module\t#Reactions_in_module\t0%_present\t<50%_present\t>=50%-<75%_present" \
             "\t>=75%-<100%_present\t100%_present\tKO_cutoffs\t" \
             "Taxonomy_(Class)_0%_present\tTaxonomy_(Class)_<50%_present\tTaxonomy_(Class)_>=50%-<75%_present" \
             "\tTaxonomy_(Class)_>=75%-<100%_present\tTaxonomy_(Class)_100%_present\n"
    sf.write(header)
    for mod, counts in modpres_calc.items():
        modname = modnamedict[mod]
        fullmod = mod + "|" + modname
        zerocount = counts["0%"]
        modkosnum = konumdict[mod]
        modstepsnum = stepsnumdict[mod]
        modcuts = cutoff_dict[mod]
        modtax = mod_tax_count[mod]
        tw = fullmod + "\t" + str(modkosnum) + "\t" + str(modstepsnum) + "\t" + str(counts["0%"]) + "\t" + str(counts["<50%"]) + "\t" + \
             str(counts["<75%"]) + "\t" + str(counts[">=75%"]) + "\t" + str(counts["100%"]) + "\t" + \
             modcuts + "\t" + modtax["0%"] + "\t" + modtax["<50%"] + "\t" + modtax["<75%"] + "\t" + modtax[">=75%"] + "\t" + modtax["100%"] + "\n"
        sf.write(tw)

print("Module completeness in percent calculated")
