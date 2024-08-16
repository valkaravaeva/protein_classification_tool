#this is a script to create the table for arcog hits
import sys
import os
from Bio import SeqIO

wrkdir = sys.argv[1] #path to the directory of arcog database files
all_arcog = sys.argv[2] #your filtered arcog output file
taxfile = sys.argv[3] #taxonomy.tsv file
savefile = sys.argv[4] #path to save your file
arcog2acc_ar18 = os.path.join(wrkdir,"tmp.ar18/ar18.ar14.02.csv")
arcogdef_ar18 = os.path.join(wrkdir,"tmp.ar18/arCOGdef.tab")
arcog2acc_asgard = os.path.join(wrkdir,"asgard20/asCOGs.2023-05.csv")
arcogdef_asgard = os.path.join(wrkdir,"asgard20/asCOGs.2020-10.def.tab")
funclass = os.path.join(wrkdir,"funclass.tab")

print("Start creating accession - genome map")
currdir = os.path.dirname(taxfile)
accgendict = {}
with open(taxfile) as gf:
    for l in gf:
        currgenome = l.strip().split("\t")[0]
        faaname = currgenome + ".faa"
        genomefaa = os.path.join(currdir,faaname)
        with open(genomefaa) as pf:
            for rec in SeqIO.parse(pf,"fasta"):
                acc = rec.id
                if acc not in accgendict:
                    accgendict[acc] = [currgenome]
                else:
                    if currgenome not in accgendict[acc]:
                        accgendict[acc].append(currgenome)
print("Finished creating accession - genome map")

arcoghits = {}
with open(all_arcog) as of:
    for line in of:
        sp = line.strip().split("\t")
        acc = sp[0]
        if acc not in accgendict: #i ran arcog on the old set that was not filtered for quality, so i am excluding those hits
            continue
        arcogacc = sp[1]
        arcoghits[acc] = arcogacc

taxdict = {}
with open(taxfile) as to:
    for line in to:
        sp = line.strip().split("\t")
        gen = sp[0]
        tax = sp[1:]
        taxdict[gen] = tax

funclassdict = {}
with open(funclass) as fuf:
    for line in fuf:
        fsp = line.strip().split("\t")
        fc = fsp[0]
        class_ = "_".join(fsp[-1].split(" ")).replace(",","").replace(":","")
        funclassdict[fc] = class_

defdict = {}
with open(arcogdef_ar18, "rb") as def18:
    lines = [x.decode(errors='replace') for x in def18.readlines()]
    for line in lines:
        sp = line.strip().split("\t")
        arcog = sp[0]
        if len(sp[1]) == 1:
            fc = sp[1] + "_" + funclassdict[sp[1]]
        else:
            fcs = list(sp[1])
            fcl = []
            for f in fcs:
                ff = f + "_" + funclassdict[f]
                fcl.append(ff)
            fc = "|".join(fcl)
        genename = sp[2]
        if genename == "-" or genename == "":
            genedef = "_".join(sp[3].split(" ")).replace(",","").replace(":","")
        else:
            genedef = "_".join(sp[3].split(" ")).replace(",","").replace(":","") + "_" + genename
        defdict[arcog] = [genedef,fc]

with open(arcogdef_asgard,"rb") as defasg:
    lines = [x.decode(errors='replace') for x in defasg.readlines()]
    for line in lines:
        sp = line.strip().split("\t")
        ascog = sp[0]
        arcog = sp[1]
        if sp[2] == "-":
            fc = sp[2]
        elif len(sp[2]) == 1:
            fc = sp[2] + "_" + funclassdict[sp[2]]
        else:
            fcs = list(sp[2])
            fcl = []
            for f in fcs:
                ff = f + "_" + funclassdict[f]
                fcl.append(ff)
            fc = "|".join(fcl)
        genename = sp[3]
        if genename == "-" or genename == "":
            genedef = "_".join(sp[4].split(" ")).replace(",","").replace(":","")
        else:
            genedef = "_".join(sp[4].split(" ")).replace(",","").replace(":","") + "_" + genename
        defdict[ascog] = [arcog,genedef,fc]

modeldict = {}
with open(arcog2acc_ar18) as ar18:
    for line in ar18:
        sp = line.strip().split(",")
        acc = sp[2]
        try:
            model = sp[6]
        except IndexError:
            continue
        if acc not in modeldict:
            modeldict[acc] = [model]
        else:
            if model not in modeldict[acc]:
                modeldict[acc].append(model)

with open(arcog2acc_asgard) as asg:
    for line in asg:
        sp = line.strip().split(",")
        acc = sp[2]
        try:
            model = sp[6]
        except IndexError:
            continue
        if acc not in modeldict:
            modeldict[acc] = [model]
        else:
            if model not in modeldict[acc]:
                modeldict[acc].append(model)
print("dicts prepared")

with open(savefile,"w+") as svf:
    header = "Accession\tArcog_hit\tArcog_model\tArcog_definition\tArcog_functional_category\tGenome\tNCBI_taxid\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tStrain\n"
    svf.write(header)
    for acc,arcogacc in arcoghits.items():
        gens = accgendict[acc]
        doms = []
        phyls = []
        classes = []
        orders = []
        families = []
        genera = []
        speciess = []
        strains = []
        taxids = []
        for gen in gens:
            gentax = taxdict[gen]
            dom = gentax[0]
            if dom not in doms:
                doms.append(dom)
            phyl = gentax[1]
            if phyl not in phyls:
                phyls.append(phyl)
            class_ = gentax[2]
            if class_ not in classes:
                classes.append(class_)
            ord = gentax[3]
            if ord not in orders:
                orders.append(ord)
            fam = gentax[4]
            if fam not in families:
                families.append(fam)
            genus = gentax[5]
            if genus not in genera:
                genera.append(genus)
            species = gentax[6]
            if species not in speciess:
                speciess.append(species)
            strain = gentax[7]
            if strain not in strains:
                strains.append(strain)
        if arcogacc in modeldict: #some arcogs have no model
            arcogmodels = modeldict[arcogacc] #get the models for each hit; each arcog can have >1 model
            arcogdefs = []
            arcog_fcs = []
            for m in arcogmodels: #get def and fc for each model, ascogs also have corresponding arcogs in the def; some arcogs/ascogs have no definition
                if m in defdict:
                    if len(defdict[m]) == 2:
                        mdef = defdict[m][0]
                        mfc = defdict[m][1]
                        arcogdefs.append(mdef)
                        if mfc not in arcog_fcs:
                            arcog_fcs.append(mfc)
                    elif len(defdict[m]) == 3:
                        mdef = "_".join(defdict[m][:2])
                        mfc = defdict[m][2]
                        arcogdefs.append(mdef)
                        if mfc not in arcog_fcs:
                            arcog_fcs.append(mfc)
                else:
                    mdef = "no_definition_for_this_model"
                    mfc = "no_functional_category_for_this_model"
                    arcogdefs.append(mdef)
                    if mfc not in arcog_fcs:
                        arcog_fcs.append(mfc)
            arcogmodel = "|".join(arcogmodels)
            arcogdef = "|".join(arcogdefs)
            arcogfc = "|".join(arcog_fcs)
        else:
            arcogmodel = "no_model_for_this_arcog_sequence"
            arcogdef = "no_definition_for_this_model"
            arcogfc = "no_functional_category_for_this_model"
        twl = [acc,arcogacc,arcogmodel,arcogdef,arcogfc,"|".join(gens),"|".join(doms),"|".join(phyls),"|".join(classes),"|".join(orders),
               "|".join(families),"|".join(genera),"|".join(speciess),"|".join(strains)]
        tw = "\t".join(twl) + "\n"
        svf.write(tw)
