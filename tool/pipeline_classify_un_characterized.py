##classify sequences into "characterized" and "uncharacterized" as described in the publication

import sys
fulltable = sys.argv[1] #path to the full table containing all annotations
unchtigrs = sys.argv[2] #path to the manually curated list of uncharacterized TIGRFAMs
pthrexclfile = sys.argv[3] #path to the manually curated list of uncharacterized PANTHERS
accgentable = sys.argv[4] #path to a file that contains mappings of accessions to the genomes in a tabular format
dbfile = sys.argv[5] #path to the genome database file
save_unch = sys.argv[6] #path to save uncharacterized
save_ch = sys.argv[7] #path to save characterized
save_summary = sys.argv[8] #path to save summary counts

pthr_excl = set()
with open(pthrexclfile) as pef:
    for line in pef:
        pthr_excl.add(line.strip())

unch_tigrdict = {}
with open(unchtigrs) as utf:
    for line in utf:
        usp = line.strip().split("\t")
        utigr = usp[0].split(".")[0]
        unch_tigrdict[utigr] = 1

print("Unch tigrs and panthers loaded")

characterized = {}
uncharacterized = {}
with open(fulltable) as ff:
    next(ff)
    for line in ff:
        fsp = line.strip().split("\t")
        anno = fsp[1:]
        acc = fsp[0]
        ko = fsp[1]
        kegg_pathway = fsp[2]
        tigr = fsp[7]
        panther = fsp[9]
        hamap = fsp[12]
        if ko != "No_KO": #ko present
            if "99997_Function_unknown" in kegg_pathway:
                uncharacterized[acc] = anno #ko is uncharacterized
            elif "99996_General_function_prediction_only" in kegg_pathway:
                uncharacterized[acc] = anno  # ko is uncharacterized
            else:
                characterized[acc] = anno #ko is characterized
        else: #ko absent
            if panther != "-": #panther present
                pthracc = panther.split("|")[0]
                if pthracc not in pthr_excl:
                    characterized[acc] = anno #ko is characterized
                else:
                    uncharacterized[acc] = anno #ko is uncharacterized
            else: #ko,panther absent
                if tigr != "-": #tigr present
                    if tigr not in unch_tigrdict: #if tigr not uncharacterized
                        characterized[acc] = anno
                    else:
                        uncharacterized[acc] = anno
                else: #ko,panther,tigr absent
                    if hamap != "-": #hamap present
                        characterized[acc] = anno
                    else: #hamap,ko,panther,tigr absent
                        uncharacterized[acc] = anno

print("Classification loaded")

##check if there are sequences in genomes that had no annotations at all, add them to uncharacterized
domaindict = {}
phylumdict = {}
classdict = {}
orderdict = {}
familydict = {}
genusdict = {}
speciesdict = {}
straindict = {}
with open(dbfile) as dbf:
    for line in dbf:
        dbsp = line.strip().split("\t")
        genome = dbsp[0]
        dom = dbsp[1]
        phyl = dbsp[2].replace(" ","_")
        class_ = dbsp[3].replace(" ","_")
        ord = dbsp[4].replace(" ","_")
        fam = dbsp[5].replace(" ","_")
        genus = dbsp[6].replace(" ","_")
        species = dbsp[7].replace(" ","_")
        strain = dbsp[8].replace(" ","_")
        domaindict[genome] = dom
        phylumdict[genome] = phyl
        classdict[genome] = class_
        orderdict[genome] = ord
        familydict[genome] = fam
        genusdict[genome] = genus
        speciesdict[genome] = species
        straindict[genome] = strain

print("genomes loaded")

allaccs = {}
with open(accgentable) as af:
    for line in af:
        asp = line.strip().split("\t")
        gen = asp[0]
        acc = asp[1]
        if acc not in allaccs:
            allaccs[acc] = [gen]
        else:
            if gen not in allaccs[acc]:
                allaccs[acc].append(gen)

print("all accs loaded")

present = {}
archlen = 0
with open(fulltable) as ff:
    for line in ff:
        fsp = line.strip().split("\t")
        fsplen = len(fsp)
        facc = fsp[0]
        fgen = fsp[19]
        archlen = fsplen - 2
        present[facc] = 1

for k,v in allaccs.items():
    if k not in present:
        currgens = v
        curr_doms = []
        curr_phyls = []
        curr_classes = []
        curr_orders = []
        curr_fams = []
        curr_genera = []
        curr_species = []
        curr_strains = []
        for g in currgens:
            gdom = domaindict[g]
            if gdom not in curr_doms:
                curr_doms.append(gdom)
            gphyl = phylumdict[g]
            if gphyl not in curr_phyls:
                curr_phyls.append(gphyl)
            gclass = classdict[g]
            if gclass not in curr_classes:
                curr_classes.append(gclass)
            gord = orderdict[g]
            if gord not in curr_orders:
                curr_orders.append(gord)
            gfam = familydict[g]
            if gfam not in curr_fams:
                curr_fams.append(gfam)
            ggen = genusdict[g]
            if ggen not in curr_genera:
                curr_genera.append(ggen)
            gspec = speciesdict[g]
            if gspec not in curr_species:
                curr_species.append(gspec)
            gstr = straindict[g]
            if gstr not in curr_strains:
                curr_strains.append(gstr)
        vl = ["-"]*18 + ["|".join(currgens),"|".join(curr_doms),"|".join(curr_phyls),"|".join(curr_classes),"|".join(curr_orders),"|".join(curr_fams),"|".join(curr_genera),"|".join(curr_species),"|".join(curr_strains)]
        uncharacterized[k] = vl

with open(save_summary,"w+") as sumf:
    unchsv = "Uncharacterized: " + str(len(uncharacterized.keys())) + "\n"
    chsv = "Characterized: " + str(len(characterized.keys())) + "\n"
    sumsv = "U+C: " + str(len(uncharacterized.keys()) + len(characterized.keys())) + "\n"
    sumf.write(unchsv)
    sumf.write(chsv)
    sumf.write(sumsv)

with open(save_ch,"w+") as svf1:
    for k,v in characterized.items():
        tw = k + "\t" + "\t".join(v) + "\n"
        svf1.write(tw)

with open(save_unch,"w+") as svf2:
    for k,v in uncharacterized.items():
        tw = k + "\t" + "\t".join(v) + "\n"
        svf2.write(tw)
