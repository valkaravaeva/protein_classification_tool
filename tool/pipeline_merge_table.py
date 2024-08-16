#this is a script to create the interpro/kegg table for the new archaea genomes
import sys
import os
hmmer = sys.argv[1] #merged filtered hmmer file
interpro = sys.argv[2] #merged parsed interpro file
accgen = sys.argv[3] #accession - genome map
currdir = os.path.abspath("./")
taxfile = sys.argv[4] #taxonomy file
kopath = sys.argv[5] #"/lisc/scratch/ecogenomics/karavaeva/kegg/ko00001_level_C.tsv" #special kegg path file
savefile = sys.argv[6] #path to savefile

konamedict = {}
kopathwaydict = {}
with open(kopath) as kp:
    for line in kp:
        psp = line.strip().split("\t")
        ko = psp[0]
        koname = psp[1]
        path = psp[-1]
        kopathwaydict[ko] = path
        konamedict[ko] = koname

accgendict = {}
with open(accgen) as agn:
    for line in agn:
        sp = line.strip().split("\t")
        gen = sp[0]
        acc = sp[1]
        if acc not in accgendict:
            accgendict[acc] = [gen]
        else:
            if gen not in accgendict[acc]:
                accgendict[acc].append(gen)

taxdict = {}
with open(taxfile) as tn:
    for line in tn:
        sp = line.strip().split("\t")
        gen = sp[0]
        tax = [x.replace(" ", "_") for x in sp[1:]]
        taxdict[gen] = tax

kostdict = {}
with open(hmmer) as stf:
    for line in stf:
        sp = line.strip().split("\t")
        acc = sp[0]
        ko = sp[1]
        koname = konamedict[ko]
        kostdict[acc] = ko + "|" + koname

##Accession	KO_std	KO_20%	KEGG_pathway_std	KEGG_pathway_20%	SUPERFAMILY	Interpro	PFAM	CDD	NCBIfam	Gene3D	PANTHER	ProSitePatterns	FunFam	Hamap	PIRSProSiteProfiles	Coils	MobiDBLite	SMART	PRINTS	genome	domain	phylum	class	order	family	genus	species	strain	taxid
with open(savefile,"w+") as svf:
    header = "Accession\tKO_std\tKEGG_pathway_std\tSUPERFAMILY\tInterpro\tPFAM\tCDD\tNCBIfam\t" \
             "Gene3D\tPANTHER\tProSitePatterns\tFunFam\tHamap\tPIRSF\tPIRSProSiteProfiles\tCoils\tMobiDBLite\tSMART\tPRINTS\t" \
             "genome\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies\tstrain\n"
    svf.write(header)
    with open(interpro) as intf:
        next(intf)
        for line in intf:
            insp = line.strip().split("\t")
            acc = insp[0]
            interproanno = "\t".join(insp[1:])
            if acc in kostdict:
                kostd = kostdict[acc]
                onlykostd = kostd.split("|")[0]
                kostdpath = kopathwaydict[onlykostd]
            else:
                kostd = "No_KO"
                kostdpath = "No_KO"
            gens = accgendict[acc]
            phyls = []
            classes = []
            orders = []
            families = []
            genera = []
            speciess = []
            strains = []
            doms = []
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
            taxonomy = "|".join(gens) + "\t" + "\t".join(doms) + "\t" + "|".join(phyls) + "\t" + "|".join(classes) + "\t" + "|".join(orders) + "\t" + "|".join(families) + "\t" + "|".join(genera) + "\t" + "|".join(speciess) + "\t" + "|".join(strains)
            tw = acc + "\t" + kostd + "\t" + kostdpath + "\t" + interproanno + "\t" + taxonomy + "\n"
            svf.write(tw)
