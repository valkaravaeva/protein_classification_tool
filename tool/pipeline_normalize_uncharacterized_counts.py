#normalize the counts of uncharacterized proteins
# #unch / #unique_cds -> per genome; take mean/median per taxon
import os
from statistics import median,mean
import sys

dbfile = sys.argv[1] #path to taxonomy file
save_unch = sys.argv[2] #path to uncharacterized file

gen2cds = {}
currdir = os.path.dirname(dbfile)
unique_taxa = {"phylum":{},"class":{},"order":{},"family":{},"genus": {},"species":{},"strain":{},"genome":{}} # level:{taxname:[genomes]}
with open(dbfile) as dbf:
    for line in dbf:
        dbsp = line.strip().split("\t")
        gen = dbsp[0]
        currgenfile = gen + ".faa"
        currfaa = os.path.join(currdir,currgenfile)
        cds = 0
        with open(currfaa) as caf:
            for lfa in caf:
                if lfa.startswith(">"):
                    cds += 1
        gen2cds[gen] = cds
        dom = dbsp[1].replace(" ","_")
        phylum = dbsp[2].replace(" ","_")
        class_ = dbsp[3].replace(" ","_")
        order = dbsp[4].replace(" ","_")
        family = dbsp[5].replace(" ","_")
        genus = dbsp[6].replace(" ","_")
        species = dbsp[7].replace(" ","_")
        strain = dbsp[8].replace(" ","_")
        if species == "":
            if "Unclassified" in dbsp[6]:
                species = dbsp[6].replace(" ","_")
            else:
                species = "Unclassified_" + dbsp[6].replace(" ","_")
        if phylum not in unique_taxa["phylum"]:
            unique_taxa["phylum"][phylum] = [gen]
        else:
            if gen not in unique_taxa["phylum"][phylum]:
                unique_taxa["phylum"][phylum].append(gen)
        if class_ not in unique_taxa["class"]:
            unique_taxa["class"][class_] = [gen]
        else:
            if gen not in unique_taxa["class"][class_]:
                unique_taxa["class"][class_].append(gen)
        if order not in unique_taxa["order"]:
            unique_taxa["order"][order] = [gen]
        else:
            if gen not in unique_taxa["order"][order]:
                unique_taxa["order"][order].append(gen)
        if family not in unique_taxa["family"]:
            unique_taxa["family"][family] = [gen]
        else:
            if gen not in unique_taxa["family"][family]:
                unique_taxa["family"][family].append(gen)
        if genus not in unique_taxa["genus"]:
            unique_taxa["genus"][genus] = [gen]
        else:
            if gen not in unique_taxa["genus"][genus]:
                unique_taxa["genus"][genus].append(gen)
        if species not in unique_taxa["species"]:
            unique_taxa["species"][species] = [gen]
        else:
            if gen not in unique_taxa["species"][species]:
                unique_taxa["species"][species].append(gen)
        if strain not in unique_taxa["strain"]:
            unique_taxa["strain"][strain] = [gen]
        else:
            if gen not in unique_taxa["strain"][strain]:
                unique_taxa["strain"][strain].append(gen)
        if gen not in unique_taxa["genome"]:
            unique_taxa["genome"][gen] = [gen]
        else:
            if gen not in unique_taxa["genome"][gen]:
                unique_taxa["genome"][gen].append(gen)

gen2unch_counts = {}
with open(save_unch) as uf:
    for line in uf:
        usp = line.strip().split("\t")
        gens = usp[19].split("|")
        for gen in gens:
            if gen not in gen2unch_counts:
                gen2unch_counts[gen] = 1
            else:
                gen2unch_counts[gen] += 1

norm_gen2unch = {}
for k,v in gen2unch_counts.items():
    curr_cds = gen2cds[k]
    normvalue = v / curr_cds
    norm_gen2unch[k] = normvalue

for lvl,taxa in unique_taxa.items():
    curr_savefile = "uncharacterized_counts_normalized_per_" + lvl + ".tsv"
    with open(curr_savefile,"w+") as svf:
        header = lvl + "\tmedian_count\tmean_count\tmax_count\tmin_count\tnum_genomes\n"
        svf.write(header)
        for tax,gens in taxa.items():
            numgens = len(gens)
            taxcounts = []
            for gen in gens:
                curr_counts = norm_gen2unch[gen]
                taxcounts.append(curr_counts)
            curr_median = median(taxcounts)
            curr_mean = mean(taxcounts)
            curr_max = max(taxcounts)
            curr_min = min(taxcounts)
            tw = tax + "\t" + "{:.4f}".format(curr_median) + "\t" + "{:.4f}".format(curr_mean) + "\t" + "{:.4f}".format(curr_max) + "\t" + "{:.4f}".format(curr_min) + "\t" + str(numgens) + "\n"
            svf.write(tw)
