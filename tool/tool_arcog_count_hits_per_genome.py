# this is a script to count arcog hits per genome
import sys
import statistics
import os
from Bio import SeqIO

arcogfile = sys.argv[1] #arcog map file
taxfile = sys.argv[2] #taxonomy.tsv
savefile = sys.argv[3] #"arcog_hits_counts_per_genome.tsv"
savecounts = sys.argv[4] #"arcog_hits_mean_median_per_phylum.tsv"

count_accgendict = {}
taxdict = {}
currdir = os.path.dirname(taxfile)
with open(taxfile) as to:
    next(to)
    for line in to:
        sp = line.strip().split("\t")
        gen = sp[0]
        faaname = gen + ".faa"
        genomefaa = os.path.join(currdir, faaname)
        with open(genomefaa) as pf:
            for rec in SeqIO.parse(pf, "fasta"):
                acc = rec.id
                if gen not in count_accgendict:
                    count_accgendict[gen] = set()
                    count_accgendict[gen].add(acc)
                else:
                    count_accgendict[gen].add(acc)
        tax = sp[1:]
        taxdict[gen] = tax

counthitdict = {}
counthit_only_classified = {}
unch_arcog = {"S_Function_unknown", "no_functional_category_for_this_model", "4_POORLY_CHARACTERIZED",
              "R_General_function_prediction_only"}
with open(arcogfile) as agf:
    next(agf)
    for line in agf:
        sp = line.strip().split("\t")
        acc = sp[0]
        gens = sp[5].split("|")
        fcs = sp[4].split("|")
        for gen in gens:
            if gen not in counthitdict:
                counthitdict[gen] = set()
                counthitdict[gen].add(acc)
            else:
                counthitdict[gen].add(acc)
        if len(fcs) == 1:
            if fcs[0] not in unch_arcog:
                for gen in gens:
                    if gen not in counthit_only_classified:
                        counthit_only_classified[gen] = set()
                        counthit_only_classified[gen].add(acc)
                    else:
                        counthit_only_classified[gen].add(acc)
            else:
                continue
        else:
            if set(fcs).issubset(unch_arcog):
                continue
            else:
                for gen in gens:
                    if gen not in counthit_only_classified:
                        counthit_only_classified[gen] = set()
                        counthit_only_classified[gen].add(acc)
                    else:
                        counthit_only_classified[gen].add(acc)

uncl_per_phylum = {}
with open(savefile, "w+") as svf:
    header = "Genome\tTotal_sequences\t#_sequences_with_arcog_hits\t%_sequences_with_arcog_hits\t" \
             "#_sequences_with_classified_arcog_hits\t%_sequences_with_classified_arcog_hits\tPercent_unclassified\t" \
             "Domain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tStrain\n"
    svf.write(header)
    for gen, totalaccs in count_accgendict.items():
        tax = "\t".join(taxdict[gen])
        phyl = taxdict[gen][1].replace(" ", "_")
        total = len(totalaccs)
        genhits = len(counthitdict[gen])
        percent = (genhits / total) * 100
        classhits = len(counthit_only_classified[gen])
        perc_class = (classhits / total) * 100
        unchhits_perc = 100 - perc_class
        if phyl not in uncl_per_phylum:
            uncl_per_phylum[phyl] = [unchhits_perc]
        else:
            uncl_per_phylum[phyl].append(unchhits_perc)
        twl = [gen, total, genhits, percent, classhits, perc_class, unchhits_perc]
        tw = "\t".join([str(x) for x in twl]) + "\t" + tax + "\n"
        svf.write(tw)

# calc mean and median percentage of unclassified by arcog hits
with open(savecounts, "w+") as sf:
    header = "Phylum\tMean\tMedian\tMax\tMin\n"
    sf.write(header)
    for phyl, values in uncl_per_phylum.items():
        currmean = statistics.mean(values)
        currmedian = statistics.median(values)
        currmax = max(values)
        currmin = min(values)
        print(phyl, currmean, currmedian, currmax, currmin)
        twl = [phyl, currmean, currmedian, currmax, currmin]
        tw = "\t".join([str(x) for x in twl]) + "\n"
        sf.write(tw)
