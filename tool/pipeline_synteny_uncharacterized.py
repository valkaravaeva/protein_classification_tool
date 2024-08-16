#this is a script to do synteny analysis per uncharacterized protein, not continuous; exclude old genomes that were filtered out by quality
import sys
from pathlib import Path
import os

taxfile = sys.argv[1] #taxonomy file
unchfile = sys.argv[2] #files with uncharacterized hits
fulltable = sys.argv[3] #full table
savefile_merged = sys.argv[4] #"synteny_merged_table.tsv"
savefile_indexed = sys.argv[5] #"synteny_indexed_table.tsv"
savefile_neighbors = sys.argv[6] #"synteny_neighbors.tsv"
savefile_count_blocks = sys.argv[7] #"synteny_count_uncharacterized_blocks.tsv"
save_counts_overten = sys.argv[8] #"synteny_count_uncharacterized_blocks_more_than_ten.tsv"
save_counts_per_genome = sys.argv[9] #"synteny_count_uncharacterized_blocks_per_genome.tsv"

syntdir = os.path.dirname(taxfile)
pathlist_feat = Path(syntdir).glob('*_feature_table.txt') #convert each path to str!
pathlist_gff = Path(syntdir).glob('*.gff') #convert each path to str!

#create a merged table with all data for all genomes
#marked unchs; keep only necessary info -- add annotations later once the neighbors are extracted
unchhitsdict = {}
with open(unchfile) as unchf:
    for line in unchf:
        sp = line.strip().split("\t")
        acc = sp[0]
        unchhitsdict[acc] = sp
print("unch hits loaded")

filtered_old_genomes = set()
taxdict = {}
with open(taxfile) as tf:
    for line in tf:
        tsp = line.strip().split("\t")
        gen = tsp[0]
        class_ = tsp[3]
        taxdict[gen] = class_
        filtered_old_genomes.add(gen)

#unch?\tacc\tgenome\tchromosome\tstart\tend\tstrand\tdescription
with open(savefile_merged, "w+") as svf:
    header = "Uncharacterized\tAccession\tGenome\tChromosome\tStart\tEnd\tStrand\tDescription\n"
    svf.write(header)
    for file in pathlist_feat:
        curr_file = str(file)
        curr_genome = os.path.basename(curr_file).split("_feature_table.txt")[0]
        if curr_genome not in filtered_old_genomes:
            continue
        print("Processing: ", curr_genome)
        with open(curr_file) as cf:
            next(cf)
            for line in cf:
                sp = line.strip().split("\t")
                feature = sp[0]
                if feature != "gene":
                    chrom = sp[6]
                    start = sp[7]
                    end = sp[8]
                    strand = sp[9]
                    acc = sp[11]
                    if acc == "":
                        if sp[10] != "":
                            acc = sp[10]
                        else:
                            if "pseudo" in sp[-1]:
                                acc = "pseudogene" + "_" + start + "|" + end
                            elif "partial" in sp[-1]:
                                acc = "No_accession_" + sp[-1] + "_" + start + "|" + end
                            else:
                                acc = "No_accession" + "_" + start + "|" + end
                    descr = sp[13].replace(" ", "_").replace(",","")
                    if acc in unchhitsdict:
                        unchstat = "uncharacterized_hit"
                    else:
                        unchstat = "-"
                    twl = [unchstat,acc,curr_genome,chrom,start,end,strand,descr]
                    tw = "\t".join(twl) + "\n"
                    svf.write(tw)
    for file in pathlist_gff:
        curr_file = str(file)
        curr_genome = os.path.basename(curr_file).split(".gff")[0]
        print("Processing: ", curr_genome)
        with open(curr_file) as cf:
            next(cf)
            for line in cf:
                sp = line.strip().split("\t")
                chrom = sp[0]
                typeentry = sp[2]
                start = sp[3]
                end = sp[4]
                strand = sp[6]
                if typeentry == "CDS":
                    rest = sp[-1]
                    acc = rest.split(";")[0].split("=")[1]
                    if "product" in rest:
                        try:
                            descr = rest.split(";")[2].split("=")[1]
                        except:
                            descr = rest.split(";")[1].split("=")[1]
                    else:
                        descr = "no_description"
                else:
                    acc = typeentry
                    descr = typeentry
                if acc in unchhitsdict:
                    unchstat = "uncharacterized_hit"
                else:
                    unchstat = "-"
                twl = [unchstat,acc,curr_genome,chrom,start,end,strand,descr]
                tw = "\t".join(twl) + "\n"
                svf.write(tw)

##index the proteins, create merged accessions
idxdict = {}
annodict = {}
with open(savefile_merged) as mf:
    next(mf)
    for line in mf:
        msp = line.strip().split("\t")
        acc = msp[1]
        gen = msp[2]
        mergedacc = acc + "_" + gen
        annodict[mergedacc] = msp
        chrom = msp[3]
        if gen not in idxdict:
            idxdict[gen] = {chrom:[mergedacc]}
        else:
            if chrom not in idxdict[gen]:
                idxdict[gen][chrom] = [mergedacc]
            else:
                idxdict[gen][chrom].append(mergedacc)

with open(savefile_indexed,"w+") as idxf:
    header = "Uncharacterized\tAccession\tGenome\tMerged_accession\tChromosome\tIndex\tStart\tEnd\tStrand\tDescription\n"
    idxf.write(header)
    for genome,chrom2idx in idxdict.items():
        print("Processing: ", genome)
        for chr,accessions in chrom2idx.items():
            for idx,acc in enumerate(accessions):
                curr_anno = annodict[acc]
                twl = curr_anno[:3] + [acc] + [curr_anno[3],str(idx)] + curr_anno[4:]
                tw = "\t".join(twl) + "\n"
                idxf.write(tw)

# get neighbors and count blocks
total_genomes = len(taxdict.keys())
print("Total genomes: ", total_genomes)

hitdict = {} #hit: [gen,chr]
syntdict = {} #genome: {chr:{idx:merged_acc}}
with open(savefile_indexed) as syntf:
    next(syntf)
    for line in syntf:
        sp = line.strip().split("\t")
        hit = sp[0]
        prot = sp[1]
        gen = sp[2]
        merged_acc = sp[3]
        chr = sp[4]
        idx = int(sp[5])
        start = sp[6]
        start_acc = merged_acc + "|" + start
        if hit != "-":
            hitdict[start_acc] = [gen,chr,idx]
        if gen not in syntdict:
            syntdict[gen] = {}
            syntdict[gen][chr] = {idx:start_acc}
        else:
            if chr not in syntdict[gen]:
                syntdict[gen][chr] = {idx:start_acc}
            else:
                syntdict[gen][chr][idx] = start_acc

print("synteny loaded")

koorpantherdict = {}
with open(fulltable) as fuf:
    next(fuf)
    for line in fuf:
        fsp = line.strip().split("\t")
        acc = fsp[0]
        gens = fsp[19].split("|")
        ko = fsp[1].replace(";","").replace(",","")
        pthr = fsp[9].replace(";","").replace(",","")
        pfam = fsp[5].replace(";","").replace(",","")
        if ko != "No_KO":
            anno = ko + ":" + pfam
        else:
            anno = pthr + ":" + pfam
        for gen in gens:
            mergedacc = acc + "_" + gen
            if anno == "-:-":
                continue
            koorpantherdict[mergedacc] = anno

print("annotation loaded")

synteny = {}
synteny_blocks = {}
for hit,positions in hitdict.items():
    curr_genome = positions[0]
    curr_chromosome = positions[1]
    curr_idx = int(positions[2])
    curr_chr_neighborhood = syntdict[curr_genome][curr_chromosome]
    try:
        before2 = curr_chr_neighborhood[curr_idx-2]
        synteny[before2] = 1
    except:
        before2 = "End_chromosome/contig"
    try:
        before1 = curr_chr_neighborhood[curr_idx-1]
        synteny[before1] = 1
    except:
        before1 = "End_chromosome/contig"
    synteny[hit] = 1
    try:
        after1 = curr_chr_neighborhood[curr_idx+1]
        synteny[after1] = 1
    except:
        after1 = "End_chromosome/contig"
    try:
        after2 = curr_chr_neighborhood[curr_idx+2]
        synteny[after2] = 1
    except:
        after2 = "End_chromosome/contig"
    curr_block = [before2,before1,hit,after1,after2]
    if curr_genome not in synteny_blocks:
        synteny_blocks[curr_genome] = [curr_block]
    else:
        synteny_blocks[curr_genome].append(curr_block)

print("neighbors loaded")

unchs = set()
pseudoorrna = {}
with open(savefile_neighbors,"w+") as svf:
    with open(savefile_indexed) as sf:
        next(sf)
        for line in sf:
            ssp = line.strip().split("\t")
            unch = ssp[0]
            acc = ssp[1]
            gen = ssp[2]
            mergedacc = acc + "_" + gen
            start = ssp[6]
            start_acc = mergedacc + "|" + start
            if unch == "-":
                if mergedacc in koorpantherdict:
                    anno = koorpantherdict[mergedacc]
                else:
                    if "pseudogene" in mergedacc:
                        anno = "NO_KO|PANTHER_AND_PFAM_" + ssp[-1]
                        pseudoorrna[mergedacc] = ssp[-1]
                    else:
                        anno = "NO_KO|PANTHER_AND_PFAM_" + ssp[-1]
                        pseudoorrna[mergedacc] = ssp[-1]
            else:
                unchs.add(mergedacc)
                if mergedacc in koorpantherdict:
                    anno = "Uncharacterized_(" + koorpantherdict[mergedacc] + ")"
                else:
                    anno = "Uncharacterized_(NO_KO|PANTHER_AND_PFAM)"
            if start_acc in synteny:
                tw = gen + "\t" + anno + "\t" + line.strip() + "\n"
                svf.write(tw)

countblocks = {}
blocks_per_genome = {}
countblocks_per_genome = {}
for genome, blocks in synteny_blocks.items():
    for block in blocks:
        annoblock = []
        for pr in block:
            prm = pr.split("|")[0]
            if "End_chromosome" in prm:
                continue
            if len(pr.split("|")) > 2:
                prm = pr.split("|")[0] + "|" + pr.split("|")[1]
            if prm in koorpantherdict:
                pranno = koorpantherdict[prm]
            else:
                if prm in pseudoorrna:
                    pranno = "NO_KO|PANTHER_AND_PFAM_" + pseudoorrna[prm]
                else:
                    pranno = "NO_KO|PANTHER_AND_PFAM"
            if prm in unchs:
                fullanno = "Uncharacterized_(" + pranno + ")"
            else:
                fullanno = pranno
            annoblock.append(fullanno)
        jb = ";".join(annoblock)
        if jb not in countblocks:
            countblocks[jb] = 1
        else:
            countblocks[jb] += 1
        if jb not in blocks_per_genome:
            blocks_per_genome[jb] = [genome]
            countblocks_per_genome[jb] = 1
        else:
            countblocks_per_genome[jb] += 1
            if genome not in blocks_per_genome[jb]:
                blocks_per_genome[jb].append(genome)

sorted_countblocks = dict(sorted(countblocks.items(), key=lambda item: item[1], reverse=True))

from collections import Counter
with open(savefile_count_blocks,"w+") as cbf:
    header = "Synteny_block\tTotal_occurences\t#_genomes\t%_genomes\tClasses_of_genomes_(#_genomes)\n"
    cbf.write(header)
    for k,v in sorted_countblocks.items():
        curr_genomes = blocks_per_genome[k]
        curr_gencount = len(curr_genomes)
        curr_classes = []
        for gen in curr_genomes:
            curr_class = taxdict[gen]
            curr_classes.append(curr_class)
        countclass = dict(Counter(curr_classes))
        sorted_countclass = dict(sorted(countclass.items(), key=lambda item: item[1], reverse=True))
        twlclasses = []
        for kk,vv in sorted_countclass.items():
            kv = kk + "_(" + str(vv) + ")"
            twlclasses.append(kv)
        twclasses = ";".join(twlclasses)
        curr_genpercent = (curr_gencount / total_genomes) * 100
        tw = k + "\t" + str(v) + "\t" + str(curr_gencount) + "\t" + "{:.2f}".format(round(curr_genpercent, 2)) + "\t" + twclasses + "\n"
        cbf.write(tw)

with open(save_counts_overten,"w+") as tenf:
    header = "Synteny_block\tTotal_occurences\t#_genomes\t%_genomes\tClasses_of_genomes_(#_genomes)\n"
    tenf.write(header)
    for k,v in sorted_countblocks.items():
        if v >=10:
            curr_genomes = blocks_per_genome[k]
            curr_gencount = len(curr_genomes)
            curr_classes = []
            for gen in curr_genomes:
                curr_class = taxdict[gen]
                curr_classes.append(curr_class)
            countclass = dict(Counter(curr_classes))
            sorted_countclass = dict(sorted(countclass.items(), key=lambda item: item[1], reverse=True))
            twlclasses = []
            for kk, vv in sorted_countclass.items():
                kv = kk + "_(" + str(vv) + ")"
                twlclasses.append(kv)
            twclasses = ";".join(twlclasses)
            curr_genpercent = (curr_gencount / total_genomes) * 100
            tw = k + "\t" + str(v) + "\t" + str(curr_gencount) + "\t" + "{:.2f}".format(
                round(curr_genpercent, 2)) + "\t" + twclasses + "\n"
            tenf.write(tw)

with open(save_counts_per_genome,"w+") as sgf:
    for k,v in blocks_per_genome.items():
        curr_count = countblocks_per_genome[k]
        tw = k + "\t" + str(curr_count) + "\t" + ";".join(v) + "\n"
        sgf.write(tw)
