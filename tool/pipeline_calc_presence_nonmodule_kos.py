##create a table with KOs that have no affiliated module -> presence/absence per genome; then group counts per KO

import sys
from collections import Counter
characthits = sys.argv[1] #char hits
taxonomy = sys.argv[2] #taxonomy file
uniquelist = sys.argv[3] # nonmodule_kos_list.txt #list of kos with no module affiliation: KO|KO_name per line)
savepergenome = sys.argv[4] # nonmodule_ko_presence_per_genome.tsv
savefilecounts = sys.argv[5] # nonmodule_ko_presence_percentage.tsv

uniquekos = {}
with open(uniquelist) as uf:
    for line in uf:
        uniquekos[line.strip()] = 1

genome2class = {}
with open(taxonomy) as tf:
    for line in tf:
        tsp = line.strip().split("\t")
        gen = tsp[0]
        class_ = tsp[3]
        genome2class[gen] = class_

allgenomes = {}
unique_nonmodule_kos = {}
with open(characthits) as chf:
    for line in chf:
        chsp = line.strip().split("\t")
        ko = chsp[1].split("|")
        if ko == ["No_KO"]:
            continue
        kos = []
        kosnames = []
        gens = chsp[19].split("|")
        for gen in gens:
            allgenomes[gen] = 1
            for idx,item in enumerate(ko):
                if idx % 2:
                    kosnames.append(item)
                else:
                    kos.append(item)
            zipko = list(zip(kos,kosnames))
            mod = chsp[-1]
            if mod == "-":
                for zk in zipko:
                    jzk = "|".join(zk)
                    if gen not in unique_nonmodule_kos:
                        unique_nonmodule_kos[gen] = {jzk:1}
                    else:
                        if jzk not in unique_nonmodule_kos[gen]:
                            unique_nonmodule_kos[gen][jzk] = 1
                        else:
                            unique_nonmodule_kos[gen][jzk] += 1

counts_group = {}
with open(savepergenome, "w+") as svf:
    for gen in allgenomes:
        try:
            curr_hits = unique_nonmodule_kos[gen]
        except KeyError:
            curr_hits = []
        for ko in uniquekos:
            if ko not in counts_group:
                counts_group[ko] = []
            if ko in curr_hits:
                counts_group[ko].append(gen)
                pres = curr_hits[ko]
                tw = gen + "\t" + ko + "\t" + str(pres) + "\n"
                svf.write(tw)
            else:
                pres = "0"
                tw = gen + "\t" + ko + "\t" + pres + "\n"
                svf.write(tw)

count_tax = {}
for k,vs in counts_group.items():
    count_tax[k] = []
    for v in vs:
        curr_tax = genome2class[v]
        count_tax[k].append(curr_tax)

totalgenomes = len(genome2class.keys())
with open(savefilecounts,"w+") as sf:
    header = "Nonmodule_KO\tPresent_in_#_genomes\t_Absent_in_#_genomes\tTaxonomy_(class)\n"
    sf.write(header)
    for k,vs in counts_group.items():
        curr_tax = Counter(count_tax[k])
        sorted_currtax = {k: v for k, v in sorted(curr_tax.items(), key=lambda pair: pair[1], reverse=True)}
        jointaxlist = []
        for t,vv in sorted_currtax.items():
            jointaxlist.append(":".join([t,str(vv)]))
        jtax = "|".join(jointaxlist)
        curr_counts = len(vs)
        absentgens = totalgenomes - curr_counts
        tw = k + "\t" + str(curr_counts) + "\t" + str(absentgens) + "\t" + jtax + "\n"
        sf.write(tw)

print("Nonmodule KOs presence calculated")
