#this is a script to prepare matrices for heatmaps for each module (manual modules) that has hits --> prepare total heatmap for now
import sys
import os
pergenome_presence = sys.argv[1]#path to "modules_presence_per_genome.tsv"
savedir = os.path.dirname(pergenome_presence)
totalpresence = os.path.join(savedir,"module_completeness_per_genome_with_taxonomy.tsv")
taxonomy = sys.argv[2] #path to taxonomy.tsv
modstats = sys.argv[3] #path to "module_completeness_percent.tsv"

taxdict = {}
count_taxa = {}
ranked_orders = set()
with open(taxonomy) as tf:
    for line in tf:
        tsp = line.strip().split("\t")
        genome = tsp[0]
        phyl = tsp[2]
        taxdict[genome] = phyl
        ranked_orders.add(phyl)

with open(totalpresence,"w+") as svf:
    header = "Genome\tModule\tCompleteness\tTaxon\n"
    svf.write(header)
    with open(pergenome_presence) as pf:
        for line in pf:
            if line.startswith("Module"):
                continue
            psp = line.strip().split("\t")
            mod = psp[0]
            gen = psp[1]
            comp = psp[2]
            twl = [gen,mod,comp,taxdict[gen]]
            tw = "\t".join(twl) + "\n"
            svf.write(tw)

##save each table per module, transform into wide format, counting genomes
##if more than >=75% of module present -> present
##if 2+ steps in module: if at least 2 steps present, and they constitute >= 25% of module (i.e. <=8 steps in total) -> partial; else, absent
##if 2 steps in module if 1 step present -> partial, if both -> present
##else: absent

mod_steps = {}
with open(modstats) as ms:
    next(ms)
    for line in ms:
        msp = line.strip().split("\t")
        m = msp[0].split("|")[0]
        step = msp[2]
        mod_steps[m] = int(step)

permodule = {}
with open(totalpresence) as tf:
    next(tf)
    for line in tf:
        tsp = line.strip().split("\t")
        mod = tsp[1]
        if mod not in permodule:
            permodule[mod] = [tsp]
        else:
            permodule[mod].append(tsp)

for mod,entries in permodule.items():
    currname = "matrix_module_completeness_per_taxon_" + mod + ".tsv"
    curr_savefile = os.path.join(savedir,currname)
    curr_counts = {}
    for e in entries:
        gen = e[0]
        comp = float(e[2])
        order = e[3]
        mod = e[1]
        msteps = mod_steps[mod]
        if order not in curr_counts:
            curr_counts[order] = {"Absent":0,"Partial":0,"Present":0}
        if msteps <= 2:
            if comp == 0:
                curr_counts[order]["Absent"] += 1
            elif comp >= 0.75:
                curr_counts[order]["Present"] += 1
            else:
                curr_counts[order]["Partial"] += 1
        elif msteps > 2:
            if comp == 0:
                curr_counts[order]["Absent"] += 1
            elif comp >= 0.75:
                curr_counts[order]["Present"] += 1
            else:
                if msteps <=8:
                    if comp >= 0.25:
                        curr_counts[order]["Partial"] += 1
                    else:
                        curr_counts[order]["Absent"] += 1
                elif msteps == 9:
                    if comp >= 0.33:
                        curr_counts[order]["Partial"] += 1
                    else:
                        curr_counts[order]["Absent"] += 1
                elif msteps == 11:
                    if comp >= 0.36:
                        curr_counts[order]["Partial"] += 1
                    else:
                        curr_counts[order]["Absent"] += 1
                elif msteps == 13:
                    if comp >= 0.3:
                        curr_counts[order]["Partial"] += 1
                    else:
                        curr_counts[order]["Absent"] += 1
    curr_counts_percent = {}
    for curr_ord,counts in curr_counts.items():
        if curr_ord not in curr_counts_percent:
            curr_counts_percent[curr_ord] = {"Absent":0,"Partial":0,"Present":0}
        curr_order_count = count_taxa[curr_ord]
        for frac,value in counts.items():
            curr_perc = value / curr_order_count * 100
            curr_counts_percent[curr_ord][frac] = curr_perc
    with open(curr_savefile,"w+") as sf:
        header = "Order\tPresent\tPartial\tAbsent\n"
        sf.write(header)
        for rank in ranked_orders:
            curr_c = curr_counts_percent[rank]
            twl = [rank,str(curr_c["Present"]),str(curr_c["Partial"]),str(curr_c["Absent"])]
            tw = "\t".join(twl) + "\n"
            sf.write(tw)

##put marker modules into one matrix, as they are one step "modules", "present" and "partial" modules are saved separately
marker_modules = ["M00XX1","M00XX3","M00XX4","M00175","M00X14","M00X15","M00X16","M00X17","M00X18",
                  "M00XX6","M00XX7","M00XX8","M00X10","M00X11","M00X12","M00X13","DsrAB","Qmo"]

sorted_taxa = []
whichcompletenesses = ["PRESENT","PARTIAL"]
for whichcompleteness in whichcompletenesses:
    savename = "matrix_markers" + whichcompleteness + ".tsv"
    save_marker_matrix = os.path.join(savedir,savename)
    save_marker_dict = {} #tax:values
    for m in marker_modules:
        print("Marker: ", m)
        currname = "matrix_manual_module_completeness_per_collapsed_taxon_" + m + ".tsv"
        currfile = os.path.join(savedir,currname)
        with open(currfile) as cf:
            next(cf)
            for line in cf:
                csp = line.strip().split("\t")
                tax = csp[0]
                sorted_taxa.append(tax)
                if whichcompleteness == "PRESENT":
                    present = csp[1]
                else:
                    present = csp[2]
                if tax not in save_marker_dict:
                    save_marker_dict[tax] = []
                    save_marker_dict[tax].append(present)
                else:
                    save_marker_dict[tax].append(present)
    used = set()
    with open(save_marker_matrix,"w+") as sfm:
        header = "Taxon\tArsenate_reductase\tbd_oxidase\tcyt_c_oxidase\tNitrogenase_Nif|Vnf\tNapAB\tNirK|NirS\tNorBC\tNosZ\tNarGH\tH2:polysulfide_oxidoreductase\tSoeAB\tSOR\tSqr\tSreABC\tTqoAD\tTtrAB\tDsrAB\tQmo\n"
        sfm.write(header)
        for t in sorted_taxa:
            currpres = save_marker_dict[t]
            tw = t + "\t" + "\t".join(currpres) + "\n"
            if tw not in used:
                used.add(tw)
                sfm.write(tw)
