import sys
import os
annofile = sys.argv[1]  # path to file from 01_archaea35k_merge_all_annotations.py
currdir = os.path.dirname(annofile)
savefile_matrix = sys.argv[2]  # where to save presence absence matrix output to
savefile_counts = sys.argv[3]  # path to save the file with counts
summary_save_counts = sys.argv[4]  # path to save the file with the summary of the counts

# calculate presence/absence
with open(savefile_matrix, "w+") as svf:
    with open(annofile) as af:
        for line in af:
            if line.startswith("Accession"):
                header = "Accession\tKO_std\tKEGG_pathway_std\tSUPERFAMILY\tInterpro\tPFAM\tCDD\tNCBIfam\t" \
                         "Gene3D\tPANTHER\tProSitePatterns\tFunFam\tHamap\tPIRSF\tPIRSProSiteProfiles\tCoils\tMobiDBLite\tSMART\tPRINTS\t" \
                         "genome\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies\tstrain\n"
                svf.write(header)
                continue
            asp = line.strip().split("\t")
            acc = asp[0]
            curr_presence = [0 for x in range(18)]
            kostd = asp[1]
            if kostd != "No_KO":
                curr_presence[0] = 1
                curr_presence[1] = 1
            supfam = asp[2]
            if supfam != "-":
                curr_presence[2] = 1
            ipr = asp[3]
            if ipr != "-|-":
                curr_presence[3] = 1
            pfam = asp[4]
            if pfam != "-":
                curr_presence[4] = 1
            cdd = asp[5]
            if cdd != "-":
                curr_presence[5] = 1
            tigr = asp[6]
            if tigr != "-":
                curr_presence[6] = 1
            gene3d = asp[7]
            if gene3d != "-":
                curr_presence[7] = 1
            panther = asp[8]
            if panther != "-":
                curr_presence[8] = 1
            prositepatterns = asp[9]
            if prositepatterns != "-":
                curr_presence[9] = 1
            funfam = asp[10]
            if funfam != "-":
                curr_presence[10] = 1
            hamap = asp[11]
            if hamap != "-":
                curr_presence[11] = 1
            pirsf = asp[12]
            if pirsf != "-":
                curr_presence[12] = 1
            prositeprofiles = asp[13]
            if prositeprofiles != "-":
                curr_presence[13] = 1
            coils = asp[14]
            if coils != "-":
                curr_presence[14] = 1
            mobidblite = asp[15]
            if mobidblite != "-":
                curr_presence[15] = 1
            smart = asp[16]
            if smart != "-":
                curr_presence[16] = 1
            prints = asp[17]
            if prints != "-":
                curr_presence[17] = 1
            curr_tw = acc + "\t" + "\t".join([str(x) for x in curr_presence]) + "\n"
            svf.write(curr_tw)

# count presence of annotations
annoidxdict = {"KOstd": 1, "KEGG_pathway_std": 2, "SUPERFAMILY": 3, "Interpro": 4, "PFAM": 5, "CDD": 6, "NCBIfam": 7,
               "Gene3D": 8, "PANTHER": 9, "ProSitePatterns": 10, "FunFam": 11, "Hamap": 12, "PIRSF": 13, "ProSiteProfiles": 14, "Coils": 15,
               "MobiDBLite": 16, "SMART": 17, "PRINTS": 18}

idxannodict = {v: k for k, v in annoidxdict.items()}

total = 0
presencedict = {}
simplecount = {}
with open(savefile_matrix) as mf:
    next(mf)
    for line in mf:
        total += 1
        msp = line.strip().split("\t")
        acc = msp[0]
        presencedict[acc] = []
        curr_presence = [int(x) for x in msp[1:]]
        for idx, i in enumerate(curr_presence):
            if i == 1:
                annopresent = idxannodict[idx + 1]
                if annopresent not in simplecount:
                    simplecount[annopresent] = 1
                else:
                    simplecount[annopresent] += 1
                if "KEGG_pathway" in annopresent:
                    continue
                presencedict[acc].append(annopresent)

sorted_simplecount = {k: v for k, v in sorted(simplecount.items(), key=lambda pair: pair[1], reverse=True)}

with open(summary_save_counts, "w+") as sf:
    header = "Annotations_present\tCount_occurences\tPercentage_out_of_total_" + str(total) + "\n"
    sf.write(header)
    for k, v in sorted_simplecount.items():
        percentage = v / total * 100
        tw = k + "\t" + str(v) + "\t" + "{:.2f}".format(percentage) + "\n"
        sf.write(tw)

sorted_presencedict = {}
for k, v in presencedict.items():
    if not v:
        jv = "NO_ANNOTATIONS_PRESENT"
    else:
        jv = "|".join(sorted(v))
    sorted_presencedict[k] = jv

countdict = {}
for k, v in sorted_presencedict.items():
    if v not in countdict:
        countdict[v] = 1
    else:
        countdict[v] += 1

sorted_countdict = {k: v for k, v in sorted(countdict.items(), key=lambda pair: pair[1], reverse=True)}

with open(savefile_counts, "w+") as svf:
    header = "Annotations_present\tCount_occurences\tPercentage_out_of_total_" + str(total) + "\n"
    svf.write(header)
    for k, v in sorted_countdict.items():
        percentage = v / total * 100
        tw = k + "\t" + str(v) + "\t" + "{:.2f}".format(percentage) + "\n"
        svf.write(tw)
