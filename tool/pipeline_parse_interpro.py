###parse raw interpro output (in standard format)
import sys
import os
genomefile = sys.argv[1] #input list of genomes
savefile = sys.argv[2] #output name

currdir = os.path.dirname(genomefile)

annoaccdict = {}
iprdict = {}
dbs = []
with open(genomefile) as gf:
    for l in gf:
        currgen = l.strip().split("\t")[0] + "_interpro.tsv"
        interprofile = os.path.join(currdir,currgen)
        print("Start parsing interpro ", interprofile)
        with open(interprofile) as intf:
            for line in intf:
                isp = line.strip().split("\t")
                acc = isp[0]
                annodb = isp[3]
                annoacc = isp[4]
                if annodb not in dbs:
                    dbs.append(annodb)
                annoname = isp[5].replace(" ", "_")
                ipracc = isp[11]
                ipranno = isp[12].replace(" ", "_")
                jipr = "|".join([ipracc,ipranno])
                janno = "|".join([annoacc,annoname])

                ###interpro dict
                if acc not in iprdict:
                    iprdict[acc] = [jipr]
                else:
                    iprdict[acc].append(jipr)

                #all other annos dict
                if acc not in annoaccdict:
                    annoaccdict[acc] = {annodb:[janno]}
                else:
                    if annodb not in annoaccdict[acc]:
                        annoaccdict[acc][annodb] = [janno]
                    else:
                        annoaccdict[acc][annodb].append(janno)

##dbs: 'CDD', 'Pfam', 'Gene3D', 'PANTHER', 'SUPERFAMILY', 'ProSitePatterns', 'NCBIfam', 'FunFam', 'Hamap', 'PIRSF', 'ProSiteProfiles', 'Coils', 'MobiDBLite', 'SMART', 'PRINTS'
with open(savefile,"w+") as svf:
    header = "Accession\tSUPERFAMILY\tInterpro\tPFAM\tCDD\tNCBIfam\tGene3D\tPANTHER\tProSitePatterns\tFunFam\tHamap\tPIRSF\tProSiteProfiles\tCoils\tMobiDBLite\tSMART\tPRINTS\n"
    svf.write(header)
    for acc in annoaccdict.keys():
        curr_interpro = ";".join(iprdict[acc])
        curr_annos = annoaccdict[acc]
        curr_pfam = "-"
        curr_cdd = "-"
        curr_gene3d = "-"
        curr_panther = "-"
        curr_superfamily = "-"
        curr_prositepatterns = "-"
        curr_ncbifam = "-"
        curr_funfam = "-"
        curr_hamap = "-"
        curr_pirsf = "-"
        curr_prositeprofiles = "-"
        curr_coils = "-"
        curr_mobidblite = "-"
        curr_smart = "-"
        curr_prints = "-"
        for db,anno in curr_annos.items():
            addanno = ";".join(anno)
            if db == "CDD":
                curr_cdd = addanno
            elif db == "Pfam":
                curr_pfam = addanno
            elif db == "Gene3D":
                curr_gene3d = addanno
            elif db == "PANTHER":
                curr_panther = addanno
            elif db == "SUPERFAMILY":
                curr_superfamily = addanno
            elif db == "ProSitePatterns":
                curr_prositepatterns = addanno
            elif db == "NCBIfam":
                curr_ncbifam = addanno
            elif db == "FunFam":
                curr_funfam = addanno
            elif db == "Hamap":
                curr_hamap = addanno
            elif db == "PIRSF":
                curr_pirsf = addanno
            elif db == "ProSiteProfiles":
                curr_prositeprofiles = addanno
            elif db == "Coils":
                curr_coils = addanno
            elif db == "MobiDBLite":
                curr_mobidblite = addanno
            elif db == "SMART":
                curr_smart = addanno
            elif db == "PRINTS":
                curr_prints = addanno
        twl = [acc,curr_superfamily,curr_interpro,curr_pfam,curr_cdd,curr_ncbifam,curr_gene3d,curr_panther,curr_prositepatterns,
               curr_funfam,curr_hamap,curr_pirsf,curr_prositeprofiles,curr_coils,curr_mobidblite,curr_smart,curr_prints]
        tw = "\t".join(twl) + "\n"
        svf.write(tw)

print("Finished parsing interpro ", interprofile)