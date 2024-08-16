#prepare tables for uncharacterized boxplots: archaea and bacteria pipeline, archaea arcogs
#columns are phylum, all observations in a row
import sys
typedata = sys.argv[1] #options: arcog, pipeline
table = sys.argv[2] #path to arcog/pipeline table counts per genome
taxonomy = sys.argv[3] #path to taxonomy.tsv file
savefile = sys.argv[4] #path to save the file, boxplot_table.tsv or your custom name

if typedata == "arcog":
    ##arcogs
    with open(savefile,"w+") as svf:
        header = "Genome\tUncharacterized_percent\tTaxon\n"
        svf.write(header)
        with open(table) as tf:
            next(tf)
            for line in tf:
                tsp = line.strip().split("\t")
                gen = tsp[0]
                unch = tsp[6]
                phyl = tsp[8].replace(" ", "_")
                tw = gen + "\t" + unch + "\t" + phyl + "\n"
                svf.write(tw)
else:
    ##pipeline
    taxdict = {}
    with open(taxonomy) as tf:
        for line in tf:
            tsp = line.strip().split("\t")
            genome = tsp[0]
            phyl = tsp[2]
            taxdict[genome] = phyl

    with open(savefile,"w+") as svf:
        header = "Genome\tUncharacterized_percent\tTaxon\n"
        svf.write(header)
        with open(table) as tf:
            next(tf)
            for line in tf:
                tsp = line.strip().split("\t")
                gen = tsp[0]
                unch = tsp[1]
                phyl = taxdict[gen].replace(" ", "_")
                tw = gen + "\t" + unch + "\t" + phyl + "\n"
                svf.write(tw)
