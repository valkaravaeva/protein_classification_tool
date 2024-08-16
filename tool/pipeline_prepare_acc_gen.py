import sys
import os
from Bio import SeqIO

genomelist = sys.argv[1] #path to genome list file
savefile = sys.argv[2] #path to output file

print("Start creating accession - genome map")
currdir = os.path.dirname(genomelist)
with open(savefile,"w+") as svf:
    with open(genomelist) as gf:
        for l in gf:
            currgenome = l.strip().split("\t")[0]
            faaname = currgenome + ".faa"
            genomefaa = os.path.join(currdir,faaname)
            with open(genomefaa) as pf:
                for rec in SeqIO.parse(pf,"fasta"):
                    acc = rec.id
                    tw = currgenome + "\t" + acc + "\n"
                    svf.write(tw)

print("Finished creating accession - genome map")
