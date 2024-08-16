#this is a script to filter arcogs by best hit
import os
import sys
taxfile = sys.argv[1] #taxonomy.tsv
outfile = sys.argv[2] #savefile

besthits = {}
currdir = os.path.dirname(taxfile)
with open(taxfile) as tf:
    for line in tf:
        currgen = line.strip().split("\t")[0] + "_arcogs.tsv"
        infile = os.path.join(currdir,currgen)
        with open(infile) as ff:
            for line in ff:
                sp = line.strip().split("\t")
                q = sp[0]
                s = sp[1]
                if "gi|" in s:
                    continue
                id = float(sp[2])
                ev = float(sp[10]) #10
                if id < 25:
                    continue
                if ev > 0.0000000001:
                    continue
                bitscore = float(sp[11]) #11
                if q not in besthits:
                    besthits[q] = [s,id,ev,bitscore]
                else:
                    curr_bh = besthits[q]
                    curr_id = curr_bh[1]
                    curr_ev = curr_bh[2]
                    curr_bt = curr_bh[3]
                    if bitscore > curr_bt:
                        besthits[q] = [s, id, ev, bitscore]
                    elif bitscore == curr_bt:
                        if id > curr_id:
                            besthits[q] = [s,id,ev,bitscore]
                        elif id == curr_id:
                            if ev <= curr_ev:
                                besthits[q] = [s, id, ev, bitscore]

with open(outfile,"w+") as svf:
    for k,v in besthits.items():
        strv = [str(vv) for vv in v]
        tw = k + "\t" + "\t".join(strv) + "\n"
        svf.write(tw)
