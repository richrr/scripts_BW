import sys
import os
from utils import *

# usage: python /data/rodriguesrr/scripts/python/map-JAMS-FeatTable.py /data/rodriguesrr/Marie_Vetizou/MV023M/data/JAMS_output_JM/PD1_MarieMouse_Inulin/Tables


infolder = sys.argv[1]

# need a one big file which is created by prefixing the filename (as column 1) to the contents of the FeatTable file
BIGLIST = ['Type\tID\tDescription']

files = os.listdir(infolder)
for f in files:
    if not f.startswith('.') and "_FeatTable.tsv" in f:
        print(f)
        if "_resfinder_" in f:
            continue
        elif "_LKT_" in f: # this is formated so that the LKT id (last column in this file) is used as the ID and the full line is joined to make the Description
            lines = read_file(f)
            tmp = []
            for l in lines:
                if 'Domain\tKingdom\tPhylum\tClass' not in l:
                    lktname = l.strip().split('\t')[-1]
                    lkt = l.strip().replace("\t", ";")
                    tmp.append(f+'\t'+lktname+'\t'+lkt)
            BIGLIST.extend(tmp)
            continue
        else: # these are assumed to be three column files
            lines = read_file(f)
            tmp = [ f+'\t'+l.strip() for l in lines if 'Accession\tDescription\n' not in l]
            BIGLIST.extend(tmp)


writeLIST_to_file(BIGLIST, infolder+"/FeatTable.Map.tsv")
