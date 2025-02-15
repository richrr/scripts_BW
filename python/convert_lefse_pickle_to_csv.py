import os,sys,math,pickle
import pandas as pd

# usage: 
#cd /data/MicrobiomeCore/analysis/DMAH01_20180516/2020-reanalysis/stool_mas_hc_q2d2/exported_202010261223
#python /data/rodriguesrr/scripts/python/convert_lefse_pickle_to_csv.py otu_table_w_tax.txt.lefse.in DFfile.csv

#####################################################################
# WRITE DICTIONARY TO FILE, ARG: DICTIONARY , FILENAME
######################################################################
def writeDICT_to_file(DICT, filename, delim='\t', method='w'):
    outfile = open(filename, method)
    for k, v in DICT.items():
        string = delim.join(str(e) for e in v) # convert a list (of string or float) to string
        outfile.write("%s%s%s\n" % (k , delim, string))
    outfile.close()
    print('DICT written to file %s' % (filename))
    return


# https://github.com/abalter/lefse/blob/c8f2078c3501723f24a580a42f771bf9f5316f3b/lefse.py#L38
def load_data(input_file, nnorm = False):
    with open(input_file, 'rb') as inputf:
        inp = pickle.load(inputf)
    if nnorm: return inp['feats'],inp['cls'],inp['class_sl'],inp['subclass_sl'],inp['class_hierarchy'],inp['norm']
    else: return inp['feats'],inp['cls'],inp['class_sl'],inp['subclass_sl'],inp['class_hierarchy']

inpFile = sys.argv[1]
outFile = sys.argv[2]
# from https://github.com/abalter/lefse/blob/master/run_lefse.py
feats,cls,class_sl,subclass_sl,class_hierarchy = load_data(inpFile)

#print(feats)
#print(cls)

writeDICT_to_file(cls, outFile, delim=',', method='w')
writeDICT_to_file(feats, outFile, delim=',', method='a')

print("OK")
