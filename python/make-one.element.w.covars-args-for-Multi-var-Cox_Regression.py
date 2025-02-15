#ml python/3.10
import os
import sys
from utilsPy3 import *

#usage:
# cd /data/MicrobiomeCore/analysis/ADAD01_20240215/survival_analysis/multivar-lkt-diag.yr-swap.HSA.or.not
#python /data/rodriguesrr/scripts/python/make-one.element.w.covars-args-for-Multi-var-Cox_Regression.py ./ ../multivar-infile-lkt-diag.yr-swap.HSA.or.not.txt-file-for-Multi-var-Cox_Regression.txt X83765acd445533175c2385605d4f8581__D_0__Bacteria..D_1__Firmicutes..D_2__Clostridia..D_3__Clostridiales..D_4__Ruminococcaceae..D_5__Butyricicoccus..D_6__uncultured.bacterium,dc9ddb6b3c51c6926b26dc6d65ff66c3__D_0__Bacteria..D_1__Bacteroidetes..D_2__Bacteroidia..D_3__Bacteroidales..D_4__Bacteroidaceae..D_5__Bacteroides,X08bd43071a06abf56df235fca6655bea__D_0__Bacteria..D_1__Firmicutes..D_2__Clostridia..D_3__Clostridiales..D_4__Lachnospiraceae..D_5__.Ruminococcus..torques.group,Cancer_HSA_2_or_not_1,Age_at_diagnosis_years 91.66,1570,441.5,11.9,1.5 Cancer_HSA_2_or_not_1,Age_at_diagnosis_years


base_cmd = "Rscript /data/rodriguesrr/scripts/R/Evaluate_cutpoints/Multi-var-Cox_Regression.R"

odir = sys.argv[1]

infile = sys.argv[2]

#comma_sep_strings
lkts = sys.argv[3]
lkt_list = lkts.split(",")

cutps = sys.argv[4]
cutp_list = cutps.split(",")

covars = sys.argv[5]
covars_list = covars.split(",")


cmd_list = []
lkt_cutp_dict = dict(zip(lkt_list, cutp_list)) # make dictionary from two lists
for lkt, cutp in lkt_cutp_dict.items():
    if lkt not in covars_list:
        #print(lkt + " "  + cutp)
        #the order of input doesn't matter to the MV Cox_Regression
        #puts lkt first and then covariates.
        #arg_lkt = [ lkt ]
        #arg_lkt.extend(covars_list)
        #puts covariates first and then lkt.
        arg_lkt = list(covars_list)
        arg_lkt.append(lkt)
        #print(arg_lkt)
        arg_cutp = [ lkt_cutp_dict[i] for i in arg_lkt ]
        #print(arg_cutp)
        full_cmd = " ".join([ base_cmd , odir , infile , ",".join(arg_lkt) , ",".join(arg_cutp) ])
        cmd_list.append(full_cmd)

writeLIST_to_file(cmd_list, odir+'/run-cmds.txt')
print("Run cmd file.")
