import os
import sys
from utilsPy3 import *
from glob import glob
import subprocess, shlex

# python3 /data/rodriguesrr/scripts/python/parse-cutpoint-results.py lkt-pfs

PATH = sys.argv[1]
files = glob(PATH + '/**/*Results.csv', recursive=True)
#print(files)

# default
def pre_formatted(PATH, files):
  header = '\t'.join(["OS type", "Dir", "Timestamp", "LKT", "Categories", "Method", "Results", "Cutpoint", "Biomarker < Cutpoint", "Biomarker > Cutpoint", "HR", "CI", "P-value"])
  reslist = [header]
  for f in files:
    print(f)
    conts = readItemList(f, 2, sep=',')
    conts = [c.replace('"','') for c in conts]
    conts[0] = f.replace('/', '\t')
    conts = '\t'.join(conts)
    print(conts)
    reslist.append(conts)

  writeLIST_to_file(reslist, PATH+'-parsed-result.tsv')
  print("Finished parsing formatted results.")




def get_val(f, pattern, rank):
    # https://stackoverflow.com/questions/24844707/escape-whitespaces-in-linux-path-and-file-names
    # https://stackoverflow.com/questions/7956865/python-subprocess-grep
    conts = subprocess.check_output("grep -m1 -A1 '"+ pattern +"' "+shlex.quote(f), shell=True)  #
    #print(conts)
    conts = conts.decode()  #  convert bytes to string https://stackoverflow.com/questions/606191/convert-bytes-to-a-string
    print(conts)
    contents = conts.split('\n')
    #print(contents)
    values = contents[1].split() # https://stackoverflow.com/questions/2492415/how-can-i-split-by-1-or-more-occurrences-of-a-delimiter-in-python
    #print(values)
    result = values[rank]
    
    # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # we get symbol for significant pvals and the actual pvalue when not significant
    # an alternative could be to use the lines with "wald test" or others to get actual pvalues, 
    # but they are rounded to two digits after decimal, so will increase FDR calulcations
    if(result in [".", "*", "**", "***"]): # get actual pvalue instead of symbol
        result = values[rank-1]
    #print(result)
    return(result)

    
# typically output of Multi-var-Cox_Regression.R or use-specific-cutpoints-from-one-cohort-for-another-cohort.R
def not_formatted(PATH, files): 
  header = '\t'.join(["Dir", "Timestamp", "LKT", "Results", "HR", "P-value"])
  reslist = [header]
  for f in files:
    print(f)
    
    HR = get_val(f, 'exp(coef)', 2)
    #print(HR)

    pval = get_val(f, 'Pr(>|z|)', -1)
    #print(pval)    
    
    conts = [f.replace('/', '\t'), HR, pval]
    conts = '\t'.join(conts)
    print(conts)
    
    reslist.append(conts)

  writeLIST_to_file(reslist, PATH+'-parsed-result.tsv')
  print("Finished parsing unformatted results.")



#pre_formatted(PATH, files)
not_formatted(PATH, files)

print("Done!")
