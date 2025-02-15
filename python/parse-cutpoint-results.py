import os
import sys
from utilsPy3 import *
from glob import glob
import subprocess, shlex

import csv


### old way: python3 /data/rodriguesrr/scripts/python/parse-cutpoint-results.py lkt-pfs
# python3 /data/rodriguesrr/scripts/python/parse-cutpoint-results.py lkt-pfs default
# python3 /data/rodriguesrr/scripts/python/parse-cutpoint-results.py analysis-results Multi-var-Cox_Regression


PATH = sys.argv[1]
files = glob(PATH + '/**/*Results.csv', recursive=True)
#print(files)

"""
# typically output of Results-formatted.csv from Multi-var-Cox_Regression.R
def pre_formatted_Multi_var_Cox_Regression(PATH, files):
  header = ','.join(["File","ID","HR","CI","P-value"])
  reslist = [header]
  for f in files:
    res = read_file(f)
    resout = [ f + ',' + i.rstrip() for i in res[1:] ] # ignore the header line and prefix filename to each line
    reslist.extend(resout)
  writeLIST_to_file(reslist, PATH+'-Multi-var-Cox_Regression-parsed-result.csv')
  print("Finished parsing formatted results from Multi_var_Cox_Regression.")
"""


def parse_Input_stats(f):
  print(f)
  if '/2/cutp/' in f:
    conts = readItemList(f, 1, sep=',')
    #print(conts)
    conts = [c.replace('"','') for c in conts]
    #print(conts)
    
    outp = '\t'.join(conts[1:]) # excludes the first word "Result"
    return(outp)
        


def parse_Results_stats(f):
  print(f)
  if '/2/cutp/' in f:
    hr_data = []
    p_values = {}
    in_p_value_section = False  # Flag to track when we reach the p-value section

    with open(f, 'r') as file:
        reader = csv.reader(file)
        headers = next(reader)  # Read header row

        # Read HR table rows until we reach the "P-value:" section
        for row in reader: # each row is a list, where elements are split based on the delimiter.
            if not row or "P-value:" in row[0]:  # Check for empty lines or P-value section start
                in_p_value_section = True
                break
            hr_data.append('\t'.join(row))

        # Process remaining lines for p-values
        for line in file:
            line = line.strip()
            if line and "p=" in line:
                key, value = line.split("=")
                p_values[key.strip()] = float(value)
                if key == 'Score_(logrank)_test_p':
                    hr_data.append(value)
                    

    outp = '\t'.join(hr_data)
    return outp


def new_default(PATH):
  files = glob(PATH + '/**/*-stats.csv', recursive=True)

  header = '\t'.join(["Dir", "Timestamp", "LKT", "Categories", "Method", "Cutpoint", "Biomarker < Cutpoint", "Biomarker >= Cutpoint", "Biomarker.Category", "HR", "HR.low", "HR.high", "Pr(>|z|)", "Score_(logrank)_test P-value"])  
  reslist = [header]
  
  INP_DICT = {}
  RES_DICT = {}
  # use dict to map/merge the correct stats
  for f in files:
    #print(f)
    if "Input-stats.csv" in f: # Input-stats.csv (similar format to the old Results.csv files)
      outp = parse_Input_stats(f)
      INP_DICT[f.replace('/Input-stats.csv', '')] = outp
    elif "Results-stats.csv" in f: # Results-stats.csv (dataframe followed by the lines of p-values)
      outp = parse_Results_stats(f)
      RES_DICT[f.replace('/Results-stats.csv', '')] = outp
  
  # Concatenate values for matching keys, then covert to list
  merged_dict = {key: INP_DICT[key] + '\t' + RES_DICT[key] for key in INP_DICT if key in RES_DICT}
  result_list = [key.replace('/', '\t') + "\t" + value for key, value in merged_dict.items()]

  reslist.extend(result_list)

  writeLIST_to_file(reslist, PATH+'-parsed-result.tsv')
  print("Finished parsing formatted results.")



def parse_categorical_Input_stats(f):
  print(f)

  lines = read_file(f)
  dict = list_to_dict_all_vals(lines, delim=',', joiner='\t', string="current", key_column=0)

  # Convert values to the desired format # '8 ( 16% )'
  #formatted_dict = {
  #  k: v if k == 'biomarker' else f"{v.split(',')[0]} ( {round(float(v.split(',')[1]))}% )"
  #  for k, v in dict.items()
  #}

  #print(dict)
  return(dict)


def parse_categorical_Results_stats(f):
  print(f)

  hr_data = []
  p_values = {}
  in_p_value_section = False  # Flag to track when we reach the p-value section

  with open(f, 'r') as file:
      reader = csv.reader(file)
      headers = next(reader)  # Read header row

      # Read HR table rows until we reach the "P-value:" section
      for row in reader: # each row is a list, where elements are split based on the delimiter.
          if not row or "P-value:" in row[0]:  # Check for empty lines or P-value section start
              in_p_value_section = True
              break
          hr_data.append(','.join(row))

      # Process remaining lines for p-values
      for line in file:
          line = line.strip()
          if line and "p=" in line:
              key, value = line.split("=")
              p_values[key.strip()] = float(value)
              if key == 'Score_(logrank)_test_p':
                  hr_data = [val + ',' + value for val in hr_data] # join at the end of every row in the list   #hr_data.append(value)
                  
  
  dict = list_to_dict_all_vals(hr_data, delim=',', joiner='\t', string="current", key_column=0)
  #print(dict)
  return(dict)




def merge_and_format_data(dict1, dict2):
    result = []

    for group, biomarkers in dict1.items():
        if 'biomarker' not in biomarkers:
            continue

        # Extract available biomarkers and reference biomarker
        biomarkers_set = sorted(set(biomarkers.keys()) - {'biomarker'})  # Sorted for consistency
        existing_keys = set(dict2.get(group, {}).keys())
        reference_biomarker = next((b for b in biomarkers_set if f'x{b}' not in existing_keys), None)


        # Skip if no reference biomarker is found
        if not reference_biomarker:
            continue

        for biomarker in biomarkers_set:
            if biomarker == reference_biomarker:
                continue  # Skip reference biomarker itself

            # Get Freq and Percent values
            freq1, percent1 = map(float, biomarkers[biomarker].split(','))
            freq2, percent2 = map(float, biomarkers[reference_biomarker].split(',')) 

            # Get HR values from dict2 or set to NaN if not available
            hr_values = dict2.get(group, {}).get(f'x{biomarker}', 'NaN,NaN,NaN,NaN,NaN').split(',')
            hr, hr_low, hr_high, p_value, score_value = hr_values

            # Extract Dir, Timestamp, and Category from the group key
            dir_value, timestamp, category = group.split('/')[:3]

            # Append formatted row
            result.append([dir_value, timestamp, category, biomarker, reference_biomarker, freq1, freq2, percent1, percent2, hr, hr_low, hr_high, p_value, score_value])

    return result


def output_csv(file_path, data, header):
    # Writing the header and data to CSV
    with open(file_path, mode='w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow(header)
        writer.writerows(data)



def new_categorical(PATH):
  files = glob(PATH + '/**/*-stats.csv', recursive=True)
  header = ["Dir", "Timestamp", "LKT",  "biomarker1", "biomarker2", "Freq1", "Freq2", "Percent1", "Percent2", "HR", "HR.low", "HR.high", "Pr(>|z|)", "Score_(logrank)_test P-value"]
  
  INP_DICT = {}
  RES_DICT = {}
  # use dict to map/merge the correct stats
  for f in files:
    #print(f)
    if "Input-stats.csv" in f: # Input-stats.csv (similar format to the old Results.csv files)
      outp = parse_categorical_Input_stats(f)
      INP_DICT[f.replace('/Input-stats.csv', '')] = outp  
    elif "Results-stats.csv" in f: # Results-stats.csv (dataframe followed by the lines of p-values)
      outp = parse_categorical_Results_stats(f)
      RES_DICT[f.replace('/Results-stats.csv', '')] = outp
  

  # Merge, format, and output CSV
  formatted_data = merge_and_format_data(INP_DICT, RES_DICT)
  output_csv(PATH+'-parsed-result.tsv', formatted_data, header)
  print("Finished parsing formatted results.") 



def new_Multi_var_Cox_Regression(PATH):
  files = glob(PATH + '/**/*-stats.csv', recursive=True)
  header = [','.join(["Dir", "Timestamp", "LKT",  "Filename", "Category", "HR", "HR.low", "HR.high", "Pr(>|z|)", "Score_(logrank)_test P-value"])]
  
  #INP_DICT = {}
  RES_DICT = {}
  # use dict to map/merge the correct stats
  for f in files:
    print(f)
    #if "Input-stats.csv" in f: # Input-stats.csv (similar format to the old Results.csv files)
    #  outp = parse_categorical_Input_stats(f)
    #  INP_DICT[f.replace('/Input-stats.csv', '')] = outp  
    if "Results-stats.csv" in f: # Results-stats.csv (dataframe followed by the lines of p-values)
      outp = parse_categorical_Results_stats(f)
      RES_DICT[f.replace('/', ',')] = outp
  

  # output CSV
  writeLIST_to_file(header, PATH+'-parsed-result.csv')
  with open(PATH+'-parsed-result.csv', "a", newline="") as file:
    writer = csv.writer(file)
    for row, values in RES_DICT.items():
        for key, value in values.items():
            writer.writerow(row.split(',') + key.split(',') + value.split(','))
  print("Finished parsing formatted results.") 



"""
# default, e.g., output from Evaluate_cutpoints/EvaluateCutpoints_multiple-analyses_analysis.R
def pre_formatted(PATH, files):
  header = '\t'.join(["Dir", "Timestamp", "LKT", "Categories", "Method", "Results", "Cutpoint", "Biomarker < Cutpoint", "Biomarker >= Cutpoint", "HR", "CI", "P-value"]) # "OS type", 
  reslist = [header]
  for f in files:
    #print(f)
    if '/2/cutp/' in f:
        conts = readItemList(f, 2, sep=',')
        #print(conts)
        conts = [c.replace('"','') for c in conts]
        #print(conts)
        conts[0] = f.replace('/', '\t')
        conts = '\t'.join(conts)
        #print(conts)
        reslist.append(conts)

  writeLIST_to_file(reslist, PATH+'-parsed-result.tsv')
  print("Finished parsing formatted results.")


# default, e.g., output from Evaluate_cutpoints/Cox_Regression_Plot_Kaplan_Meyer_for_categories_multiple-analyses.R
def pre_formatted_categorical(PATH, files):
  header = '\t'.join(["Dir", "Timestamp", "LKT", "Results", "biomarker1", "biomarker2", "Freq1", "Freq2", "Percent1", "Percent2", "HR", "CI", "P-value"]) #"OS type", 
  reslist = [header]
  for f in files:
    #print(f)
    conts = readItemList(f, 1, sep=',')
    #print(conts)
    if(len(conts)==10): # only two categories
        conts = [c.replace('"','') for c in conts]
        #print(conts)
        conts[0] = f.replace('/', '\t')
        conts = '\t'.join(conts)
        #print(conts)
        reslist.append(conts)

  writeLIST_to_file(reslist, PATH+'-parsed-result.tsv')
  print("Finished parsing formatted categorical results. Results with more than 2 categories are skipped.")



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

    
# typically output of Results.csv from Multi-var-Cox_Regression.R 
# or use-specific-cutpoints-from-one-cohort-for-another-cohort.R
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
"""

if (len(sys.argv) == 2 or sys.argv[2] == 'new.default'): # either don't provide sys.argv[2] or specify as default
    new_default(PATH)
elif (sys.argv[2] == 'new.categorical'):
    new_categorical(PATH)
elif (sys.argv[2] == 'new.Multi_var_Cox_Regression'):
    new_Multi_var_Cox_Regression(PATH)
    """
elif (sys.argv[2] == 'default'):
    pre_formatted(PATH, files)
elif (sys.argv[2] == 'categorical'):
    pre_formatted_categorical(PATH, files)
elif (sys.argv[2] == 'Multi-var-Cox_Regression'):
    files = glob(PATH + '/**/*Results-formatted.csv', recursive=True)
    pre_formatted_Multi_var_Cox_Regression(PATH, files)
elif (sys.argv[2] == 'unformatted'):
    not_formatted(PATH, files)
    """
else:
    sys.exit("Wrong arguments.")

print("Done!")
