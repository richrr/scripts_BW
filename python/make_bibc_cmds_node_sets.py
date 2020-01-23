import os
import sys
from utilsPy3 import *
import pandas as pd

# args:
# node file of bibc: 2 column tab delimited with node name and node type
# node set file: 3 column tab delimited with:
    ## name of the node set
    ## second column takes entries as "fixed" , "type" or "node"
    ## third column is "group" id to match (for fixed and type) or comma separated list of node ids (for node)
    ## each line is a set of nodes that will be the target group, everything else is "Ignore"
# edge file of bibc: 3 column tab delimited with p1, p2, and combined coeff.
'''
ml python/3.6
cd /data/rodriguesrr/Marie_Vetizou/MV023M/analysis/transnet-FF-FM/all/corr/net/v2/cyto/separateFishers/bibcs/edges-w-mo
python /data/rodriguesrr/scripts/python/make_bibc_cmds_node_sets.py nodes.tsv nodeset.txt edges.tsv
ml unload python/3.6
'''

########### edit this description ############
# this code takes a node file and makes a new node for each type mentioned and makes
# the remaining nodes in that type to ignore. mostly useful to make per pheno node files
### see if/how this code is different than the make-bibc-cmds-per-node.R ###
###############################################

node_file = sys.argv[1]
#out_node_file = node_file.split('/')[-1]
print(node_file)

#type = sys.argv[2]
#print type

nodeset_file = sys.argv[2]
print(nodeset_file)

edge_file = sys.argv[3]
print(edge_file)




data = pd.read_csv(node_file, sep="\t", header=0)
#print(data)


def create_dir(dirName):
  if not os.path.exists(dirName):
    os.mkdir(dirName)
    print("Directory " , dirName ,  " Created ")
  else:
    print("Directory " , dirName ,  " already exists")


nodesets = read_file(nodeset_file)

as_is = list() # these lines will not be changed in the output files
nodes_in_as_is = list()
cmds_list = list()


for line in nodesets:
    #nodes compiled so far
    nodes_seen_so_far = list()

    line = line.strip()
    if not line:
        continue
    print(line)
    filen, note, categ = line.split('\t')
    #print(filen)
    #print(note)
    #print(categ)
    if note == "fixed":
        tmp = data[data['type'] == categ]
        tmp.loc[:, 'type'] = 'otu'
        #print(tmp)

        # https://datatofish.com/convert-pandas-dataframe-to-list/
        #as_is.extend([tmp.columns.values.tolist()] + tmp.values.tolist())
        as_is.extend([['name', 'Class']] + tmp.values.tolist()) # the output column headers for the node file should be "name" and "Class"
        #print(as_is)

        nodes_in_as_is.extend(list(tmp['node']))
        continue

    grp = ''

    if note == "type":
        if '<==>' in categ:
            categs = categ.split('<==>')
            grp = data[data['type'].isin(categs)]
        else:
            grp = data[data['type'] == categ]


    if note == "node":
        categs = categ.split(',')
        categs = [item.strip('"') for item in categs]
        grp = data[data['node'].isin(categs)]

    grp.loc[:, 'type'] = 'gene'
    print(grp)

    nodes_seen_so_far.extend(nodes_in_as_is)
    nodes_seen_so_far.extend(list(grp['node']))

    nodes_not_seen_so_far = set(list(data['node'])) - set(nodes_seen_so_far)
    #print(nodes_not_seen_so_far)

    # everything else will be set to 'ignore'
    ignr = data[data['node'].isin(nodes_not_seen_so_far)]
    ignr.loc[:, 'type'] = 'ignore'
    #print(ignr)

    outlist_node = list()
    outlist_node.extend([as_is, grp.values.tolist() + ignr.values.tolist()])

    #print(outlist_node[[:]])

    flat_list = [item for sublist in outlist_node for item in sublist]
    #print(flat_list)

    out_list = [','.join(sublist) for sublist in flat_list]
    #print(out_list)

    filen = ''.join(e for e in filen if e.isalnum())
    create_dir(filen)
    new_file_name = filen+'/'+filen+".csv"
    writeLIST_to_file(out_list, new_file_name)

    ###### bibc commands ###############
    cmds1 = ' '.join(['Rscript /nfs3/PHARM/Morgun_Lab/richrr/scripts/R/prep-for-betw-central.R', edge_file, new_file_name, new_file_name])
    cmds2 = 'SGE_Batch -c "/nfs3/PHARM/Morgun_Lab/richrr/scripts/R/betweennessadapter node_type microbe degs ./' + new_file_name + '.graphml > bibc-' + filen + '-results.txt" -m 50G -F 50G -q samwise -M rodrrich@oregonstate.edu -r log_' + filen
    cmds_list.append(cmds1)
    cmds_list.append(cmds2)


writeLIST_to_file(cmds_list, node_file + '.bibc-cmds.txt')


sys.exit(0)
