#-*- coding: utf-8 -*-

# Biowulf: ml python/3.8
# python /data/rodriguesrr/scripts/python/NN_bibc/import_network_data.py --input net.csv --source partner1 --target partner2
"""
Created on Mon Apr 29 12:01:46 2019

Written in Python v3.5.3

@author: Nolan K Newman <newmanno@oregonstate.edu>

This script takes a network file consisting of node correlations and outputs a pickled file of a networkx graph, as well as separate calculations if the user specifies the --dev arg
See README file for file input format

Developer notes:
	3/12/21 - Users can now specify which columns the source and target nodes are in. Also, if a user wants to calculate deviation from expected they must now supply the --dev command

"""

###################################################################################
 ##################### import network and make networkx object ###################
###################################################################################

import argparse
import csv
import numpy as np 
import pickle
import networkx as nx
from collections import Counter
from math import factorial

# Create a dictionary function to easily add keys and values to dictionary
class dictionary(dict):
    def __init__(self):
        self = dict()
        
    # Add key/value pair
    def add(self, key, value):
        self[key] = value

corr_dict = dictionary()
fc = {}

parser = argparse.ArgumentParser(description='Example: python import_network_data.py --input <network file> --source partner1 --target partner2 \n\n See README.md for more info\n\n')
parser.add_argument("--input", help = 'Network file (see README.md for example)')
parser.add_argument("--source", default = 'partner1', help = 'Name of the source node column for the edge')
parser.add_argument("--target", default = 'partner2', help = 'Name of the target node column for the edge')
parser.add_argument("--dev", action = 'store_true', help = 'Do you want to calculate deviation from expected? (Usually do not use this arg. Requires larger input file format)')
parser.add_argument("--num_groups", help = 'Usually do not use this arg. Number of groups correlations were initially performed in, required if --dev is supplied')

args = parser.parse_args()

net_file = args.input
net_file_trimmed = net_file[:-4] # trim the ".csv" or ".txt" from the input file string   

if args.dev:
    groups = int(args.num_groups)

# Counters for puc calculation
puc_compliant = 0
puc_noncompliant = 0

row_count = 0

# import specified file into python
with open(net_file) as csvfile:
    file = csv.reader(csvfile, delimiter = ',')
    for row in file:
            
        # Take the index of the source and target node in the header of the file
        if row_count == 0: 
            p1 = int(row.index(args.source))
            p2 = int(row.index(args.target)) + 1
        
        nodes = row[p1:p2]
        print(nodes)
        list_to_tuple = tuple(nodes)
        corr_dict.add(list_to_tuple,row[3:len(row)])
        
        if args.dev:
            fc_node1_column = 11 + groups
            fc_node2_column = 12 + groups
            
            # Find FC direction of each node
            fc[nodes[0]] = row[fc_node1_column].strip()
            fc[nodes[1]] = row[fc_node2_column].strip()
            
            # Is each edge PUC-compliant?
            puc_col = 14 + groups
            
            if row[puc_col].strip() == str(1):
                puc_compliant += 1
            elif row[puc_col].strip() == str(-1):
                puc_noncompliant += 1
            
        row_count += 1

    csvfile.close()

# This removes the header of the data frame essentially
print(corr_dict)
del corr_dict[args.source, args.target]

if args.dev:
    del fc['partner1']
    del fc['partner2']

# function to get unique values
def unique(list1):
    x = np.array(list1)
    return np.unique(x)

# add all nodes to a list
node_list = []
for key,value in corr_dict.items():
    node_list.append(key[0])
    node_list.append(key[1]) 
        
# find all nodes and store them in a new list
unique_nodes = list(unique(node_list))

G = nx.Graph() 


# Find the absolute value of each rho edge (to be added for node strength)
#print("Node1\tnode2\tabs(Edge weight)")
for key,value in corr_dict.items():
    #print(key[0] + "<==>" + key[1])
    #print(value[groups + 1])
    if args.dev:
        G.add_edge(key[0],key[1],weight = abs(float(value[groups + 1])))
    else:
        G.add_edge(key[0],key[1])


# Save to pickle
pickle_out = open(net_file_trimmed + ".pickle", "wb")

### If --dev is supplied as an argument, calculate the devation from expected
if args.dev:    
    # Get a dictionary of all correlation directions and count all positive and negative edges in observed network
    rho_column = 7 + groups
    pos_corr = 0 # counter for the number of positive edges
    neg_corr = 0 # counter for the number of negative edges
    
    for key,value in corr_dict.items():
        try:    
            if str(value[rho_column].strip()) == '1':
                pos_corr += 1
            elif str(value[rho_column].strip()) == '-1':
                neg_corr += 1
        except:
            print("ERROR: an incorrect value was supplied for correlation directions. Aborting.")         
        
    nedges = pos_corr + neg_corr       
    
    # Count the number of positive and negative nodes       
    nodedir = Counter(fc.values())
    pos_nodes = nodedir['1']
    neg_nodes = nodedir['-1']         
    total_nodes = int(pos_nodes) + int(neg_nodes)
    obs_edge_node_ratio = nedges / total_nodes
    
    obs_posneg_node_ratio = int(pos_nodes) / int(neg_nodes)
    obs_negpos_node_ratio = int(neg_nodes) / int(pos_nodes)
    
    # Find the ratio of positive:negative edges (and vice versa) in the observed graph
    if int(neg_corr) != 0:
        obs_posneg_ratio = int(pos_corr)/int(neg_corr)
    else:
        obs_posneg_ratio = 1.0
        
    if int(pos_corr) != 0:
        obs_negpos_ratio = int(neg_corr)/int(pos_corr)
    else:
        obs_negpos_ratio = 1.0
    
    # Find the number of edges in a full graph      
    expec_pos = int(factorial(pos_nodes)/(2 * factorial(pos_nodes - 2)) + factorial(neg_nodes)/(2 * factorial(neg_nodes - 2)))           
    expec_neg = pos_nodes * neg_nodes
    expec_total = expec_pos + expec_neg          
    expec_edge_node_ratio = expec_total / total_nodes
    
    # Find the ratio of positive:negative edges (and vice versa) in a full graph
    ideal_ratio_posneg = expec_pos/expec_neg
    ideal_ratio_negpos = expec_neg/expec_pos
    
    # Find the ratio of positive:negative edges (and vice versa) in a complete graph
    ideal_ratio_posneg = expec_pos/expec_neg
    ideal_ratio_negpos = expec_neg/expec_pos
    
    #Calculate the non-normalized deviation from the expected (full) graph
    dev_posneg = obs_posneg_ratio/ideal_ratio_posneg
    dev_negpos = obs_negpos_ratio/ideal_ratio_negpos
                          
    #Calculate the normalized deviation from the expected (full) graph
    dev_norm_posneg = (obs_posneg_ratio - ideal_ratio_posneg) / ideal_ratio_posneg
    dev_norm_negpos = (obs_negpos_ratio - ideal_ratio_negpos) / ideal_ratio_negpos
    
    # calculate the normalized deviation of the edge:node (density) from the full graph
    dens_dev = (abs(obs_edge_node_ratio - expec_edge_node_ratio)) / expec_edge_node_ratio
    
    # Calculate PUC (the proportion of edges that do not follow the expected direction)
    puc = puc_noncompliant / nx.number_of_edges(G)
    
    dev_dict = {}
    dev_dict['OBSERVED_number_nodes_positive'] = pos_nodes
    dev_dict['OBSERVED_number_nodes_negative'] = neg_nodes
    dev_dict['OBSERVED_number_edges_positive'] = pos_corr
    dev_dict['OBSERVED_number_edges_negative'] = neg_corr
    dev_dict['OBSERVED_ratio_pos_to_neg_nodes'] = round(obs_posneg_node_ratio, 2)
    dev_dict['OBSERVED_ratio_neg_to_pos_nodes'] = round(obs_negpos_node_ratio, 2)
    dev_dict['OBSERVED_edge_node_ratio'] = round(obs_edge_node_ratio, 2)
    dev_dict['OBSERVED_ratio_pos_to_neg_edges'] = round(obs_posneg_ratio, 2)
    dev_dict['OBSERVED_ratio_neg_to_pos_edges'] = round(obs_negpos_ratio, 2)
    dev_dict['IDEAL_total_number_edges_full_graph'] = expec_total
    dev_dict['IDEAL_number_positive_edges_full_graph'] = expec_pos
    dev_dict['IDEAL_number_negative_edges_full_graph'] = expec_neg
    dev_dict['IDEAL_density_full_graph'] = round(expec_edge_node_ratio, 2)
    dev_dict['IDEAL_ratio_pos_to_neg_edges'] = round(ideal_ratio_posneg, 2)
    dev_dict['IDEAL_ratio_neg_to_pos_edges'] = round(ideal_ratio_negpos, 2)
    
    dev_dict['DEVIATION_posneg_deviation_nonnormalized'] = round(dev_posneg, 2)
    dev_dict['DEVIATION_negpos_deviation_nonnormalized'] = round(dev_negpos, 2)
    dev_dict['DEVIATION_posneg_deviation_normalized'] = round(dev_norm_posneg, 2)
    dev_dict['DEVIATION_negpos_deviation_normalized'] = round(dev_norm_negpos, 2)
    dev_dict['DEVIATION_density_deviation'] = round(dens_dev, 2)
    dev_dict['PUC'] = round(puc, 2)
    dev_dict['PUC_noncompliant_edge_number'] = round(puc_noncompliant, 2)
    dev_dict['PUC_compliant_edge_number'] = round(puc_compliant, 2)

    pickle_list = [G, dev_dict] # output both tht network and deviation dictionary
    pickle.dump(pickle_list, pickle_out)
    pickle_out.close()

### otherwise just output the network in the pickle
else:
    pickle.dump(G, pickle_out)
    pickle_out.close()

print("There are " + str(G.number_of_nodes()) + " nodes in the network and " + str(G.number_of_edges()) + " edges in the network.")

if args.dev:
    print("Successfully saved all items to pickle. Now run calc_network_properties.py")
else:
    print("Successfully saved network to pickle.")
