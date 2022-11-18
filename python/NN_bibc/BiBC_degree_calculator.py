# -*- coding: utf-8 -*-

# Biowulf: ml python/3.8
# python /data/rodriguesrr/scripts/python/NN_bibc/BiBC_degree_calculator.py net.pickle --bibc_groups node_types --bibc_calc_type bibc --node_map map.csv --node_groups gene micro

"""
Created on Thu Aug 20 10:08:49 2020

@author: Nolan
"""

"""
Author: Nolan K Newman <newmanno@oregonstate.edu>
Last updated: 7/28/20

Written in Python v3.5.3

Description:
Takes as input the pickled network file created from import_network_data.py and calculates just the degree and BiBC for each node, which are two of the most critical node properties in finding a causal node of a disease being modeled by a network

Example usage:
python BiBC_degree_calculator.py <pickled file> --bibc --bibc_groups node_types --bibc_calc_type rbc --node_map <node map csv> --node_groups gene pheno

"""

import pickle
import networkx as nx
import numpy as np 
import re
import argparse
#import community  # commented by RR
#from community import community_louvain
from statistics import mean, median
from sortedcontainers import SortedDict
import csv
from collections import defaultdict
from networkx import all_shortest_paths
from collections import OrderedDict
from math import factorial

####### Get user input ########
parser = argparse.ArgumentParser(description='Example: python BiBC_degree_calculator.py <pickled network file> --bibc_groups node_types --bibc_calc_type bibc --node_map <node map csv> --node_groups <type1> <type2> --log')

# Required args
# pickle output from import_network_data2.py
parser.add_argument('pickle', help = 'The pickle file created with import_network_data.py')

# Flags and optional arguments

parser.add_argument("--bibc_groups", choices = ['node_types', 'modularity'], help = 'What to compute BiBC on, either distinct groups or on the two most modular regions of the network')
parser.add_argument("--bibc_calc_type", choices = ['rbc', 'bibc'], help = 'Would you like to normalize based on amount of nodes in each group (rbc) or not (bibc)?')
parser.add_argument("--node_map", help = 'Required if node_types is specified for --bibc_groups. CSV of nodes and their types (i.e. otu, pheno, gene, etc.)')
parser.add_argument("--node_groups", nargs = 2, help = 'Required if node_types is specified for --bibc_groups. Its the two groups of nodes to calculate BiBC/RBC on')
parser.add_argument("--log", action = 'store_true', help = 'Optional log file for progress of BiBC calculation')

parser.add_argument("--outputfileprefix", help = 'Output filename prefix')


args = parser.parse_args()

if args.bibc_groups == "node_types":
    bibc_choice = "node_types"
    node_input_file = args.node_map
    node_type1 = args.node_groups[0]
    node_type2 = args.node_groups[1]
elif args.bibc_groups == "modularity":
    bibc_choice = "modularity"

bibc_calc_type = args.bibc_calc_type
outputfilename = args.outputfileprefix if args.outputfileprefix is not None else node_input_file + "_"

if __name__ == '__main__':

    # Unpack the pickle
    p = open(args.pickle, "rb")
    p = pickle.load(p)
    
    
    if len(p) == 2:
        G = p[0] # Graph stored in first position of pickle
    else:
        G = p
        
    # try:
    #     G = p[0] # Graph stored in first position of pickle
    #     dev_dict = p[1] # Deviation dictionary stored in second position of pickle
    # except: 
    #     print("File being imported did not have deviations calculated")
    #     G = p

    # Function that takes as input the node_type list from the user and creates a dictionary of the type for each node. Then, for the nodes
    # that are in the netowrk input file, it assigns types to them based on the typing from the node_type file. It then outputs 2 dictionaries
    # of nodes, one that includes phenotypes and one that contains OTUs. These then get passed to the restricted_betweenness_centrality function,
    # where rbc is calculated for each node.
    def assign_node_type(node_list_file, gc_nodes, type1, type2):

        otu_and_pheno_dict = {}

        # Empty lists to hold all the OTUs in the node_type input file from the command line
        type1_list = []
        type2_list = []            

        node_type_dict = {}

        # Add all node-type pairs from the input file into the node_type_dict
        with open(node_list_file) as node_file:
            node_file = csv.reader(node_file, delimiter = ',')
            
            for row in node_file:
                node_type_dict[row[0].strip(' ')] = row[1].strip(' ')

        # Search the previously created dictionary and, for each 'otu' value in the second column of the input file, assign 
        # the corresponding key to otu_list, then do the same thing for each 'pheno' value and its corresponding list 
        for key,value in node_type_dict.items():
            try:
                if re.search(type1, value):
                    type1_list.append(key)
                elif re.match(type2, value):
                    type2_list.append(key)
            except:
                print("Unexpected value in the 'type' column of node_type input file.")

        # From the otu_list/pheno_list, only take the nodes that are present in the giant component. This is what the intersect1d function does. 
        # This is because I don't want to generate a new node type file for every network and this way I can keep using the same one. 

        type1_for_dict = np.intersect1d(type1_list, list(gc_nodes)) # list() added by RR
        print("\nCommon nodes between node_type1 input and giant component:")
        print(type1_for_dict)        

        type2_for_dict = np.intersect1d(type2_list, list(gc_nodes)) # list() added by RR
        print("\nCommon nodes between node_type2 input and giant component:")
        print(type2_for_dict)        

        # Add the nodes from node_type1 and node_type2 that are exclusive to the network to their respective dictionaries, 
        # then return the dictionaries and use them to call the restricted_betweenness_centrality function
        otu_and_pheno_dict['Type1'] = type1_for_dict
        otu_and_pheno_dict['Type2'] = type2_for_dict

        return(otu_and_pheno_dict)

    # Function that finds which nodes belong to the two most modular portions of the giant component, then returns
    # those nodes as a dictionary. These then get passed to the restricted_betweenness_centrality function below.
    def bibc_mod(nodes_from_gc):
        # Split the gc into 'clusters', clustering by the modularity
        part = community_louvain.best_partition(nodes_from_gc)
        
        # The previous method returns a dictionary, where keys are nodes and
        # values are the cluster they belong in. We want to find the two largest 
        # clusters, so we first find which unique clusters there are, then 
        # parse through those and assign nodes
        unique_clusters = set(part.values())

        # Add the clusters to a dictionary, where keys are clusters and values are nodes in that cluster
        mod_dict = {}

        for val in unique_clusters:
            nodes_in_cluster = []
            for key,value in part.items():
                if value == val:
                    nodes_in_cluster.append(key)
    
            mod_dict[val] = nodes_in_cluster

        # Sort the mod_dict by the length (size) of the clusters
        sorted_clust = sorted(mod_dict, key = lambda k: len(mod_dict[k]), reverse = True)
    
        # make lists of nodes in the largest clusters. sorted_clust[0] is the key of the largest cluster in mod_dict 
        large_mod = mod_dict[sorted_clust[0]]
        second_large_mod = mod_dict[sorted_clust[1]]

        # Add the lists to a dictionary, which then get returned to be passed to restricted_betweenness_centrality()
        return_mod = {}
        return_mod['mod1'] = large_mod
        return_mod['mod2'] = second_large_mod
        
        print("\n Network partitions")
        print("Largest GC partition (", str(len(return_mod['mod1'])), " nodes):")
        print(return_mod['mod1'])
        print("\n")
        print("Second largest GC partition (", str(len(return_mod['mod2'])), " nodes):")    
        print(return_mod['mod2'])
        print("\n")

        return(return_mod)        
    
    # Calculates BiBC (more correctly called restricted betweenness centrality) for each node
    def restricted_betweenness_centrality(G,nodes_0,nodes_1,type):
        '''
        Restricted betweenness centrality that only computes centrality using paths with sources
        in nodes_0 and targets in nodes_1 (or vice versa, which double counts).

        Returns three dictionaries of centralities, one for nodes_0, one for nodes_1,
        and one for the other nodes in G not in either set.

        If one thinks carefully about normalization:
        -centrality values for nodes in group 0 should be divided by (N0-1)N1
        -group 1: N0(N1-1)
        -others: N0*N1

        (Should just be able to normalize by len(nodes_0)*len(nodes_1))
        '''
        
        flatten = lambda l: [item for sublist in l for item in sublist]
        rbc = defaultdict(int)
        bibc_counter = 1
        g1_node_num = len(nodes_0)
        
        for s in nodes_0:
            pct_done = bibc_counter/g1_node_num * 100


            if args.log:
                with open(LOGFILENAME, "w") as logfile:
                    logfile.write("BiBC calculation status: %.1f percent completed" % pct_done)
            else:
                print("BiBC calculation status: %.1f percent completed" % pct_done)
                    
            for t in nodes_1:
                # might want to avoid putting the whole thing in memory
                # betweenness centrality does not count the endpoints (v not in s,t)
                paths_st = [x for x in list(all_shortest_paths(G,s,t)) if len(x) > 2]
                n_paths = len(paths_st)
                nodes_to_update = flatten([p[1:-1] for p in paths_st])
                for n in nodes_to_update:
                    rbc[n] += 1/n_paths
        
            bibc_counter += 1
        
        # split the dictionary in three
        rbc_n0 = {}.fromkeys(nodes_0)
        rbc_n1 = {}.fromkeys(nodes_1)
        rbc_other = {}
        
        for n in G.nodes():
            if n in nodes_0:
                rbc_n0[n] = rbc[n]
            elif n in nodes_1:
                rbc_n1[n] = rbc[n]
            else:
                rbc_other[n] = rbc[n]
    
        # If the user specifies rbc as the bibc calculation type, then normalize each node
        if type == "rbc":
            # Normalize each node - code suggested by Kevin and added by Nolan
            for i in rbc_n0:
                rbc_n0[i] = rbc_n0[i]/((len(nodes_0)-1) * len(nodes_1))
    
            for i in rbc_n1:
                rbc_n1[i] = rbc_n1[i]/((len(nodes_1)-1) * len(nodes_0))
    
            for i in rbc_other:
                rbc_other[i] = rbc_other[i]/(len(nodes_0) * len(nodes_1))

        elif type == "bibc":
            print("Normalization was not conducted on the set of nodes in the rbc function. BiBC was calculated instead.")

        # Return the list of BiBC of each node from each of the two groups
        return rbc_n0,rbc_n1,rbc_other


    ################################################################################
    ######################## Calculate network properties ##########################
    ################################################################################

    network_name = args.pickle[:-7]
    OFILENAME = outputfilename + network_name + "_degree_BiBC_" + node_type1 + "_" + node_type2 + ".tsv"
    LOGFILENAME = outputfilename + network_name + "_degree_BiBC_" + node_type1 + "_" + node_type2 + "_LOG.txt"
    TMPFILE = "tmp.txt"

    with open(TMPFILE, "w") as f:

        subg = sorted(nx.connected_components(G), key=len, reverse=True) # added by RR
        #subg = sorted(nx.connected_component_subgraphs(G), key = len, reverse = True) # worked in older nx version
        #print(subg)

        ###    Empty cell   ###
        f.write("ID")
        f.write("\t")    

        ###    Node list    ###
        node_names = list(G.nodes)
        print("G.nodes printed here")
        print(G.nodes)
        node_names_sort = sorted(node_names)    
        [f.write(i + "\t") for i in node_names_sort]
        f.write("\n")   
        
        ###   Node degree   ###
        print("Finding each node's degree...")
        node_degrees = ""
        for i in node_names_sort:
            node_degrees = node_degrees + str(G.degree(i)) + "\t"  
        f.write("Node_degrees\t" + node_degrees + "\n")    
        
        ### Bi-BC ###
        # Get the nodes from /just/ the giant component because the rbc function won't work on networks with multiple subgraphs 
        gc_nodes = max(nx.connected_components(G), key=len) # added by RR
        gc = G.subgraph(gc_nodes)
        #print(gc.nodes)
        #print(list(gc.edges))
        #gc = max(nx.connected_component_subgraphs(G), key=len) # worked in older nx version
        #gc_nodes = gc.nodes() # gc was already made earlier in the mean geodesic pathlength function
        
        
        # Make an empty list and string to add the output to. These will be updated right away if there are multiple giant components, 
        # or will be updated at the end if there is only one gc 
        otu_pheno_value_list = []        
        otu_pheno_value_str = ""
        
        # The following try-except statement is intended to be used ONLY if there is an issue in calculating rbc, such as there not being 
        # enough nodes (too restrictive of cutoffs), there only being one modular region of the giant component (if the suer specifies 
        # modularity for bibc parameter), etc. It should (ideally) not be used to catch other errors.
        print("Finding each node's BiBC...")

        # First check if there are multiple connected components and if so then...
        if len(subg) > 1:
            print("Multiple components were found in the graph")

            # Check if the two giant comps are the same size. If so then do not bother calculating rbc
            if len(subg[0]) == len(subg[1]):
                print("There are multiple giant components. BiBC of each node = NA.")

                # Make a list that is the length of the number of nodes in the whole network and write "NA" for the BiBC of each 
                # node, then add those to a string and write the string to the output file
                for i in range(len(node_names_sort)):
                    otu_pheno_value_list.append("NA")
            
                for i in otu_pheno_value_list:
                    otu_pheno_value_str = otu_pheno_value_str + i + "\t"

                f.write("BiBC\t" + otu_pheno_value_str + "\n")
            
            # Otherwise, run restricted_betweenness_centrality function, passing it only the group of nodes the user wants (the arg to "--bibc")        
            elif (len(subg[0]) != len(subg[1])):
                print("BiBC being calculated for giant component.")

                # If the user wants to only calculate bibc on node_types.....
                if bibc_choice == "node_types":
                    # Pass the gc nodes to the function that will assign the correct node types to each of the nodes    
                    otu_pheno_types = assign_node_type(node_input_file, gc_nodes, node_type1, node_type2)
    
                    # Calculate rbc using the above function that Kevin wrote, which takes the outputs of assign_node_type
                    rbc = restricted_betweenness_centrality(gc, otu_pheno_types['Type1'], otu_pheno_types['Type2'], bibc_calc_type)
                    print(rbc)                

                # Otherwise, if they wish to use modularity as the BiBC parameter...
                elif bibc_choice == "modularity":
                    nodes_in_gc_for_bibc_mod = max(nx.connected_component_subgraphs(G), key=len)
                    nodes_from_bibc_mod = bibc_mod(nodes_in_gc_for_bibc_mod)
                    rbc = restricted_betweenness_centrality(gc, nodes_from_bibc_mod['mod1'], nodes_from_bibc_mod['mod2'], bibc_calc_type)
            
                print(rbc)

                # Combine the rbc function output into one single dictionary
                merged_rbc = {**rbc[0], **rbc[1], **rbc[2]}
    
                # Loop through the list of sorted node names and for each one create a new listing in bibc_dict_w_NAs that describes
                # 1) whether the node is present in the network or not and 2) what the BiBC is of that node
                bibc_dict_w_NAs = {}
                for i in node_names_sort:
                    # Check if node is in giant comp by checking to see if it was output from the rbc function (which only accepts the giant comp as input)
                    # If it is present, then just create a key of that node, using the calculated BiBC value
                    if i in merged_rbc:
                        bibc_dict_w_NAs[i] = merged_rbc[i]
                    # If it is not present, then still create a key of the current node, but assign it NA since it was not in the giant comp. 
                    else:
                        bibc_dict_w_NAs[i] = "NA"
    
                # Order the previously made dictionary by key name so it can be input into the properties file
                ordered_bibc = OrderedDict(sorted(bibc_dict_w_NAs.items()))
                    
                for key,value in ordered_bibc.items():
                    otu_pheno_value_list.append(value) # otu_pheno_value_list created prior to checking if there are multiple giant components
            
                # Loop through the list containing just the values, inj order of node names and add each to the otu_pheno_value_str
                for i in otu_pheno_value_list:
                    otu_pheno_value_str = otu_pheno_value_str + str(i) + "\t"
        
                f.write("BiBC\t" + otu_pheno_value_str + "\n")


        # If there is only one giant comp, then still run restricted_betweenness_centrality function. Note that this is just copied and pasted from the above else if
        # statement, since this one only works if there is only one giant component and the other is for checking if there is more than one giant component.        
        else:
            print("There is only one component. BiBC being calculated on the entire graph.")

            # If the user wants to only calculate bibc on node_types.....
            if bibc_choice == "node_types":
                # Pass the gc nodes to the function that will assign the correct node types to each of the nodes    
                otu_pheno_types = assign_node_type(node_input_file, gc_nodes, node_type1, node_type2)
        
                # Calculate rbc using the above function that Kevin wrote, which takes the outputs of assign_node_type
                rbc = restricted_betweenness_centrality(gc, otu_pheno_types['Type1'], otu_pheno_types['Type2'], bibc_calc_type)
            
            # Otherwise, if they wish to use modularity as the BiBC parameter...
            elif bibc_choice == "modularity":
                nodes_in_gc_for_bibc_mod = max(nx.connected_component_subgraphs(G), key=len)
                nodes_from_bibc_mod = bibc_mod(nodes_in_gc_for_bibc_mod)

                #nodes_from_bibc_mod = bibc_mod(gc_nodes)
                rbc = restricted_betweenness_centrality(gc, nodes_from_bibc_mod['mod1'], nodes_from_bibc_mod['mod2'], bibc_calc_type)

            # Combine the rbc function output into one single dictionary
            merged_rbc = {**rbc[0], **rbc[1], **rbc[2]}

            # Loop through the list of sorted node names and for each one create a new listing in bibc_dict_w_NAs that describes
            # 1) whether the node is present in the network or not and 2) what the BiBC is of that node
            bibc_dict_w_NAs = {}
            
            for i in node_names_sort:
                # Check if node is in giant comp by checking to see if it was output from the rbc function (which only accepts the giant comp as input)
                # If it is present, then just create a key of that node, using the calculated BiBC value
                if i in merged_rbc:
                    bibc_dict_w_NAs[i] = merged_rbc[i]
                # If it is not present, then still create a key of the current node, but assign it NA since it was not in the giant comp. 
                else:
                    bibc_dict_w_NAs[i] = "NA"

            # Order the previously made dictionary by key name so it can be input into the properties file
            ordered_bibc = OrderedDict(sorted(bibc_dict_w_NAs.items()))
                
            for key,value in ordered_bibc.items():
                otu_pheno_value_list.append(value) # otu_pheno_value_list created prior to checking if there are multiple giant components
    
            # Loop through the list containing just the values, inj order of node names and add each to the otu_pheno_value_str
            for i in otu_pheno_value_list:
                otu_pheno_value_str = otu_pheno_value_str + str(i) + "\t"
    
            f.write("BiBC\t" + otu_pheno_value_str + "\n")

        # Overwrite any previous file with same name instead of appending    
        f.truncate()
                
    f.close()

print("\nNetwork and node properties have been calculated. Check the pickle_properties.txt file.\n")


# tranpose added by RR (Aug 9 2022)
# https://stackoverflow.com/questions/1729824/an-efficient-way-to-transpose-a-file-in-bash
#for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip())):
#    print(' '.join(c))

# https://stackoverflow.com/questions/4869189/how-to-transpose-a-dataset-in-a-csv-file
import csv
#from itertools import izip # not needed in python 3; zip doesn't like files in bytes mode so removed b
a = zip(*csv.reader(open(TMPFILE, "r"), delimiter="\t"))
csv.writer(open(OFILENAME, "w"), delimiter="\t").writerows(a)
