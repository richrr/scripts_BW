# aug 20 2021
import_network_data.py
Takes a csv network file (like the output from your network creation pipeline). You tell it which columns the source and target nodes are in and it will create a networkx-format network, saved as a .pickle file for easy reading into the next script

Example command:
python import_network_data.py --input <network file> --source partner1 --target partner2



BiBC_degree_calculator.py
Takes the .pickle file and calculates BiBC and degree for each node in the network. 
Needs a mapping file without a header to tell it which node belongs to which data type, i.e,:
gene1,gene
gene2,gene
gene3,gene
ASV1,micro
ASV2,micro
ASV5,micro

Example:
python BiBC_degree_calculator.py <pickled network file> --bibc_groups node_types --bibc_calc_type bibc --node_map <node map csv> --node_groups gene micro --log
