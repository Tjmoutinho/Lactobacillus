# -*- coding: utf-8 -*-

import cobra
import pickle
import sys

with open(sys.argv[1], 'r') as file:
    for line in file:
        genome_id = str(line)

print(genome_id)

# file_name = "../metabolic_output/%s.data" % (genome_id) 
# pickle.dump(total_dataset_dict, open(file_name, "wb"))

# pickle.dump(genome_id, open("output.data", "wb"))
        
model = cobra.io.read_sbml_model('../gap_models/'+ genome_id +'.xml')

