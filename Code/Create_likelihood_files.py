import probanno
import cobra
import pickle
import copy
import json
import glob

# Open all useful files
# genome_ids_list
model_paths = glob.glob('/scratch/tjm4k/Lactobacillus/models/*.xml')
genome_ids = [x.replace("/scratch/tjm4k/Lactobacillus/models/","").replace(".xml","") for x in model_paths]
# Existing Prob files
existing_files = glob.glob('../likelihoods/*.probs')
existing_files = [x.replace("../likelihoods/","").replace(".probs","") for x in existing_files]

for genome_id in genome_ids:
	if not genome_id in existing_files:
		reaction_probabilities = probanno.generate_reaction_probabilities('/scratch/tjm4k/Lactobacillus/fastas/'+ genome_id +'.faa', '/scratch/tjm4k/Lactobacillus/Data/GramPositive.json', genome_id = genome_id)
		pickle.dump(reaction_probabilities, open('/scratch/tjm4k/Lactobacillus/likelihoods/'+ genome_id +'.probs', "wb"))
