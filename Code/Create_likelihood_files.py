import probanno
import cobra
import pickle
import copy
import json
import glob

# Open all useful files
# model_paths = glob.glob('/scratch/tjm4k/Lactobacillus/models/*.xml')
# genome_ids = [x.replace("/scratch/tjm4k/Lactobacillus/models/","").replace(".xml","") for x in model_paths]

# Read in additional bifido genomes
# count = 0
# with open('../Data/Bifido_genomes.csv') as csvfile:
#     genome_ids = []
#     for line in csvfile:
#         if count == 0:
#             count = 1
#             continue
#         genome_ids.append(line.strip().split(',')[0])

# Existing Prob files
existing_files = glob.glob('../likelihoods_V4/*.probs')
existing_files = [x.replace("../likelihoods_V4/","").replace(".probs","") for x in existing_files]

genome_ids = pickle.load(open("../Data/Filtered_lacto_genome_ids.data", "rb"))

for genome_id in genome_ids[0:2]:
	if not genome_id in existing_files:
		try:
# 			reaction_probabilities = probanno.generate_reaction_probabilities('/scratch/tjm4k/Lactobacillus/fastas/'+ genome_id +'.faa', '/scratch/tjm4k/Lactobacillus/Data/GramPositive.json', genome_id = genome_id)
			reaction_probabilities = probanno.generate_reaction_probabilities('/scratch/tjm4k/Lactobacillus/fastas/'+ genome_id +'.faa', '/scratch/tjm4k/Lactobacillus/Data/universal_w_OK_rxns.json', genome_id = genome_id)
			pickle.dump(reaction_probabilities, open('/scratch/tjm4k/Lactobacillus/likelihoods_V4/'+ genome_id +'.probs', "wb"))
		except:
			pass
