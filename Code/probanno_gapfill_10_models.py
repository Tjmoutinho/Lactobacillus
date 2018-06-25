import probanno
import cobra
import pickle
import copy

model_ids = ['1002365.5','220668.9','712961.3','33959.180','568703.30','1597.20','1596.54','557433.4','1215914.3','1587.42']

for model_id in model_ids:
	# Open all useful files
	fasta_path = "/scratch/tjm4k/Lactobacillus/fastas/%s.faa" % (model_id)
	model_path = "/scratch/tjm4k/Lactobacillus/gap_models/%s.xml" % (model_id)
	reaction_probabilities = probanno.generate_reaction_probabilities(fasta_path, '/scratch/tjm4k/Lactobacillus/Data/GramPositive.json', genome_id=model_id)
	# pickle.dump(reaction_probabilities, open("/scratch/tjm4k/Lactobacillus/Data/rxn_probs.dic", "wb"))
	universal = cobra.io.load_json_model("/scratch/tjm4k/Lactobacillus/Data/GramPosUni.json")
	model = cobra.io.read_sbml_model(model_path)
	model.solver = 'gurobi'

	# Create Universal Model with only reactions that have genetic evidence
	# rxns_w_probs = []
	# for rxn in universal.reactions:
	#     try:
	#         if rxn_probs[rxn.id] > 0:
	#             rxns_w_probs.append(rxn.id)
	#         elif rxn_probs[rxn.id] == 0:
	#             pass
	#     except:
	#         pass
	# Create list of reactions without probs from the universal model
	# uni_rxns = set([rxn.id for rxn in universal.reactions])
	# rxns_to_remove = list(uni_rxns - set(rxns_w_probs))
	# Remove Reactions from Universal Model
	# uni_slim = copy.deepcopy(universal)
	# uni_slim.remove_reactions(rxns_to_remove)

	# Try to gapfill 10 times for each model
	output_path = "/scratch/tjm4k/Lactobacillus/Data/%s_solutions.list" % model_id
	for i in range(10):
		try:
			rxn_list = probanno.probabilistic_gapfill(model, universal, reaction_probabilities)
			# rxn_list = probanno.probabilistic_gapfill(model, uni_slim, reaction_probabilities)
			pickle.dump(rxn_list, open(output_path, "wb"))
			print("%(model_id)s passed after %(i)s loops.") % {"model_id" : model_id, "i" : i}
			break
		except:
			pass
