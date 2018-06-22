import probanno
import cobra
import pickle
import copy

# Open all useful files
reaction_probabilities = probanno.generate_reaction_probabilities('/scratch/tjm4k/Lactobacillus/fastas/1002365.5.faa', '/scratch/tjm4k/Lactobacillus/Data/GramPositive.json', genome_id='1002365.5')
pickle.dump(reaction_probabilities, open("/scratch/tjm4k/Lactobacillus/Data/rxn_probs.dic", "wb"))
universal = cobra.io.load_json_model("/scratch/tjm4k/Lactobacillus/Data/GramPosUni.json")
model = cobra.io.read_sbml_model('/scratch/tjm4k/Lactobacillus/gap_models/1002365.5.xml')
model.solver = 'gurobi'

# Create Universal Model with only reactions that have genetic evidence
rxns_w_probs = []
for rxn in universal.reactions:
    try:
        if rxn_probs[rxn.id] > 0:
            rxns_w_probs.append(rxn.id)
        elif rxn_probs[rxn.id] == 0:
            pass
    except:
        pass
# Create list of reactions without probs from the universal model
uni_rxns = set([rxn.id for rxn in universal.reactions])
rxns_to_remove = list(uni_rxns - set(rxns_w_probs))
# Remove Reactions from Universal Model
uni_slim = copy.deepcopy(universal)
uni_slim.remove_reactions(rxns_to_remove)
uni_slim

# Try to gapfill 100 times
for i in range(100):
	try:
		rxn_list = probanno.probabilistic_gapfill(model, uni_slim, reaction_probabilities)
		pickle.dump(rxn_list, open("/scratch/tjm4k/Lactobacillus/Data/probanno_test.list", "wb"))
		print("Passed after %s loops.") % (i)
		break
	except:
		pass
