import probanno
import cobra
import pickle

reaction_probabilities = probanno.generate_reaction_probabilities('/scratch/tjm4k/Lactobacillus/fastas/1002365.5.faa', '/scratch/tjm4k/Lactobacillus/Data/GramPositive.json', genome_id='1002365.5')
pickle.dump(reaction_probabilities, open("/scratch/tjm4k/Lactobacillus/Data/rxn_probs.dic", "wb"))
# universal_model = cobra.io.load_json_model("/scratch/tjm4k/Lactobacillus/Data/GramPosUni.json")
# model = cobra.io.read_sbml_model('/scratch/tjm4k/Lactobacillus/gap_models/1002365.5.xml')

# for i in range(100):
# 	try:
# 		rxn_list = probanno.probabilistic_gapfill(model, universal_model, reaction_probabilities)
# 		pickle.dump(rxn_list, open("/scratch/tjm4k/Lactobacillus/Data/probanno_test.list", "wb"))
# 		print("Passed after %s loops.") % (i)
# 		break
# 	except:
# 		pass