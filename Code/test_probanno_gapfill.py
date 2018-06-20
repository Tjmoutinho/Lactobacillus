import probanno

reaction_probabilities = probanno.generate_reaction_probabilities('/scratch/tjm4k/Lactobacillus/fastas/1002365.5.faa', '/scratch/tjm4k/Lactobacillus/Data/GramPositive.json', genome_id='1002365.5')
universal_model = cobra.io.read_json_model("/scratch/tjm4k/Lactobacillus/Data/GramPosUni.json")
model = cobra.io.read_sbml_model('/scratch/tjm4k/Lactobacillus/gap_models/1002365.5.xml')
rxn_list = probabilistic_gapfill(model, universal_model, reaction_probabilities, clean_exchange_rxns=True, default_penalties=None, dm_rxns=False, ex_rxns=False, **solver_parameters)
