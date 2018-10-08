# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import absolute_import

import cobra
import cobra.test
import numpy as np
import csv
import glob
import pickle
import pandas as pd
import math
import copy
import time
import random
import time
import sys

from copy import deepcopy
from collections import defaultdict
from cobra.flux_analysis import sample
from cobra.core.solution import get_solution
from cobra.flux_analysis.sampling import OptGPSampler
from cobra.manipulation.delete import *
from cobra.medium import find_boundary_types
from cobra.flux_analysis import pfba

from warnings import warn
from itertools import chain
from optlang.symbolics import Zero
from cobra.util import solver as sutil
from cobra.core.solution import get_solution

import logging
LOGGER = logging.getLogger(__name__)

def set_media(model, media, universal, verbose=False):

    # Find and close all exchange reactions in the model
    model_rxns = [rxn.id for rxn in model.reactions]
    for rxn in model_rxns:
        if rxn.startswith('EX_') and rxn.endswith('_e'):
            model.reactions.get_by_id(rxn).lower_bound = 0.0

    # Check for existence of exchange reactions for the media metabolites in the model
    for metabolite in media:
        met = metabolite[1]+'_e'
        if 'EX_'+met in model_rxns:
            model.reactions.get_by_id('EX_'+met).lower_bound = -1000.
        else:
            # Create exchange reaction and add to model
            if verbose:
                print("added exchange rxn for " + met)
            new_exchange = cobra.Reaction('EX_'+met)
            new_exchange.name = met + ' exchange'
            met_obj = universal.metabolites.get_by_id(met)
            new_exchange.add_metabolites({met_obj:-1})
            new_exchange.lower_bound = -1000.
            new_exchange.upper_bound = 1000.
            model.add_reaction(new_exchange)
            model.repair()

# Basal Synthetic Media
bsm = [
    ['H+','cpd00067'],
    ['H2O','cpd00001'],
    ['CO2','cpd00011'],
    ['O2','cpd00007'],
    ['N2','cpd00528'], 
#     ['H2','cpd11640'], # Only with no O2
    
    ['K+','cpd00205'],
    ['Na+','cpd00971'],
    ['Mg','cpd00254'],
    ['Mn2+','cpd00030'],
    ['Fe2+','cpd10515'], # Iron ion in heme
    ['Ca2+','cpd00063'], # Calcium pantothenate;cpd19112
    
#     ['Vitamin B12r','cpd00423'], # C62H91CoN13O14P : cobalamin;cpd03424;cpd00730 : not present in any exchange reactions
    ['Cobinamide','cpd03422'], #EXs : related to cobalamin (B12) Added to ensure cells have access to B12
#     ['BIOT','cpd00104'], # C10H15N2O3S : biotin B7
    ['PAN','cpd00644'], # C9H16NO5 : Pantothenate B5
#     ['Folate','cpd00393'], # C19H17N7O6 : B9
#     ['Niacin','cpd00218'], # C6H4NO2 : B3
#     ['Pyridoxal','cpd00215'], # C8H9NO3 : B6
#     ['Riboflavin','cpd00220'], # C17H19N4O6 : B2
#     ['thiamin','cpd00305'], # C12H17N4OS : B1
    
    ['Thioglycolate','cpd01415'], # C2H3O2S : not present in any exchange reactions
    ['Acetate','cpd00029'], # C2H3O2 : not present in any exchange reactions
    ['Citrate','cpd00137'], # C6H5O7 : Consider removing. 
#     ['Polysorbate 60','cpd24450'], # C35H68O10 : Almost tween 80 : not present in any reactions
#     ['Ethyl acetate','cpd00633'], # C4H8O2 : not present in any exchange reactions, only present in one reaction at all
    
    ['ABEE','cpd00443'] # C7H6NO2 : aminobenzoate : not present in any exchange reactions
]

# M9 default carbon, nitrogen, phosphorous, and sulfur sources
M9_sources = [
#     ['D-Glucose','cpd00027'],
    ['NH3','cpd00013'], # this is actually NH4 : ammonium
    ['Phosphate','cpd00009'],
    ['Sulfate','cpd00048']
]

# DNA/RNA related metabolites
rna_bases = [
    ['Adenosine','cpd00182'], #EXs : In BSM (as adenine)
    ['Cytosine','cpd00307'], #EXs : 
    ['Guanosine','cpd00311'], #EXs : In BSM (as Guanine)
    ['Thymidine','cpd00184'], #EXs : In BSM
    ['Uridine','cpd00249'], #EXs : In BSM (as uracil)
]

# Check to see if these metabolites are used in pathways? Should I add some of these to media? 
# Yes for ATP, and GTP. (TTP, CTP as well?)

# Amino Acid related metabolites
products = [
    ['D_Alanine','cpd00117'], #EXs : 
    ['D_Glutamate','cpd00186'], #EXs : 
    ['D_Methionine','cpd00637'], #EXs : 
    ['D_Serine','cpd00550'], #EXs : 
    ['Glycine','cpd00033'], #EXs : 1
    ['L_Alanine','cpd00035'], #EXs : 2
    ['L_Arginine','cpd00051'], #EXs : 3
    ['L_Asparagine','cpd00132'], #EXs : 4
    ['L_Aspartate','cpd00041'], #EXs : 5
    ['L_Cysteine','cpd00084'], #EXs : 7
    ['L_Glutamate','cpd00023'], #EXs : 8
    ['L_Glutamine','cpd00053'], #EXs : 9
    ['L_Histidine','cpd00119'], #EXs : 10
    ['L_Isoleucine','cpd00322'], #EXs : 11
    ['L_Leucine','cpd00107'], #EXs : 12
    ['L_Lysine','cpd00039'], #EXs : 13
    ['L_Methionine','cpd00060'], #EXs : 14
    ['L_Phenylalanine','cpd00066'], #EXs : 15
    ['L_Proline','cpd00129'], #EXs : 16
    ['L_Serine','cpd00054'], #EXs : 17
    ['L_Threonine','cpd00161'], #EXs : 18
    ['L_Tryptophan','cpd00065'], #EXs : 19
    ['L_Tyrosine','cpd00069'], #EXs : 20
    ['L_Valine','cpd00156'], #EXs : 21
    
    ['Choline','cpd00098'],
    ['Citrulline','cpd00274'],
#     ['Cytosine','cpd00307'],
    ['Hypoxanthine','cpd00226'],
    ['Inosine','cpd00246'],
    ['Tyramine','cpd00374'],
    
    ['D_Lactate','cpd00221'],
    ['L_Lactate','cpd00159'],
#     ['Acetate','cpd00029'],
    ['Butyrate','cpd00211'],
    ['Formate','cpd00047'],
    ['Fumarate','cpd00106'],
    ['Propionate','cpd00141'],
#     ['isobutyrate','cpd01711'], 
#     ['Valerate','cpd00597'],
#     ['Isovaleric acid','cpd05178'], # Not in universal?
    
    ['Chorismate','cpd00216'],
#     ['Deoxycholate','cpd02733'],
#     ['Hexanoate','cpd01113'],
    ['Succinate','cpd00036'],
    ['Urocanate','cpd00581'],
    
    ['Adenine','cpd00128'],
    ['AMP','cpd00018'],
    ['UMP','cpd00091'],
    ['Uracil','cpd00092'],
    
    ['B1_Thiamin','cpd00305'], # C12H17N4OS : B1
    ['B2_Riboflavin','cpd00220'], # C17H19N4O6 : B2
    ['B3_Niacin','cpd00218'], # C6H4NO2 : B3
    ['B6_Pyridoxal','cpd00215'], # C8H9NO3 : B6
    ['B7_BIOT','cpd00104'], # C10H15N2O3S : biotin B7
    ['B9_Folate','cpd00393'], # C19H17N7O6 : B9
#     ['B12r','cpd00423'], # C62H91CoN13O14P : cobalamin;cpd03424;cpd00730 : not present in any exchange reactions
    
    ['Ethanol','cpd00363'],
    ['GABA','cpd00281'],
    ['H2O2','cpd00025'],
    ['TMAO','cpd00811'], # (CH3)3NO
    ['3_Hydroxypropanal','cpd00714']
]

carbon_sources = [
    ['D_Glucose','cpd00027'],
    ['Galactose','cpd00108'], #EXs : 
    ['D_Mannose','cpd00138'], #EXs : related to mucin
#     ['D_Hexose','cpd00548']
]

def add_pfba_likely(model, likelihoods, objective=None, fraction_of_optimum=1.0):
    if objective is not None:
        model.objective = objective
    if model.solver.objective.name == '_pfba_objective':
        raise ValueError('The model already has a pFBA objective.')
    sutil.fix_objective_as_constraint(model, fraction=fraction_of_optimum)
    reaction_variables = ((rxn.forward_variable, rxn.reverse_variable)
                          for rxn in model.reactions)
    variables = chain(*reaction_variables)
    dict1 = {}
    fail_report = []
    model_reactions = [rxn.id.split('_')[0] for rxn in model.reactions if rxn.id.startswith('rxn')]
    for v in variables:
        if set([str(v.name.split('_')[0])]).issubset(set(model_reactions)) and str(v.name.split('_')[0]).startswith('rxn'):
            rxn_id = (v.name.split('_')[0] + '_c')
            try:
                dict1[v] = max([0.0, 1.0 - likelihoods[rxn_id]])
            except:
                try:
                    dict1[v] = 1.0
                except:
                    print('FAILED')
                    pass
                pass
            
        elif str(v.name.split('_')[0]).startswith('DM'):
            dict1[v] = 1.0
        else:
            fail_report.append(1)
    model.objective = model.problem.Objective(Zero, direction='min', sloppy=True, name="_pfba_objective")
    model.objective.set_linear_coefficients(dict1)
    
# Actually prune all unused metabolites and reactions (innate function does not work)
def removeUnused(model):
    removed_cpd = set()
    removed_rxn = set()
    unused_current_cpd = 1
    unused_current_rxn = 1
    
    while unused_current_cpd != 0 or unused_current_rxn != 0:
        unused_cpd = prune_unused_metabolites(model)
        removed_cpd |= set(unused_cpd)
        unused_rxn = prune_unused_reactions(model)
        removed_rxn |= set(unused_rxn)
        
        unused_current_cpd = len(unused_cpd)
        unused_current_rxn = len(unused_rxn)
    
    print('Pruned ' + str(len(removed_cpd)) + ' metabolites from model')
    print('Pruned ' + str(len(removed_rxn)) + ' reactions from model')
        
    return(list(removed_cpd), list(removed_rxn))

# Load in models

t = time.time()

sys.stdout.write('Loading in models...')

import sys
with open(sys.argv[1], 'r') as file:
    for line in file:
        genome_id = str(line)

universal = cobra.io.load_json_model("../Data/GramPosUni.json")
universal_orig = cobra.io.load_json_model("../Data/GramPosUni.json")

# genome_id = '220668.9'
# model = cobra.io.read_sbml_model('../gap_models/'+ genome_id +'.xml')
likelihoods = pickle.load(open('../likelihoods/'+ genome_id +'.probs'))

sys.stdout.write('Adding Water...')

# Ensure free diffusion of water
universal.reactions.get_by_id('rxn05319_c').name = "Water transport"
universal.reactions.get_by_id('rxn05319_c').bounds = (-1000., 1000.)

# print(len(universal.reactions))

sys.stdout.write('Remove no-like rxns...')
rxns_to_remove = []
for rxn in universal.reactions:
    try:
        likelihoods[str(rxn.id)]
    except:
        rxns_to_remove.append(rxn)
        pass

for rxn in rxns_to_remove:
    universal.reactions.get_by_id(rxn.id).remove_from_model(remove_orphans = True)
unused_c, unused_r = removeUnused(universal)

sys.stdout.write('Set-up Universal...')

# Add demand for all metabolites in Universal model to stop blocked reactions
all_mets = []
for met in universal.metabolites:
    if (met.id.endswith('_c')):
        universal.add_boundary(met, type='demand')

print(str(round(time.time() - t)) + 'seconds to complete')

# Run through each amino acid to check for production
global_time = time.time()
counter = 0

total_dataset_dict = {}
carb_idx = 0

for carbon in carbon_sources:
    # Create and set specific Media List
    media_list = bsm + M9_sources + rna_bases + [carbon]
    set_media(universal, media_list, universal_orig, verbose=False) 

    product_idx = 0
    for product_list in products: # ADJUST THIS RANGE
        t = time.time()
        sys.stdout.write('\n'+ 'Loop' + str(counter) + ' ')
        product = product_list[1]+'_c'
        product_name = product_list[0]

        with universal as temp_universal:
            # Optimize with filled pathway
            sys.stdout.write('pFBA...')
#             solution = pfba(new_model, objective = demand)
            demand = universal.reactions.get_by_id('DM_'+product_list[1]+'_c')
            temp_universal.objective = demand
            add_pfba_likely(temp_universal, likelihoods, objective=demand)
            solution = temp_universal.optimize()
            sys.stdout.write(str(round(temp_universal.slim_optimize())) + '...')

            # Find reactions that carry flux
            df = solution.fluxes.to_frame()
            active = df.loc[(abs(df['fluxes'])) > 0.1]

            # Acquire likelihood scores for reactions that carry flux
            flux_rxns = []
            like_list = []
            for rxn in list(active.index):
                if rxn.startswith('rxn'):
                    try:
                        flux_rxns.append([str(rxn),likelihoods[str(rxn)]])
                        like_list.append(likelihoods[str(rxn)])
                    except:
                        pass
            avg_like = np.mean(like_list)

            sys.stdout.write('Ave likelihood of: ' + product + ' is ' + str(avg_like) + '...')

            counter += 1

            report_dict = {}
            report_dict['Model_ID'] = genome_id
            report_dict['Carbon'] = carbon
            report_dict['objective'] = product_name
            report_dict['avg_path_like'] = avg_like
            report_dict['reactions_w_flux'] = flux_rxns
            report_dict['active_rxns'] = active

            report_dict_ID = genome_id + ':' + str(carb_idx) + '.' + str(product_idx)
            total_dataset_dict[report_dict_ID] = report_dict
            product_idx += 1 # Keep track to which product is being maximized

    carb_idx += 1
file_name = "../metabolic_output_V3/%s.data" % (genome_id) # CHANGE BACK TO CORRECT NAME
pickle.dump(total_dataset_dict, open(file_name, "wb"))

print('\n'+ str((time.time() - global_time)/60) + 'mins to complete!')


