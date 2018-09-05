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
    
    ['Vitamin B12r','cpd00423'], # C62H91CoN13O14P : cobalamin;cpd03424;cpd00730 : not present in any exchange reactions
    ['Cobinamide','cpd03422'], #EXs : related to cobalamin (B12) Added to ensure cells have access to B12
    ['BIOT','cpd00104'], # C10H15N2O3S : biotin B7
    ['PAN','cpd00644'], # C9H16NO5 : Pantothenate B5
    ['Folate','cpd00393'], # C19H17N7O6 : B9
    ['Niacin','cpd00218'], # C6H4NO2 : B3
    ['Pyridoxal','cpd00215'], # C8H9NO3 : B6
    ['Riboflavin','cpd00220'], # C17H19N4O6 : B2
    ['thiamin','cpd00305'], # C12H17N4OS : B1
    
#     ['Phosphate','cpd00009'], # HO4P : In M9 Defaults
    
    ['Thioglycolate','cpd01415'], # C2H3O2S : not present in any exchange reactions
#     ['Sulfate','cpd00048'], # O4S : In M9 Defaults
    
    ['Acetate','cpd00029'], # C2H3O2 : not present in any exchange reactions
    ['Citrate','cpd00137'], # C6H5O7 : Consider removing. 
#     ['Polysorbate 60','cpd24450'], # C35H68O10 : Almost tween 80 : not present in any reactions
#     ['Ethyl acetate','cpd00633'], # C4H8O2 : not present in any exchange reactions, only present in one reaction at all
    
    ['ABEE','cpd00443'] # C7H6NO2 : aminobenzoate : not present in any exchange reactions
]

# Potentially add to BSM (from M9 media)
M9_ions = [
    ['Cl-','cpd00099'],
    ['Co2+','cpd00149'],
    ['Cu2+','cpd00058'],
    ['Fe3','cpd10516'],
#     ['Sodium molybdate','cpd11145'], # This doesn't connect to anything
    ['Ni2+','cpd00244'],
    ['Selenate','cpd03396'],
    ['Selenite','cpd03387'],
    ['Zn2+','cpd00034']
]

# M9 default carbon, nitrogen, phosphorous, and sulfur sources
M9_sources = [
#     ['D-Glucose','cpd00027'],
#     ['NH3','cpd00013'], # this is actually NH4 : ammonium
    ['Phosphate','cpd00009'],
    ['Sulfate','cpd00048']
]

# Vitamins
vit_k = [
#     ['BIOT','cpd00104'], #EXs : Biotin
#     ['Cobinamide','cpd03422'], #EXs : related to cobalamin (B12)
#     ['Folate','cpd00393'], #EXs : 
    ['Menaquinone 7','cpd11606'], #EXs : Vitamine K2 : Add when there is no O2
#     ['Niacin','cpd00218'], #EXs : 
#     ['PAN','cpd00644'], #EXs : Pantothenate
#     ['Pyridoxal','cpd00215'], #EXs : 
#     ['Riboflavin','cpd00220'], #EXs : 
#     ['Thiamin','cpd00305'] #EXs : 
]

# For aerobic simulations, O2 was added with a lower bound of −20 and to 0 for anaerobic simulations.

# DNA/RNA related metabolites
rna_bases = [
#     ['35ccmp','cpd00696'], #EXs : 
#     ['AMP','cpd00018'], #EXs : 
    ['Adenosine','cpd00182'], #EXs : In BSM (as adenine)
#     ['Adenosine 3-5-bisphosphate','cpd00045'], #EXs : 
    ['Cytosine','cpd00307'], #EXs : 
#     ['Deoxyadenosine','cpd00438'], #EXs : 
#     ['Deoxycytidine','cpd00654'], #EXs : 
#     ['Deoxyguanosine','cpd00277'], #EXs : In BSM
#     ['Deoxyinosine','cpd03279'], #EXs : 
#     ['Deoxyuridine','cpd00412'], #EXs : 
#     ['GMP','cpd00126'], #EXs : 
#     ['GTP','cpd00038'], #EXs : 
    ['Guanosine','cpd00311'], #EXs : In BSM (as Guanine)
#     ['Inosine','cpd00246'], #EXs : 
#     ['HYXN','cpd00226'], #EXs : Hypoxanthine
#     ['Nicotinamide ribonucleotide','cpd00355'], #EXs : 
#     ['TTP','cpd00357'], #EXs : Deoxythymidine triphosphate
    ['Thymidine','cpd00184'], #EXs : In BSM
#     ['Thyminose','cpd01242'], #EXs : deoxyribose
#     ['Uracil','cpd00092'], #EXs : 
    ['Uridine','cpd00249'], #EXs : In BSM (as uracil)
#     ['XAN','cpd00309'], #EXs : Xanthine
#     ['Xanthosine','cpd01217'], #EXs : 
#     ['dATP','cpd00115'], #EXs : 
#     ['dGTP','cpd00241'], #EXs : 
#     ['dTMP','cpd00298'] #EXs : 
]

# Amino Acid related metabolites
products = [
    ['D-Alanine','cpd00117'], #EXs : 
    ['D-Glutamate','cpd00186'], #EXs : 
    ['D-Methionine','cpd00637'], #EXs : 
    ['D-Serine','cpd00550'], #EXs : 
    ['Glycine','cpd00033'], #EXs : 1
    ['L-Alanine','cpd00035'], #EXs : 2
    ['L-Arginine','cpd00051'], #EXs : 3
    ['L-Asparagine','cpd00132'], #EXs : 4
    ['L-Aspartate','cpd00041'], #EXs : 5
    ['L-Cysteine','cpd00084'], #EXs : 7
    ['L-Glutamate','cpd00023'], #EXs : 8
    ['L-Glutamine','cpd00053'], #EXs : 9
    ['L-Histidine','cpd00119'], #EXs : 10
    ['L-Isoleucine','cpd00322'], #EXs : 11
    ['L-Leucine','cpd00107'], #EXs : 12
    ['L-Lysine','cpd00039'], #EXs : 13
    ['L-Methionine','cpd00060'], #EXs : 14
    ['L-Phenylalanine','cpd00066'], #EXs : 15
    ['L-Proline','cpd00129'], #EXs : 16
    ['L-Serine','cpd00054'], #EXs : 17
    ['L-Threonine','cpd00161'], #EXs : 18
    ['L-Tryptophan','cpd00065'], #EXs : 19
    ['L-Tyrosine','cpd00069'], #EXs : 20
    ['L-Valine','cpd00156'], #EXs : 21
    ['H2O2','cpd00025'],
    ['L-Lactate','cpd00159'],
    ['Acetate','cpd00029'],
    ['Butyrate','cpd00211'],
#     ['isobutyrate','cpd01711'], 
    ['GABA','cpd00281'],
    ['ethanol','cpd00363'],
    ['Propionate','cpd00141'],
    ['formate','cpd00047'],
    ['Valerate','cpd00597'],
#     ['Isovaleric acid','cpd05178'], # Not in universal?
    ['TMAO','cpd00811'], # (CH3)3NO
#     ['Indole-3-(carb)aldehyde','cpd05401'],
#     ['Acetaldehyde','cpd00071'],
    ['Deoxycholate','cpd02733'],
    ['Chorismate','cpd00216'],
    ['Hexanoate','cpd01113']
]

carbon_sources = [
    ['D-Glucose','cpd00027'],
    ['4-Hydroxybenzoate','cpd00136'], #EXs : found in coconuts
    ['2-keto-3-deoxygluconate','cpd00176'], #EXs : degraded pectin product
    ['Amylotriose','cpd01262'], #EXs : 
#     ['CELB','cpd00158'], #EXs : Cellobiose
    ['D-Fructose','cpd00082'], #EXs : 
    ['D-Mannitol','cpd00314'], #EXs : sweetener the is poorly absorbed in the gut
    ['D-Mannose','cpd00138'], #EXs : related to mucin
    ['Ribose','cpd00105'], #EXs : 
    ['Dextrin','cpd11594'], #EXs : 
    ['Dulcose','cpd01171'], #EXs : Galactitol
    ['GLCN','cpd00222'], #EXs : Gluconate 
    ['GLUM','cpd00276'], #EXs : Glucosamine
    ['Galactose','cpd00108'], #EXs : 
    ['L-Arabinose','cpd00224'], #EXs : 
    ['L-Inositol','cpd00121'], #EXs : 
    ['L-Lactate','cpd00159'], #EXs : 
    ['L-Malate','cpd00130'], #EXs : 
    ['Glycerol','cpd00100'], #EXs : 
#     ['LACT','cpd00208'], #EXs : lactose
    ['Maltohexaose','cpd01329'], #EXs : 
#     ['Maltose','cpd00179'], #EXs : 
    ['Melibiose','cpd03198'], #EXs : 
    ['Palmitate','cpd00214'], #EXs : 
    ['Propionate','cpd00141'], #EXs : 
#     ['Salicin','cpd01030'], #EXs : 
    ['Sorbitol','cpd00588'], #EXs : 
    ['Stachyose','cpd01133'], #EXs : 
    ['Succinate','cpd00036'], #EXs : 
#     ['Sucrose','cpd00076'], #EXs : 
#     ['TRHL','cpd00794'], #EXs : Trehalose
    ['Ursin','cpd03696'], #EXs : Arbutin
    ['Xylose','cpd00154'], #EXs : 
    ['hexadecenoate','cpd15237'] #EXs : 
]

# Nitrogen Sources
nitrogen_sources = [
    ['NH3','cpd00013'], #EXs : 
    ['Allantoin','cpd01092'], #EXs : degradation product of purines
    ['BET','cpd00540'], #EXs : Betaine
    ['Choline','cpd00098'], #EXs : Found in milk
    ['Nitrate','cpd00209'], #EXs : 
    ['Nitrite','cpd00075'], #EXs : 
    ['Spermidine','cpd00264'], #EXs : 
    ['Urea','cpd00073'], #EXs : 
    ['crotonobetaine','cpd08305'] #EXs : 
]

# pFBA gapfiller
def pfba_gapfill(model, reaction_bag, likelihoods, obj=None, obj_lb=10., obj_constraint=False,
                 iters=1, tasks=None, task_lb=0.05, 
                 add_exchanges=True, extracellular='e'):

    start_time = time.time()
    # Save some basic network info for downstream membership testing
    orig_rxn_ids = set([str(x.id) for x in model.reactions])
    orig_cpd_ids = set([str(y.id) for y in model.metabolites])

    # Get model objective reaction ID
    if obj == None:
        obj = get_objective(model)
    else:
        obj = obj
    
    # Modify universal reaction bag
    new_rxn_ids = set()
    with reaction_bag as universal:

        # Remove overlapping reactions from universal bag, and reset objective if needed
        orig_rxns = list(copy.deepcopy(model.reactions))
            
        # Add pFBA to universal model and add model reactions
        
        add_pfba_likely(universal, likelihoods)
        
        universal.add_reactions(orig_rxns)
        
        # If previous objective not set as constraint, set minimum lower bound
        if obj_constraint == False: 
            universal.reactions.get_by_id(obj).lower_bound = obj_lb
        
        # Run FBA and save solution
        solution = universal.optimize()
#         print([bound.id for bound in universal.boundary if bound.lower_bound != 0.0 and bound.upper_bound != 0.0 and bound.id.startswith('EX')])
        
        # Identify which reactions carry flux in solution
        rxns = list(solution.fluxes.index)
        fluxes = list(solution.fluxes)
        for flux in range(0, len(fluxes)):
            if abs(fluxes[flux]) > 1e-6:
                new_rxn_ids |= set([rxns[flux]])
        
    # Screen new reaction IDs
    if obj in new_rxn_ids: new_rxn_ids.remove(obj)
    for rxn in orig_rxn_ids:
        try:
            new_rxn_ids.remove(rxn)
        except:
            continue
    
    # Get reactions and metabolites to be added to the model
    new_rxns = copy.deepcopy([reaction_bag.reactions.get_by_id(rxn) for rxn in new_rxn_ids])
    new_cpd_ids = set()
    for rxn in new_rxns: new_cpd_ids |= set([str(x.id) for x in list(rxn.metabolites)])
    new_cpd_ids = new_cpd_ids.difference(orig_cpd_ids)
    new_cpds = copy.deepcopy([reaction_bag.metabolites.get_by_id(cpd) for cpd in new_cpd_ids])
    
    # Copy model and gapfill
    new_model = copy.deepcopy(model)
    new_model.add_metabolites(new_cpds)
    new_model.add_reactions(new_rxns)
    
    duration = int(round(time.time() - start_time))
    print('Took ' + str(duration) + ' seconds to gapfill ' + str(len(new_rxn_ids)) + \
          ' reactions and ' + str(len(new_cpd_ids)) + ' metabolites.') 
    
    return {'NewModel':new_model, 'gaps':new_rxn_ids, 'mets':new_cpd_ids}

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

# Using pFBA with new media components (RNA bases + thymidine...) 
# Remove reaction likelihoods of zero from model
# Add demands for all metabolites in model to avoid any reactions being blocked
# Weighted pFBA gapfill solution

t = time.time()

sys.stdout.write('Loading in models...')

universal = cobra.io.load_json_model("../Data/GramPosUni.json")
genome_id = '220668.9'
model = cobra.io.read_sbml_model('../gap_models/'+ genome_id +'.xml')
likelihoods = pickle.load(open('../likelihoods/'+ genome_id +'.probs'))

sys.stdout.write('Adding Water...')

# Ensure free diffusion of water
model.reactions.get_by_id('rxn05319_c').name = "Water transport"
model.reactions.get_by_id('rxn05319_c').bounds = (-1000., 1000.)

sys.stdout.write('Set-up Universal...')

# Add demand for all metabolites in Universal model to stop blocked reactions
all_mets = []
for met in universal.metabolites:
    if (met.id.endswith('_c')):
        universal.add_boundary(met, type='demand')
# universal.solver = 'gurobi'

### Set Up Model: remove low likelihood reactions
sys.stdout.write('Set-up Model...')
low_like_model = []
for rxn in model.reactions:
    if rxn.id.startswith('rxn'):
        try:
            if likelihoods[rxn.id] <= 0.1:
                low_like_model.append(rxn.id)
        except:
            pass
model_rxns_to_remove = [model.reactions.get_by_id(rxn) for rxn in low_like_model]
model.remove_reactions(model_rxns_to_remove)
# model.solver = 'gurobi'

# Remove model reactions from universal
orig_rxn_ids = set([str(x.id) for x in model.reactions])
univ_rxn_ids = set([str(z.id) for z in universal.reactions])
overlap_rxn_ids = univ_rxn_ids.intersection(orig_rxn_ids)
sys.stdout.write('removing reactions from universal...')
for rxn in overlap_rxn_ids: 
    universal.reactions.get_by_id(rxn).remove_from_model()

print(str(round(time.time() - t)) + 'seconds to complete')

# Run through each amino acid to check for production
global_time = time.time()
counter = 0

total_dataset_dict = {}
carb_idx = 0

for carbon in carbon_sources[0:1]:
    nit_idx = 0
    for nitrogen in nitrogen_sources[0:1]:
        # Create and set specific Media List
        media_list = bsm + M9_sources + rna_bases + [nitrogen] + [carbon]
        set_media(model, media_list, universal, verbose=False)
        
        product_idx = 0
        for product_list in products:
            t = time.time()
            sys.stdout.write('\n'+ 'Loop' + str(counter) + ' ')
            product = product_list[1]+'_c'
            product_name = product_list[0]
            
            with model as temp_model, universal as temp_universal:
                
                try:
                    metabolite = temp_model.metabolites.get_by_id(product)
                    demand = temp_model.add_boundary(metabolite, type='demand')
                    temp_model.objective = demand
                except:
                    if set(product).issubset(set([met.id for met in model.metabolites])) == 0:
                        temp_model.add_metabolites(copy.deepcopy(temp_universal.metabolites.get_by_id(product)))
                        metabolite = temp_model.metabolites.get_by_id(product)
                        demand = temp_model.add_boundary(metabolite, type='demand')
                        temp_model.objective = demand
                
                temp_universal.reactions.get_by_id(demand.id).remove_from_model()
                
                sys.stdout.write('Gapfilling ' + product + '...')
                dont_continue = 0

                try:
                    new_gapfill_data = pfba_gapfill(temp_model, temp_universal, likelihoods, obj=demand.id, obj_lb=10., obj_constraint=False, iters=1, add_exchanges=False)
                except:
                    try:
                        sys.stdout.write('Restart Gapfilling...')
                        new_gapfill_data = pfba_gapfill(temp_model, temp_universal, likelihoods, obj=demand.id, obj_lb=10., obj_constraint=False, iters=1, add_exchanges=False)
                    except:
                        sys.stdout.write('Failed Gapfilling...')
                        dont_continue = 1
                        pass
                
                gaps_to_fill = new_gapfill_data['gaps']
                new_model = new_gapfill_data['NewModel']
                mets_added = new_gapfill_data['mets']
                
                if dont_continue == 0:
                    # Optimize with filled pathway
                    sys.stdout.write('pFBA...')
                    solution = pfba(new_model, objective = demand)
                    sys.stdout.write(str(round(new_model.slim_optimize())) + '...')
                    
                    if new_model.slim_optimize() <= 1.0:
                        dont_continue = 1 
                    
                    if dont_continue == 0:
                    
                        # Find reactions that carry flux
                        df = solution.fluxes.to_frame()
                        active = df.loc[(abs(df['fluxes'])) > 0.1]

                        demand_list = []
                        for rxn_id in active.index:
                            if rxn_id.startswith('DM') and rxn_id != demand.id:
                                demand_list.append(rxn_id)
                                print('demands added')

                        # Acquire likelihood scores for reactions that carry flux
                        flux_rxns = []
                        like_list = []
                        gap_flux_rxns = []
                        gap_like_list = []
                        path_flux_rxns = []
                        path_like_list = []
                        for rxn in list(active.index):
                            if rxn in gaps_to_fill and rxn.startswith('rxn'):
                                try:
                                    gap_flux_rxns.append([str(rxn),likelihoods[str(rxn)]])
                                    gap_like_list.append(likelihoods[str(rxn)])
                                except:
                                    pass
                            if rxn not in gaps_to_fill and rxn.startswith('rxn'):
                                try:
                                    path_flux_rxns.append([str(rxn),likelihoods[str(rxn)]])
                                    path_like_list.append(likelihoods[str(rxn)])
                                except:
                                    pass
                            if rxn.startswith('rxn'):
                                try:
                                    flux_rxns.append([str(rxn),likelihoods[str(rxn)]])
                                    like_list.append(likelihoods[str(rxn)])
                                except:
                                    pass
                        avg_like = np.mean(like_list)

                        opt_before = model.slim_optimize()
                        if opt_before > 0.0:
                            gap_avg_like = []
                        else:
                            gap_avg_like = np.mean(gap_like_list)
                        path_avg_like = np.mean(path_like_list)
                        sys.stdout.write('Ave likelihood of: ' + product + ' is ' + str(avg_like) + '...')

                        counter += 1

                        report_dict = {}
                        report_dict['Model_ID'] = genome_id
                        report_dict['Carbon'] = carbon
                        report_dict['Nitrogen'] = nitrogen
                        report_dict['objective'] = product_name
                        report_dict['opt_before'] = opt_before
                        report_dict['opt_after'] = new_model.slim_optimize()
                        report_dict['avg_path_like'] = avg_like
                        report_dict['gap_avg_like'] = gap_avg_like
                        report_dict['path_avg_like'] = path_avg_like
                        report_dict['gaps_filled'] = gaps_to_fill
                        report_dict['mets_added'] = mets_added
                        report_dict['reactions_w_flux'] = flux_rxns
                        report_dict['gaps_w_flux'] = gap_flux_rxns
                        report_dict['path_w_flux'] = path_flux_rxns
                        report_dict['active_rxns'] = active
                        report_dict['demands'] = demand_list

                        report_dict_ID = genome_id + ':' + str(carb_idx) + '.' + str(nit_idx) + '.' + str(product_idx)
                        total_dataset_dict[report_dict_ID] = report_dict
                        product_idx += 1 # Keep track to which product is being maximized
                    
#                     elapsed = time.time() - t
#                     sys.stdout.write('Run time: ' + str(elapsed/60) + " [mins]...")
                    
                if dont_continue == 1:
                    counter += 1
                    
                    report_dict = {}
                    report_dict['Model_ID'] = genome_id
                    report_dict['Carbon'] = carbon
                    report_dict['Nitrogen'] = nitrogen
                    report_dict['objective'] = product_name
                    report_dict['opt_before'] = "Failed to gapfill"
                    report_dict['opt_after'] = "Failed to gapfill"
                    report_dict['avg_path_like'] = "Failed to gapfill"
                    report_dict['gap_avg_like'] = "Failed to gapfill"
                    report_dict['path_avg_like'] = "Failed to gapfill"
                    report_dict['gaps_filled'] = "Failed to gapfill"
                    report_dict['mets_added'] = "Failed to gapfill"
                    report_dict['reactions_w_flux'] = "Failed to gapfill"
                    report_dict['gaps_w_flux'] = "Failed to gapfill"
                    report_dict['path_w_flux'] = "Failed to gapfill"
                    report_dict['active_rxns'] = "Failed to gapfill"
                    report_dict['demands'] = "Failed to gapfill"
                    
                    report_dict_ID = genome_id + ':' + str(carb_idx) + '.' + str(nit_idx) + '.' + str(product_idx)
                    total_dataset_dict[report_dict_ID] = report_dict
                    product_idx += 1 #Keep track to which product is being maximized
                    
#                     elapsed = time.time() - t
#                     sys.stdout.write('Run time: ' + str(elapsed/60) + " [mins]...")
        nit_idx += 1
    carb_idx += 1
file_name = "../metabolic_output/%s.data" % (genome_id) 
pickle.dump(total_dataset_dict, open(file_name, "wb"))

print('\n'+ str((time.time() - global_time)/60) + 'mins to complete!')
