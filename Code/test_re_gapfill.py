import cobra
print cobra.__version__
import cobra.flux_analysis
import glob
import time
import pickle
import copy

# Identify gapfilled reactions
def findGapFilled(model):
    gapfilled = []
    for index in model.reactions:
        if not index.id in ['bio1']:
            if len(list(index.genes)) == 0:
                if not index in model.boundary:
                    gapfilled.append(index.id)

#     if len(gapfilled) > 0:
#         print(str(len(gapfilled)) + ' reactions not associated with genes')
    
    return gapfilled

# Remove gapfilled rxns from the PATRIC model
def pruneGaps(model, orphans, gap_object):
    model_temp = copy.deepcopy(model)
    gaps_total = set()
    for index in gap_object:
        gaps = [x for x in index['reactions']]
        gaps_total |= set(gaps)

    model_gaps = []
    for x in gaps_total:
        try:
            model_gaps.append(model_temp.reactions.get_by_id(x).id)
        except:
            pass

    remove = [x for x in model_gaps if x in orphans]    

    print('Reactions removed: ' + str(len(remove)))
    model_temp.remove_reactions(remove)

    return(model_temp)

universal = cobra.io.load_json_model("../Data/GramPosUni.json")

# Test Gapfill using just reactions removed
genome_id = '220668.9'

file_name = "../models/%s.xml" % (genome_id)
model = cobra.io.read_sbml_model(file_name)
orphans = findGapFilled(model)

pickle_path = "../gapfilled/%s.gf" % (genome_id)
gap_object = pickle.load(open(pickle_path, "rb"))
pruned_model = pruneGaps(model, orphans, gap_object)

loaded_model_path = "../gap_models/%s.xml" % (genome_id)
loaded_model = cobra.io.read_sbml_model(loaded_model_path)

solution = model.optimize()
print(solution)

solution = pruned_model.optimize()
print(solution)

solution = loaded_model.optimize()
print(solution)

gaps_total = set()
for index in gap_object:
    gaps = [x for x in index['reactions']]
    gaps_total |= set(gaps)

print(len(gaps_total))

model_gaps = []
for x in gaps_total:
    try:
        model_gaps.append(model.reactions.get_by_id(x).id)
    except:
        pass

remove = [x for x in model_gaps if x in orphans]    

print('Reactions removed: ' + str(len(remove)))

# Create new mini-universal with just the removed reactions to test to be sure that gapfill is working
mini_uni = cobra.Model("mini_universal_reactions")
for i in remove:
    reaction = model.reactions.get_by_id(i)
    mini_uni.add_reaction(reaction.copy())
    model.remove_reactions([reaction])

model.solver = 'gurobi'

solution = gapfill(model, mini_uni, demand_reactions=False)
for reaction in solution[0]:
    print(reaction.id)