from subprocess import run

from datadispatch.access import select

from scripts import *
from scripts import A_CellSeg

#=-=#=-=#=- Logging Run -=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#

# ----- Log Git Address ----- # 
# TODO:
# safety check: make sure that the git branch is up to date so that we can 
# cross reference our records with the code it was run on

# ----- Enter in Database ----- # 
# TODO:
# safety check: make sure that the git branch is up to date so that we can 
# cross reference our records with the code it was run on



#=-=#=-=#=- Precalculating Hashes -=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#

# ----- Generate "Config Hashes" ----- # 
hashes = {}
all_config = order_snakefood()
for stepname in all_config.sections():
    h = hash_parameters(
        stepname
    )  
    hashes[stepname] = h


# ----- Get experiment Names and IDs ----- # 
where = all_config['Run Info']['where'].strip('"')
results = select('Experiment', where=where)

ids = {res.meta.MERFISHExperimentID: res.name for res in results}
print(ids)
hashes['EXPS'] = ids

# ----- Generate "Aggregate Hash" ----- # 
# this hash exists to represent the combonation of all experiments selected by the where statement
exp_names = '$'.join(sorted(ids.values()))
hashes['AGGREGATE'] = hash_strli(exp_names)

#=-=#=-=#=- Generating Targets -=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#
# A1 Cellpose: 
# Resegmentation of raw imaging data using cellpose:
A1_Cellpose_out = '{id}_{A1}.checkpoint'.format(**{'id': '{eid}', 'A1':hashes['A1_Cellpose']})
A1_Cellpose_target = [
    '{id}_{A1}.checkpoint'.format(**{'id': id, 'A1':hashes['A1_Cellpose']})
    for id in ids.keys()
]

# 


###############################################################################
#=============================================================================#
###############################################################################


rule A1_Cellpose:
    input:
        A1_Cellpose_target
# rule _cook_snakefood:
#     input: 
#     output: 
#     run:
        
# rule example:
#     input: 
#         # Example:
#         # {section_name1}_{section_name2}_...{}.something
#     output: 
#         # {section_name1}_{section_name2}_...{}.something

#     run: 
#         # subprocess(script include all parameters) -> print out and capture any output not saved to file
#         # do merbot interaction here

# rule A_test1:
#     # input: 
#     output: 
#     run: 

# rule B_test2:
#     input: 
#     output: 
#     run: 

# rule C_test3:
#     input: 
#     output: 
#     run: 

rule A1:
    output:
        A1_Cellpose_out

    threads:64
    run:
        name=ids[wildcards.eid]
        A_CellSeg.Cellpose(name, output)