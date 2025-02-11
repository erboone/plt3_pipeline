
import datadispatch
import mftools

import subprocess as sub
from datetime import datetime
import os,sys
from string import Template

import datadispatch.access as db
from datadispatch.orm import ParamLog

from scripts import _1_
from scripts import _2_
from scripts.meta import *

_output = f'{os.getenv("HOME")}/pipeline/_output'

###############################################################################
#===================== Hashing config and Logging  ===========================#
###############################################################################

# ----- Log Git Address ----- # 
# TODO: Add the other major supporting libraries (mftools) + safety check
# safety check: make sure that the git branch is up to date so that we can 
# cross reference our records with the code it was run on
if not branch_up_to_date():
    pass
    #raise RuntimeError("The branch is not up to date; run 'git status'")
commit = get_commit_hash('TODO: this parameter does nothing')

# ----- Copy snakemake file ----- # 
datetime_str = str(datetime.now()).\
                split('.')[0].\
                replace(' ', '_').\
                replace('-', '.')
sub.run(['cp', 'snakefood.ini', f'{_output}/.leftovers/{datetime_str}.ini'])


#=-=#=-=#=- Precalculating Hashes -=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#

# ----- Generate "Config Hashes" ----- # 
hashes:dict[str, any] = {}
all_config = order_snakefood()
check_RunInfo(all_config)

for stepname in all_config.sections():
    h = hash_parameters(                # in scripts.meta
        stepname
    )  
    hashes[stepname] = h


# ----- Get experiment Names and IDs ----- # 
where = all_config['Run Info']['where'].strip('"')
results = db.select('Experiment', where=where)

ids = {res.meta.MERFISHExperimentID: res.name for res in results}
print(ids)
hashes['EXPS'] = ids

# ----- Generate "Aggregate Hash" ----- # 
# This hash exists to represent the combonation of all experiments selected by the where statement
exp_names = '$'.join(sorted(ids.values()))
hashes['AGGREGATE'] = hash_strli(exp_names) # in scripts.meta

# ----- Add Misc Hashes ----- # 
hashes['COMMIT'] = commit
hashes['SNAKEMAKE_CALL'] = datetime_str
hashes['DESCRIPTION'] = all_config['Run Info']['description'].strip('"')

#=-=#=-=#=- Enter all info in Database -=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#
log_parameters(hashes, all_config) # in scripts.meta


###############################################################################
#======================== Generating Targets  ================================#
###############################################################################
# NOTE: General pattern is 
### ex.) step description -> Variable example ###
# 1.) Define for step -> A1_T:Template
# 2.) Fill hashes (Template.safe_substitute) -> A1_H:String
# 3.) Fill formats 
# 4.) Full wildcards (String.format) -> A1_{StepName}_targets:String

# for easy alteration; base_T = base template
base_T = '{_o}/{_step}/{_file}'

# --- A1 Segmentation:  ------------------------------------------------------- 
# Resegmentation of raw imaging data using cellpose and saving resulting cell
# by gene table in the form of a scanpy H5 AnnData object (.h5ad)
t = {
    '_o': _output,
    '_step': 'A1_Segmentation',
    '_file':'${exp_id}_${a11}.${form}',
}
h = {
    'a11':hashes['A11_Cellpose']
}
w = ('exp_id', [id for id in ids.keys()])

A11_Segmentation_target, A11_target = assemble_target(
    template=t,
    hashes=h,
    format=S_CHECKPOINT,
    wildcards=w
)
A12_SaveRawScanpy_target, A12_target = assemble_target(
    template=t,
    hashes=h,
    format=S_H5AD,
    wildcards=w
)   


# --- A2 Annotation:  ---------------------------------------------------------
# Automated cell label transfer from reference data
t = {
    '_o': _output,
    '_step': 'A2_Annotation',
    '_file':'${agg}_${a11}_${a21}.${form}',
}
h = {
    'a11':hashes['A11_Cellpose'],
    'a21':hashes['A21_HarmAnnotation'],
    'agg':hashes['AGGREGATE']
}

A21_Annotation_target, A21_target = assemble_target(
    template=t,
    hashes=h,
    format=S_H5AD,
    wildcards=()
)


print('A11_Segmentation_target:' ,A11_Segmentation_target)
print('A11_target:', A11_target)
print('A12_SaveRawScanpy_target:', A12_SaveRawScanpy_target)
print('A12_target:', A12_target)
print('A21_Annotation_target:', A21_Annotation_target)
print('A21_target:', A21_target)

###############################################################################
#=========================== Snakemake Rules  ================================#
###############################################################################

# rule all:
#     input:
#         A11_Segmentation_target,
#         A12_SaveRawScanpy_target,
#         A21_Annotation_target


rule A11_Segmentation:
    input:
        A11_Segmentation_target
    run:
        print('DONE:', A11_Segmentation_target)
        print('updating database parameter log with success...')
        mark_success(datetime_str)


rule A12_SaveRawScanpy:
    input:
        A12_SaveRawScanpy_target
    run:
        print(A12_SaveRawScanpy_target)
        print('updating database parameter log with success...')
        mark_success(datetime_str)



rule A21_Annotation:
    input:
        A21_Annotation_target
    run:
        print('updating database parameter log with success...')
        mark_success(datetime_str)


###############################################################################
# ---- Behind the scenes: rules to produce targets, call at your own risk --- #
###############################################################################

rule A11:
    output:
        A11_target
    threads:64
    run:
        name=ids[wildcards.exp_id]
        print('NAME HERE', name)
        _1_._1_Cellpose(name, input, output, hashes, commit)

rule A12:
    input:
        A11_target
    output:
        A12_target
    run:
        name=ids[wildcards.exp_id]
        print('NAME HERE', name)
        _1_._2_SaveRawScanpy(name, input, output, hashes, commit)


rule A21:
    input:
        A12_SaveRawScanpy_target
    output:
        A21_target
    run:
        _2_._1_Annotation(input, output, hashes, commit)


# rule B1:
#     input:
#         *A12_SaveRawScanpy_target
#     run:
#         #some script here
#         print('reached')