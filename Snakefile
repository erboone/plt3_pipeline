import warnings

with warnings.catch_warnings(action="ignore"):
    import datadispatch
    import mftools

    import subprocess as sub
    from datetime import datetime
    import os, sys
    from string import Template

    import datadispatch.access as db
    from datadispatch.orm import ParamLog

    import scripts.MERSCOPE._1_Segmentation as _1_
    import scripts.MERSCOPE._2_Preproc as _2_
    import scripts.MERSCOPE._3_Annotation as _3_
    from scripts.meta import *

_output = f'{os.getenv("HOME")}/pipeline/_output'

###############################################################################
# ===================== Hashing config and Logging  ===========================#
###############################################################################

# ----- Log Git Address ----- #
# TODO: Add the other major supporting libraries (mftools) + safety check
# safety check: make sure that the git branch is up to date so that we can
# cross reference our records with the code it was run on
if not branch_up_to_date():
    pass
    # raise RuntimeError("The branch is not up to date; run 'git status'")
commit = get_commit_hash("TODO: this parameter does nothing")

# ----- Copy snakemake file ----- #
datetime_str = str(datetime.now()).split(".")[0].replace(" ", "_").replace("-", ".")
sub.run(["cp", "snakefood.ini", f"{_output}/.leftovers/{datetime_str}.ini"])


# =-=#=-=#=- Precalculating Hashes -=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#

# ----- Generate "Config Hashes" ----- #
hashes: dict[str, any] = {}
all_config = order_snakefood()
check_RunInfo(all_config)

for stepname in all_config.sections():
    h = hash_parameters(stepname)  # in scripts.meta
    hashes[stepname] = h


# ----- Get experiment Names and IDs ----- #
where = all_config["Run Info"]["where"].strip('"')
results = db.select("Experiment", where=where)

print(results)
ids = {f"{res.BICANID}": res.name for res in results}
print(ids)
hashes["EXPS"] = ids

# ----- Generate "Aggregate Hash" ----- #
# This hash exists to represent the combonation of all experiments selected by the where statement
exp_names = "$".join(sorted(ids.values()))
hashes["AGGREGATE"] = hash_strli(exp_names)  # in scripts.meta

# ----- Add Misc Hashes ----- #
hashes["COMMIT"] = commit
hashes["SNAKEMAKE_CALL"] = datetime_str
hashes["DESCRIPTION"] = all_config["Run Info"]["description"].strip('"')

# =-=#=-=#=- Enter all info in Database -=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#
log_parameters(hashes, all_config)  # in scripts.meta


###############################################################################
# ======================== Generating Targets  ================================#
###############################################################################
# NOTE: General pattern is
### ex.) step description -> Variable example ###
# 1.) Define for step -> A1_T:Template
# 2.) Fill hashes (Template.safe_substitute) -> A1_H:String
# 3.) Fill formats
# 4.) Full wildcards (String.format) -> A1_{StepName}_targets:String

# for easy alteration; base_T = base template
base_T = "{_o}/{_step}/{_file}"

# --- A Segmentation:  -------------------------------------------------------
# Resegmentation of raw imaging data using cellpose and saving resulting cell
# by gene table in the form of a scanpy H5 AnnData object (.h5ad)
t = {
    "_o": _output,
    "_step": "A1_Segmentation",
    "_file": "${exp_id}_${a1}.${form}",
}
h = {"a1": hashes["A1_Cellpose"]}
w = ("exp_id", [id for id in ids.keys()])

A1_Cellpose_target, A1_target = assemble_target(
    template=t, hashes=h, format=S_CHECKPOINT, wildcards=w
)

# --- A Preproc:  -------------------------------------------------------
# Intermediate step after cells have been segmented and cell_by_gene data clarified.

t = {
    "_o": _output,
    "_step": "A2_Preproc",
    "_file": "${exp_id}_${a1}.${form}",
}
h = {"a1": hashes["A1_Cellpose"]}
w = ("exp_id", [id for id in ids.keys()])

A2_SaveRawScanpy_target, A2_target = assemble_target(
    template=t, hashes=h, format=S_H5AD, wildcards=w
)

t = {
    "_o": _output,
    "_step": "A2_Preproc",
    "_file": "${exp_id}_${a1}_${a3}.filt.${form}",
}
h = {"a1": hashes["A1_Cellpose"], "a3": hashes["Filter"]}
w = ("exp_id", [id for id in ids.keys()])

A3_Filtering_target, A3_target = assemble_target(
    template=t, hashes=h, format=S_H5AD, wildcards=w
)

# --- A Annotation:  ---------------------------------------------------------
# Automated cell label transfer from reference data

t = {
    "_o": _output,
    "_step": "A3_Annotation",
    "_file": "${agg}_${a1}_${a4}.${form}",
}
h = {
    "a1": hashes["A1_Cellpose"],
    "a4": hashes["A4_HarmAnnotation"],
    "agg": hashes["AGGREGATE"],
}

A4_Annotation_target, A4_target = assemble_target(
    template=t, hashes=h, format=S_H5AD, wildcards=()
)

# --- B Preproc:  ---------------------------------------------------------

t = {
    "_o": _output,
    "_step": "B2_Preproc",
    "_file": "${exp_id}.${form}",
}
h = {}
w = ("exp_id", [id for id in ids.keys()])

B1_SaveRawScanpy_target, B1_target = assemble_target(
    template=t, hashes=h, format=S_H5AD, wildcards=w
)

t = {
    "_o": _output,
    "_step": "B2_Preproc",
    "_file": "${exp_id}_${b2}.filt.${form}",
}
h = {"b2": hashes["Filter"]}
w = ("exp_id", [id for id in ids.keys()])

B2_Filtering_target, B2_target = assemble_target(
    template=t, hashes=h, format=S_H5AD, wildcards=w
)

t = {
    "_o": _output,
    "_step": "B2_Preproc",
    "_file": "${exp_id}_${b2}.QC.${form}",
}
h = {"b2": hashes["Filter"]}
w = ("exp_id", [id for id in ids.keys()])

BQC1_PostprocQC_target, BQC1_target = assemble_target(
    template=t, hashes=h, format=S_IMAGE, wildcards=w
)

t = {
    "_o": _output,
    "_step": "B3_Annotation",
    "_file": "${agg}.${b2}.${b3}.${form}",
}
h = {
    "b2": hashes["Filter"],
    "b3": hashes["A4_Annotation"],
    "agg": hashes["AGGREGATE"]
}

A4_Annotation_target, A4_target = assemble_target(
    template=t, hashes=h, format=S_H5AD, wildcards=()
)


print()
print("-" * 80)
print()
print("A1_Cellpose_target:", A1_Cellpose_target)
print("A1_target:", A1_target)
print("A2_SaveRawScanpy_target:", A2_SaveRawScanpy_target)
print("A2_target:", A2_target)
print("A3_Filtering_target:", A3_Filtering_target) 
print("A3_target:", A3_target)
print("A4_Annotation_target:", A4_Annotation_target)
print("A4_target:", A4_target)
print()
print("B1_SaveRawScanpy_target:", B1_SaveRawScanpy_target)
print("B1_target:", B1_target)
print("B2_Filterin_target:", B2_Filtering_target)
print("B2_target:", B2_target)
print("BQC1_PostprocQC_target:", BQC1_PostprocQC_target)
print("BQC1_target:", BQC1_target)
print("B3_Annotation_target:", B3_Annotation_target)
print("B3_target:", B3_target)

print()
print("-" * 80)
print()
###############################################################################
# ====================== Snakemake Rules - Path A  ============================#
###############################################################################

# rule all:
#     input:
#         A1_Cellpose_target,
#         A2_SaveRawScanpy_target,
#         A4_Annotation_target


rule A1_Segmentation:
    input:
        A1_Cellpose_target,
    run:
        print("DONE:", A1_Cellpose_target)
        print("updating database parameter log with success...")
        mark_success(datetime_str)


rule A2_SaveRawScanpy:
    input:
        A2_SaveRawScanpy_target,
    run:
        print(A2_SaveRawScanpy_target)
        print("updating database parameter log with success...")
        mark_success(datetime_str)


rule A3_Preproc:
    input:
        A3_Filtering_target,
    run:
        print("updating database parameter log with success...")
        mark_success(datetime_str)


rule A4_Annotation:
    input:
        A4_Annotation_target,
    run:
        print("updating database parameter log with success...")
        mark_success(datetime_str)


# ---- Behind the scenes: rules to produce targets, call at your own risk --- #


rule A1:
    output:
        A1_target,
    threads: 64
    run:
        name = ids[wildcards.exp_id]
        print("NAME HERE", name, "Output", output)
        _1_.A1_Cellpose(name, input, output, hashes, commit)


rule A2:
    input:
        A1_target,
    output:
        A2_target,
    run:
        name = ids[wildcards.exp_id]
        print("NAME HERE", name)
        _2_.A2_SaveRawScanpy(name, input, output, hashes, commit)


rule A3:
    input:
        A2_target,
    output:
        A3_target,
    run:
        name = ids[wildcards.exp_id]
        print("NAME HERE", name)
        _2_.A3_Preprocessing(name, input, output, hashes, commit)


rule A4:
    input:
        A3_Filtering_target,
    output:
        A4_Annotation_target,
    run:
        _3_.A4_HarmAnnotation(input, output, hashes, commit)


# rule B1:
#     input:
#         *A2_SaveRawScanpy_target
#     run:
#         #some script here
#         print('reached')

###############################################################################
# ====================== Snakemake Rules - Path B  ============================#
###############################################################################


rule B1_SaveRawScanpy:
    input:
        B1_SaveRawScanpy_target,
    run:
        print("updating database parameter log with success...")
        mark_success(datetime_str)


rule B2_Filtering:
    input:
        B2_Filtering_target,
    run:
        print("updating database parameter log with success...")
        mark_success(datetime_str)


rule BQC1_PostprocQC:
    input:
        BQC1_PostprocQC_target
    run:
        print("updating database parameter log with success...")
        mark_success(datetime_str)

rule B3_Annotation:
    input:
        BQC1_PostprocQC_target
    run:
        print("updating database parameter log with success...")
        mark_success(datetime_str)


rule B1:
    output:
        B1_target,
    threads: 64
    wildcard_constraints:
        exp_id="[A-z0-9]{3,5}",
    run:
        name = ids[wildcards.exp_id]
        print("NAME HERE", name, "Output", output)
        _2_.B1_SaveRawScanpy(name, input, output, hashes, commit)


rule B2:
    input:
        B1_target,
    output:
        B2_target,
    run:
        name = ids[wildcards.exp_id]
        print("NAME HERE", name)
        _2_.B2_Preprocessing(name, input, output, hashes, commit)


rule BQC1:
    input:
        B2_target,
    output:
        BQC1_target,
    run:
        name = ids[wildcards.exp_id]
        _2_.QC_1_postsegqc(input, output, hashes, commit)

rule B3:
