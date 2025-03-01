import os, sys
import json
from configparser import ConfigParser
from pathlib import Path
from itertools import chain
import warnings

import numpy as np
from skimage.measure import regionprops # This has been moved to plotting
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.axes import Axes 
from cellpose import models
import pandas as pd
from scipy.stats import pearsonr
import scanpy as sc
import scrublet as scr

from mftools.fileio import ImageDataset, MerfishAnalysis
from mftools.plotting import fov_show
from mftools.segmentation import CellSegmentation
from mftools.barcodes import assign_to_cells, link_cell_ids, create_cell_by_gene_table
from mftools.cellgene import create_scanpy_object
from mftools.mfexperiment import MerscopeExperiment
from mftools import MerData as md
from datadispatch.access import select
from datadispatch.orm import Experiment

from ..meta import order_snakefood
from ..meta import OUTPUT

# Load config 
save_conf = order_snakefood()
preproc_conf = order_snakefood()



################################################################################
#   Save raw scanpy
################################################################################
"""_summary_
"""


def A2_SaveRawScanpy(
        exp_name, 
        input, 
        output, 
        hashes,
        commit):
    # -------------------- Load and create bare scanpy object --------------- #
    EXPERIMENT_NAME = exp_name
    result = select('Experiment', where=f"Experiment.name == '{EXPERIMENT_NAME}'")
    if len(result) != 1:
        raise RuntimeError('Experiment search critera did not find unique entry.') 
    
    res = result[0]

    cpconf = order_snakefood('A11_Cellpose')
    alt_save = cpconf['alt_path']
    if alt_save:
        alt_paths = {'cellpose': Path(alt_save), 'masks': Path(alt_save) / 'masks'}
    else:
        alt_paths = {}

    e = MerscopeExperiment(res.root.path, res.name, alt_paths=alt_paths)

    mdata = e.create_scanpy_object()
    # ------------------------ Decorate object ------------------------------ #

    mdata.log = ('step_name', '_1_Segmentation')
    mdata.log = ('step_hash', hashes['A11_Cellpose']) 
    # ^^^ This is unusual, normally try to only use a steps hash in it's own 
    # step. However, this object does not exist in the Cellpose step sooo...
    mdata.log = ('function', 'A2_SaveRawScanpy')
    mdata.log = ('commit_id', commit)
    mdata.log = ('normalization', None)


    mdata.obs[md.BATCH_KEY] = pd.Series(res.meta.BICANExperimentID)
    mdata.obs[md.CTYPE_KEY] = pd.Series(None)

    # -------------------------- Save object -------------------------------- #

    mdata.safe_write(Path(str(output)))

def B1_SaveRawScanpy(
        exp_name, 
        input, 
        output, 
        hashes,
        commit):
    # -------------------- Load and create bare scanpy object --------------- #
    EXPERIMENT_NAME = exp_name
    result = select('Experiment', where=f"Experiment.name == '{EXPERIMENT_NAME}'")
    if len(result) != 1:
        raise RuntimeError('Experiment search critera did not find unique entry.') 
    
    res = result[0]

    cpconf = order_snakefood('A11_Cellpose')
    alt_save = cpconf['alt_path']
    if alt_save:
        alt_paths = {'cellpose': Path(alt_save), 'masks': Path(alt_save) / 'masks'}
    else:
        # This is the only difference between A2 and B1
        alt_paths = {
            'cellpose': f"{res.root}/*output*/{res.name}/", 
            'masks': None
        }

    e = MerscopeExperiment(res.root.path, res.name, alt_paths=alt_paths)

    mdata = e.create_scanpy_object()
    # ------------------------ Decorate object ------------------------------ #

    mdata.log = ('step_name', '_1_Segmentation')
    mdata.log = ('step_hash', None) 
    # ^^^ This is unusual, normally try to only use a steps hash in it's own 
    # step. However, this object does not exist in the Cellpose step sooo...
    mdata.log = ('function', 'B2_SaveRawScanpy')
    mdata.log = ('commit_id', commit)
    mdata.log = ('normalization', None)


    mdata.obs[md.BATCH_KEY] = pd.Series(res.meta.BICANExperimentID)
    mdata.obs[md.CTYPE_KEY] = pd.Series(None)

    # -------------------------- Save object -------------------------------- #

    mdata.safe_write(Path(str(output)))

################################################################################
#   Pre-processing
################################################################################
"""_summary_
"""

REMOVE_DOUBLETS = preproc_conf['remove_doublets']
MIN_COUNTS = preproc_conf['min_counts']
MIN_GENES = preproc_conf['min_genes']
VOLUME_KEY = 'volume'

def _filter(mdata:md):
    print(f'Filtering: prefilter')
    print(f'X is {mdata.X.dtype}\tmin:max::{mdata.X.min()}:{mdata.X.max()}')
    print()
    mdata.layers['counts'] = mdata.X.copy()
    # Basic filtering
    ncells = mdata.shape[0]
    sc.pp.filter_cells(mdata, min_counts=MIN_COUNTS)
    diff = ncells - mdata.shape[0]; print(f"{diff}({diff/ncells:.2f})cells removed by min counts"); ncells = mdata.shape[0]
    sc.pp.filter_cells(mdata, min_genes=MIN_GENES)
    diff = ncells - mdata.shape[0]; print(f"{diff}({diff/ncells:.2f})cells removed by min genes"); ncells = mdata.shape[0]

    # Scrublet
    print(type(REMOVE_DOUBLETS), int(REMOVE_DOUBLETS))
    if REMOVE_DOUBLETS:
        scrub = scr.Scrublet(mdata.X)
        dblt_score, dblt_pred = scrub.scrub_doublets(log_transform=True)
        mdata = mdata[dblt_pred]
        diff = ncells - mdata.shape[0]; print(f"{diff}({diff/ncells:.2f})cells removed by scrublet"); ncells = mdata.shape[0]

    # Volume
    if VOLUME_KEY in mdata.obs.columns:
        mdata = mdata[
            (mdata.obs[VOLUME_KEY] > 100) &
            (mdata.obs[VOLUME_KEY] < (mdata.obs[VOLUME_KEY].median() * 3))
        ]
        diff = ncells - mdata.shape[0]; print(f"{diff}({diff/ncells:.2f})cells removed by volume"); ncells = mdata.shape[0]


    # Transformation
    sc.pp.normalize_total(mdata)
    sc.pp.log1p(mdata)

    return mdata

def A3_Preporcessing(
        exp_name, 
        input, 
        output, 
        hashes,
        commit):
    mdata = md(sc.read_h5ad(input))
    filt_mdata = _filter(mdata)

    filt_mdata.log = ('step_name', 'Preprocessing')
    filt_mdata.log = ('step_hash', hashes['Preprocessing'])
    filt_mdata.log = ('function', 'A3_Preprocessing')
    filt_mdata.log = ('commit_id', commit)
    filt_mdata.log = ('normalization', "sc.pp.normalize_total(_default params_), sc.pplog1p()")

    filt_mdata.safe_write(output)

def B2_Preprocessing(
        exp_name, 
        input, 
        output, 
        hashes,
        commit):
    
    mdata = md(sc.read_h5ad(input))
    filt_mdata = _filter(mdata)

    filt_mdata.log = ('step_name', 'Preprocessing')
    filt_mdata.log = ('step_hash', hashes['Preprocessing'])
    filt_mdata.log = ('function', 'A3_Preprocessing')
    filt_mdata.log = ('commit_id', commit)
    filt_mdata.log = ('normalization', "sc.pp.normalize_total(_default params_), sc.pplog1p()")

    filt_mdata.safe_write(output)

################################################################################
#   Quality control 
################################################################################
"""_summary_
"""

# ----- Helper functions ------------------------------------------------------

def _bulk_correlation(means, bulkref, name, output):
    common = means.index.intersection(bulkref.index)    
    counts1 = means.loc[common]
    counts2 = bulkref.loc[common]
    corr = pearsonr(np.log1p(counts1), np.log1p(counts2))
    
    plt.figure(figsize=(5,5),dpi=120)
    plt.scatter(counts1, counts2, s=1)
    plt.xscale("log")
    plt.yscale("log")
    plt.title(f"{name} r={corr[0]:.03f}")

    return corr

# ----- Quality control -------------------------------------------------------

def QC_1_postsegqc(
        input,
        outfile,
        hashes,
        commit):
    mdata:md = md(adata=sc.read_h5ad(input))
    print("reached: we made it!")
    print(sc.pp.calculate_qc_metrics(mdata, percent_top=None)[0]['total_counts'].median())
    qcconf = order_snakefood('QC')
    bulkref_path = qcconf['bulkref_path']

    try:
        mdata.check_requirements()
    except:
        warnings.warn("Not all MerData requirements are met; continue with caution")

    # TODO: implement a way to load reference data, while checking that it 
    # has been properly sanitized
    bulkref = pd.read_csv(bulkref_path, sep='\t', skiprows=2).reset_index()
    bulkref = bulkref.drop_duplicates(subset='Description')
    bulkref = bulkref.drop(columns=['id', 'Name']).set_index('Description').mean(axis=1)
    bulkref.name = 'bulkref'

    means = sc.get.obs_df(mdata, list(mdata.var_names)).mean(axis=0)
    print(means)
    print(bulkref)
