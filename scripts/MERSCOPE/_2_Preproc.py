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
preproc_conf = order_snakefood('Filter')


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

    cpconf = order_snakefood('A1_Cellpose')
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
        res = result.pop() 
    else:
        res = result[0]

    cpconf = order_snakefood('A1_Cellpose')
    alt_save = cpconf['alt_path']
    if alt_save:
        alt_paths = {'cellpose': Path(alt_save), 'masks': Path(alt_save) / 'masks'}
    else:
        # This is the only difference between A2 and B1
        alt_paths = {
            'cellpose': f"{res.rootdir}/*output*/{res.name}/", 
            'masks': f"{res.rootdir}/*output*/{res.name}/"
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

    mdata.obs[md.BATCH_KEY] = pd.Series(res.BICANID, dtype='category')
    mdata.obs[md.CTYPE_KEY] = pd.Series(None)
    mdata.obs['fov_y'] = pd.Series(None)
    mdata.obs['fov_x'] = pd.Series(None)
    mdata.obs['global_x'] = pd.Series(None)
    mdata.obs['global_y'] = pd.Series(None)

    # -------------------------- Save object -------------------------------- #
    # TODO: Unsafe write here, figure out a way to use safe write without the 
    # metadata produced by the cellpose pipeline
    print(Path(str(output)))
    mdata.safe_write(Path(str(output)))

################################################################################
#   Pre-processing
################################################################################
"""_summary_
"""

REMOVE_DOUBLETS = preproc_conf.getboolean('remove_doublets')
MIN_COUNTS = int(preproc_conf['min_counts'])
MIN_GENES = int(preproc_conf['min_genes'])
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
    if REMOVE_DOUBLETS:
        scrub = scr.Scrublet(mdata.X)
        try:
            dblt_score, dblt_pred = scrub.scrub_doublets(log_transform=True)
        except ValueError:
            warnings.warn("Using fewer components for PCA: LOOK INTO THIS; why so few?.")
            dblt_score, dblt_pred = scrub.scrub_doublets(log_transform=True, n_prin_comps=2)

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

    return md(adata=mdata)

def A3_Preprocessing(
        exp_name, 
        input, 
        output, 
        hashes,
        commit):
    mdata = md(sc.read_h5ad(Path(str(input))))
    filt_mdata = _filter(mdata)

    filt_mdata.log = ('step_name', 'Preprocessing')
    filt_mdata.log = ('step_hash', hashes['Filter'])
    filt_mdata.log = ('function', 'A3_Preprocessing')
    filt_mdata.log = ('commit_id', commit)
    filt_mdata.log = ('normalization', "sc.pp.normalize_total(_default params_), sc.pplog1p()")
    mdata.history.append(mdata.uns.pop(mdata.LOG_KEY))

    filt_mdata.safe_write(Path(str(output)))

def B2_Preprocessing(
        exp_name, 
        input, 
        output, 
        hashes,
        commit):
    
    mdata = md(adata=sc.read_h5ad(Path(str(input))))
    mdata = _filter(mdata)

    mdata.log = ('step_name', 'Preprocessing')
    mdata.log = ('step_hash', hashes['Filter'])
    mdata.log = ('function', 'A3_Preprocessing')
    mdata.log = ('commit_id', commit)
    mdata.log = ('normalization', "sc.pp.normalize_total(_default params_), sc.pplog1p()")

    mdata.safe_write(Path(str(output)))

################################################################################
#   Quality control 
################################################################################
"""_summary_
"""

# ----- Helper functions ------------------------------------------------------

def _bulk_correlation(means, bulkref, name, ax:plt.Axes):
    print(means)
    print(bulkref)
    common = means.index.astype(str).intersection(bulkref.index.astype(str))
    counts1 = means.loc[common]
    counts2 = bulkref.loc[common]
    corr = pearsonr(np.log1p(counts1), np.log1p(counts2))
    
    ax.scatter(counts1, counts2, s=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(f"{name} r={corr[0]:.03f}")
    ax.get_figure().savefig(f'_testfig1.{name}.png')
    return ax, corr.statistic

# ----- Quality control -------------------------------------------------------

def QC_1_postsegqc(
        input,
        output,
        hashes,
        commit):
    mdata:md = md(adata=sc.read_h5ad(Path(str(input))))
    print("reached: we made it!")
    print(sc.pp.calculate_qc_metrics(mdata, percent_top=None)[0]['total_counts'].median())
    qcconf = order_snakefood('QC1_postsegqc')
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
    # bulkref = bulkref.apply(np.log1p)
    bulkref.name = 'bulkref'

    mdata.X = mdata.layers['counts']
    
    statistics = []
    regions = list(mdata.obs['region'].cat.categories)
    for reg in regions:
        fig, axs = plt.subplots(2, 2, figsize=(8, 8))

        reg_mdata = mdata[mdata.obs['region'] == reg]


        qc_obs, qc_var = sc.pp.calculate_qc_metrics(reg_mdata, percent_top=None)
        
        n_cells = qc_obs.shape[0]
        med_trx_p_cell = qc_obs['total_counts'].median()
        med_trx_p_cell_p_gene = med_trx_p_cell / qc_var.shape[0]
        total_trans_postfilt = qc_obs['total_counts'].sum()
        med_gene_p_cell = qc_obs['n_genes_by_counts'].median()
        
        sc.pp.normalize_total(reg_mdata, target_sum=1e6, inplace=True)
        # sc.pp.log1p(reg_mdata)
        means = sc.get.obs_df(reg_mdata, list(reg_mdata.var_names)).mean(axis=0)
    
        corr_ax, corr_stat = _bulk_correlation(means, bulkref, reg, axs[0][0])
        
        sc.pp.pca(reg_mdata)
        sc.pp.neighbors(reg_mdata)
        sc.tl.umap(reg_mdata)
        sc.tl.leiden(reg_mdata)
        sc.pl.embedding(reg_mdata, basis='umap', color='leiden', ax=axs[1][0])
        sc.pl.embedding(reg_mdata, basis='spatial', color='total_counts', ax=axs[1][1])

        statistics.append([reg, corr_stat, med_trx_p_cell, med_gene_p_cell])

        fig.savefig(Path(str(output)))
        plt.clf()   

    print('region', 'pearsonr', 'median t.p.c.', 'median g.p.c.', sep='\t')
    for l in statistics:
        print('\t'.join([str(i) for i in l]))
