import os, sys
import json
from configparser import ConfigParser
from pathlib import Path
from itertools import chain
import warnings
from glob import glob
import re

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
import anndata as an
# import scrublet as scr

from mftools.fileio import ImageDataset, MerfishAnalysis
from mftools.plotting import fov_show
from mftools.segmentation import CellSegmentation
from mftools.barcodes import assign_to_cells, link_cell_ids, create_cell_by_gene_table
from mftools.cellgene import create_scanpy_object
from mftools.mfexperiment import MerscopeExperiment
from mftools import MerData as md
from datadispatch.access import select
from datadispatch.orm import Experiment
from mftools.plotting import Correlation

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
        alt_paths = {'cellpose': str(Path(alt_save)), 'masks': str(Path(alt_save) / 'masks')}
    elif glob(f"{res.rootdir}/*output*/"):
        # This is the only difference between A2 and B1
        alt_paths = {
            'cellpose': f"{res.rootdir}/*output*/{res.name}/", 
            'masks': f"{res.rootdir}/*output*/{res.name}/"
        }
    else:
        # This is the only difference between A2 and B1
        alt_paths = {
            'cellpose': f"{res.rootdir}/{res.name}/", 
            'masks': f"{res.rootdir}/{res.name}/",
            'data': f"/mnt/merfish15/MERSCOPE/*data*/202407221121_20240722M176BICANRen22_VMSC10002/",
            'settings': f"/mnt/merfish15/MERSCOPE/*data*/202407221121_20240722M176BICANRen22_VMSC10002/settings"
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

    try:
        mdata.obs[md.BATCH_KEY] = pd.Series(res.BICANID, dtype='category')
    except:
        mdata.obs[md.BATCH_KEY] = pd.Series(res.region, dtype='category')
    mdata.obs[md.CTYPE_KEY] = pd.Series("unassigned", dtype='category')
    mdata.obs['fov_y'] = pd.Series(None, dtype=int)
    mdata.obs['fov_x'] = pd.Series(None, dtype=int)
    mdata.obs['global_x'] = pd.Series(None, dtype=float)
    mdata.obs['global_y'] = pd.Series(None, dtype=float)

    # -------------------------- Save object -------------------------------- #
    # TODO: Unsafe write here, figure out a way to use safe write without the 
    # metadata produced by the cellpose pipeline
    print(Path(str(output)))
    mdata.safe_write(Path(str(output)))

def C1_SaveRawScanpy( 
        input, 
        output, 
        hashes,
        commit):
    
    def parse_slide_quan(phrase):
        sample = re.search("[A-Z]{3}[0-9]{4}.*\$", phrase)
        _region_str = sample.group(0)
        section_code = re.search('[EQ]0[0-9]', _region_str)
        _region_str = _region_str[:section_code.span(0)[0]] + _region_str[section_code.span(0)[1]:]
        reg = _region_str[7:10]
        return [sample.group()[:7], reg.lower().strip('$'), 'quan']

    def parse_slide_ecker(phrase):
        sample = re.search('[A-Z]{3}-{,1}[0-9]{4}', phrase)
        _region_str = phrase[:sample.span(0)[0]] + phrase[sample.span(0)[1]:]
        middle_part = re.search('4[Xx]1.*0[0-9]{1}', _region_str)
        _, reg, code = middle_part.group(0).split('-')[:3]
        return [sample.group().replace('-',''), reg.lower(), 'ecker']


    mdata = sc.read_h5ad("/data/_MERSCOPE/BICAN_BG_proseg_mapped_filt.h5ad")

    words = mdata.obs[['region', 'slide']].drop_duplicates(ignore_index=True)
    for tup in words.itertuples():
        _, r, s = (tup[0], tup[1], tup[2])
        phrase =     r+'$'+s
        sub_mdata = mdata[(mdata.obs['region'] == r) & (mdata.obs['slide'] == s)].copy()
        try:
            donor, br_region, rep = parse_slide_quan(phrase)
        except AttributeError:
            donor, br_region, rep = parse_slide_ecker(phrase)
        
        sub_mdata.obs['donor'] = donor
        sub_mdata.obs['br_region'] = br_region
        sub_mdata.obs['rep'] = rep
        
        sub_mdata.write(filename=str(Path(str(output)).parent / f"{donor}_{br_region}_{rep}.h5ad"))

    with open(str(output), 'w') as f:
        f.write('This is just a flag')


################################################################################
#   Pre-processing
################################################################################
"""_summary_
"""

REMOVE_DOUBLETS = preproc_conf['remove_doublets']
MIN_COUNTS = int(preproc_conf['min_counts'])
MIN_GENES = int(preproc_conf['min_genes'])
VOLUME_KEY = 'volume'
BANKSY_KEY = '_BANKSY'

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
    # if REMOVE_DOUBLETS:
    #     scrub = scr.Scrublet(mdata.X)
    #     try:
    #         dblt_score, dblt_pred = scrub.scrub_doublets(log_transform=True)
    #     except ValueError:
    #         warnings.warn("Using fewer components for PCA: LOOK INTO THIS; why so few?.")
    #         dblt_score, dblt_pred = scrub.scrub_doublets(log_transform=True, n_prin_comps=2)

    #     mdata = mdata[dblt_pred]
    #     diff = ncells - mdata.shape[0]; print(f"{diff}({diff/ncells:.2f})cells removed by scrublet"); ncells = mdata.shape[0]

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

def _remove_alpha(text):
    result = ""
    for char in text:
        if not char.isalpha():  # Check if the character is NOT an alphabet
            result += char
    return result

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
    mdata.log = ('function', 'B2_Preprocessing')
    mdata.log = ('commit_id', commit)
    mdata.log = ('normalization', "sc.pp.normalize_total(_default params_), sc.pplog1p()")

    # Add banksy here
    from banksy_utils.filter_utils import filter_cells
    from banksy.initialize_banksy import initialize_banksy
    from banksy.embed_banksy import generate_banksy_matrix
    from banksy_utils.umap_pca import pca_umap
    from banksy.cluster_methods import run_Leiden_partition
    from banksy.main import median_dist_to_nearest_neighbour

    # set params
    # ==========
    plot_graph_weights = True
    k_geom = 15 # only for fixed type
    max_m = 1 # azumithal transform up to kth order
    nbr_weight_decay = "scaled_gaussian" # can also be "reciprocal", "uniform" or "ranked"
    coord_keys = ('center_x', 'center_y', 'X_spatial')
    WM_GENES = order_snakefood("Filter")['wm_genes']
    # Find median distance to closest neighbours, the median distance will be `sigma`

    # The following are the main hyperparameters for BANKSY
    resolutions = [0.2] # clustering resolution for UMAP
    pca_dims = [20] # Dimensionality in which PCA reduces to
    lambda_list = [0.8] # list of lambda parameters

    banksy_adata = mdata.copy()
    nbrs = median_dist_to_nearest_neighbour(banksy_adata, key=coord_keys[2])
    banksy_dict = initialize_banksy(
        banksy_adata,
        coord_keys,
        k_geom,
        nbr_weight_decay=nbr_weight_decay,
        max_m=max_m,
        plt_edge_hist=True,
        plt_nbr_weights=True,
        plt_agf_angles=False, # takes long time to plot
        plt_theta=True,
    )
    banksy_dict, banksy_matrix = generate_banksy_matrix(banksy_adata, banksy_dict, lambda_list, max_m)
    pca_umap(banksy_dict,
            pca_dims = pca_dims,
            add_umap = True,
            plt_remaining_var = False,
            )
    
    results_df, max_num_labels = run_Leiden_partition(
        banksy_dict,
        resolutions,
        num_nn = 50,
        num_iterations = -1,
        match_labels = True,
    )

    
    labels = results_df.loc['scaled_gaussian_pc20_nc0.80_r0.20', 'labels']
    banksy_adata.obs[BANKSY_KEY] = labels.dense
    banksy_adata.obs[BANKSY_KEY] = banksy_adata.obs[BANKSY_KEY].astype('string') + '_banksy'
    sc.pp.scale(banksy_adata)
    expression = sc.get.obs_df(banksy_adata, keys=banksy_adata.var_names.to_list() + [BANKSY_KEY])
    # wm_genes = expression.groupby(by=BANKSY_KEY).mean()[banksy_adata.var_names.intersection(WM_GENES)]
    # wm_clust = wm_genes.idxmax().mode()
    # banksy_adata.obs.loc[banksy_adata.obs[BANKSY_KEY].isin(wm_clust.values), BANKSY_KEY] = 'white_matter'
    
    wm_markers = banksy_adata.var_names.intersection(WM_GENES)
    wm_genes = expression.groupby(by=BANKSY_KEY).mean()[wm_markers]
    wm_clust = wm_genes.idxmax().mode()
    mdata.obs[BANKSY_KEY] = banksy_adata.obs[BANKSY_KEY].astype('category')
    dend = sc.tl.dendrogram(mdata, groupby=BANKSY_KEY, inplace=False)
    wm_ind = int(wm_clust[0].split('_')[0])
    all_wm_clust = [f"{i}_banksy" for i,x in enumerate(dend['correlation_matrix'][wm_ind] > .7) if x]
    print(all_wm_clust)
    mdata.obs[BANKSY_KEY] = mdata.obs[BANKSY_KEY].cat.add_categories('white_matter')
    mdata.obs.loc[mdata.obs[BANKSY_KEY].isin(all_wm_clust), BANKSY_KEY] = 'white_matter'
    mdata.obs[BANKSY_KEY] = mdata.obs[BANKSY_KEY].cat.remove_unused_categories()    
    mdata.obs[BANKSY_KEY] = mdata.obs[BANKSY_KEY].astype('object')
    print(mdata.obs[BANKSY_KEY])
    _nametemp = "{date}_BICAN_4x1-{reg}-Q-{num}-{dev}-{samp}.h5ad"
    for region_str in mdata.obs['region'].unique():
        subadata = mdata[mdata.obs['region'] == region_str].copy()
        print(subadata.obs.dtypes)
        date, _, dev = exp_name.split('_')
        sample = re.search('[A-Z]{3}[0-9]{4}', exp_name)
        _region_str = region_str[:sample.span(0)[0]] + region_str[sample.span(0)[1]:]
        section_code = re.search('[EQ]0[0-9]', _region_str)
        _region_str = _region_str[:section_code.span(0)[0]] + _region_str[section_code.span(0)[1]:]
        reg = _region_str

        # name = _nametemp.format([
        #     date,
        #     reg,
        #     "".join([char for char in section_code.group(0) if not char.isalpha()]), 
        #     dev,
        #     sample.group(0)
        # ])
        name = f"{date}_BICAN_4x1-{reg}-Q-{''.join([char for char in section_code.group(0) if not char.isalpha()])}-{dev}-{sample.group(0)}.h5ad"
        print(name)
        subadata.write(name)

    mdata.safe_write(Path(str(output)))

def C2_Preprocessing(
        input, 
        output, 
        hashes,
        commit):
    
    dir = Path(str(input))
    paths = pd.Series([g for g in glob(str(dir.parent /'*.h5ad'))], name='paths')
    names = pd.Series([str(Path(g).name)[:-5] for g in glob(str(dir.parent /'*.h5ad'))], name='names')
    
    dirdf = pd.DataFrame()
    dirdf['path'] = paths
    dirdf['name'] = names

    sep = lambda x: x.split("_")[-3:]
    df = pd.DataFrame(dirdf['name'].map(sep).to_list())
    dirdf[['donor', 'br_region', 'rep']] = df
    print(dirdf)
    dirdf = dirdf.set_index(['donor', 'br_region', 'rep'])

    for br_reg in ['pu']:#dirdf.index.get_level_values(1).tolist():
        qpaths = list(dirdf.loc[:, br_reg, 'quan']['path'].values)
        epaths = list(dirdf.loc[:, br_reg, 'ecker']['path'].values)

        qdata = an.concat([sc.read_h5ad(p) for p in qpaths])
        edata = an.concat([sc.read_h5ad(p) for p in epaths])
        edata.obsm['X_spatial'][0] += 6000

        mdata = an.concat([qdata, edata])

        # Add banksy here
        from banksy_utils.filter_utils import filter_cells
        from banksy.initialize_banksy import initialize_banksy
        from banksy.embed_banksy import generate_banksy_matrix
        from banksy_utils.umap_pca import pca_umap
        from banksy.cluster_methods import run_Leiden_partition
        from banksy.main import median_dist_to_nearest_neighbour

        # set params
        # ==========
        plot_graph_weights = True
        k_geom = 15 # only for fixed type
        max_m = 1 # azumithal transform up to kth order
        nbr_weight_decay = "scaled_gaussian" # can also be "reciprocal", "uniform" or "ranked"
        coord_keys = ('center_x', 'center_y', 'X_spatial')
        WM_GENES = order_snakefood("Filter")['wm_genes']
        # Find median distance to closest neighbours, the median distance will be `sigma`

        # The following are the main hyperparameters for BANKSY
        resolutions = [0.2] # clustering resolution for UMAP
        pca_dims = [20] # Dimensionality in which PCA reduces to
        lambda_list = [0.8] # list of lambda parameters

        banksy_adata = mdata.copy()
        nbrs = median_dist_to_nearest_neighbour(banksy_adata, key=coord_keys[2])
        banksy_dict = initialize_banksy(
            banksy_adata,
            coord_keys,
            k_geom,
            nbr_weight_decay=nbr_weight_decay,
            max_m=max_m,
            plt_edge_hist=True,
            plt_nbr_weights=True,
            plt_agf_angles=False, # takes long time to plot
            plt_theta=True,
        )
        banksy_dict, banksy_matrix = generate_banksy_matrix(banksy_adata, banksy_dict, lambda_list, max_m)
        pca_umap(banksy_dict,
                pca_dims = pca_dims,
                add_umap = True,
                plt_remaining_var = False,
                )
        
        results_df, max_num_labels = run_Leiden_partition(
            banksy_dict,
            resolutions,
            num_nn = 50,
            num_iterations = -1,
            match_labels = True,
        )

        
        labels = results_df.loc['scaled_gaussian_pc20_nc0.80_r0.20', 'labels']
        banksy_adata.obs[BANKSY_KEY] = labels.dense
        banksy_adata.obs[BANKSY_KEY] = banksy_adata.obs[BANKSY_KEY].astype('string') + '_banksy'
        sc.pp.scale(banksy_adata)
        expression = sc.get.obs_df(banksy_adata, keys=banksy_adata.var_names.to_list() + [BANKSY_KEY])
        # wm_genes = expression.groupby(by=BANKSY_KEY).mean()[banksy_adata.var_names.intersection(WM_GENES)]
        # wm_clust = wm_genes.idxmax().mode()
        # banksy_adata.obs.loc[banksy_adata.obs[BANKSY_KEY].isin(wm_clust.values), BANKSY_KEY] = 'white_matter'
        
        wm_markers = banksy_adata.var_names.intersection(WM_GENES)
        wm_genes = expression.groupby(by=BANKSY_KEY).mean()[wm_markers]
        wm_clust = wm_genes.idxmax().mode()
        mdata.obs[BANKSY_KEY] = banksy_adata.obs[BANKSY_KEY].astype('category')
        dend = sc.tl.dendrogram(mdata, groupby=BANKSY_KEY, inplace=False)
        wm_ind = int(wm_clust[0].split('_')[0])
        all_wm_clust = [f"{i}_banksy" for i,x in enumerate(dend['correlation_matrix'][wm_ind] > .6) if x]
        print(all_wm_clust)
        mdata.obs[BANKSY_KEY] = mdata.obs[BANKSY_KEY].cat.add_categories('white_matter')
        mdata.obs.loc[mdata.obs[BANKSY_KEY].isin(all_wm_clust), BANKSY_KEY] = 'white_matter'
        mdata.obs[BANKSY_KEY] = mdata.obs[BANKSY_KEY].cat.remove_unused_categories()    
        mdata.obs[BANKSY_KEY] = mdata.obs[BANKSY_KEY].astype('object')
        print(mdata.obs[BANKSY_KEY])

        mdata.write(Path(str(output)).parent / f"{br_reg}.h5ad")

    with open(str(output), 'w') as f:
        f.write('This is just a flag')






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
    qc2conf = order_snakefood('QC2_postannoqc')
    CONTROL_GENES = qc2conf
    

    try:
        mdata.check_requirements()
    except:
        warnings.warn("Not all MerData requirements are met; continue with caution")

    
    # TODO: implement a way to load reference data, while checking that it 
    # has been properly sanitized
    bulkref = pd.read_csv(bulkref_path, sep='\t', skiprows=2).reset_index()
    print(bulkref)
    bulkref = bulkref.drop_duplicates(subset='Description')
    bulkref = bulkref.drop(columns=['id', 'Name']).set_index('Description').mean(axis=1)
    # bulkref = bulkref.apply(np.log1p)
    bulkref.name = 'bulkref'

    mdata.X = mdata.layers['counts']
    statistics = []
    regions = list(mdata.obs['region'].cat.categories)
    fig, axs = plt.subplots(2, 2, figsize=(8, 8))
    for reg in regions:
        reg_mdata = mdata[(mdata.obs['region'] == reg) & (mdata.obs[BANKSY_KEY] != 'white_matter')].copy()

        qc_obs, qc_var = sc.pp.calculate_qc_metrics(reg_mdata, percent_top=None)
        
        n_cells = qc_obs.shape[0]
        med_trx_p_cell = qc_obs['total_counts'].median()
        med_trx_p_cell_p_gene = med_trx_p_cell / qc_var.shape[0]
        total_trans_postfilt = qc_obs['total_counts'].sum()
        med_gene_p_cell = qc_obs['n_genes_by_counts'].median()
        
        statistics.append([reg, n_cells, med_trx_p_cell, med_gene_p_cell])
        # Replicate correlation
        eckerpath = "/mnt/_MERSCOPE/ecker/"
        sample = reg[:7]
        regname = reg[7:]
        efiles = glob(eckerpath + '*.h5ad')
        _wc = ".*"
        print(f"Searching for ecker sample with patterns: {sample} {regname}...", end='')
        found = False
        for f in efiles:
            name = Path(f).name
            if re.search(_wc + _wc.join(sample) + _wc, name) and \
                re.search(_wc + regname + _wc, name):
                ecker_adata = sc.read_h5ad(f)
                found = True
                break
        if found:
            print("found: name")
            corr, corrfig = Correlation(ecker_adata, reg_mdata, ax1=f"E-{regname}-{sample}", ax2=f"Qu-{regname}-{sample}", 
            highlight=CONTROL_GENES, logscale=True)
            corrfig.savefig(Path(str(output)).parent / f"{regname}-{sample}.png")

        else:
            print("not found; skipping.")
            continue

    if 'rep' in mdata.obs.columns:        
        sub_mdata_q = mdata[mdata.obs['rep'] == 'quan']
        sub_mdata_e = mdata[mdata.obs['rep'] == 'ecker']
        res, fig = Correlation(sub_mdata_q, sub_mdata_e, ax1='quan', ax2='ecker', highlight=CONTROL_GENES)
        fig.savefig(Path(str(output)).parent / f"{reg}.repcorr.png")

    if 'br_region' in mdata.obs.columns:    
        from itertools import combinations
        for br_reg1, br_reg2 in  combinations(mdata.obs['br_region'].unique().tolist(), 2):
            sub_mdata_q = mdata[mdata.obs['br_region'] == br_reg1]
            sub_mdata_e = mdata[mdata.obs['rep'] == br_reg2]
            res, fig = Correlation(sub_mdata_q, sub_mdata_e, ax1='quan', ax2='ecker', highlight=CONTROL_GENES)
            fig.savefig(Path(str(output)).parent / f"{reg}.regcorr.png")


    sc.pp.normalize_total(mdata, target_sum=1e6, inplace=True)
    # sc.pp.log1p(reg_mdata)
    means = sc.get.obs_df(mdata, list(mdata.var_names)).mean(axis=0)

    corr_ax, corr_stat = _bulk_correlation(means, bulkref, reg, axs[0][0])

    sc.pp.pca(mdata)
    sc.pp.neighbors(mdata)
    sc.tl.umap(mdata)
    sc.tl.leiden(mdata)
    sc.pl.embedding(mdata, basis='umap', color='leiden', ax=axs[1][0])
    sc.pl.embedding(mdata, basis='spatial', color='total_counts', ax=axs[1][1])
    sc.pl.embedding(mdata, basis='spatial', color=BANKSY_KEY, ax=axs[0][1])


    fig.savefig(Path(str(output)).parent / 'QC.png')

    pd.DataFrame(statistics, columns=['region', 'n_cells', 'median t.p.c.', 'median g.p.c.']).to_csv(Path(str(output)).parent / 'stats.csv')