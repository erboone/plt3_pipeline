import warnings as warn
from pathlib import Path
from datetime import datetime
import os
from glob import glob
import re

import scanpy as sc
import pandas as pd
import anndata as an
import numpy as np
import scrublet as scr
from sklearn.neighbors import KNeighborsClassifier
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import matplotlib
import seaborn as sns

from datadispatch.access import select
from mftools.scanpy_wrap import MerData as md
from mftools.mfexperiment import MerscopeExperiment
from mftools.plotting import Correlation, celltype_corr

from ..meta import order_snakefood
from ..meta import OUTPUT
from ..meta import read_ref_table
from . import _helper as help

SOURCE_KEY = 'source'
_PCA_HARM = 'X_pca_harmony'
_NEIGH_HARM = 'neighbors_harmony'

def _check_sanitized(refdata:an):
    #TODO: check that the reference adata object has been properly sanitized.
    # refdata.obs[CELLTYPE_KEY]
    try:
        refdata.obs[md.CTYPE_KEY]
    except KeyError as e:
        raise KeyError("Make sure that the reference data object being loaded has been sanitized.") from e
    
def get_ref_path(adata):
    for region_str in adata.obs['region'].unique():
        sample = re.search('[A-Z]{3}[0-9]{4}', region_str)
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
        return [read_ref_table(reg)]



def A4_HarmAnnotation(input, output, hashes, commit):
    # TODO: make all of these capable of taking a list or not
    STEPNAME = 'A4_HarmAnnotation'
    anconf = order_snakefood(STEPNAME)
    filtconf = order_snakefood('Filter')
    CELLTYPE_KEY = anconf['celltype_key']
    if not isinstance(CELLTYPE_KEY, list):
        CELLTYPE_KEY = [CELLTYPE_KEY]
    MIN_COUNTS = int(filtconf['min_counts'])
    MIN_GENES = int(filtconf['min_genes'])
    OUTPATH = Path(str(output))
    print(OUTPUT)
    os.makedirs(OUTPATH.parent, exist_ok=True)

    # QC_PATH = Path()
    try:
        REMOVE_DOUBLETS = filtconf['remove_doublets']
    except:
        warn.warn("Cannot parse text in 'anconf['remove_doublets']' as a boolean; defaulting to True")
        print("Cannot parse text in 'anconf['remove_doublets']' as a boolean; defaulting to True")
        REMOVE_DOUBLETS = True

    print('\n --- Preparing AnnData objects --- ')
    input_h5ad_paths = input

    # MERSCOPE data
    #TODO: Pull this into a function
    adatas = []
    for path in input_h5ad_paths:
        print(f'{path} Loading...')
        adata = sc.read_h5ad(path)
        if 'counts' in adata.layers:
            adata.X = adata.layers['counts']
        adatas.append(adata)
    merdata = sc.concat(adatas)
    print("Combined adata X min:max", f'{merdata.X.min()}:{merdata.X.max()}')

    ref_h5ad_path = get_ref_path(merdata)
    print('Using ref path:', ref_h5ad_path)

    # Reference data
    #TODO: Pull this into a function
    # TODO: Save concatenated reference?
    adatas = []
    for path in ref_h5ad_path:
        print(f'{path} Loading...')
        adata = sc.read_h5ad(path)
        adatas.append(adata)
    refdata = sc.concat(adatas)
    print(*refdata.obs.columns)
    # refdata = refdata[~refdata.obs[CELLTYPE_KEY].isna()]
    print("Refdata:", refdata.shape)
    for ctkey in CELLTYPE_KEY:
        refdata.obs[ctkey] = refdata.obs[ctkey].cat.add_categories('unlabeled')
        refdata.obs.loc[refdata.obs[ctkey].isna(), ctkey] = 'unlabeled'
    try:
        refdata.uns.pop(f'{md.CTYPE_KEY}_colors')
    except KeyError:
        pass

    # clust_table = pd.read_csv('/mnt/merfish18/BICAN/Reference_Data/Linnarsson/media-2.csv').dropna(how='all')
    # clust_table = clust_table[['Cluster', 'Supercluster', 'Class auto-annotation', 'Neurotransmitter auto-annotation']]
    # clust_table['Cluster'] = clust_table['Cluster'].astype('int')
    # clust_table.set_index('Cluster', inplace=True)
    # refdata.obs = refdata.obs.join(clust_table, on='Clusters', lsuffix='$')

    # Subset both to common genes only
    common_genes = merdata.var_names.intersection(refdata.var_names)
    merdata = merdata[:, common_genes].copy()
    refdata = refdata[:, common_genes].copy()

    ref_counts = pd.DataFrame(refdata.obs.groupby(by=CELLTYPE_KEY[0]).size(), columns=['reference_counts'])
    # normalize and scale data
    for desig, a in zip(['merFISH data', 'Reference data'],[merdata, refdata]):
        # print(f'\n{desig}:')
        # print(f'X is {a.X.dtype}\tmin:max::{a.X.min()}:{a.X.max()}')
        # print()
        print(f'Normalizing {desig} data...')
        # Basic filtering
        ncells = a.shape[0]
        sc.pp.filter_cells(a, min_counts=MIN_COUNTS)
        diff = ncells - a.shape[0]; print(f"{diff}({diff/ncells:.2f})cells removed by min counts"); ncells = a.shape[0]
        sc.pp.filter_cells(a, min_genes=MIN_GENES)
        diff = ncells - a.shape[0]; print(f"{diff}({diff/ncells:.2f})cells removed by min genes"); ncells = a.shape[0]

        # # Scrublet
        # print(type(REMOVE_DOUBLETS), int(REMOVE_DOUBLETS))
        # if REMOVE_DOUBLETS:
        #     scrub = scr.Scrublet(a.X)
        #     dblt_score, dblt_pred = scrub.scrub_doublets(log_transform=True)
        #     a = a[dblt_pred]
        #     diff = ncells - a.shape[0]; print(f"{diff}({diff/ncells:.2f})cells removed by scrublet"); ncells = a.shape[0]

        # # Volume
        # if desig == 'merFISH data':
        #     a = a[
        #         (a.obs['volume'] > 100) &
        #         (a.obs['volume'] < (a.obs['volume'].median() * 3))
        #     ]
        #     diff = ncells - a.shape[0]; print(f"{diff}({diff/ncells:.2f})cells removed by volume"); ncells = a.shape[0]


        # Transformation 
        sc.pp.normalize_total(a)
        sc.pp.log1p(a)
        sc.pp.scale(a)

    print()

    plt.style.use('dark_background')
    regadata = []
    try:
        merdata.obs['region']
    except KeyError:
        # if 'region' not deliniated, region name is assigned based on hash of
        # shape and min/max of X this implies that if the hash of a region is 
        # the same as another, it is almost certainly from the h5ad file. 
        
        code = str(hash(f'{merdata.shape[0]}{merdata.shape[1]}{merdata.X.min()}{merdata.X.max()}'))[-8]
        merdata.obs['region'] = f'unnamed_{code}'
        print(f'regions not specified, using randomized hash: unnamed_{code}')

    for region in merdata.obs['region'].unique():
        
        # Combining merFISH and reference data sets
        merdata_view = merdata[merdata.obs['region'] == region]
        print(f'Concatenating...merFISH:{merdata_view.shape} ref:{refdata.shape}')
        adata = an.concat(
            [merdata_view, refdata], 
            label=SOURCE_KEY, 
            keys=['merdata', 'refdata'],
            join='outer'
        )
        if 'original_cell_id' in merdata_view.obs.columns:
            adata.obs['original_cell_id'] = adata.obs['original_cell_id'].replace(np.nan, -1)
            adata.obs['original_cell_id'] = adata.obs['original_cell_id'].astype(int)
            adata.obs.loc[merdata_view.obs_names, 'original_cell_id'] = merdata_view.obs['original_cell_id']

        adata.X = np.nan_to_num(adata.X, nan=0, copy=True)

        # Set up plots for integration QC
        fig, axs = plt.subplots(1, 3, figsize=(9, 3))
        print("Preintegration UMAP", end=' -- ')
        # Preintegration umap
        print("PCA", end=' ')
        adata.obsm.pop('X_pca')
        sc.pp.pca(adata)
        print("-> Neighbors", end=' ')
        sc.pp.neighbors(adata, n_neighbors=30, metric='euclidean')
        print("-> UMAP", end=' ')
        sc.tl.umap(adata, min_dist=1)
        sc.pl.embedding(adata, basis='umap', color=SOURCE_KEY, ax=axs[0])
        
        # Integration
        print('\n --- Harmony Integration --- ')
        sc.external.pp.harmony_integrate(adata, key=SOURCE_KEY, basis='X_pca', adjusted_basis=_PCA_HARM, max_iter_harmony=50)
        sc.pp.neighbors(adata, n_neighbors=30, metric='euclidean', use_rep=_PCA_HARM, key_added=_NEIGH_HARM)
        sc.tl.umap(adata, min_dist=.3, neighbors_key=_NEIGH_HARM)
        sc.pl.embedding(adata, basis='umap', color=SOURCE_KEY, ax=axs[1])

        newadata = adata[adata.obs[SOURCE_KEY] == "merdata"]
        traindata = adata[adata.obs[SOURCE_KEY] == "refdata"] 
        
        nn = KNeighborsClassifier(n_jobs=16)
        nn.fit(traindata.obsm["X_umap"], traindata.obs[CELLTYPE_KEY[0]])
        pred = nn.predict(newadata.obsm["X_umap"])
        newadata.obs[md.CTYPE_KEY] = pred

        if len(CELLTYPE_KEY) > 1:
            for key in CELLTYPE_KEY[1:]:
                nn = KNeighborsClassifier(n_jobs=16)
                nn.fit(traindata.obsm["X_umap"], traindata.obs[key])
                pred = nn.predict(newadata.obsm["X_umap"])
                newadata.obs[key] = pred
        
        regadata.append(newadata)
        # # Postintegration
        # sc.pp.neighbors(newadata, n_neighbors=30, metric='correlation')
        # sc.tl.umap(newadata, min_dist=0.1)
        # sc.pl.embedding(newadata, basis='umap', color='region', ax=axs[2])
        
        plt.savefig(OUTPATH.parent / f'{region}_integration.png')

    

    recomb = an.concat(regadata)
    recomb.obs[[md.CTYPE_KEY] + CELLTYPE_KEY[1:]].to_csv(OUTPATH.stem + ".anno.csv")
    recomb.uns['reference_counts'] = ref_counts
    # TODO: add Merdata meta annotation and safewrite

    recomb.write(OUTPATH)

B3_HarmAnnotation = A4_HarmAnnotation

# TODO: move this an unify
def embedding_highlight(
        adata:sc.AnnData,
        highlight:str | list[str],
        basis:str,
        color:str,
        save:str=None,
        scale = 1,
        face='black',
        dpi:int=600,
        kwargs={},
        na_color=(.96, .96, .96)
        ):
    
    sc.set_figure_params(dpi_save=dpi, frameon=False)
    if face == 'black':
        plt.rc('text', color='w')
    elif face == 'white':
        plt.rc('text', color='black')

    # Mimics default size calculation of sc.pl.embedding()
    s = scale * (100000 / adata.shape[0])

    ax = sc.pl.embedding(
        adata=adata,
        basis=basis,
        show=False,
        na_color=na_color,
        **kwargs
    
    )
    ax.set_facecolor(face)
    ax.get_figure().set_facecolor(face)
    
    ax = sc.pl.embedding(
        adata=adata[adata.obs[color].isin(highlight)],
        basis=basis,
        color=color,
        ax=ax,
        size=s,
        show=False,
        **kwargs
    )

    return ax.get_figure()


###############################################################################
###############################################################################

def QC_2_postanno(
        input,
        output,
        hashes,
        commit):
    INPATH = Path(str(input))
    OUTPATH = Path(str(output))
    print(INPATH)
    print(OUTPATH)

    qcconf = order_snakefood('QC2_postannoqc')
    BULKREF_PATH = qcconf['bulkref_path']
    HIGHLIGHT = qcconf['spatial_highlight'] 
    GENES_OF_INT = qcconf['genes_of_interest']
    CONTROL_GENES = qcconf['control_genes']

    spatial_key = order_snakefood('A4_HarmAnnotation')['celltype_key'][0]
    MASKS = qcconf['masks']
    

    if qcconf["spatial_ctkey"] is not None:
        SPATCELLTYPE_KEY = qcconf['spatial_ctkey']
    else:
        SPATCELLTYPE_KEY = md.CTYPE_KEY
    OUTPATH_DIR = OUTPATH.parent
    os.makedirs(OUTPATH_DIR, exist_ok=True)


    # Bulk reference ETL
    bulkref = pd.read_csv(BULKREF_PATH, sep='\t', skiprows=2).reset_index()
    bulkref = bulkref.drop_duplicates(subset='Description')
    bulkref = bulkref.drop(columns=['id', 'Name']).set_index('Description').mean(axis=1)
    # bulkref = bulkref.apply(np.log1p)
    bulkref.name = 'bulkref'

    adata:an.AnnData = sc.read_h5ad(INPATH)
    adata.obs[md.CTYPE_KEY] = adata.obs[md.CTYPE_KEY].astype('category')
    sc.pl.embedding(adata, basis='umap', color=md.CTYPE_KEY, show=False)
    stat_df = pd.DataFrame(columns=[
        'sample',
        'premask.n.cells',
        'premask.blk.corr',
        'n.cells',
        'avg.trx.p.cell',
        'med.trx.p.cell',
        'med.genes.p.cell',
        'blk.corr'
    ])
    ctrlstat_df = pd.DataFrame(columns=[
        'sample',
        'total.cells',
        'med.ctlx.p.cell',
        'avg.ctlx.p.cell',
        'med.ctlg.p.cell',
        'ctl.blk.corr',
        'std.ctl.cells.p.fov',
        'avg.ctrlg.rank'
    ])

    stat_df.set_index('sample', inplace=True)
    ctrlstat_df.set_index('sample', inplace=True)


    # add mask columns
    cell_masks = pd.Series(name='mask')
    mask_paths = [Path(g) for g in glob(MASKS)]
    for mp in mask_paths:
        label = mp.name.split('_')[1]
        _cells = pd.Index(pd.read_csv(mp).values.flatten())
        cell_masks = cell_masks.reindex(cell_masks.index.append(_cells))
        cell_masks.loc[_cells] = label

    adata.obs['mask'] = cell_masks.astype('category')

    print(cell_masks)
    # print(adata.obs['mask'].groupby(by='mask').size())
        
    for reg in adata.obs['region'].unique():
        _stats = []
        _sub_adata = adata[(adata.obs['region'] == reg)]
        _stats.append(_sub_adata.shape[0]) # for stats: premask.n.cell
        _stats.append(Correlation(_sub_adata, bulkref, ax1=f"{reg}.unmasked", ax2="bulk", highlight=CONTROL_GENES)[0])  # for stats: premask.blk.corr
        sub_adata = _sub_adata[(_sub_adata.obs['mask'].isna())].copy()
        
        sub_adata.X = sub_adata.layers['counts']
        _stats.append(sub_adata.shape[0]) # For stats: n.cells
        obs_qc, var_qc = sc.pp.calculate_qc_metrics(sub_adata, percent_top=None)
        print(obs_qc, var_qc)
        _stats.append(obs_qc['total_counts'].mean()) # for stats: med.trx.p.cell
        _stats.append(obs_qc['total_counts'].median()) # for stats: med.trx.p.cell
        _stats.append(obs_qc['n_genes_by_counts'].median()) # For stats: med.genes.p.cell


        # Bulk correlation
        plt.clf()
        plt.style.use('default')
        bcorr, fig = Correlation(sub_adata, bulkref, ax1=reg, ax2="bulk", highlight=CONTROL_GENES, logscale=True)
        _stats.append(bcorr) # blk.corr
        fig.savefig(OUTPATH_DIR / f'{reg}_5bulkcorr.png', bbox_inches='tight')
        plt.clf()

        plt.style.use('default')
        # Top markers
        nums = sub_adata.obs.groupby(by=md.CTYPE_KEY).size()
        drop = nums[nums < 2].index
        temp = sub_adata[~sub_adata.obs[md.CTYPE_KEY].isin(drop)]
        sc.tl.rank_genes_groups(temp, groupby=md.CTYPE_KEY, n_genes=3)
        fig = sc.pl.rank_genes_groups_dotplot(temp, standard_scale='var', return_fig=True)
        fig.savefig(OUTPATH_DIR / f'{reg}_1.1mark.png', bbox_inches='tight')
        
        # Minimal marker table
        sc.pl.dotplot(temp, groupby=md.CTYPE_KEY, var_names=GENES_OF_INT, standard_scale='var', dendrogram=True, 
                      title=f'{reg} marker genes', return_fig=True) 
        fig.savefig(OUTPATH_DIR / f'{reg}_1.2mark.png', bbox_inches='tight')
        
        # Cell counts
        fig, ax = plt.subplots(1, 1, figsize=(5, 7))
        a = sub_adata.copy()# [pu_adata.obs['batch'] == batch]
        counts_df = pd.DataFrame(a.obs.groupby(by=md.CTYPE_KEY).size(), columns=['counts']).reset_index()
        # print(counts_df)
        abs_values = counts_df['counts']
        tot = sum(abs_values.astype('int'))
        rel_values = abs_values.apply(lambda x: x/tot) * 100
        lbls = [f' {p[0]} ({p[1]:.0f}%)' for p in zip(abs_values, rel_values)]
        ax = sns.barplot(counts_df,
                            x='counts', y=md.CTYPE_KEY, gap=.05, ax=ax)
        ax.bar_label(container=ax.containers[0], labels=lbls)
        ax.set_title(f"{reg} counts")
        ax.get_figure().savefig(OUTPATH_DIR / f'{reg}_2.1cellcts.png', bbox_inches='tight')

        # Cell counts (compare to ref)
        fig, ax1 = plt.subplots(1, 1, figsize=(5, 7))
        counts_df.set_index(md.CTYPE_KEY, inplace=True)
        counts_df['ref_counts'] = a.uns['reference_counts']['reference_counts']
        counts_df['norm_counts'] = counts_df['counts'] / counts_df['counts'].sum()
        counts_df['norm_ref_counts'] = counts_df['ref_counts'] / counts_df['ref_counts'].sum()
        counts_df = counts_df.drop(['counts', 'ref_counts'], axis=1)
        counts_df = counts_df.stack().reset_index()
        counts_df.columns = [md.CTYPE_KEY, 'source', 'counts']
        # print(counts_df)
        # print(counts_df)
        # abs_values = counts_df['counts']
        # tot = sum(abs_values.astype('int'))
        # rel_values = abs_values.apply(lambda x: x/tot) * 100
        # lbls = [f' {p[0]} ({p[1]:.0f}%)' for p in zip(abs_values, rel_values)]
        ax = sns.catplot(counts_df, x='counts', y='_CELLTYPE', hue='source', kind='bar')
        plt.savefig(OUTPATH_DIR / f'{reg}_2.2cellcts.png', bbox_inches='tight')
        plt.clf()

        # Cell counts (of reference)
        fig, ax = plt.subplots(1, 1, figsize=(5, 7))
        a = sub_adata.copy()# [pu_adata.obs['batch'] == batch]
        counts_df = a.uns['reference_counts'].reset_index()
        counts_df.columns = [md.CTYPE_KEY, 'reference_counts']
        abs_values = counts_df['reference_counts']
        tot = sum(abs_values.astype('int'))
        rel_values = abs_values.apply(lambda x: x/tot) * 100
        lbls = [f' {p[0]} ({p[1]:.0f}%)' for p in zip(abs_values, rel_values)]
        ax = sns.barplot(counts_df,
                            x='reference_counts', y=md.CTYPE_KEY, gap=.05, ax=ax, color='tab:orange')
        ax.bar_label(container=ax.containers[0], labels=lbls)
        ax.set_title("Reference counts")
        ax.get_figure().savefig(OUTPATH_DIR / f'{reg}_2.3cellcts.png', bbox_inches='tight')
        plt.clf()


        # spatial
        plt.style.use('dark_background')
        # fig = sc.pl.embedding(sub_adata, basis='spatial', color=md.CTYPE_KEY, return_fig=True, color_map='Set1')
        # fig.savefig(OUTPATH_DIR / f'{reg}_3spat.png')
        if isinstance(HIGHLIGHT, dict):
            H_PALETTE= HIGHLIGHT
            HIGHLIGHT=HIGHLIGHT.keys()
        fig = embedding_highlight(sub_adata, highlight=HIGHLIGHT, basis='spatial', color=SPATCELLTYPE_KEY,
                                   scale=1, dpi=1000, na_color=(.30, .30, .30),
                                   kwargs={'palette':H_PALETTE})
        fig.savefig(OUTPATH_DIR / f'{reg}_3.1spat.png', bbox_inches='tight') # bbox_extra_artists=(lgd,text),
        plt.style.use('dark_background')
        # fig = sc.pl.embedding(sub_adata, basis='spatial', color=md.CTYPE_KEY, return_fig=True, color_map='Set1')
        # fig.savefig(OUTPATH_DIR / f'{reg}_3spat.png')
        if isinstance(HIGHLIGHT, dict):
            H_PALETTE= HIGHLIGHT
            HIGHLIGHT=HIGHLIGHT.keys()
        fig = sc.pl.embedding(sub_adata, basis='spatial', color=md.CTYPE_KEY, return_fig=True)
        fig.savefig(OUTPATH_DIR / f'{reg}_3.2spat.png', bbox_inches='tight') # bbox_extra_artists=(lgd,text),

        # umap
        fig = sc.pl.embedding(sub_adata, basis='umap', color=md.CTYPE_KEY, return_fig=True)
        fig.savefig(OUTPATH_DIR / f'{reg}_4umap.png', bbox_inches='tight')
        plt.clf()
        # plt.rcParams.update(plt.rcParamsDefault)

        # cell-type correlations
        matplotlib.rcParams.update(matplotlib.rcParamsDefault)
        anconf = order_snakefood('A4_HarmAnnotation')
        ref_h5ad_path = get_ref_path(sub_adata)
        adatas = []
        for path in ref_h5ad_path:
            print(f'{path} Loading...')
            _adata = sc.read_h5ad(path)
            adatas.append(_adata)
        refdata = sc.concat(adatas)
        ax = celltype_corr(sub_adata, refdata, (md.CTYPE_KEY, spatial_key))
        ax.get_figure().savefig(OUTPATH_DIR / f'{reg}_5ctcorr.png', bbox_inches='tight')
    

        # Control genes spatial distribution and stats 
        plt.style.use('dark_background')
        fig = sc.pl.embedding(sub_adata, basis='spatial', color=CONTROL_GENES, return_fig=True, vmax=1)
        fig.savefig(OUTPATH_DIR / f'{reg}_6ctrlspat.png', bbox_inches='tight')
        plt.clf()
        
        _ctrlstats = []
        cont_adata = sub_adata[:, CONTROL_GENES]
        obs_qc, var_qc = sc.pp.calculate_qc_metrics(cont_adata, percent_top=None)
         
        _ctrlstats.append((obs_qc['total_counts'] > 1).sum()) # for stats: med.ctlx.p.cell
        _ctrlstats.append(obs_qc['total_counts'].median()) # for stats: med.ctlx.p.cell
        _ctrlstats.append(obs_qc['total_counts'].mean()) # for stats: avg.ctlx.p.cell
        _ctrlstats.append(obs_qc['n_genes_by_counts'].median()) # for stats: med.ctlg.p.cell  
        _ctrlstats.append(Correlation(sub_adata, bulkref, ax1=reg, ax2="bulk", highlight=CONTROL_GENES, logscale=True)[0]) # for stats: ctl.blk.corr  
        _ctrlstats.append(cont_adata.obs[obs_qc['total_counts'] > 0].groupby(by='fov').size().std())
        
        norm_cont_adata = sc.pp.normalize_total(sub_adata, copy=True)
        ctgrp = sc.get.obs_df(norm_cont_adata, keys=(CONTROL_GENES + [md.CTYPE_KEY]))
        ctgrp_avgexp = pd.DataFrame(ctgrp.groupby(by=md.CTYPE_KEY).mean())
        ctgrp_avgexp = ctgrp_avgexp.stack().reset_index()
        ctgrp_avgexp.columns = ['celltype', 'gene', 'avg.count']
        fig = sns.catplot(ctgrp_avgexp, x='gene', y='avg.count', hue='celltype', kind='bar', aspect=1.4)
        for ax in fig.axes.flat: ax.tick_params("x", rotation=90); ax.set_ylim(0, 1.5)
        fig.savefig(OUTPATH_DIR / f'{reg}_tempname.png', bbox_inches='tight')

        # Gene counts analysis, abs, detection
        plt.style.use('default')
        fig, ax = plt.subplots(1, 1, figsize=(20, 5))
        gene_ex_data = pd.DataFrame(sub_adata.X.sum(axis=0), index=sub_adata.var_names.astype(str), columns=['counts'])
        gene_ex_data['counts'] = np.log10(gene_ex_data['counts'].values)
        gene_ex_data['cat'] = "non-control"
        gene_ex_data.loc[CONTROL_GENES, 'cat'] = "control"
        palette = {"non-control":'tab:blue', "control":'tab:red'}
        gene_ex_data['xticks'] = gene_ex_data.index
        print(gene_ex_data)
        ctrl_gene_index = gene_ex_data.index.difference(CONTROL_GENES)
        gene_ex_data.loc[ctrl_gene_index, 'xticks'] = ""
        gene_ex_data.sort_values(by='counts', ascending=False, inplace=True)
        gene_ex_data['rank'] = range(1, gene_ex_data.shape[0]+1)
        _ctrlstats.append(gene_ex_data.loc[ctrl_gene_index, 'rank'].mean())
        gene_ex_data.reset_index(inplace=True, names='genes')
        ax = sns.barplot(gene_ex_data, x='genes', y='counts', hue='cat',
                          log_scale=False, ax=ax)#, palette=palette, aspect=3, )
        ax.set_xticklabels(gene_ex_data['xticks'])
        ax.tick_params("x", rotation=90)
        ax.get_figure().savefig(OUTPATH_DIR / f'{reg}_7absgene.png', bbox_inches='tight')
        plt.clf()
        
        stat_df.loc[reg] = _stats
        print(ctrlstat_df.columns)
        print(_ctrlstats)
        ctrlstat_df.loc[reg] = _ctrlstats
        plt.close("all")
        
    _sub_adata = adata.copy()
    reg = "combined"
    _stats.append(_sub_adata.shape[0]) # for stats: premask.n.cell
    _stats.append(Correlation(_sub_adata, bulkref, ax1=f"{reg}.unmasked", ax2="bulk", highlight=CONTROL_GENES)[0])  # for stats: premask.blk.corr
    sub_adata = _sub_adata[(_sub_adata.obs['mask'].isna())].copy()

    sub_adata.X = sub_adata.layers['counts']
    _stats.append(sub_adata.shape[0]) # For stats: n.cells
    obs_qc, var_qc = sc.pp.calculate_qc_metrics(sub_adata, percent_top=None)
    print(obs_qc, var_qc)
    _stats.append(obs_qc['total_counts'].mean()) # for stats: med.trx.p.cell
    _stats.append(obs_qc['total_counts'].median()) # for stats: med.trx.p.cell
    _stats.append(obs_qc['n_genes_by_counts'].median()) # For stats: med.genes.p.cell


    # Bulk correlation
    plt.clf()
    plt.style.use('default')
    bcorr, fig = Correlation(sub_adata, bulkref, ax1=reg, ax2="bulk", highlight=CONTROL_GENES, logscale=True)
    _stats.append(bcorr) # blk.corr
    fig.savefig(OUTPATH_DIR / f'{reg}_5bulkcorr.png', bbox_inches='tight')
    plt.clf()

    plt.style.use('default')
    # Top markers
    nums = sub_adata.obs.groupby(by=md.CTYPE_KEY).size()
    drop = nums[nums < 2].index
    temp = sub_adata[~sub_adata.obs[md.CTYPE_KEY].isin(drop)]
    sc.tl.rank_genes_groups(temp, groupby=md.CTYPE_KEY, n_genes=3)
    fig = sc.pl.rank_genes_groups_dotplot(temp, standard_scale='var', return_fig=True)
    fig.savefig(OUTPATH_DIR / f'{reg}_1.1mark.png', bbox_inches='tight')

    # Minimal marker table
    sc.pl.dotplot(temp, groupby=md.CTYPE_KEY, var_names=GENES_OF_INT, standard_scale='var', dendrogram=True, 
            title=f'{reg} marker genes', return_fig=True) 
    fig.savefig(OUTPATH_DIR / f'{reg}_1.2mark.png', bbox_inches='tight')

    # Cell counts
    fig, ax = plt.subplots(1, 1, figsize=(5, 7))
    a = sub_adata.copy()# [pu_adata.obs['batch'] == batch]
    counts_df = pd.DataFrame(a.obs.groupby(by=md.CTYPE_KEY).size(), columns=['counts']).reset_index()
    # print(counts_df)
    abs_values = counts_df['counts']
    tot = sum(abs_values.astype('int'))
    rel_values = abs_values.apply(lambda x: x/tot) * 100
    lbls = [f' {p[0]} ({p[1]:.0f}%)' for p in zip(abs_values, rel_values)]
    ax = sns.barplot(counts_df,
                x='counts', y=md.CTYPE_KEY, gap=.05, ax=ax)
    ax.bar_label(container=ax.containers[0], labels=lbls)
    ax.set_title(f"{reg} counts")
    ax.get_figure().savefig(OUTPATH_DIR / f'{reg}_2.1cellcts.png', bbox_inches='tight')

    # Cell counts (compare to ref)
    fig, ax1 = plt.subplots(1, 1, figsize=(5, 7))
    counts_df.set_index(md.CTYPE_KEY, inplace=True)
    counts_df['ref_counts'] = a.uns['reference_counts']['reference_counts']
    counts_df['norm_counts'] = counts_df['counts'] / counts_df['counts'].sum()
    counts_df['norm_ref_counts'] = counts_df['ref_counts'] / counts_df['ref_counts'].sum()
    counts_df = counts_df.drop(['counts', 'ref_counts'], axis=1)
    counts_df = counts_df.stack().reset_index()
    counts_df.columns = [md.CTYPE_KEY, 'source', 'counts']
    # print(counts_df)
    # print(counts_df)
    # abs_values = counts_df['counts']
    # tot = sum(abs_values.astype('int'))
    # rel_values = abs_values.apply(lambda x: x/tot) * 100
    # lbls = [f' {p[0]} ({p[1]:.0f}%)' for p in zip(abs_values, rel_values)]
    ax = sns.catplot(counts_df, x='counts', y='_CELLTYPE', hue='source', kind='bar')
    plt.savefig(OUTPATH_DIR / f'{reg}_2.2cellcts.png', bbox_inches='tight')
    plt.clf()

    # Cell counts (of reference)
    fig, ax = plt.subplots(1, 1, figsize=(5, 7))
    a = sub_adata.copy()# [pu_adata.obs['batch'] == batch]
    counts_df = a.uns['reference_counts'].reset_index()
    counts_df.columns = [md.CTYPE_KEY, 'reference_counts']
    abs_values = counts_df['reference_counts']
    tot = sum(abs_values.astype('int'))
    rel_values = abs_values.apply(lambda x: x/tot) * 100
    lbls = [f' {p[0]} ({p[1]:.0f}%)' for p in zip(abs_values, rel_values)]
    ax = sns.barplot(counts_df,
                x='reference_counts', y=md.CTYPE_KEY, gap=.05, ax=ax, color='tab:orange')
    ax.bar_label(container=ax.containers[0], labels=lbls)
    ax.set_title("Reference counts")
    ax.get_figure().savefig(OUTPATH_DIR / f'{reg}_2.3cellcts.png', bbox_inches='tight')
    plt.clf()


    # spatial
    plt.style.use('dark_background')
    # fig = sc.pl.embedding(sub_adata, basis='spatial', color=md.CTYPE_KEY, return_fig=True, color_map='Set1')
    # fig.savefig(OUTPATH_DIR / f'{reg}_3spat.png')
    if isinstance(HIGHLIGHT, dict):
        H_PALETTE= HIGHLIGHT
        HIGHLIGHT=HIGHLIGHT.keys()
    fig = embedding_highlight(sub_adata, highlight=HIGHLIGHT, basis='spatial', color=SPATCELLTYPE_KEY,
                        scale=1, dpi=1000, na_color=(.30, .30, .30),
                        kwargs={'palette':H_PALETTE})
    fig.savefig(OUTPATH_DIR / f'{reg}_3.1spat.png', bbox_inches='tight') # bbox_extra_artists=(lgd,text),
    plt.style.use('dark_background')
    # fig = sc.pl.embedding(sub_adata, basis='spatial', color=md.CTYPE_KEY, return_fig=True, color_map='Set1')
    # fig.savefig(OUTPATH_DIR / f'{reg}_3spat.png')
    if isinstance(HIGHLIGHT, dict):
        H_PALETTE= HIGHLIGHT
        HIGHLIGHT=HIGHLIGHT.keys()
    fig = sc.pl.embedding(sub_adata, basis='spatial', color=md.CTYPE_KEY, return_fig=True)
    fig.savefig(OUTPATH_DIR / f'{reg}_3.2spat.png', bbox_inches='tight') # bbox_extra_artists=(lgd,text),

    # umap
    sub_adata.obsm.pop('X_pca')
    sc.pp.pca(sub_adata)
    sc.pp.neighbors(sub_adata, n_neighbors=30, metric='euclidean')
    sc.pp.umap()
    fig = sc.pl.embedding(sub_adata, basis='umap', color=md.CTYPE_KEY, return_fig=True)
    fig.savefig(OUTPATH_DIR / f'{reg}_4umap.png', bbox_inches='tight')
    plt.clf()
    # plt.rcParams.update(plt.rcParamsDefault)

    # cell-type correlations
    matplotlib.rcParams.update(matplotlib.rcParamsDefault)
    anconf = order_snakefood('A4_HarmAnnotation')
    ref_h5ad_path = get_ref_path(sub_adata)
    adatas = []
    for path in ref_h5ad_path:
        print(f'{path} Loading...')
        _adata = sc.read_h5ad(path)
        adatas.append(_adata)
    refdata = sc.concat(adatas)
    if refdata.X.min() > -.01 and refdata.X.max() < 13:
        refdata.X = np.exp(refdata.X.toarray()) - 1
    ax = celltype_corr(sub_adata, refdata, (md.CTYPE_KEY, spatial_key))
    ax.get_figure().savefig(OUTPATH_DIR / f'{reg}_5ctcorr.png', bbox_inches='tight')


    # Control genes spatial distribution and stats 
    plt.style.use('dark_background')
    fig = sc.pl.embedding(sub_adata, basis='spatial', color=CONTROL_GENES, return_fig=True, vmax=1)
    fig.savefig(OUTPATH_DIR / f'{reg}_6ctrlspat.png', bbox_inches='tight')
    plt.clf()

    _ctrlstats = []
    cont_adata = sub_adata[:, CONTROL_GENES]
    obs_qc, var_qc = sc.pp.calculate_qc_metrics(cont_adata, percent_top=None)

    _ctrlstats.append((obs_qc['total_counts'] > 1).sum()) # for stats: med.ctlx.p.cell
    _ctrlstats.append(obs_qc['total_counts'].median()) # for stats: med.ctlx.p.cell
    _ctrlstats.append(obs_qc['total_counts'].mean()) # for stats: avg.ctlx.p.cell
    _ctrlstats.append(obs_qc['n_genes_by_counts'].median()) # for stats: med.ctlg.p.cell  
    _ctrlstats.append(Correlation(sub_adata, bulkref, ax1=reg, ax2="bulk", highlight=CONTROL_GENES, logscale=True)[0]) # for stats: ctl.blk.corr  
    _ctrlstats.append(cont_adata.obs[obs_qc['total_counts'] > 0].groupby(by='fov').size().std())

    norm_cont_adata = sc.pp.normalize_total(sub_adata, copy=True)
    ctgrp = sc.get.obs_df(norm_cont_adata, keys=(CONTROL_GENES + [md.CTYPE_KEY]))
    ctgrp_avgexp = pd.DataFrame(ctgrp.groupby(by=md.CTYPE_KEY).mean())
    ctgrp_avgexp = ctgrp_avgexp.stack().reset_index()
    ctgrp_avgexp.columns = ['celltype', 'gene', 'avg.count']
    fig = sns.catplot(ctgrp_avgexp, x='gene', y='avg.count', hue='celltype', kind='bar', aspect=1.4)
    for ax in fig.axes.flat: ax.tick_params("x", rotation=90); ax.set_ylim(0, 1.5)
    fig.savefig(OUTPATH_DIR / f'{reg}_tempname.png', bbox_inches='tight')

    # Gene counts analysis, abs, detection
    plt.style.use('default')
    fig, ax = plt.subplots(1, 1, figsize=(20, 5))
    gene_ex_data = pd.DataFrame(sub_adata.X.sum(axis=0), index=sub_adata.var_names.astype(str), columns=['counts'])
    gene_ex_data['counts'] = np.log10(gene_ex_data['counts'].values)
    gene_ex_data['cat'] = "non-control"
    gene_ex_data.loc[CONTROL_GENES, 'cat'] = "control"
    palette = {"non-control":'tab:blue', "control":'tab:red'}
    gene_ex_data['xticks'] = gene_ex_data.index
    print(gene_ex_data)
    ctrl_gene_index = gene_ex_data.index.difference(CONTROL_GENES)
    gene_ex_data.loc[ctrl_gene_index, 'xticks'] = ""
    gene_ex_data.sort_values(by='counts', ascending=False, inplace=True)
    gene_ex_data['rank'] = range(1, gene_ex_data.shape[0]+1)
    _ctrlstats.append(gene_ex_data.loc[ctrl_gene_index, 'rank'].mean())
    gene_ex_data.reset_index(inplace=True, names='genes')
    ax = sns.barplot(gene_ex_data, x='genes', y='counts', hue='cat',
                log_scale=False, ax=ax)#, palette=palette, aspect=3, )
    ax.set_xticklabels(gene_ex_data['xticks'])
    ax.tick_params("x", rotation=90)
    ax.get_figure().savefig(OUTPATH_DIR / f'{reg}_7absgene.png', bbox_inches='tight')
    plt.clf()

    stat_df.loc[reg] = _stats
    print(ctrlstat_df.columns)
    print(_ctrlstats)
    ctrlstat_df.loc[reg] = _ctrlstats
    plt.close("all")

    stat_df.to_csv(OUTPATH, sep='\t')
    ctrlstat_df.to_csv(OUTPATH.parent / 'ctrlstats.csv', sep='\t')

if __name__ == "__main__":
    l = ['test', 'this', 'list']
    print(f"Experiment.name in '{expnames!s}'")
    