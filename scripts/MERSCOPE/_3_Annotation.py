import warnings as warn
from pathlib import Path
from datetime import datetime
import os

import scanpy as sc
import pandas as pd
import anndata as an
import numpy as np
import scrublet as scr
from sklearn.neighbors import KNeighborsClassifier
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import seaborn as sns

from datadispatch.access import select
from mftools.scanpy_wrap import MerData as md
from mftools.mfexperiment import MerscopeExperiment
from mftools.plotting import celltype_corr

from ..meta import order_snakefood
from ..meta import OUTPUT

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
    os.makedirs(OUTPATH.parent, exist_ok=True)

    # QC_PATH = Path()
    try:
        REMOVE_DOUBLETS = filtconf['remove_doublets']
    except:
        warn.warn("Cannot parse text in 'anconf['remove_doublets']' as a boolean; defaulting to True")
        print("Cannot parse text in 'anconf['remove_doublets']' as a boolean; defaulting to True")
        REMOVE_DOUBLETS = True

    ref_h5ad_path = anconf['ref_path']
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

        np.nan_to_num(adata.X, nan=0, copy=False)

        # Set up plots for integration QC
        fig, axs = plt.subplots(1, 3, figsize=(9, 3))
        print("Preintegration UMAP", end=' ')
        # Preintegration umap
        print("PCA", end=' ')
        sc.pp.pca(adata)
        print("->Neighbors", end=' ')
        sc.pp.neighbors(adata, n_neighbors=30, metric='euclidean')
        print("->UMAP", end=' ')
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
    
    recomb.uns['reference_counts'] = pd.DataFrame(a.obs.groupby(by=CELLTYPE_KEY[0]).size(), columns=['reference_counts'])
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
    if qcconf["spatial_ctkey"] is not None:
        SPATCELLTYPE_KEY = qcconf['spatial_ctkey']
    else:
        SPATCELLTYPE_KEY = md.CTYPE_KEY
    OUTPATH_DIR = OUTPATH.parent
    os.makedirs(OUTPATH_DIR, exist_ok=True)

    adata:an.AnnData = sc.read_h5ad(INPATH)
    adata.obs[md.CTYPE_KEY] = adata.obs[md.CTYPE_KEY].astype('category')
    sc.pl.embedding(adata, basis='umap', color=md.CTYPE_KEY, show=False)
    stat_df = pd.DataFrame(columns=['sample', 'n.cells', 'med.trx.p.cell', 'med.genes.p.cell'])
    stat_df.set_index('sample', inplace=True)

    for reg in adata.obs['region'].unique():
        sub_adata = adata[adata.obs['region'] == reg]

        sub_adata.X = sub_adata.layers['counts']
        n_cell = sub_adata.shape[0]
        obs_qc, var_qc = sc.pp.calculate_qc_metrics(sub_adata, percent_top=None)
        print(obs_qc, var_qc)
        mtpc = obs_qc['total_counts'].median()
        obs_qc.to_csv('test.csv')
        print()
        print("Celltype correlation")

        mgpc = obs_qc['n_genes_by_counts'].median()

        stat_df.loc[reg] = [n_cell, mtpc, mgpc]

        plt.style.use('default')
        # Top markers
        nums = sub_adata.obs.groupby(by=md.CTYPE_KEY).size()
        drop = nums[nums < 2].index
        temp = sub_adata[~sub_adata.obs[md.CTYPE_KEY].isin(drop)]
        sc.tl.rank_genes_groups(temp, groupby=md.CTYPE_KEY, n_genes=3)
        fig = sc.pl.rank_genes_groups_dotplot(temp, standard_scale='var', return_fig=True)
        fig.savefig(OUTPATH_DIR / f'{reg}_1.1mark.png', bbox_inches='tight')
        
        # Minimal marker table
        sc.pl.dotplot(temp, groupby=md.CTYPE_KEY, var_names=GENES_OF_INT, standard_scale='var', dendrogram=True, return_fig=True) 
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
        ax.get_figure().savefig(OUTPATH_DIR / f'{reg}_2cellcts.png', bbox_inches='tight')

        # Cell counts (compare to ref)
        fig, ax1 = plt.subplots(1, 1, figsize=(5, 7))
        ax2 = ax1.twiny()
        counts_df.set_index(md.CTYPE_KEY, inplace=True)
        counts_df['ref_counts'] = a.uns['reference_counts']['reference_counts']
        tempcopy = counts_df.copy()
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
        sns.catplot(counts_df, x='counts', y='_CELLTYPE', hue='source', kind='bar')
        plt.savefig(OUTPATH_DIR / f'{reg}_2.2cellcts.png', bbox_inches='tight')
        # fig = sns.catplot(counts_df, kind='bar',
        #                     x='counts', y=CTKEY, hue='source', gap=.05, ax=ax, log_scale=True)
        # ax1.get_figure().savefig()
        plt.clf()
        tempcopy.to_csv(OUTPATH_DIR / f'{reg}_2.2cellcts.csv', sep="\t")
        test = pd.DataFrame((tempcopy['counts'] / tempcopy['counts'].sum()) * 100)
        test = test.join((tempcopy['ref_counts'] / tempcopy['ref_counts'].sum()) * 100)
        ax = sns.heatmap(test)
        # ax.title.set_text(f"Percent of sample")
        # ax.get_figure().savefig(f"/home/erboone/HubMap_multOrg/plotdump/{reg}_2.3cellcts.png", bbox_inches='tight')


        # spatial
        plt.style.use('dark_background')
        # fig = sc.pl.embedding(sub_adata, basis='spatial', color=md.CTYPE_KEY, return_fig=True, color_map='Set1')
        # fig.savefig(OUTPATH_DIR / f'{reg}_3spat.png')
        fig = embedding_highlight(sub_adata, highlight=HIGHLIGHT, basis='spatial', color=SPATCELLTYPE_KEY, scale=1, face='black', dpi=1000, na_color=(.30, .30, .30))
        fig.savefig(OUTPATH_DIR / f'{reg}_3spat.png', bbox_inches='tight') # bbox_extra_artists=(lgd,text),

        # umap
        fig = sc.pl.embedding(sub_adata, basis='umap', color=md.CTYPE_KEY, return_fig=True)
        fig.savefig(OUTPATH_DIR / f'{reg}_4umap.png', bbox_inches='tight')
        # plt.rcParams.update(plt.rcParamsDefault)

        # cell-type correlations
        # STEPNAME = 'A4_HarmAnnotation'
        # anconf = order_snakefood(STEPNAME)
        # ref_h5ad_path = anconf['ref_path']
        # adatas = []
        # for path in ref_h5ad_path:
        #     print(f'{path} Loading...')
        #     adata = sc.read_h5ad(path)
        #     adatas.append(adata)
        # refdata = sc.concat(adatas)
        # fig = celltype_corr(adata, refdata, md.CTYPE_KEY)
        # fig.savefig(OUTPATH_DIR / f'{reg}_5ctcorr.png', bbox_inches='tight')


    stat_df.to_csv(OUTPATH, sep='\t')


if __name__ == "__main__":
    l = ['test', 'this', 'list']
    print(f"Experiment.name in '{expnames!s}'")
    