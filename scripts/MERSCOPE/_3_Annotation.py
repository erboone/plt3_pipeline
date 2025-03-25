import warnings as warn
from pathlib import Path
from datetime import datetime
import scanpy as sc
import pandas as pd
import anndata as an
import numpy as np
import scrublet as scr
from sklearn.neighbors import KNeighborsClassifier
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from datadispatch.access import select
from mftools.scanpy_wrap import MerData as md
from mftools.mfexperiment import MerscopeExperiment

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
    STEPNAME = 'A4_HarmAnnotation'
    anconf = order_snakefood(STEPNAME)
    filtconf = order_snakefood('Filter')
    CELLTYPE_KEY = str(anconf['celltype_key'])
    MIN_COUNTS = int(filtconf['min_counts'])
    MIN_GENES = int(filtconf['min_genes'])
    OUTPATH = Path(str(output))
    # QC_PATH = Path()
    try:
        REMOVE_DOUBLETS = eval(filtconf['remove_doublets'])
    except:
        warn.warn("Cannot parse text in 'anconf['remove_doublets']' as a boolean; defaulting to True")
        print("Cannot parse text in 'anconf['remove_doublets']' as a boolean; defaulting to True")
        REMOVE_DOUBLETS = True

    ref_h5ad_path = str(anconf['ref_path']).split(',')
    print('\n --- Preparing AnnData objects --- ')
    input_h5ad_paths = input
    print(input)
    # MERSCOPE data
    adatas = []
    for path in input_h5ad_paths:
        adata = sc.read_h5ad(path)
        adatas.append(adata)
    merdata = sc.concat(adatas)
    print(merdata.X.min(), merdata.X.max())

    # Reference data
    adatas = []
    for path in ref_h5ad_path:
        adata = sc.read_h5ad(path)
        adatas.append(adata)
    refdata = sc.concat(adatas)
    clust_table = pd.read_csv(Path(path).parent.parent / 'media-2.csv').dropna(how='all')
    clust_table = clust_table[['Cluster', 'Supercluster', 'Class auto-annotation', 'Neurotransmitter auto-annotation']]
    clust_table['Cluster'] = clust_table['Cluster'].astype('int')
    clust_table.set_index('Cluster', inplace=True)
    refdata.obs = refdata.obs.join(clust_table, on='Clusters')

    # Subset both to common genes only
    common_genes = merdata.var_names.intersection(refdata.var_names)
    merdata = merdata[:, common_genes].copy()
    refdata = refdata[:, common_genes].copy()

    # normalize and scale data
    for desig, a in zip(['merFISH data', 'Reference data'],[merdata, refdata]):
        # print(f'\n{desig}:')
        # print(f'X is {a.X.dtype}\tmin:max::{a.X.min()}:{a.X.max()}')
        # print()

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
    for region in merdata.obs['region'].unique():

        # Combining merFISH and reference data setss
        adata = an.concat(
            [merdata[merdata.obs['region'] == region], refdata], 
            label=SOURCE_KEY, 
            keys=['merdata', 'refdata'],
            join='outer'
        )

        # Set up plots for integration QC
        fig, axs = plt.subplots(1, 3)
        # Preintegration umap
        sc.pp.pca(adata)
        sc.pp.neighbors(adata, n_neighbors=30, metric='euclidean')
        sc.tl.umap(adata, min_dist=.3)
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
        nn.fit(traindata.obsm["X_umap"], traindata.obs[CELLTYPE_KEY])
        pred = nn.predict(newadata.obsm["X_umap"])
        newadata.obs[md.CTYPE_KEY] = pred

        newadata.write(OUTPATH.parent / f'{region}_{datetime.now()}.h5ad')
        
        # Postintegration
        sc.pp.neighbors(newadata, n_neighbors=30, metric='correlation')
        sc.tl.umap(newadata, min_dist=0.1)
        sc.pl.embedding(newadata, basis='umap', color='region', ax=axs[2])
        plt.savefig('z_EMPTY.png')
    
    with open(OUTPATH, 'w') as f:
        f.write(F"FINISHED: {datetime.now()}\n\n")


B3_HarmAnnotation = A4_HarmAnnotation

    
if __name__ == "__main__":
    l = ['test', 'this', 'list']
    print(f"Experiment.name in '{expnames!s}'")
    