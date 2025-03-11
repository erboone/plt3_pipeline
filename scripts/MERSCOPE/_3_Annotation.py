import warnings as warn

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

from ..meta import order_snakefood
from ..meta import OUTPUT

SOURCE_KEY = 'source'

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
    MIN_COUNTS = int(filtconf['min_counts'])
    MIN_GENES = int(filtconf['min_genes'])
    # QC_PATH = Path()
    try:
        REMOVE_DOUBLETS = eval(anconf['remove_doublets'])
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
    nonnormalized_merdata = merdata.copy()
    # Reference data
    adatas = []
    for path in ref_h5ad_path:
        adata = sc.read_h5ad(path)
        adatas.append(adata)
    refdata = sc.concat(adatas)

    for desig, a in zip(['merFISH data', 'Reference data'],[merdata, refdata]):
        print(f'\n{desig}:')
        print(f'X is {a.X.dtype}\tmin:max::{a.X.min()}:{a.X.max()}')
        print()

        # Basic filtering
        ncells = a.shape[0]
        sc.pp.filter_cells(a, min_counts=MIN_COUNTS)
        diff = ncells - a.shape[0]; print(f"{diff}({diff/ncells:.2f})cells removed by min counts"); ncells = a.shape[0]
        sc.pp.filter_cells(a, min_genes=MIN_GENES)
        diff = ncells - a.shape[0]; print(f"{diff}({diff/ncells:.2f})cells removed by min genes"); ncells = a.shape[0]

        # Scrublet
        print(type(REMOVE_DOUBLETS), int(REMOVE_DOUBLETS))
        if REMOVE_DOUBLETS:
            scrub = scr.Scrublet(a.X)
            dblt_score, dblt_pred = scrub.scrub_doublets(log_transform=True)
            a = a[dblt_pred]
            diff = ncells - a.shape[0]; print(f"{diff}({diff/ncells:.2f})cells removed by scrublet"); ncells = a.shape[0]

        # Volume
        if desig == 'merFISH data':
            a = a[
                (a.obs['volume'] > 100) &
                (a.obs['volume'] < (a.obs['volume'].median() * 3))
            ]
            diff = ncells - a.shape[0]; print(f"{diff}({diff/ncells:.2f})cells removed by volume"); ncells = a.shape[0]


        # Transformation
        sc.pp.normalize_total(a)
        sc.pp.log1p(a)

    # Combining merFISH and reference data setss
    print()
    adata = an.concat([merdata, refdata], label=SOURCE_KEY, keys=['merdata', 'refdata'])
    fig, axs = plt.subplots(1, 2)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata, n_neighbors=30, metric='correlation')
    sc.tl.umap(adata, min_dist=0.1)
    sc.pl.embedding(adata, basis='umap', color=SOURCE_KEY, ax=axs[0])
    
    print('\n --- Harmony Integration --- ')
    sc.external.pp.harmony_integrate(adata, key=SOURCE_KEY, basis='X_pca', adjusted_basis='X_pca_harmony', max_iter_harmony=50)
    sc.pp.neighbors(adata, n_neighbors=30, metric='correlation')
    sc.tl.umap(adata, min_dist=0.1)
    sc.pl.embedding(adata, basis='umap', color=SOURCE_KEY, ax=axs[1])
    plt.savefig(f'_output/quality_control/{output}.png')

    traindata = adata[adata.obs[SOURCE_KEY] == "refdata"] 
    nn = KNeighborsClassifier(n_jobs=16)
    nn.fit(traindata.obsm["X_umap_harmony"], traindata.obs["seurat_clusters"])
    pred = nn.predict(adata[adata.obs[SOURCE_KEY] == "0"].obsm["X_umap_harmony"])
    nonnormalized_merdata.obs.loc[adata.obs[SOURCE_KEY] == "0", md.CTYPE_KEY] = pred

    nonnormalized_merdata.write('test.h5ad')

    
if __name__ == "__main__":
    l = ['test', 'this', 'list']
    print(f"Experiment.name in '{expnames!s}'")
    