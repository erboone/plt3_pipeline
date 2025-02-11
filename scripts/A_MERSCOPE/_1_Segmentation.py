#!/usr/bin/env python
# coding: utf-8

import os, sys
from configparser import ConfigParser
from pathlib import Path
from itertools import chain

import numpy as np
from skimage.measure import regionprops # This has been moved to plotting
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.axes import Axes 
from cellpose import models
import pandas as pd

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

def _1_Cellpose(
        exp_name,
        input,
        outfile,
        hashes,
        commit,
        fov_range:tuple|None=None
        ) -> None:
    # ----- Load Config ----- #
    # TODO: change this to work with config files generated from merbot
    cpconf = order_snakefood('A11_Cellpose')

    # ----- Parse Args ------ #
    # TODO: upgrade this
    # EXPERIMENT_NAME = sys.argv[1]
    EXPERIMENT_NAME = exp_name
    # OUTPUT = sys.argv[2]
    OUTPUT = outfile
    PATHS = {
        'checkpoint': OUTPUT
    }

    print('SEARCH THIS STRING', EXPERIMENT_NAME, OUTPUT)
    # PATH PARAMS 
    # CHANGE ONLY IF NEEDED: these should be standardized in the future and creation of these dirs automated from
    # globalized config file
    mer_rawdata_dir = "data"
    mer_output_dir = "output"
    cellpose_dir = f"./cellpose_{EXPERIMENT_NAME}" # TODO: for testing, change later
    masks_dir = "masks" 

    # Script params
    ZSLICE = int(cpconf['z_slice'])
    CHANNEL = cpconf['channel']
    MODEL= ""
    MODEL_PATH = '/home/eboone/CellSegmentation/cp_models/CP_20240227_000004_staged'
    ALT_SAVE = cpconf['alt_path']
    TEST_FOVS = list(map(int, cpconf.get('test_fovs').split(',')))

    # Path assembly TODO: remove this when we can trust MERBOT to generate these for
    # image_dataset_path = f"{MERSCOPE_DIR}/{MER_RAWDATA_DIR}/{EXPERIMENT_NAME}/" # Used to load imageset
    # expiriment_out_path = f"{MERSCOPE_DIR}/{MER_OUTPUT_DIR}/{EXPERIMENT_NAME}/" # Used to find barcodes 
    # cellpose_out_path = f"./{CELLPOSE_DIR}/" # used to save final cellpose output
    # masks_out_path = f"{cellpose_out_path}{MASKS_DIR}/" # used to save masks

    # print(f"Looking for images in {image_dataset_path}")
    # print(f"Looking for barcodes in {expiriment_out_path}")
    # print(f"Writing segmenatation output {cellpose_out_path}; {masks_out_path}")


    # # Check directory structure
    # needed_dirs = [image_dataset_path, expiriment_out_path, cellpose_out_path, masks_out_path]
    # print(needed_dirs)
    # for path in needed_dirs:
    #     os.makedirs(path, exist_ok=True)

    #     if os.path.exists(path):
    #         if not os.access(path, os.W_OK) or not os.access(path, os.R_OK):
    #             raise RuntimeError(f"You do not have read/write permissions in directory \"{path}\"")
    #     else:
    #         raise RuntimeError(f"Attempted to create \"{path}\" but failed. Check that you have permission.")
        

    # Use file path like format ("/mnt/merfish12/MERSCOPE/merfish_raw_data/202401261424_20240126M134UWA7648CX22PuS6_VMSC10102
    result = select('Experiment', where=f"Experiment.name == '{EXPERIMENT_NAME}'")

    if len(result) != 1:
        raise RuntimeError('Experiment search critera did not find unique entry.') 


    
    if ALT_SAVE:
        alt_paths = {'cellpose': ALT_SAVE, 'masks': f"{ALT_SAVE}/masks"}
    else:
        alt_paths = {}

    experiment = MerscopeExperiment(result[0].root.path, EXPERIMENT_NAME, 
                alt_paths=alt_paths,
                seg_kwargs={
                    'zslice': ZSLICE, 
                    'channel': CHANNEL},
                img_kwargs={}
            )
    e = experiment

    fig:Figure; axs:list[list[Axes]]
    # fig, axs = plt.subplots(2, len(TEST_FOVS))

    # TODO: come back to this; implement plotting and saving a couple of test segmentations
    Path(e.files['cellpose']).parent.mkdir(parents=True, exist_ok=True)
    for i, fov in enumerate(TEST_FOVS):
        fig, ax = plt.subplots(1, 1)
        fov_show(seg=e.seg, imgs=e.imgs, fov=fov, show=False, ax=ax)
        # TODO: include implementation to create the qc directory if not exist
        fig.savefig(f"{e.files['cellpose']}/testseg_{fov}")

    from datetime import datetime
    print("start metadata:", datetime.now())
    # either load or create meatadata, then save if not already
    metadata = e.seg.metadata
    e.save_cell_metadata(metadata)
    print("end metadata:", datetime.now())


    # load barcodes, assign to cells, save to disk
    print("start barcodes:", datetime.now())
    barcodes = e.load_barcode_table()
    assign_to_cells(barcodes, e.seg)
    link_cell_ids(barcodes, e.seg.linked_cells)
    e.save_barcode_table(barcodes)
    print("end barcodes:", datetime.now())

    barcodes = e.load_barcode_table()
    cbgtab = create_cell_by_gene_table(barcodes)
    cbgtab.index = cbgtab.index.astype(int)
    e.save_cell_by_gene_table(cbgtab)
    # TODO: Create and write the scanpy object somewhere.

    with open(outfile, 'w') as f:
        f.write('FINISHED: ', datetime.now(), '\n\n')
        
        for k, p in PATHS.items():
            f.write(f'{k}: {p}')

def _2_SaveRawScanpy(
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
    alt_save = cpconf['alt_save']
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
    mdata.log = ('function', '_2_SaveRawScanpy')
    mdata.log = ('commit_id', commit)
    mdata.log = ('normalization', None)


    mdata.obs[md.BATCH_KEY] = pd.Series(res.meta.BICANExperimentID)
    mdata.obs[md.CTYPE_KEY] = pd.Series(None)

    # -------------------------- Save object -------------------------------- #

    mdata.write(Path(str(output)))

