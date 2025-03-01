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

def A1_Cellpose(
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
    CHANNEL = cpconf['channel'].split(',')
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

    # TODO: Load parameters into class easily
    # with open(f"{e.files['settings']}/microscope_parameters.json") as f:
    #     microscope_params = json.load(f)

    fig:Figure; axs:list[list[Axes]]
    # fig, axs = plt.subplots(2, len(TEST_FOVS))
    
    for i, fov in enumerate(TEST_FOVS):
        fig, ax = plt.subplots(1, 1)
        fov_show(seg=e.seg, imgs=e.imgs, fov=fov, show=False, ax=ax)
        # TODO: include implementation to create the qc directory if not exist
        fig.savefig(f"{e.files['cellpose']}/testseg_{fov}")

    from datetime import datetime
    print("start metadata:", datetime.now())
    # either load or create meatadata, then save if not already
    metadata:pd.DataFrame = e.seg.metadata
    if 'region' not in metadata.columns: # TODO: move this into the logic of 'create_metadata'
        old_metadata = pd.read_csv(Path(e.files['output']) / 'cell_metadata.csv')
        fov_to_reg = old_metadata['fov', 'region'].reset_index(drop=True).reindex('fov').drop_duplicates()
        metadata = metadata.join(fov_to_reg, on='fov')
    e.save_cell_metadata(metadata)
    print("end metadata:", datetime.now())

    # load barcodes, assign to cells, save to disk
    print("start barcodes:", datetime.now())
    # TODO: The fundimental issue here is that we explicitly want to ask for the cannontical barcodes
    # that come off the machine, then modify. We do not have the ability to do that here easily
    cannon_barcodes = e.load_barcode_table()
    new_barcodes = assign_to_cells(cannon_barcodes, e.seg, e.settings['image_dimensions'][0])
                                #    flip_x=microscope_params['flip_horizontal'],
                                #    flip_y=microscope_params['flip_vertical'])
    link_cell_ids(new_barcodes, e.seg.linked_cells)
    e.save_barcode_table(new_barcodes)
    print("end barcodes:", datetime.now())

    cbgtab = create_cell_by_gene_table(new_barcodes)
    cbgtab.index = cbgtab.index.astype(int)
    e.save_cell_by_gene_table(cbgtab)
    # TODO: Create and write the scanpy object somewhere.

    # Replot test FOV's with transcripts assigned
    for i, fov in enumerate(TEST_FOVS):
        fig, ax = plt.subplots(1, 1)
        fov_show(seg=e.seg, imgs=e.imgs, fov=fov, show=False, ax=ax, plot_transcripts=True,
                 output=MerfishAnalysis(e.files['cellpose']))
        # TODO: include implementation to create the qc directory if not exist
        fig.savefig(f"{e.files['cellpose']}/testseg_{fov}")
    
    with open(str(outfile), 'w') as f:
        f.write(F"FINISHED: {datetime.now()}\n\n")
        
        for k, p in PATHS.items():
            f.write(f'{k}: {p}')



