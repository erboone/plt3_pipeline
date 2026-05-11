    
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

def _4_ControlQC(input, output, hashes, commit)
    exp_name = input
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
