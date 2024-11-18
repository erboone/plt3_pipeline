#!/usr/bin/env python

import sys

import pandas as pd
from datadispatch.orm import Experiment
from mftools.mfexperiment import MerscopeExperiment

ROOT = sys.argv[1]
NAME = sys.argv[2]

e = MerscopeExperiment(ROOT, NAME)
mfadata = e.create_scanpy_object()

# TODO: improve loading of reference data
ref = pd.read_csv('/mnt/merfish18/BICAN/Reference_Data/bulkRNA/gn_reads.csv')
