import scanpy as sc
import pandas as pd
import anndata as an
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from datadispatch.access import select

def _1_Annotation(input, output):
    print('\n --- Beginning Annotation --- ')
    input_h5ad_files = list(input)
    print("Annotating") 

    print(type(output), output)

    
if __name__ == "__main__":
    l = ['test', 'this', 'list']
    print(f"Experiment.name in '{expnames!s}'")
    