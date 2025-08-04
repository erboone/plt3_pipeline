from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _bulk_correlation(means, bulkref, name, ax:plt.Axes, highlight=None):
    print(means)
    print(bulkref)
    common = means.index.astype(str).intersection(bulkref.index.astype(str))
    counts1 = means.loc[common]
    counts2 = bulkref.loc[common]
    corr = pearsonr(np.log1p(counts1), np.log1p(counts2))
    c_ser = pd.Series(data='tab:blue', index=counts1.var_names)
    c_ser.loc[highlight] = 'tab:red'

    ax.scatter(counts1, counts2, s=1, c=c_ser.values)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(f"{name} r={corr[0]:.03f}")
    return ax, corr.statistic
 