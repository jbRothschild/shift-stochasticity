import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
from fig_variables import stochasticity, K, variability
import fig_variables as fv

from matplotlib.colors import LogNorm #color in log the maps

plt.style.use('parameters.mplstyle')  # particularIMporting

mte = []
"""
for i, delta in enumerate(variability):
    mte.append(ss.mteSum1D(cap[K], stochasticity, cap[K], ss.maxPop(cap[K], stochasticity, delta), delta, "sum1d"))
np.save("../data/heat_MTE_K100_log.npy", mte)
"""
MTE2 = np.load("../data/heat_MTE_K100_log.npy")
if np.isinf(np.log10((np.asarray(MTE2)).min())) or np.isnan(np.log10((np.asarray(MTE2)).min())):
    minimum = 0
else:
    minimum = int(np.log10((np.asarray(MTE2)).min()))
if np.isinf(np.log10((np.asarray(MTE2)).max())) or np.isnan(np.log10((np.asarray(MTE2)).max())):
    maximum = int(np.log10(np.finfo(np.float64).max))
else:
    maximum = int(np.log10((np.asarray(MTE2)).max()))

fv.plot_heatmap(np.asarray(MTE2), stochasticity, variability, 'Figure2', 'mean time to extinction')
