import latex_plots as lp
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns #niceeee colors
import PopDyn_variation.shiftStochasticity as ss #All our models
from matplotlib.colors import LogNorm #color in log the maps
from matplotlib.ticker import LogLocator #log axes
from mpl_toolkits.axes_grid.inset_locator import inset_axes #axes within axes. How deep can we go????
import matplotlib.patches as patches #curved arrows

from fig_variables import *

cwd = os.getcwd()

fig, ax  = lp.newfig(0.6)
"""
mte = []
for i, delta in enumerate(variability):
    mte.append(ss.mteSum1D(cap[K], stochasticity, cap[K], ss.maxPop(cap[K], stochasticity, delta), delta, "sum1d"))
np.save(cwd + "/Data/heat_MTE_K100_log.npy", mte)
"""
MTE2 = np.load(cwd + "/Data/heat_MTE_K100_log.npy")
if np.isinf(np.log10((np.asarray(MTE2)).min())) or np.isnan(np.log10((np.asarray(MTE2)).min())):
    minimum = 0
else:
    minimum = int(np.log10((np.asarray(MTE2)).min()))
if np.isinf(np.log10((np.asarray(MTE2)).max())) or np.isnan(np.log10((np.asarray(MTE2)).max())):
    maximum = int(np.log10(np.finfo(np.float64).max))
else:
    maximum = int(np.log10((np.asarray(MTE2)).max()))

cax = ax.contourf(stochasticity, variability, np.asarray(MTE2), cmap=plt.cm.inferno, norm=LogNorm(), levels=np.logspace(minimum, maximum, maximum))
cbar = fig.colorbar(cax, ticks=[10**minimum, 10**int((maximum-minimum)/3), 10**int((maximum-minimum)*2/3), 10**maximum])
for c in ax.collections:
    c.set_edgecolor("face")
cbar.ax.set_ylabel(r'$\tau_e$')
ax.set_xlabel(r"$q$")
ax.set_ylabel(r"$\delta$")
ax.set_yscale("log")
ax.minorticks_off()
lp.savefig("Figure2")
plt.close(fig)
