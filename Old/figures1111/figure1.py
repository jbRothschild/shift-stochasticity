import latex_plots as lp
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns #niceeee colors
import shiftStochasticity as ss #All our models
from matplotlib.colors import LogNorm #color in log the maps
from matplotlib.ticker import LogLocator #log axes
from mpl_toolkits.axes_grid.inset_locator import inset_axes #axes within axes. How deep can we go????
import matplotlib.patches as patches #curved arrows

from fig_variables import *

cwd = os.getcwd()

#__________________________A____________________________________

fig, ax  = lp.newfig(0.6)
"""
pdf = []
for i, stochas in enumerate(stochasticity_i):
    pdf.append(ss.statDistributionAlgo(cap[K], stochas, cap[K], MAX, variability[var], 10))
    print "1A"
np.save(cwd + "/Data/pdf_stoch.npy", pdf)
"""
PDF = np.load(cwd + "/Data/pdf_stoch.npy")

for i, stochas in enumerate(stochasticity_i):
    ax.plot(range(1,MAX+1), PDF[i], color=colors_gradient[i])

ax.set_xlabel("Population")
ax.set_ylabel("Probability density function")
ax.set_yscale("log")
ax.set_ylim(10**(-30), 10**(0))
ax.get_yaxis().set_major_locator(LogLocator(numticks=5))
ax.minorticks_off()
ax.set_xlim(1, 200)

lp.savefig("Figure1-A")
plt.close(fig)

#__________________________B____________________________________

fig, ax  = lp.newfig(0.6)
"""
pdf = []
for i, varia in enumerate(variability_i):
    pdf.append(ss.statDistributionAlgo(cap[K], stochasticity[sto], cap[K], MAX, varia, 10))
    print "1B"
np.save(cwd + "/Data/pdf_varia.npy", pdf)
"""
PDF = np.load(cwd + "/Data/pdf_varia.npy")

for i, varia in enumerate(variability_i):
    ax.plot(range(1,MAX+1), PDF[i], color=colors_gradient[i])

ax.set_xlabel("Population")
ax.set_ylabel("Probability density function")
ax.set_yscale("log")
ax.set_ylim(10**(-30), 10**(0))
ax.get_yaxis().set_major_locator(LogLocator(numticks=5))
ax.minorticks_off()
ax.set_xlim(1, 200)

lp.savefig("Figure1-B")
plt.close(fig)
