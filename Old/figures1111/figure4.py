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

maximum = ss.maxPop(cap[K], stochasticity[sto], variability[var])
population = np.linspace(1,maximum,maximum)

PDF = []
#ORDER IS VERY IMPORTANT HERE
#FP QSD
PDF.append(ss.pdfFP_full_normalized(population, stochasticity[sto], cap[K], variability[var]))
#FP GAUSSIAN
PDF.append(ss.pdfFP_gaussian(population, stochasticity[sto], cap[K], variability[var]))
#WKB REALSPACE
PDF.append(ss.distributionWKB_RS(population, cap[K], stochasticity[sto], cap[K], maximum, variability[var]))
#QSD ALGORITHM
PDF.append(ss.statDistributionAlgo(cap[K], stochasticity[sto], cap[K], maximum, variability[var], 10))

for i in range(0,len(PDF)):
    ax.plot(population, PDF[i], color=colors_techniques[i], label = techniques[i], linestyle=lines[i])

ax.set_xlabel("Population")
ax.set_ylabel("Probability density function")
ax.set_yscale("log")
ax.set_ylim(10**(-20), 10**(0))
ax.get_yaxis().set_major_locator(LogLocator(numticks=5))
ax.minorticks_off()
ax.set_xlim(1, 200)
ax.legend(loc = 'lower center')

lp.savefig("Figure4")
plt.close(fig)
