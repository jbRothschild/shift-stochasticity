import numpy as np
import sys, os
sys.path.append(os.path.abspath('..'))
import mte.shiftStochasticity as ss
from fig_variables import stochasticity, K, variability
import fig_variables as fv
from fig_variables import *
import matplotlib.pyplot as plt
plt.style.use('parameters.mplstyle')  # particularIMporting
import seaborn as sns #niceeee colors
from matplotlib.ticker import LogLocator #log axes

#__________________________A____________________________________

figname = 'Fig1A'
plt.figure()
ax = plt.gca()

pdf = []
"""
for i, stochas in enumerate(stochasticity_i):
    pdf.append(ss.statDistributionAlgo(ss.maxPop(cap[K], stochas, variability[var]), stochas, cap[K], variability[var], 10))
    print("1A")
np.save("../data/pdf_stoch.npy", pdf)
"""
PDF = np.load("../data/pdf_stoch.npy")

for i, stochas in enumerate(stochasticity_i):
    ax.plot(range(1,len(PDF[i])+1), PDF[i], color=colors_gradient2[i])

ax.set_xlabel("Population")
ax.set_ylabel("Probability density function")
ax.set_yscale("log")
ax.set_ylim(10**(-30), 10**(0))
ax.get_yaxis().set_major_locator(LogLocator(numticks=5))
ax.minorticks_off()
ax.set_xlim(1, 200)

plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')
plt.close()


#__________________________B____________________________________

figname = 'Fig1B'
plt.figure()
ax = plt.gca()
pdf = []
"""
for i, varia in enumerate(variability_i):
    pdf.append(ss.statDistributionAlgo(ss.maxPop(cap[K], stochasticity[sto], varia), stochasticity[sto], cap[K], varia, 10))
    print("1B")
np.save("../data/pdf_varia.npy", pdf)
"""

PDF = np.load("../data/pdf_varia.npy")

for i, varia in enumerate(variability_i):
    ax.plot(range(1,len(PDF[i])+1), PDF[i], color=colors_gradient2[i])

ax.set_xlabel("Population")
ax.set_ylabel("Probability density function")
ax.set_yscale("log")
ax.set_ylim(10**(-30), 10**(0))
ax.get_yaxis().set_major_locator(LogLocator(numticks=5))
ax.minorticks_off()
ax.set_xlim(1, 200)

plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')
plt.close()
