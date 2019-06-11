import numpy as np
import sys, os
sys.path.append(os.path.abspath('..'))
import mte.shiftStochasticity as ss
from fig_variables import stochasticity, K, variability, stochasticity_i, variability_i, cap
from fig_variables import *
import fig_variables as fv
import matplotlib.pyplot as plt
plt.style.use('parameters.mplstyle')  # particularIMporting

from matplotlib.ticker import LogLocator #log axes


maximum = ss.maxPop(cap[K], stochasticity[sto], variability[var])
population = np.linspace(1,maximum,maximum)

PDF = []
#ORDER IS VERY IMPORTANT HERE
#FP QSD
PDF.append(ss.pdfFP_full_normalized(population, stochasticity[sto], cap[K], variability[var]))
#FP GAUSSIAN
PDF.append(ss.pdfFP_gaussian(population, stochasticity[sto], cap[K], variability[var]))
#WKB REALSPACE
PDF.append(ss.distributionWKB_RS(population, stochasticity[sto], cap[K], variability[var]))
#QSD ALGORITHM
PDF.append(ss.statDistributionAlgo(population, stochasticity[sto], cap[K], variability[var]))



figname = 'Fig4'
plt.figure()
ax = plt.gca()

for i in range(0,len(PDF)):
    ax.plot(population, PDF[i], color=colors_techniques[i], label = techniques[i], linestyle=lines[i], lw=3)

ax.set_xlabel("Population")
ax.set_ylabel("Probability density function")
ax.set_yscale("log")
ax.set_ylim(10**(-20), 10**(0))
ax.get_yaxis().set_major_locator(LogLocator(numticks=5))
ax.minorticks_off()
ax.set_xlim(1, 200)
ax.legend(loc = 'lower center')

plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')
plt.close()
