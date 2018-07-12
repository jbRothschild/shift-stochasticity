import latex_plots as lp #Necessary for the latex default formatting
from fig_variables import *

from matplotlib.ticker import LogLocator #log axes

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
