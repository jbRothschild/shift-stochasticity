import latex_plots as lp
from fig_variables import *
import sys
sys.path.insert(0, '../mte')
import comparaison as comp
cimport shiftStochasticity as ss

from matplotlib.colors import LogNorm #color in log the maps

exac_mean = np.load("../data/heat_exact_mean.npy")
exac_var = np.load("../data/heat_exact_var.npy")
exac_var[exac_var < 0] = 0.01 #funky stuff at mean and var
exac_mean[exac_mean < 0] = 0.01

exact_mean = exac_mean.transpose() #made a mistake and for now the saved file is weird. need to change.
exact_var = exac_var.transpose()

#FP QSD
#PDF.append(ss.pdfFP_full_normalized(population, stochasticity[sto], cap[K], variability[var]))
#FP GAUSSIAN
#PDF.append(ss.pdfFP_gaussian(population, stochasticity[sto], cap[K], variability[var]))
#WKB REALSPACE
#PDF.append(ss.distributionWKB_RS(population, cap[K], stochasticity[sto], cap[K], maximum, variability[var]))
#QSD ALGORITHM
#PDF.append(ss.statDistributionAlgo(population, stochasticity[sto], cap[K], variability[var])

#============================FP QSD=============================

FPQSD_mean, FPQSD_var = comp.mean_var(ss.pdfFP_full_normalized, stochasticity, variability, capacity)
np.save("../data/heat_FPQSD_mean.npy", FPQSD_mean)
np.save("../data/heat_FPQSD_var.npy", FPQSD_vsr)

FPQSD_mean = np.load("../data/heat_FPQSD_mean.npy")
FPQSD_vsr = np.load("../data/heat_FPQSD_vsr.npy")

FPQSD_vsr[FPQSD_vsr < 0] = 0.01 #Funky stuff happens, sometimes it's negative???
FPQSD_mean[FPQSD_mean < 0] = 0.01

#------------------------------ A ------------------------------

fig, ax  = lp.newfig(0.6)

cax = ax.contourf(stochasticity, variability, np.asarray(FPQSD_mean/exact_mean), 100, cmap=plt.cm.inferno)#, norm=LogNorm(), levels=np.logspace(minimum, maximum, maximum))
cbar = fig.colorbar(cax, ticks=[np.min(FPQSD_mean/exact_mean), 1.0, np.max(FPQSD_mean/exact_mean)])
for c in ax.collections:
    c.set_edgecolor("face")
cbar.ax.set_ylabel(r'ratio population FPQSD_mean/exact mean')
ax.set_xlabel(r"$q$")
ax.set_ylabel(r"$\delta$")
ax.set_yscale("log")
ax.minorticks_off()
lp.savefig("Figure7-A")
plt.close(fig)

#------------------------------ B ------------------------------

fig, ax  = lp.newfig(0.6)

cax = ax.contourf(stochasticity, variability, np.asarray(FPQSD_var), 100, cmap=plt.cm.inferno)#, levels=100)#, norm=LogNorm(), levels=np.logspace(minimum, maximum, maximum))
cbar = fig.colorbar(cax, ticks=[np.min(FPQSD_var/exact_var), 1.0, np.max(FPQSD_var)])
for c in ax.collections:
    c.set_edgecolor("face")
cbar.ax.set_ylabel(r'Population FPQSD_variance')
ax.set_xlabel(r"$q$")
ax.set_ylabel(r"$\delta$")
ax.set_yscale("log")
ax.minorticks_off()
lp.savefig("Figure7-B")
plt.close(fig)
