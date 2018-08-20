#import latex_plots as lp
from fig_variables import *
import sys
sys.path.insert(0, '../mte')
import compare_mean_var as comp
import shiftStochasticity as ss

from matplotlib.colors import LogNorm #color in log the maps

exac_mean = np.load("../data/heat_exact_mean.npy")
exac_var = np.load("../data/heat_exact_var.npy")
exac_var[exac_var < 0] = 0.01 #funky stuff at mean and var
exac_mean[exac_mean < 0] = 0.01

exact_mean = exac_mean.transpose() #made a mistake and for now the saved file is weird. need to change.
exact_var = exac_var.transpose()

#============================FP QSD=============================

FPQSD_mean, FPQSD_var = comp.mean_var(ss.pdfFP_full_normalized, stochasticity, variability, capacity)
np.save("../data/heat_FPQSD_mean.npy", FPQSD_mean)
np.save("../data/heat_FPQSD_var.npy", FPQSD_var)

FPQSD_mean = np.load("../data/heat_FPQSD_mean.npy")
FPQSD_var = np.load("../data/heat_FPQSD_var.npy")

FPQSD_var[FPQSD_var < 0] = 0.01 #Funky stuff happens, sometimes it's negative???
FPQSD_mean[FPQSD_mean < 0] = 0.01

#------------------------------ A ------------------------------
"""
fig, ax  = lp.newfig(0.6)

cax = ax.contourf(stochasticity, variability, np.asarray(FPQSD_mean/exact_mean), 100, cmap=plt.cm.inferno)#, norm=LogNorm(), levels=np.logspace(minimum, maximum, maximum))
cbar = fig.colorbar(cax, ticks=[np.min(FPQSD_mean/exact_mean), 1.0, np.max(FPQSD_mean/exact_mean)]) #Should we useplt.cm.RdBu
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
"""
#============================FP GAUSSIAN=============================

GAUSS_mean, GAUSS_var = comp.mean_var(ss.pdfFP_gaussian, stochasticity, variability, capacity)
np.save("../data/heat_GAUSS_mean.npy", GAUSS_mean)
np.save("../data/heat_GAUSS_var.npy", GAUSS_var)

GAUSS_mean = np.load("../data/heat_GAUSS_mean.npy")
GAUSS_var = np.load("../data/heat_GAUSS_var.npy")

GAUSS_var[GAUSS_var < 0] = 0.01 #Funky stuff happens, sometimes it's negative???
GAUSS_mean[GAUSS_mean < 0] = 0.01


#============================WKB RS=============================

WKBRS_mean, WKBRS_var = comp.mean_var(ss.distributionWKB_RS, stochasticity, variability, capacity)
np.save("../data/heat_WKBRS_mean.npy", WKBRS_mean)
np.save("../data/heat_WKBRS_var.npy", WKBRS_var)

WKB_mean = np.load("../data/heat_WKBRS_mean.npy")
WKBRS_var = np.load("../data/heat_WKBRS_var.npy")

WKBRS_var[WKBRS_var < 0] = 0.01 #Funky stuff happens, sometimes it's negative???
WKBRS_mean[WKBRS_mean < 0] = 0.01


#============================QSD ALGO=============================

QSD_mean, QSD_var = comp.mean_var(ss.statDistributionAlgo, stochasticity, variability, capacity)
np.save("../data/heat_QSD_mean.npy", QSD_mean)
np.save("../data/heat_QSD_var.npy", QSD_var)

QSD_mean = np.load("../data/heat_QSD_mean.npy")
QSD_var = np.load("../data/heat_QSD_var.npy")

QSD_var[QSD_var < 0] = 0.01 #Funky stuff happens, sometimes it's negative???
QSD_mean[QSD_mean < 0] = 0.01
