import latex_plots as lp
from fig_variables import *
import sys
sys.path.insert(0, '../mte')
import compare_mean_var as comp
import shiftStochasticity as ss

from matplotlib.colors import LogNorm #color in log the maps

exact_mean = np.load("../data/heat_exact_mean.npy")
exact_var = np.load("../data/heat_exact_var.npy")
#exac_var[exac_var < 0] = 0.01 #funky stuff at mean and var
#exac_mean[exac_mean < 0] = 0.01


#============================FP QSD=============================
"""
FPQSD_mean, FPQSD_var = comp.mean_var(ss.pdfFP_full_normalized, stochasticity, variability, capacity)
np.save("../data/heat_FPQSD_mean.npy", FPQSD_mean)
np.save("../data/heat_FPQSD_var.npy", FPQSD_var)
"""
FPQSD_mean = np.load("../data/heat_FPQSD_mean.npy")
FPQSD_var = np.load("../data/heat_FPQSD_var.npy")

#FPQSD_var[FPQSD_var < 0.0] = 0.01 #Funky stuff happens, sometimes it's negative???
#FPQSD_mean[FPQSD_mean < 0.0] = 0.01

#------------------------------ A ------------------------------

fig, ax  = lp.newfig(0.6)

cax = ax.contourf(stochasticity, variability, np.asarray(FPQSD_mean/exact_mean), 100, cmap=plt.cm.inferno)#, norm=LogNorm(), levels=np.logspace(minimum, maximum, maximum))
cbar = fig.colorbar(cax, ticks=[np.min(FPQSD_mean/exact_mean), 1.0, np.max(FPQSD_mean/exact_mean)]) #Should we useplt.cm.RdBu
for c in ax.collections:
    c.set_edgecolor("face")
cbar.ax.set_ylabel(r'ratio population FPQSD mean/exact mean')
ax.set_xlabel(r"$q$")
ax.set_ylabel(r"$\delta$")
ax.set_yscale("log")
ax.minorticks_off()
lp.savefig("Figure7-A")
plt.close(fig)

#------------------------------ B ------------------------------

fig, ax  = lp.newfig(0.6)

cax = ax.contourf(stochasticity, variability, np.asarray(FPQSD_var/exact_var), 100, cmap=plt.cm.inferno)#, levels=100)#, norm=LogNorm(), levels=np.logspace(minimum, maximum, maximum))
cbar = fig.colorbar(cax, ticks=[np.min(FPQSD_var/exact_var), 1.0, np.max(FPQSD_var/exact_var)])
for c in ax.collections:
    c.set_edgecolor("face")
cbar.ax.set_ylabel(r'Population FPQSD variance/exact variance')
ax.set_xlabel(r"$q$")
ax.set_ylabel(r"$\delta$")
ax.set_yscale("log")
ax.minorticks_off()
lp.savefig("Figure7-B")
plt.close(fig)

#============================FP GAUSSIAN=============================
"""
GAUSS_mean, GAUSS_var = comp.mean_var(ss.pdfFP_gaussian, stochasticity, variability, capacity)
np.save("../data/heat_GAUSS_mean.npy", GAUSS_mean)
np.save("../data/heat_GAUSS_var.npy", GAUSS_var)
"""
GAUSS_mean = np.load("../data/heat_GAUSS_mean.npy")
GAUSS_var = np.load("../data/heat_GAUSS_var.npy")

#GAUSS_var[GAUSS_var < 0.] = 0.01 #Funky stuff happens, sometimes it's negative???
#GAUSS_mean[GAUSS_mean < 0.] = 0.01

#------------------------------ C ------------------------------

fig, ax  = lp.newfig(0.6)

cax = ax.contourf(stochasticity, variability, np.asarray(GAUSS_mean/exact_mean), 100, cmap=plt.cm.inferno)#, norm=LogNorm(), levels=np.logspace(minimum, maximum, maximum))
cbar = fig.colorbar(cax, ticks=[np.min(GAUSS_mean/exact_mean), 1.0, np.max(GAUSS_mean/exact_mean)]) #Should we useplt.cm.RdBu
for c in ax.collections:
    c.set_edgecolor("face")
cbar.ax.set_ylabel(r'ratio population GAUSS mean/exact mean')
ax.set_xlabel(r"$q$")
ax.set_ylabel(r"$\delta$")
ax.set_yscale("log")
ax.minorticks_off()
lp.savefig("Figure7-C")
plt.close(fig)

#------------------------------ D ------------------------------

fig, ax  = lp.newfig(0.6)

cax = ax.contourf(stochasticity, variability, np.asarray(GAUSS_var/exact_var), 100, cmap=plt.cm.inferno)#, levels=100)#, norm=LogNorm(), levels=np.logspace(minimum, maximum, maximum))
cbar = fig.colorbar(cax, ticks=[np.min(GAUSS_var/exact_var), 1.0, np.max(GAUSS_var/exact_var)])
for c in ax.collections:
    c.set_edgecolor("face")
cbar.ax.set_ylabel(r'Population GAUSS variance/exact variance')
ax.set_xlabel(r"$q$")
ax.set_ylabel(r"$\delta$")
ax.set_yscale("log")
ax.minorticks_off()
lp.savefig("Figure7-D")
plt.close(fig)

#============================WKB RS=============================
"""
WKBRS_mean, WKBRS_var = comp.mean_var(ss.distributionWKB_RS, stochasticity, variability, capacity)
np.save("../data/heat_WKBRS_mean.npy", WKBRS_mean)
np.save("../data/heat_WKBRS_var.npy", WKBRS_var)
"""
WKBRS_mean = np.load("../data/heat_WKBRS_mean.npy")
WKBRS_var = np.load("../data/heat_WKBRS_var.npy")

#WKBRS_var[WKBRS_var < 0] = 0.01 #Funky stuff happens, sometimes it's negative???
#WKBRS_mean[WKBRS_mean < 0] = 0.01

#------------------------------ E ------------------------------

fig, ax  = lp.newfig(0.6)

cax = ax.contourf(stochasticity, variability, np.asarray(WKBRS_mean/exact_mean), 100, cmap=plt.cm.inferno)#, norm=LogNorm(), levels=np.logspace(minimum, maximum, maximum))
cbar = fig.colorbar(cax, ticks=[np.min(WKBRS_mean/exact_mean), 1.0, np.max(WKBRS_mean/exact_mean)]) #Should we useplt.cm.RdBu
for c in ax.collections:
    c.set_edgecolor("face")
cbar.ax.set_ylabel(r'ratio population WKBRS mean/exact mean')
ax.set_xlabel(r"$q$")
ax.set_ylabel(r"$\delta$")
ax.set_yscale("log")
ax.minorticks_off()
lp.savefig("Figure7-E")
plt.close(fig)

#------------------------------ F ------------------------------

fig, ax  = lp.newfig(0.6)

cax = ax.contourf(stochasticity, variability, np.asarray(WKBRS_var/exact_mean), 100, cmap=plt.cm.inferno)#, levels=100)#, norm=LogNorm(), levels=np.logspace(minimum, maximum, maximum))
cbar = fig.colorbar(cax, ticks=[np.min(WKBRS_var/exact_var), 1.0, np.max(WKBRS_var/exact_var)])
for c in ax.collections:
    c.set_edgecolor("face")
cbar.ax.set_ylabel(r'Population WKBRS variance/exact variance')
ax.set_xlabel(r"$q$")
ax.set_ylabel(r"$\delta$")
ax.set_yscale("log")
ax.minorticks_off()
lp.savefig("Figure7-F")
plt.close(fig)

#============================QSD ALGO=============================
"""
QSD_mean, QSD_var = comp.mean_var(ss.statDistributionAlgo, stochasticity, variability, capacity)
np.save("../data/heat_QSD_mean.npy", QSD_mean)
np.save("../data/heat_QSD_var.npy", QSD_var)
"""
QSD_mean = np.load("../data/heat_QSD_mean.npy")
QSD_var = np.load("../data/heat_QSD_var.npy")

#QSD_var[QSD_var < 0.] = 0.01 #Funky stuff happens, sometimes it's negative???
#QSD_mean[QSD_mean < 0.] = 0.01

#------------------------------ G ------------------------------

fig, ax  = lp.newfig(0.6)

cax = ax.contourf(stochasticity, variability, np.asarray(QSD_mean/exact_mean), 100, cmap=plt.cm.inferno)#, norm=LogNorm(), levels=np.logspace(minimum, maximum, maximum))
cbar = fig.colorbar(cax, ticks=[np.min(FPQSD_mean/exact_mean), 1.0, np.max(QSD_mean/exact_mean)]) #Should we useplt.cm.RdBu
for c in ax.collections:
    c.set_edgecolor("face")
cbar.ax.set_ylabel(r'ratio population QSD mean/exact mean')
ax.set_xlabel(r"$q$")
ax.set_ylabel(r"$\delta$")
ax.set_yscale("log")
ax.minorticks_off()
lp.savefig("Figure7-G")
plt.close(fig)

#------------------------------ H ------------------------------

fig, ax  = lp.newfig(0.6)

cax = ax.contourf(stochasticity, variability, np.asarray(QSD_var/exact_var), 100, cmap=plt.cm.inferno)#, levels=100)#, norm=LogNorm(), levels=np.logspace(minimum, maximum, maximum))
cbar = fig.colorbar(cax, ticks=[np.min(QSD_var/exact_var), 1.0, np.max(QSD_var/exact_mean)])
for c in ax.collections:
    c.set_edgecolor("face")
cbar.ax.set_ylabel(r'Population QSD variance/exact variance')
ax.set_xlabel(r"$q$")
ax.set_ylabel(r"$\delta$")
ax.set_yscale("log")
ax.minorticks_off()
lp.savefig("Figure7-H")
plt.close(fig)
