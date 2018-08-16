import latex_plots as lp
from fig_variables import *
import sys
sys.path.insert(0, '../mte')
import comparaison as comp

from matplotlib.colors import LogNorm #color in log the maps

#mean, var = comp.mean_var(comp.exact_solution, stochasticity, variability, capacity)
#np.save("../data/heat_exact_mean.npy", mean)
#np.save("../data/heat_exact_var.npy", var)

mean = np.load("../data/heat_exact_mean.npy")
var = np.load("../data/heat_exact_var.npy")

var[var < 0] = 0 #Funky stuff happens, sometimes it's negative???
mean[mean < 0] = 0
print np.asarray(mean).shape
print stochasticity.shape
#------------------------------ A (MEAN)------------------------------

fig, ax  = lp.newfig(0.6)

cax = ax.contourf(stochasticity, variability, np.asarray(mean).T, cmap=plt.cm.inferno, levels=100)#, norm=LogNorm(), levels=np.logspace(minimum, maximum, maximum))
cbar = fig.colorbar(cax, ticks=[np.min(mean), capacity, np.max(mean)])
for c in ax.collections:
    c.set_edgecolor("face")
cbar.ax.set_ylabel(r'Population mean')
ax.set_xlabel(r"$q$")
ax.set_ylabel(r"$\delta$")
ax.set_yscale("log")
ax.minorticks_off()
lp.savefig("Figure6-A")
plt.close(fig)

#------------------------------ B (VARIANCE)------------------------------

fig, ax  = lp.newfig(0.6)

var = np.load("../data/heat_exact_var.npy")
var[var < 0] = 0

cax = ax.contourf(stochasticity, variability, np.asarray(var), cmap=plt.cm.inferno, levels=100)#, norm=LogNorm(), levels=np.logspace(minimum, maximum, maximum))
cbar = fig.colorbar(cax, ticks=[np.min(var), (np.max(var)-np.min(var))/2., np.max(var)])
for c in ax.collections:
    c.set_edgecolor("face")
cbar.ax.set_ylabel(r'Population variance')
ax.set_xlabel(r"$q$")
ax.set_ylabel(r"$\delta$")
ax.set_yscale("log")
ax.minorticks_off()
lp.savefig("Figure6-B")
plt.close(fig)
