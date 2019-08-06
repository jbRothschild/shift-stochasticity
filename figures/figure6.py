#import latex_plots as lp
from fig_variables import *
import sys
sys.path.insert(0, '../mte')
import compare_mean_var as comp

from matplotlib.colors import LogNorm #color in log the maps

mean, var = comp.mean_var(comp.exact_solution, stochasticity, variability, capacity)
np.save("../data/heat_exact_mean.npy", mean)
np.save("../data/heat_exact_var.npy", var)

mean = np.load("../data/heat_exact_mean.npy")
var = np.load("../data/heat_exact_var.npy")


#------------------------------ A (MEAN)------------------------------

fig, ax  = lp.newfig(0.6)

cax = ax.contourf(stochasticity, variability, mean.transpose(),100, cmap=plt.cm.inferno)#, levels=25)#, norm=LogNorm(), levels=np.logspace(minimum, maximum, maximum))
cbar = fig.colorbar(cax, ticks=[np.min(mean), (np.max(mean)-np.min(mean))/2., np.max(mean)])
for c in ax.collections:
    c.set_edgecolor("face")
cbar.ax.set_ylabel(r'Population mean')
ax.set_xlabel(r"Quadratic prefactor, $q$")
ax.set_ylabel(r"Linear prefactor, $\delta$")
ax.set_yscale("log")
ax.minorticks_off()
lp.savefig("Figure6-A")
plt.close(fig)

#------------------------------ B (VARIANCE)------------------------------

fig, ax  = lp.newfig(0.6)

var = np.load("../data/heat_exact_var.npy")
var[var < 0] = 0

cax = ax.contourf(stochasticity, variability, var.transpose(),100, cmap=plt.cm.inferno)#, norm=LogNorm(), levels=np.logspace(minimum, maximum, maximum))
cbar = fig.colorbar(cax, ticks=[np.min(var), (np.max(var)-np.min(var))/2., np.max(var)])
for c in ax.collections:
    c.set_edgecolor("face")
cbar.ax.set_ylabel(r'Population variance')
ax.set_xlabel(r"Quadratic prefactor, $q$")
ax.set_ylabel(r"Linear prefactor, $\delta$")
ax.set_yscale("log")
ax.minorticks_off()
lp.savefig("Figure6-B")
plt.close(fig)
