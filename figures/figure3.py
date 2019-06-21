import numpy as np
import sys, os
sys.path.append(os.path.abspath('..'))
import mte.shiftStochasticity as ss
from fig_variables import stochasticity, K, variability, stochasticity_i, variability_i, cap
from fig_variables import *
import fig_variables as fv
import matplotlib.pyplot as plt
plt.style.use('parameters.mplstyle')  # particularIMporting

from matplotlib.colors import LogNorm #color in log the maps
from matplotlib.ticker import LogLocator #log axes
from mpl_toolkits.axes_grid.inset_locator import inset_axes #axes within axes. How deep can we go????
import matplotlib.patches as patches #curved arrows

ratio=0.37 #plot to subplot
#__________________________A____________________________________
LOAD = True
LOAD_MTE = True

if LOAD:
    mte = []
    for i, stochas in enumerate(stochasticity_i):
        mte.append(ss.mteSum1D(cap[K], stochas, cap[K], ss.maxPop(cap[K], stochas, variability), variability, "sum1d"))
        print np.size(mte)
    np.save("../data/delta_vs_mte.npy", mte)

MTE = np.load("../data/delta_vs_mte.npy")

figname = 'Fig3A'
plt.figure()
ax = plt.gca()

for i, stochas in enumerate(stochasticity_i):
    ax.plot(variability, MTE[i], color=colors_gradient[i])

ax.set_xlabel(r"$\delta$")
ax.set_ylabel(r"$\tau_e$")
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_xlim([np.min(variability),np.max(variability)])
ax.get_yaxis().set_major_locator(LogLocator(numticks=5))
ax.minorticks_off()

#plot within plot

if LOAD_MTE:
    mte = []
    for i, delta in enumerate(variability):
        mte.append(ss.mteSum1D(cap[K], stochasticity, cap[K], ss.maxPop(cap[K], stochasticity, delta), delta, "sum1d"))
    np.save("../data/heat_MTE_K100_log.npy", mte)

MTE2 = np.load("../data/heat_MTE_K100_log.npy")
if np.isinf(np.log10((np.asarray(MTE2)).min())) or np.isnan(np.log10((np.asarray(MTE2)).min())):
    minimum = 0
else:
    minimum = int(np.log10((np.asarray(MTE2)).min()))
if np.isinf(np.log10((np.asarray(MTE2)).max())) or np.isnan(np.log10((np.asarray(MTE2)).max())):
    maximum = int(np.log10(np.finfo(np.float64).max))
else:
    maximum = int(np.log10((np.asarray(MTE2)).max()))

axes_ins = inset_axes(ax,
                    width=mpl.rcParams["figure.figsize"][0]/3, # width = 30% of parent_bbox
                    height=mpl.rcParams["figure.figsize"][1]/3, # height : 1 inch
                    loc=1)
axes_ins.contourf(stochasticity, variability, np.asarray(MTE2), cmap=plt.cm.YlGnBu, norm=LogNorm(), levels=np.logspace(minimum, maximum, maximum))
for c in axes_ins.collections: #So there are no white lines appearing in 2D plot
    c.set_edgecolor("face")
axes_ins.set_yscale("log")
axes_ins.minorticks_off()
axes_ins.set_xticks([])
axes_ins.set_yticks([])
axes_ins.set_xlabel(r"$q$")
axes_ins.set_ylabel(r"$\delta$")
for i, stochas in enumerate(stochasticity_i):
    axes_ins.axvline(stochasticity_i[i], color=colors_gradient[i])

plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')
plt.close()

#__________________________B____________________________________
if LOAD:
    mte = []
    for i, varia in enumerate(variability_i):
        mte.append(ss.mteSum1D(cap[K], stochasticity, cap[K], ss.maxPop(cap[K], stochasticity, varia), varia, "sum1d"))
        print "3B"
    np.save("../data/q_vs_mte.npy", mte)

MTE = np.load("../data/q_vs_mte.npy")

figname = 'Fig3B'
plt.figure()
ax = plt.gca()

for i, varia in enumerate(variability_i):
    ax.plot(stochasticity, MTE[i], color=colors_gradient[i])

ax.set_xlabel(r"$q$")
ax.set_ylabel(r"$\tau_e$")
ax.set_yscale("log")
ax.set_xlim([np.min(stochasticity),np.max(stochasticity)])
ax.get_yaxis().set_major_locator(LogLocator(numticks=5))
ax.minorticks_off()

#plot within plot
axes_ins = inset_axes(ax,
                    width=mpl.rcParams["figure.figsize"][0]/3, # width = 30% of parent_bbox
                    height=mpl.rcParams["figure.figsize"][1]/3, # height : 1 inch
                    loc=2)

axes_ins.contourf(stochasticity, variability, np.asarray(MTE2), cmap=plt.cm.YlGnBu, norm=LogNorm(), levels=np.logspace(minimum, maximum, maximum))
for c in axes_ins.collections: #So there are no white lines appearing in 2D plot
    c.set_edgecolor("face")
axes_ins.set_yscale("log")
axes_ins.minorticks_off()
axes_ins.set_xticks([])
axes_ins.set_yticks([])
axes_ins.set_xlabel(r"$q$")
axes_ins.set_ylabel(r"$\delta$")
axes_ins.yaxis.set_label_position("right")
for i, varia in enumerate(variability_i):
    axes_ins.axhline(variability_i[i], color=colors_gradient[i])

plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')
plt.close()

"""
#__________________________A_ALT____________________________________
fig, ax  = lp.newfig(0.6)

mte = []
for i, stochas in enumerate(stochasticity_i):
    mte.append(ss.mteSum1D(cap[K], stochas, cap[K], ss.maxPop(cap[K], stochas, variability), variability, "sum1d"))
    print "3A"
np.save("../data/delta_vs_mte.npy", mte)

MTE = np.load("../data/delta_vs_mte.npy")

for i, stochas in enumerate(stochasticity_i):
    ax.plot(variability, MTE[i], color=colors_gradient2[i])

ax.set_xlabel(r"$\delta$")
ax.set_ylabel(r"$\tau_e$")
ax.set_yscale("log")
ax.set_xscale("log")
ax.get_yaxis().set_major_locator(LogLocator(numticks=5))
ax.minorticks_off()

for i, stochas in enumerate(stochasticity_i):
    axes_ins.axvline(stochasticity_i[i], color=colors_gradient2[i])

kw = dict(arrowstyle="Simple,tail_width=0.5,head_width=4,head_length=8", color="k")
a = patches.FancyArrowPatch((10**(-1./2.),10**12), (10**(1./2.),10**32),connectionstyle="arc3,rad=.5", **kw)
plt.gca().add_patch(a)

lp.savefig("Figure3-A_ALT")
plt.close(fig)

#__________________________B_ALT____________________________________
fig, ax  = lp.newfig(0.6)

mte = []
for i, varia in enumerate(variability_i):
    mte.append(ss.mteSum1D(cap[K], stochasticity, cap[K], ss.maxPop(cap[K], stochasticity, varia), varia, "sum1d"))
    print "3B"
np.save("../data/q_vs_mte.npy", mte)

MTE = np.load("../data/q_vs_mte.npy")

for i, varia in enumerate(variability_i):
    ax.plot(stochasticity, MTE[i], color=colors_gradient2[i])

ax.set_xlabel(r"$q$")
ax.set_ylabel(r"$\tau_e$")
ax.set_yscale("log")
ax.get_yaxis().set_major_locator(LogLocator(numticks=5))
ax.minorticks_off()



lp.savefig("Figure3-B_ALT")
plt.close(fig)
"""
