#import sys
#sys.path.insert(0, '/home/jrothschild/Research/')
import os
import latex_plots as lp
cwd = os.getcwd()

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns #niceeee colors
import shiftStochasticity as ss #All our models
from matplotlib.colors import LogNorm #color in log the maps
from matplotlib.ticker import LogLocator #log axes
from mpl_toolkits.axes_grid.inset_locator import inset_axes #axes within axes. How deep can we go????
import matplotlib.patches as patches #curved arrows


techniques = ['FP QSD', 'FP Gaussian', 'WKB Realspace', 'QSD Algorithm', r"$\tau[1]$", "small n approx.", "Exact Solution"]
colors_techniques = plt.cm.viridis(np.linspace(0.,1.,len(techniques))) #BuPu
lines = [':', '-', ':', '-', ':', '-', '-']
n = 10
colors_gradient = plt.cm.inferno(np.linspace(0,1,n))
colors_gradient2 = plt.cm.YlGn(np.linspace(0,1,n))


#=======================VARIABLES============================
capacity = 100.0
stochasticity = np.linspace(0.01, .99, 100)
variability = np.logspace(-1.0, 1.0, 1000)
cap = np.linspace(1.0, capacity, num=capacity)

stochasticity_i = np.linspace(0.1, 0.3, n)
#stochasticity_i = np.linspace(0.1, 0.9, n)
variability_i = np.logspace(-1.0, 0.0, n)

var = 500
sto = 5
K = 99
MAX = ss.maxPop(cap[K], stochasticity[0], variability[-1])

#=======================FIGURE 1==========================
#__________________________A____________________________________

fig, ax  = lp.newfig(0.6)
"""
pdf = []
for i, stochas in enumerate(stochasticity_i):
    pdf.append(ss.statDistributionAlgo(cap[K], stochas, cap[K], MAX, variability[var], 10))
    print "1A"
np.save(cwd + "/Data/pdf_stoch.npy", pdf)
"""
PDF = np.load(cwd + "/Data/pdf_stoch.npy")

for i, stochas in enumerate(stochasticity_i):
    ax.plot(range(1,MAX+1), PDF[i], color=colors_gradient[i])

ax.set_xlabel("Population")
ax.set_ylabel("Probability density function")
ax.set_yscale("log")
ax.set_ylim(10**(-30), 10**(0))
ax.get_yaxis().set_major_locator(LogLocator(numticks=5))
ax.minorticks_off()
ax.set_xlim(1, 200)

lp.savefig("Figure1-A")
plt.close(fig)

#__________________________B____________________________________

fig, ax  = lp.newfig(0.6)
"""
pdf = []
for i, varia in enumerate(variability_i):
    pdf.append(ss.statDistributionAlgo(cap[K], stochasticity[sto], cap[K], MAX, varia, 10))
    print "1B"
np.save(cwd + "/Data/pdf_varia.npy", pdf)
"""
PDF = np.load(cwd + "/Data/pdf_varia.npy")

for i, varia in enumerate(variability_i):
    ax.plot(range(1,MAX+1), PDF[i], color=colors_gradient[i])

ax.set_xlabel("Population")
ax.set_ylabel("Probability density function")
ax.set_yscale("log")
ax.set_ylim(10**(-30), 10**(0))
ax.get_yaxis().set_major_locator(LogLocator(numticks=5))
ax.minorticks_off()
ax.set_xlim(1, 200)

lp.savefig("Figure1-B")
plt.close(fig)

#=======================FIGURE 2==========================
fig, ax  = lp.newfig(0.6)
"""
mte = []
for i, delta in enumerate(variability):
    mte.append(ss.mteSum1D(cap[K], stochasticity, cap[K], ss.maxPop(cap[K], stochasticity, delta), delta, "sum1d"))
np.save(cwd + "/Data/heat_MTE_K100_log.npy", mte)
"""
MTE2 = np.load(cwd + "/Data/heat_MTE_K100_log.npy")
if np.isinf(np.log10((np.asarray(MTE2)).min())) or np.isnan(np.log10((np.asarray(MTE2)).min())):
    minimum = 0
else:
    minimum = int(np.log10((np.asarray(MTE2)).min()))
if np.isinf(np.log10((np.asarray(MTE2)).max())) or np.isnan(np.log10((np.asarray(MTE2)).max())):
    maximum = int(np.log10(np.finfo(np.float64).max))
else:
    maximum = int(np.log10((np.asarray(MTE2)).max()))

cax = ax.contourf(stochasticity, variability, np.asarray(MTE2), cmap=plt.cm.inferno, norm=LogNorm(), levels=np.logspace(minimum, maximum, maximum))
cbar = fig.colorbar(cax, ticks=[10**minimum, 10**int((maximum-minimum)/3), 10**int((maximum-minimum)*2/3), 10**maximum])
for c in ax.collections:
    c.set_edgecolor("face")
cbar.ax.set_ylabel(r'$\tau_e$')
ax.set_xlabel(r"$q$")
ax.set_ylabel(r"$\delta$")
ax.set_yscale("log")
ax.minorticks_off()
lp.savefig("Figure2")
plt.close(fig)

#=======================FIGURE 3==========================
ratio=0.37 #plot to subplot
#__________________________A____________________________________
fig, ax  = lp.newfig(0.6)
"""
mte = []
for i, stochas in enumerate(stochasticity_i):
    mte.append(ss.mteSum1D(cap[K], stochas, cap[K], ss.maxPop(cap[K], stochas, variability), variability, "sum1d"))
    print "3A"
np.save(cwd + "/Data/delta_vs_mte.npy", mte)
"""
MTE = np.load(cwd + "/Data/delta_vs_mte.npy")

for i, stochas in enumerate(stochasticity_i):
    ax.plot(variability, MTE[i], color=colors_gradient2[i])

ax.set_xlabel(r"$\delta$")
ax.set_ylabel(r"$\tau_e$")
ax.set_yscale("log")
ax.set_xscale("log")
ax.get_yaxis().set_major_locator(LogLocator(numticks=5))
ax.minorticks_off()

#plot within plot
axes_ins = inset_axes(ax,
                    width=ratio*lp.figsize(0.6)[0], # width = 30% of parent_bbox
                    height=ratio*lp.figsize(0.6)[1], # height : 1 inch
                    loc=1)
axes_ins.contourf(stochasticity, variability, np.asarray(MTE2), cmap=plt.cm.inferno, norm=LogNorm(), levels=np.logspace(minimum, maximum, maximum))
for c in axes_ins.collections: #So there are no white lines appearing in 2D plot
    c.set_edgecolor("face")
axes_ins.set_yscale("log")
axes_ins.minorticks_off()
axes_ins.set_xticks([])
axes_ins.set_yticks([])
axes_ins.set_xlabel(r"$q$")
axes_ins.set_ylabel(r"$\delta$")
for i, stochas in enumerate(stochasticity_i):
    axes_ins.axvline(stochasticity_i[i], color=colors_gradient2[i])

lp.savefig("Figure3-A")
plt.close(fig)

#__________________________B____________________________________
fig, ax  = lp.newfig(0.6)
"""
mte = []
for i, varia in enumerate(variability_i):
    mte.append(ss.mteSum1D(cap[K], stochasticity, cap[K], ss.maxPop(cap[K], stochasticity, varia), varia, "sum1d"))
    print "3B"
np.save(cwd + "/Data/q_vs_mte.npy", mte)
"""
MTE = np.load(cwd + "/Data/q_vs_mte.npy")

for i, varia in enumerate(variability_i):
    ax.plot(stochasticity, MTE[i], color=colors_gradient2[i])

ax.set_xlabel(r"$q$")
ax.set_ylabel(r"$\tau_e$")
ax.set_yscale("log")
ax.get_yaxis().set_major_locator(LogLocator(numticks=5))
ax.minorticks_off()

#plot within plot
axes_ins = inset_axes(ax,
                    width=ratio*lp.figsize(0.6)[0], # width = 30% of parent_bbox
                    height=ratio*lp.figsize(0.6)[1], # height : 1 inch
                    loc=2)
axes_ins.contourf(stochasticity, variability, np.asarray(MTE2), cmap=plt.cm.inferno, norm=LogNorm(), levels=np.logspace(minimum, maximum, maximum))
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
    axes_ins.axhline(variability_i[i], color=colors_gradient2[i])

lp.savefig("Figure3-B")
plt.close(fig)

#__________________________A_ALT____________________________________
fig, ax  = lp.newfig(0.6)
"""
mte = []
for i, stochas in enumerate(stochasticity_i):
    mte.append(ss.mteSum1D(cap[K], stochas, cap[K], ss.maxPop(cap[K], stochas, variability), variability, "sum1d"))
    print "3A"
np.save(cwd + "/Data/delta_vs_mte.npy", mte)
"""
MTE = np.load(cwd + "/Data/delta_vs_mte.npy")

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
"""
mte = []
for i, varia in enumerate(variability_i):
    mte.append(ss.mteSum1D(cap[K], stochasticity, cap[K], ss.maxPop(cap[K], stochasticity, varia), varia, "sum1d"))
    print "3B"
np.save(cwd + "/Data/q_vs_mte.npy", mte)
"""
MTE = np.load(cwd + "/Data/q_vs_mte.npy")

for i, varia in enumerate(variability_i):
    ax.plot(stochasticity, MTE[i], color=colors_gradient2[i])

ax.set_xlabel(r"$q$")
ax.set_ylabel(r"$\tau_e$")
ax.set_yscale("log")
ax.get_yaxis().set_major_locator(LogLocator(numticks=5))
ax.minorticks_off()



lp.savefig("Figure3-B_ALT")
plt.close(fig)


#=======================FIGURE 4==========================
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

#=======================FIGURE 5==========================
fig, ax  = lp.newfig(0.6)

MTE = []

exact_solution = ss.mteSum1D(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var], tech="sum1d")
#ORDER IS VERY IMPORTANT HERE
#FP QSD
MTE.append(ss.mte_from_FPpdf(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var]))
#FP GAUSSIAN
MTE.append(ss.mteFP_gaussian(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var]))
#WKB REALSPACE
MTE.append(ss.mteDist(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var], "WKB_RS"))
#QSD ALGORITHM
MTE.append(ss.mteDist(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var], "cuteAlgo"))
#TAU 1
MTE.append(ss.mteSum1D(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var], tech="tau1"))
#SMALL N
MTE.append(np.zeros(K+1))
#MTE.append(ss.mte_smalln_recursive_list(stochasticity[sto], cap, variability[var]))
#EXACT SOLUTION
MTE.append(ss.mteSum1D(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var], tech="sum1d"))

for i in range(0,len(MTE)):
    if i==3:
        ax.plot(cap[5:], MTE[i][5:], color=colors_techniques[i], label = techniques[i], linestyle=lines[i])
    else:
        ax.plot(cap, MTE[i], color=colors_techniques[i], label = techniques[i], linestyle=lines[i])

ax.set_xlabel("Carrying capacity, K")
ax.set_ylabel(r"$\tau_e$")
ax.set_yscale("log")
#ax.set_ylim(10**(-20), 10**(0))
ax.get_yaxis().set_major_locator(LogLocator(numticks=5))
ax.minorticks_off()
#ax.set_xlim(1, 200)
ax.legend(loc = 'upper left')

lp.savefig("Figure5")

#=======================FIGURE 5ALT==========================
fig, ax  = lp.newfig(0.6)

MTE = []

exact_solution = ss.mteSum1D(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var], tech="sum1d")

#ORDER IS VERY IMPORTANT HERE
#FP QSD
MTE.append(ss.mte_from_FPpdf(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var]))
#FP GAUSSIAN
MTE.append(ss.mteFP_gaussian(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var]))
#WKB REALSPACE
MTE.append(ss.mteDist(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var], "WKB_RS"))
#QSD ALGORITHM
MTE.append(ss.mteDist(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var], "cuteAlgo"))
#TAU 1
MTE.append(ss.mteSum1D(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var], tech="tau1"))
#SMALL N
MTE.append(np.zeros(K+1))
#MTE.append(ss.mte_smalln_recursive_list(stochasticity[sto], cap, variability[var]))

for i in range(0,len(MTE)):
    if i==3:
        ax.plot(cap, np.asarray(MTE[i][5:])/np.asarray(exact_solution[5:]), color=colors_techniques[i], label = techniques[i], linestyle=lines[i])
    else:
        ax.plot(cap, np.asarray(MTE[i])/np.asarray(exact_solution), color=colors_techniques[i], label = techniques[i], linestyle=lines[i])

ax.set_xlabel("Carrying capacity, K")
ax.set_ylabel(r"\tau_e method / \tau_e exact")
ax.set_yscale("log")
#ax.set_ylim(10**(-20), 10**(0))
ax.get_yaxis().set_major_locator(LogLocator(numticks=5))
ax.minorticks_off()
#ax.set_xlim(1, 200)
ax.legend(loc = 'upper left')

lp.savefig("Figure5_ALT")
