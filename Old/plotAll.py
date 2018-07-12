import sys
import os
sys.path.insert(0, '/home/jrothschild/Research')
cwd = os.getcwd()

import mods.prettyplot as pp
import shiftStochasticity as ss
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatter
import numpy as np

#=======================VARIABLES============================
capacity = 100.0
stochasticity = np.linspace(0.01, .99, 99)
variability = np.logspace(-2.0, 1.0, 1000)
#variability = np.logspace(0.01, 10.0, 1000)
cap = np.linspace(1.0, capacity, num=capacity)

stochasticity_i = np.linspace(0.01, 0.99, 10)
variability_i = np.logspace(-2.0, 2.0, 10)

var = 900
#var = 100
#sto = 90
sto = 10

K = 99

MTE = []
mteLegend = []
PDF = []
pdfLegend = []
#==========================PLOTS=============================
#--------------------------2D HeatMaps-----------------------
#MTE = np.load(cwd + '/Data/heat_MTE_K100.npy')
"""
for i, delta in enumerate(variability):
    MTE.append(ss.mteSum1D(cap[K], stochasticity, cap[K], ss.maxPop(cap[K], stochasticity, delta), delta, "sum1d"))
    print(i)

MTE = np.load(cwd + "/Data/heat_MTE_K100_log.npy")

if np.isinf(np.log10((np.asarray(MTE)).min())) or np.isnan(np.log10((np.asarray(MTE)).min())):
    minimum = 0
else:
    minimum = int(np.log10((np.asarray(MTE)).min()))
if np.isinf(np.log10((np.asarray(MTE)).max())) or np.isnan(np.log10((np.asarray(MTE)).max())):
    maximum = int(np.log10(np.finfo(np.float64).max))
else:
    maximum = int(np.log10((np.asarray(MTE)).max()))
#np.save(cwd + "/Data/heat_MTE_K100_log", np.asarray(MTE))

pp.plotHeat2D(stochasticity, variability, np.asarray(MTE), lvls=np.logspace(minimum, maximum, maximum), title="", xlab=r"q", ylab=r"$\delta$", zlab="$\tau_{mte}$")
plt.colorbar(ticks=[10**minimum, 10**int((maximum-minimum)/3), 10**int((maximum-minimum)*2/3), 10**maximum])
plt.yscale('log')
plt.show()

"""
#-----------------------Compare Techniques MTE--------------------
"""
MTE.append(ss.mteSum1D(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var], tech="sum1d"))
mteLegend.append("1D Sum")
print("done!")

MTE.append(ss.mteSum1D(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var], tech="tau1"))
mteLegend.append(r"$\tau[1]$")
print("done!")

#MTE.append(ss.mte1dsum_totaltau1(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var]))
#mteLegend.append("tau1")
#print("done!")

MTE.append(ss.mte_from_FPpdf(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var]))
mteLegend.append("FP QSD")
print("done!")

MTE.append(ss.mteFP_gaussian(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var]))
mteLegend.append("FP Gaussian")
print("done!")

MTE.append(ss.mteFP_quasi(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var]))
mteLegend.append("FP WKB")
print("done!")

MTE.append(ss.mteDist(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var], "WKB_RS"))
mteLegend.append("WKB Realspace")
print("done!")

#MTE.append(ss.mteDist(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var], "WKB_MS"))
#mteLegend.append("WKB Momentumspace")
#print("done!")

#MTE.append(ss.mteDist(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var], "cuteAlgo"))
#mteLegend.append("QSD Algorithm")
#print("done!")

pp.multi_plot(cap, np.asarray(MTE), mteLegend, title="q: " + str(stochasticity[sto]) + r" , $\delta$: " + str(variability[var]), xlab=r"$K$", ylab=r"$\tau_{MTE}$")
plt.yscale('log')
plt.locator_params(axis='y', numticks=7)
plt.show()

"""
#print ss.mteDist(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var], "cuteAlgo")
"""
#---------------------Compare techniques PDF--------------------
ratio = ss.mteSum1D(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var], tech="sum1d")/ss.mteSum1D(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var], tech="tau1") #depends on variable var and probably sto!
plt.plot(cap,ratio,linewidth=2) #start noticing a plateau as of K=400, depends on q and delta. Find plots and paper calcs at /Research/PopDyn_variation/All_Figures/ratio_tk_t1.png . /Research/PopDyn_variation/All_Figures/calc_ratio_tk_t1.png
#---------------------Compare techniques PDF----------------------

maximum = ss.maxPop(cap[K], stochasticity[sto], variability[var])
population = np.linspace(1,maximum,maximum)

PDF.append(ss.pdfFP_full_normalized(population, stochasticity[sto], cap[K], variability[var]))
pdfLegend.append("FP QSD")

PDF.append(ss.pdfFP_gaussian(population, stochasticity[sto], cap[K], variability[var]))
pdfLegend.append("FP Gaussian")

PDF.append(ss.statDistributionFP_quasi(population, cap[K], stochasticity[sto], cap[K], maximum, variability[var]))
pdfLegend.append("FP WKB")

PDF.append(ss.distributionWKB_RS(population, cap[K], stochasticity[sto], cap[K], maximum, variability[var]))
pdfLegend.append("WKB Realspace")

#PDF.append(ss.distributionWKB_MS(population, stochasticity[sto], cap[K], variability[var]))
#pdfLegend.append("WKB Momentumspace")

PDF.append(ss.statDistributionAlgo(cap[K], stochasticity[sto], cap[K], maximum, variability[var], 10))
pdfLegend.append("QSD Algorithm")

pp.multi_plot(population, np.asarray(PDF), pdfLegend, title="", xlab=r"Population", ylab=r"Probability")
#plt.ylim(ymax=0.1)
plt.yscale('log')
plt.xlim(0,200)
plt.ylim(10**(-48),10**(0))
plt.show()
"""

#print ss.statDistributionAlgo(cap[K], stochasticity[sto], cap[K], maximum, variability[var], 10)

#print ss.statDistributionAlgo2(cap[K], stochasticity[sto], cap[K], maximum, variability[var], 10)
#-------------------varying stochasticity-------------------------
MAX = ss.maxPop(cap[K], stochasticity[0], variability[-1])
#__________________________PDF____________________________________
"""
for stochas in stochasticity_i:
    PDF.append(ss.statDistributionAlgo(cap[K], stochas, cap[K], MAX, variability[var], 10))
    pdfLegend.append(r"$q=$"+str(stochas))
pp.multi_plot(range(1,MAX+1), np.asarray(PDF), pdfLegend, title="", xlab=r"Population", ylab=r"Probability")
plt.ylim(ymax=0.1)
plt.xlim(xmax=150)
plt.show()

"""
#__________________________MTE____________________________________
"""
for stochas in stochasticity_i:
    MTE.append(ss.mteSum1D(cap[K], stochas, cap[K], ss.maxPop(cap[K], stochas, variability), variability, "sum1d"))
    mteLegend.append(r"$q=$"+str(stochas))
    print(stochas)
pp.multi_plot(variability, np.asarray(MTE), mteLegend, title="", xlab=r"$q$", ylab=r"$\tau_{MTE}$")
plt.yscale('log')
plt.show()
"""
#-------------------varying variability-------------------------
#__________________________PDF____________________________________

for delta in variability_i:
    PDF.append(ss.statDistributionAlgo(cap[K], stochasticity[sto], cap[K], MAX, delta, 10))
    #pdfLegend.append(r"$\delta=$"+str(delta))
pp.multi_plot(range(1,MAX+1), np.asarray(PDF), title="", xlab=r"$Population$", ylab=r"Probability density")
plt.ylim(ymax=0.05)
plt.xlim(xmin=10, xmax=150)
plt.show()

#__________________________MTE____________________________________
"""
for delta in variability_i:
    MTE.append(ss.mteSum1D(cap[K], stochasticity, cap[K], ss.maxPop(cap[K], stochasticity, delta), delta, "sum1d"))
    mteLegend.append(r"$\delta=$"+str(delta))
    print(delta)
pp.multi_plot(stochasticity, np.asarray(MTE), mteLegend, title="", xlab=r"$q$", ylab=r"$\tau_{MTE}$")
plt.yscale('log')
plt.show()
"""
