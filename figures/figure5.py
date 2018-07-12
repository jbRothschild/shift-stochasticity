import latex_plots as lp #Necessary for the latex default formatting
from fig_variables import *

from matplotlib.ticker import LogLocator #log axes

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
