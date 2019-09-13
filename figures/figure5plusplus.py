import numpy as np
import sys, os
sys.path.append(os.path.abspath('..'))
import mte.shiftStochasticity as ss
from fig_variables import stochasticity, K, variability, stochasticity_i, variability_i, cap
from fig_variables import *
import fig_variables as fv
import matplotlib.pyplot as plt
plt.style.use('parameters.mplstyle')  # particularIMporting

from matplotlib.ticker import LogLocator #log axes

MTE = []; blabels=[]
cap=16;
qlist=np.linspace(0.001,1,100,endpoint=False)

exact_solution = ss.mteSum1D(cap, qlist, cap, ss.maxPop(cap, qlist, variability[var]), variability[var], tech="sum1d")
#ORDER IS VERY IMPORTANT HERE
#FP QSD
#MTE.append(ss.mte_from_FPpdf(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var])); blabels.append("FP full")
#FP full
#MTE.append(ss.mteFP_full_arrayed2(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), deltalist)); blabels.append("FP full")
#MTE.append(ss.mteFP_Laplace_arrayed2(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), deltalist)); blabels.append("FP full")
MTE.append(ss.mteFP_full222(cap, qlist, cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var])); blabels.append("FP full")
#FP full trapz integration
#!!!MTE.append(ss.mteFP_full_numeric(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var])); blabels.append("FP numeric")
#FP GAUSSIAN
MTE.append(ss.mteFP_gaussian(cap, qlist, cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var])); blabels.append("FP Gaussian")
#WKB REALSPACE
MTE.append(1.*ss.mteDist(cap, qlist, cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var], "WKB_RS")); blabels.append("WKB")
#QSD ALGORITHM
#MTE.append(1.*ss.mteDist(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), deltalist, "cuteAlgo")); blabels.append("QSD algorithm")
#TAU 1
#!!!MTE.append(ss.mteSum1D(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var], tech="tau1"))
#SMALL N
#!!!MTE.append(np.zeros(K+1))
#MTE.append(ss.mte_smalln_recursive_list(stochasticity[sto], cap, variability[var]))
#EXACT SOLUTION
dummyvar=ss.mteSum1D(cap, qlist, cap, ss.maxPop(cap, qlist, variability[var]), variability[var], tech="sum1d");
MTE.append(dummyvar); blabels.append("exact solution")

#print(MTE[0])
print("is it even printing?!?")
print(len(MTE))
#print(cap)
colourlist=["r","m","b","g","k","c"]

fnt=18;
#ax.xaxis.set_tick_params(labelsize=fnt,fontsize=fnt)
#font = {'family' : 'normal',
#        'weight' : 'bold',
#        'size'   : fnt}
#import matplotlib as mp
#mp.rc('font', **font)
#plt.rcParams.update({'font.size': fnt})
#plt.tick_params(labelsize=fnt,fontsize=fnt)

figname = 'Fig5' + '_q' + str(round(stochasticity[sto],3)) + '_d' + str(round(variability[var],3))
plt.figure()
ax = plt.gca()
#'''
pdfls = 0
for i in range(0,len(MTE)):
    ax.plot(qlist, MTE[i], color=colors_techniques[i+1], label = blabels[i], linestyle=lines[pdfls+1], lw=3)
    pdfls += 1
'''
    #    if i==3:
    if i+1==3:
#        ax.plot(cap[5:], MTE[i][5:], color=colors_techniques[i], label = techniques[pdfls], linestyle=lines[i], lw=2)
        ax.plot(deltalist, MTE[i], color=colors_techniques[i+1], label = blabels[pdfls], linestyle=lines[i+1], lw=2)
    elif i==4 or i==5:
        print("do not print t[1]")
        pdfls -= 1
    elif i==len(MTE)-1:
        ax.plot(deltalist, MTE[i], color=colors_techniques[i+3], label = blabels[i], linestyle=lines[pdfls+1], lw=3)
    else:
#        ax.plot(cap, MTE[i], color=colors_techniques[i], label = techniques[i], linestyle=lines[pdfls], lw=3)
        ax.plot(deltalist, MTE[i], color=colors_techniques[i+1], label = blabels[i], linestyle=lines[pdfls+1], lw=3)
    pdfls += 1
'''
#'''
#for i in range(0,len(MTE)):
#    ax.plot(cap[0:], MTE[i][0:], color=colourlist[i], label = blabels[i], linestyle=lines[i], lw=2)

ax.tick_params(axis='both', which='major', labelsize=fnt)########
ax.set_xlabel(r"competition mechanism $q$",fontsize=fnt)
ax.set_ylabel(r"mean time to extinction $\tau_e$",fontsize=fnt)
ax.set_yscale("log")
#ax.set_ylim(10**(-20), 10**(0))
ax.get_yaxis().set_major_locator(LogLocator(numticks=5))
ax.minorticks_off()
#ax.set_xlim(1,cap[-1])
ax.legend()#loc = 'upper right')
#plt.savefig('C:\\Users\\lenov\\Pictures\\'+figname+'-q-K16.pdf',transparent=True)#, bbox_inches='tight')
plt.show(); 

#plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
#plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')
plt.close()

"""
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
"""
