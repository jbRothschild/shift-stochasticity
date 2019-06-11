import numpy as np
import scipy.sparse as spsp
import matplotlib.pyplot as plt

import sys, os
sys.path.append(os.path.abspath('..'))
import mte.shiftStochasticity as ss
from fig_variables import stochasticity, K, variability
import fig_variables as fv
plt.style.use('parameters.mplstyle')  # particularIMporting

def Markovian(delta, stoch, cap):
    upto = min(5*cap, ss.maxPop(cap,stoch,delta))
    upperdiag = np.array([ss.deathrate(n,delta,stoch,cap) for n in range(upto+1)])#I may have mixed up upper and lower
    lowerdiag = np.array([ss.birthrate(n,delta,stoch,cap) for n in range(upto+1)])
    maindiag = -upperdiag-lowerdiag;
    data = np.array([maindiag,upperdiag,lowerdiag])
    diags = np.array([0,+1,-1])
    M = spsp.spdiags(data,diags,upto,upto)
    return M

def get_eigen(stochas, variab, cap):
    eigen0 = np.zeros((variability.shape[0], stochasticity.shape[0]))
    eigen1 = np.zeros((variability.shape[0], stochasticity.shape[0]))
    eigen2 = np.zeros((variability.shape[0], stochasticity.shape[0]))

    for i, delta in enumerate(variability):
        for j, stoch in enumerate(stochasticity):
            M = Markovian(delta, stoch, cap)
            w, v = np.linalg.eig(M.todense()) # vector!
            w = np.real(w)
            w[::-1].sort()
            eigen0[i,j] = -1./w[1]; eigen1[i,j] = -1./w[2]; eigen2[i,j] = -1./w[3]
        print 'done '+str(i)

    np.save("../data/eigen0_K100.npy", eigen0); np.save("../data/eigen1_K100.npy", eigen1); np.save("../data/eigen2_K100.npy", eigen2)
    print "\n\nReady to plot!\n\n"
    return 0

#get_eigen(stochasticity, variability, K )

eigen0 = np.load("../data/eigen0_K100.npy"); eigen1 = np.load("../data/eigen1_K100.npy"); eigen2 = np.load("../data/eigen2_K100.npy")

#fv.plot_heatmap(np.asarray(eigen0), stochasticity, variability, 'mte_eigen', 'mean time to extinction')
fv.plot_heatmap(np.asarray(eigen1), stochasticity, variability, 'correlationtime_eigen', 'correlation time')
fv.plot_heatmap(np.asarray(eigen2), stochasticity, variability, 'eigen3_eigen', 'eigenvalue 3')
