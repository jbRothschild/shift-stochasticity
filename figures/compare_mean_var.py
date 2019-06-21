import numpy as np
import sys, os
sys.path.append(os.path.abspath('..'))
import mte.shiftStochasticity as ss
from fig_variables import stochasticity, K, variability
import fig_variables as fv
import matplotlib.pyplot as plt
plt.style.use('parameters.mplstyle')  # particularIMporting

def mean_var(probability_function, stochasticity, variability, cap):

    mean = np.zeros((variability.shape[0], stochasticity.shape[0]))
    variance = np.zeros((variability.shape[0], stochasticity.shape[0]))

    for i, delta in enumerate(variability):
        for j, stoch in enumerate(stochasticity):
            maxSum = ss.maxPop(cap, stoch, delta)
            n = np.arange(1, maxSum+1).astype(float)
            pdf = probability_function(n, stoch, cap, delta)
            mean[i,j] = np.sum( np.multiply(n, pdf) )
            variance[i,j] = np.sum(np.multiply(n**2, pdf) ) - mean[i,j]**2
        print "done!" + str(i)

    np.save("../data/mean_prob_K100.npy", mean)
    np.save("../data/error_prob_K100.npy", variance)

    return 0

def exact_solution(population, stoch, cap, delta):
        q, s = ss.rates(stoch, cap, delta, np.max(population))
        q /= np.sum(q)
        return q

#mean_var(exact_solution, stochasticity, variability, K )

mean = np.load("../data/mean_prob_K100.npy")
error = np.load("../data/error_prob_K100.npy")

fv.plot_heatmap(np.asarray(mean), stochasticity, variability, 'MeanProb', 'mean population')
fv.plot_heatmap(np.asarray(error), stochasticity, variability, 'Var', 'variance population')
