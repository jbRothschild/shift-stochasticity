import numpy as np
import shiftStochasticity as ss

def mean_var(probability_function, stochasticity, variability, cap):

    mean = np.zeros((variability.shape[0], stochasticity.shape[0]))
    variance = np.zeros((variability.shape[0], stochasticity.shape[0],))

    for i, stoch in enumerate(variability):
        for j, delta in enumerate(stochasticity):
            maxSum = ss.maxPop(cap, stoch, delta)
            n = np.arange(1, maxSum+1).astype(float)
            pdf = probability_function(n, stoch, cap, delta)
            mean[i,j] = np.sum( np.multiply(n, pdf) )
            variance[i,j] = np.sum(np.multiply(n**2, pdf) ) - mean[i,j]**2

    return mean, variance

def exact_solution(population, stoch,cap, delta):
        q, s = ss.rates(stoch, cap, delta, np.max(population))
        q /= np.sum(q)
        return q
