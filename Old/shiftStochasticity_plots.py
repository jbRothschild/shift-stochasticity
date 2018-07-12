import sys
import os
sys.path.insert(0, '/home/jrothschild/Research')

#import scipy as sp
import scipy.integrate as inte
import mods.prettyplot as pp
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatter
import numpy as np
import mpmath as mp #MAB - REALLY ONLY NEED hyp3f2

#-------------Math------------------

def maxPop(capacity, stochasticity, variability):
    return (np.rint(capacity/stochasticity*(1+variability/2))).astype(np.int64)

def gaussian(x, mean, std):
    return np.sqrt(2*np.pi*std**2)**(-1) * np.exp( -( (x-mean)/(np.sqrt(2)*std) )**2 )

def birthrate(n, delta, stoch, cap):
    return (1+delta/2)*n-stoch*n**2/cap

def deathrate(n, delta, stoch, cap):
    return delta*n/2+(1-stoch)*n**2/cap

#----------1D Sum & Tau[1]-------------

def rates(stoch, cap, delta, maxSum):
    """
    Function that calculates out outputs different ratios of the product of death rates and product of birth rates.

    Args:
        stoch(int): The stochastic variable in our equations
        cap(int): The capacity
        maxSum(int): The maximal population size (determined when the birthrate is 0)
        delta: The coefficient which determines the stochasticity of the population.
    Returns:
        Q: Some weird ratio of the rates, product(birth(i-1)...birth(1)/death(i)...death(1)).
        S: Some weird ratio of the rates, product(death(i)...death(1)/birth(i)...birth(1)).
    """

    Q = np.array([deathrate(1, delta, stoch, cap)**(-1)])
    S = np.array([deathrate(1, delta, stoch, cap)*birthrate(1, delta, stoch, cap)**(-1)])

    for i in range(2,maxSum):
        Q = np.append(Q,Q[-1]*birthrate(i-1, delta, stoch, cap)*deathrate(i, delta, stoch, cap)**(-1))
        if i <= cap:
            S = np.append(S,S[-1]*deathrate(i, delta, stoch, cap)*birthrate(i, delta, stoch, cap)**(-1))

    return Q, S

def mteSum1D(fixedPoints, stoch, cap, maxSum, delta, tech="sum1d"):
    """
    Function that calculates the Mean Time Extinction from different rate equations.

    Args:
        n(array or int): The number of individuals in the population
        stoch(array or int): The stochastic variable in our equations
        cap(array or int): The capacity
        maxSum(array): The maximal population size (determined when the birthrate is 0)
        delta: The coefficient which determines the stochasticity of the population.
    Returns:
        mte: Mean Time Extinction.
    """
    mte = []

    if np.size(stoch)>1:
        for i, stoch_i in enumerate(stoch):
            q, s = rates(stoch_i, cap, delta, np.int(maxSum[i]))
            mte.append(np.sum(q))
            if tech == "sum1d":
                for j, ratio in enumerate(s):
                    mte[i] += ratio*np.sum(q[j+1:-1])

    elif np.size(cap)>1:
        for i, cap_i in enumerate(cap):
            q, s = rates(stoch, cap_i, delta, np.int(maxSum[i]))
            mte.append(np.sum(q))
            if tech == "sum1d":
                for j, ratio in enumerate(s):
                    mte[i] += ratio*np.sum(q[j+1:-1])

    elif np.size(delta)>1:
        for i, delta_i  in enumerate(delta):
            q, s = rates(stoch, cap, delta_i, np.int(maxSum[i]))
            mte.append(np.sum(q))
            if tech == "sum1d":
                for j, ratio in enumerate(s):
                    mte[i] += ratio*np.sum(q[j+1:-1])

    return np.asarray(mte)

#------------1D SUM Analytics------------------ #MAB

def mte1Dsum_tau1(fixedPoints, stoch, cap, N, delta): #MAB
    """
    Function that calculates the Mean Time for Extinction (MTE) using the (found in Mathematica) 1D sum solution for tau(1)
    n.b. this requires the hypergeometric function from mpmath

    Args:
        fixPoints(array): The fixed points of our equations. - THIS IS NOT USED
        stoch(float): The stochastic variable in our equations
        cap(array): The capacity of our population. - THIS ONLY SEEMS TO WORK WITH A FLOAT
        N: Factor that scales the FP such that n << N - really, N=K=cap - THIS IS NOT USED
        delta: The coefficient which determines the stochasticity of the population.
    Returns:
        The Mean Time for extinction
    """
    return float(2*cap*mp.hyp3f2(1,1,1-cap/stoch-delta*cap/(2*stoch),2,(2+delta*cap/2-2*stoch)/(stoch-2),stoch/(stoch-1))/(2+delta*cap-2*stoch))

def mte1dsum_totaltau1(fixedPoints, stoch, cap, N, delta):
    mte = []
    for i, x in enumerate(cap):
        print("i")
        mte.append(mte1Dsum_tau1(x, stoch, x, N[i], delta))
        print("k")
    return np.asarray(mte)

#------------small n approximation------------------ #MAB

def mtesmalln(fixedPoints, stoch, cap, N, delta): #MAB
    """
    Function that calculates the Mean Time for Extinction (MTE) assuming small n - this might work for the pdf, but should fail for the MTE

    Args:
        fixPoints(array): The fixed points of our equations. - THIS IS NOT USED
        stoch(float): The stochastic variable in our equations
        cap(array): The capacity of our population.
        N: Factor that scales the FP such that n << N - really, N=K=cap - THIS IS NOT USED
        delta: The coefficient which determines the stochasticity of the population.
    Returns:
        The Mean Time for extinction
    """
    return 2.0/(2.0+delta) #note: always less than one

#------------FP full------------------ #MAB

def FPdenom(n, stoch, cap, N, delta): #MAB
    """ denominator is integrated by hand (and confirmed with Mathematica) """
    return np.exp(-2*n/(1.0-2.0*stoch))*((-cap*(1.0+delta)-n*(1.0-2.0*stoch))**(2.*cap*(2.+delta-2.*stoch)/(1.-2.*stoch)**2))

def FPnumer(n, stoch, cap, N, delta): #MAB
    tointegrate = lambda x: FPdenom(x, stoch, cap, N, delta)*2./(birthrate(n, delta, stoch, cap)+deathrate(n, delta, stoch, cap))
    return inte.quad(tointegrate, cap, n)[0]

def mteFP_full(fixedPoints, stoch, cap, N, delta): #MAB
    """
    Function that calculates the Mean Time for Extinction (MTE) from the Fokker-Planck equation, with no further assumptions than n<<K
    n.b. THE INTEGRAL DOES NOT ALWAYS CONVERGE NICELY - USE WITH CAUTION
    """
    tointegrate = lambda x: FPnumer(x, stoch, cap, N, delta)/FPdenom(x, stoch, cap, N, delta)
    return inte.quad(tointegrate, 0., cap)[0] #[0] is the result, [1] is the error

#------------FP full except quasistationary (ie dP/dt=0)------------------ #MAB

def pdfFP_full(x, stoch, cap, delta): #MAB
#    temp = np.exp(-2*cap/(1-2*stoch))*(1/cap)*((cap*(1-2*stoch)+cap*(1+delta))**(1+(2*cap*(2+delta-2*stoch))/(1-2*stoch)**2))
#    return np.exp(-2*x/(1-2*stoch))*(1/x)*((x*(1-2*stoch)+cap*(1+delta))**(1+(2*cap*(2+delta-2*stoch))/(1-2*stoch)**2))/temp
    return np.exp(-2*x/(1-2*stoch)+2*cap/(1-2*stoch))*(cap/x)*(((x*(1-2*stoch)+cap*(1+delta))/((cap*(1-2*stoch)+cap*(1+delta))))**(1+(2*cap*(2+delta-2*stoch))/(1-2*stoch)**2))

def normalization_constant(stoch, cap, delta):
    return inte.quad(lambda x: pdfFP_full(x,stoch,cap,delta), 0, 2.*cap)[0]

def pdfFP_full_normalized(x, stoch, cap, delta):
    return pdfFP_full(x, stoch, cap, delta)/normalization_constant(stoch,cap,delta)

def mte_from_FPpdf(fixedPoints, stoch, cap, N, delta): #MAB
    """
    Function that calculates the Mean Time for Extinction (MTE) from the Fokker-Planck quasi-stationary pdf
    n.b. I was getting some error I haven't sorted out yet
    """
    return 1./(deathrate(1, delta, stoch, cap)*pdfFP_full_normalized(1., stoch, cap, delta))

#------------FP Gaussian------------------ #MAB

def pdfFP_gaussian(x, stoch, cap, delta): #MAB
    """ this is the assumed pdf that solves the quasi-stationary FP equation """
    return np.exp((-(cap-x)**2)/(cap*(2+delta-2*stoch)))/np.sqrt(np.pi*cap*(2+delta-2*stoch))

def mteFP_gaussian(fixedPoints, stoch, cap, N, delta): #MAB
    """
    Function that calculates the Mean Time for Extinction (MTE) using a Gaussian solution to the FP approximation

    Args:
        fixPoints(array): The fixed points of our equations.
        stoch(int): The stochastic variable in our equations
        cap(array): The capacity of our population.
        N: Factor that scales the FP such that n << N - really, N=K=cap
        delta: The coefficient which determines the stochasticity of the population.
    Returns:
        The Mean Time for extinction
    """
    return 2*np.sqrt(np.pi*cap*(2+delta-2*stoch))*np.exp(2*cap/(2+delta-2*stoch))/(1+delta)

#------------_ACTUALLY_ FP WKB------------------ #MAB

def actionFP(n, stoch, cap, N, delta): #MAB
    """
    Function that calculates the action according to the equation defined by the FP approximation

    Args:
        n: The number of individuals in the population
        stoch: The stochastic variable in our equations
        cap: The capacity
        N: Factor that scales the WKB such that n << N
        delta: The coefficient which determines the stochasticity of the population.
    Returns:
        action: The action from the WKB.
        action1: The second term in the exponential expansion of WKB.
        secDerivAction: The second derivative of the action from WKB. Calculated for the constant term in the distribution.
    """
    action = (1/N)*(-1/(1-2*stoch)**2)*( cap*(2+delta-2*stoch)*np.log( (cap*(1+delta)+(1-2*stoch)*n)/(cap*(2+delta-2*stoch)) ) - (1-2*stoch)*(n-cap))

    action1 = 0.0 #this is not necessarily true, but likely doesn't matter

    secDerivAction = (2*N/cap)/(2+delta-2*stoch)

    return action, action1, secDerivAction

def statDistributionFP_quasi(n, fixPoints, stoch, cap, N, delta): #MAB
    """
    Function that calculates the action according to the equation defined by the WKB approximation: Probability Distribution = constant*exp(-K*S(n)-S_1(n)) where S
    except that is uses the action from the Fokker-Planck equation

    Args:
        n(int): The number of individuals in the population
        fixPoints(array): The fixed points of our equations.
        stoch(int): The stochastic variable in our equations
        cap(array): The capacity of our population.
        N: Factor that scales the WKB such that n << N
        delta: The coefficient which determines the stochasticity of the population.
    Returns:
        distribution: The distribution as a function of capacity.
    """
    act, act1, secDerivAct = actionFP(n, stoch, cap, N, delta)
    actF, act1F, secDerivActF = actionFP(fixPoints, stoch, cap, N, delta)

    constant = np.sqrt(secDerivActF/(2*np.pi*N))*np.exp(N*actF + act1F) #constant term in the probability distribution found from WKB

    distribution = constant*np.exp(-N*act - act1)

    return distribution

def mteFP_quasi(fixedPoints, stoch, cap, N, delta): #MAB
    """
    Function that calculates the Mean Time for Extinction (MTE) using the FP action in the otherwise WKB approach

    Args:
        fixPoints(array): The fixed points of our equations.
        stoch(int): The stochastic variable in our equations
        cap(array): The capacity of our population.
        N: Factor that scales the FP such that n << N - really, N=K=cap
        delta: The coefficient which determines the stochasticity of the population.
    Returns:
        The Mean Time for extinction
    """
    return 1/(deathrate(1, delta, stoch, cap)*statDistributionFP_quasi(1, fixedPoints, stoch, cap, N, delta))

#------------WKB realspace------------------
	#I think this is more aptly named "WKB realspace" as contrasted with the gen fn WKB ("WKB momentum space")

def WKB_RS(n, stoch, cap, N, delta):
    """
    Function that calculates the action according to the equation defined by the WKB approximation: Probability Distribution = constant*exp(-N*S(n)-S_1(n)) where S is the action, N very large.

    Args:
        n: The number of individuals in the population
        stoch: The stochastic variable in our equations
        cap: The capacity
        N: Factor that scales the WKB such that n << N
        delta: The coefficient which determines the stochasticity of the population.
    Returns:
        action: The action from the WKB.
        action1: The second term in the exponential expansion of WKB.
        secDerivAction: The second derivative of the action from WKB. Calculated for the constant term in the distribution.
    """
    action = (1.0/N)*( n*np.log( deathrate(n, delta, stoch, cap)/birthrate(n, delta, stoch, cap)) + (cap*delta/(2*(1-stoch))) * np.log(cap*(deathrate(n, delta, stoch, cap))/n) + (cap*(1+delta/2)/stoch) * np.log(cap*birthrate(n, delta, stoch, cap))/n)

    action1 = (1/2)*np.log( deathrate(n, delta, stoch, cap)*birthrate(n, delta, stoch, cap) )

    secDerivAction = (n*N/cap)*( (1-stoch)/deathrate(n, delta, stoch, cap) + stoch/birthrate(n, delta, stoch, cap) )

    return action, action1, secDerivAction

def distributionWKB_RS(n, fixPoints, stoch, cap, N, delta):
    """
    Function that calculates the action according to the equation defined by the WKB approximation: Probability Distribution = constant*exp(-K*S(n)-S_1(n)) where S

    Args:
        n(int): The number of individuals in the population
        fixPoints(array): The fixed points of our equations.
        stoch(int): The stochastic variable in our equations
        cap(array): The capacity of our population.
        N: Factor that scales the WKB such that n << N
        delta: The coefficient which determines the stochasticity of the population.
    Returns:
        distribution: The distribution as a function of capacity.
    """

    act, act1, secDerivAct = WKB_RS(n, stoch, cap, N, delta)
    actF, act1F, secDerivActF = WKB_RS(fixPoints, stoch, cap, N, delta)

    constant = np.sqrt(secDerivActF/(2*np.pi*N))*np.exp(N*actF + act1F) #constant term in the probability distribution found from WKB

    distribution = constant*np.exp(-N*act - act1)

    return distribution

#------WKB momentumspace-------

def WKB_MS(stoch, cap, N, delta):
    """
    Function that calculates the action according to the equation defined by the WKB approximation: Probability Distribution = constant*exp(-N*S(n)-S_1(n)) where S is the action, N very large from the generating function of the Masters quation

    Args:
        n: The number of individuals in the population
        stoch: The stochastic variable in our equations
        cap: The capacity
        N: Factor that scales the WKB such that n << N
        delta: The coefficient which determines the stochasticity of the population.
    Returns:
        action: The action from the WKB.
        action1: The second term in the exponential expansion of WKB.
        secDerivAction: The second derivative of the action from WKB. Calculated for the constant term in the distribution.
    """
    #action =

    #action1 =

    #secDerivAction =

    return 0#action, action1, secDerivAction

def distributionWKB_MS(n, fixPoints, stoch, cap, N, delta):
    """
    Function that calculates the action according to the equation defined by the WKB approximation from generating actions: Probability Distribution = constant*exp(-K*S(n)-S_1(n)) where S is the action.

    Args:
        n(int): The number of individuals in the population
        fixPoints(array): The fixed points of our equations.
        stoch(int): The stochastic variable in our equations
        cap(array): The capacity of our population.
        N: Factor that scales the WKB such that n << N
        delta: The coefficient which determines the stochasticity of the population.
    Returns:
        distribution: The distribution as a function of capacity.
    """

    act, act1, secDerivAct = genActionWKB(n, stoch, cap, N, delta)
    actF, act1F, secDerivActF = genActionWKB(fixPoints, stoch, cap, N, delta)

    #constant =  #constant term in the probability distribution found from WKB

    #distribution = constant*np.exp(-N*act - act1)

    return 0#distribution

#------------Algorithm----------------

def statDistributionAlgo(fixedPoints, stoch, cap, N, delta, std=10):
    """
    Algorithm that calculates the quasistationary conditional probability distribution function

    Args:
        fixPoints(int): The fixed points of our equations (mean of intial guess for distribution).
        stoch(int): The stochastic variable in our equations
        cap(int): The capacity of our population.
        N(int): The maximal population size (determined when the birthrate is 0)
        delta: The coefficient which determines the stochasticity of the population.
    Returns:
        The quasistationary conditional probability distribution function
    """
    popDist = np.asarray([gaussian(n, cap, std) for n in range(0,N+1)])
    birthArray = np.asarray([birthrate(n, delta, stoch, cap) for n in range(0,N+1)])
    deathArray = np.asarray([deathrate(n, delta, stoch, cap) for n in range(0,N+1)])

    popDist = np.insert(popDist,0,0)
    birthArray = np.insert(birthArray,0,0)
    deathArray = np.insert(deathArray,0,0)
    popDist = np.append(popDist,0)
    birthArray = np.append(birthArray,0)
    deathArray = np.append(deathArray,0)

    for i in range(0,5000):
        statDist = (birthArray[:-2]*popDist[:-2] + deathArray[2:]*popDist[2:]) / (birthArray[1:-1] + deathArray[1:-1] - popDist[2]*deathArray[2])
        popDist[1:-1] = np.abs(statDist)#/max(statDist))

    return popDist[1:-2]

def quasiStatDist(fixedPoints, stoch, cap, N, delta, std=10):
    """
    Algorithm that calculates the quasistationary conditional probability distribution function at population 1 for multiple K

    Args:
        fixPoints(int): The fixed points of our equations (mean of intial guess for distribution).
        stoch(int): The stochastic variable in our equations
        cap(int): The capacity of our population.
        N(int): The maximal population size (determined when the birthrate is 0)
        delta: The coefficient which determines the stochasticity of the population.
    Returns:
        The quasistationary conditional probability distribution function
    """
    dist1 = np.asarray([])

    if np.size(stoch)>1:
        for i, stoch_i in enumerate(stoch):
            dist = statDistributionAlgo(cap, stoch_i, cap, maxPop(cap, stoch_i, delta), delta, std)
            dist1 = np.append(dist1,dist[1])
    elif np.size(cap)>1:
        for i, cap_i in enumerate(cap):
            dist = statDistributionAlgo(cap_i, stoch, cap_i, maxPop(cap_i, stoch, delta), delta, std)
            dist1 = np.append(dist1,dist[1])
    elif np.size(delta)>1:
        for i, delta_i in enumerate(delta):
            dist = statDistributionAlgo(cap, stoch, cap, maxPop(cap, stoch, delta_i), delta_i, std)
            dist1 = np.append(dist1,dist[1])

    return dist1

#--------MTE from distributions--------

def mteDist(fixedPoints, stoch, cap, N, delta, technique="WKB_RS"):
    """
    Function that calculates the Mean Time for Extinction (MTE)

    Args:
        fixPoints(array): The fixed points of our equations.
        stoch(int): The stochastic variable in our equations
        cap(array): The capacity of our population.
        N: Factor that scales the WKB such that n << N
        delta: The coefficient which determines the stochasticity of the population.
    Returns:
        The Mean Time for extinction
    """
    if technique == "WKB_RS":
        return 1/(deathrate(1, delta, stoch, cap)*distributionWKB_RS(1, fixedPoints, stoch, cap, N, delta))
    elif technique == "WKB_MS":
        return 1/(deathrate(1, delta, stoch, cap)*distributionWKB_RS(1, fixedPoints, stoch, cap, N, delta))
    elif technique == "cuteAlgo":
        statDist = quasiStatDist(fixedPoints, stoch, cap, N, delta)
        return 1/(deathrate(1, delta, stoch, cap)*statDist)
    else:
        print("No technique chosen for MTE")

def main():
    """
    Plotting the MTE as a function of the capacity given the shifting stochasticity  (q in our model) with
    birth rate = (r(1+delta/2)*n - r*stoch*n**2/cap)
    death rate = (delta/2*n + r*(1-stoch)*n**2/cap)
    """

cwd = os.getcwd()
    #=======================VARIABLES============================
capacity = 100.0
stochasticity = np.linspace(0.01, .99, 99)
variability = np.linspace(0.01, 10.0, 1000)
cap = np.linspace(1.0, capacity, num=capacity)

stochasticity_i = np.linspace(0.01, 0.99, 10)
variability_i = np.linspace(0.01, 9.99, 10)

var = 99
sto = 50
K = 99

MTE = []
mteLegend = []
PDF = []
pdfLegend = []
#==========================PLOTS=============================

#--------------------------2D HeatMaps-----------------------
"""
for i, delta in enumerate(variability):
    MTE.append(mteSum1D(cap[K], stochasticity, cap[K], maxPop(cap[K], stochasticity, delta), delta, "sum1d"))
    print(i)
np.save(cwd + "/Data/heat_MTE_K100",np.asarray(MTE))
pp.plotHeat2D(stochasticity, variability, np.asarray(MTE), lvls=np.logspace(minimum, maximum, maximum), title="Mean Time Extinction", xlab=r"q", ylab=r"$\delta$", zlab="$\tau_{mte}$")
plt.colorbar(ticks=[10**minimum, 10**int((maximum-minimum)/3), 10**int((maximum-minimum)*2/3), 10**maximum])
plt.yscale('log')
plt.show()
"""
#-----------------------Compare Techniques MTE--------------------
"""
np.save(cwd + "/Data/sum1d_MTEvsK_q50d99", mteSum1D(cap, stochasticity[sto], cap, maxPop(cap, stochasticity[sto], variability[var]), variability[var], tech="sum1d"))

np.save(cwd + "/Data/tau1_MTEvsK_q50d99", mteSum1D(cap, stochasticity[sto], cap, maxPop(cap, stochasticity[sto], variability[var]), variability[var], tech="tau1"))

#MTE.append(ss.mte1dsum_totaltau1(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var]))
#mteLegend.append("tau1")
#print("done!")

#MTE.append(ss.mte_from_FPpdf(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var]))
#mteLegend.append("FP QSD")
#print("done!")

np.save(cwd + "/Data/FPGauss_MTEvsK_q50d99", mteFP_gaussian(cap, stochasticity[sto], cap, maxPop(cap, stochasticity[sto], variability[var]), variability[var]))

np.save(cwd + "/Data/FPWKB_MTEvsK_q50d99", mteFP_quasi(cap, stochasticity[sto], cap, maxPop(cap, stochasticity[sto], variability[var]), variability[var]))

np.save(cwd + "/Data/wkbRS_MTEvsK_q50d99", mteDist(cap, stochasticity[sto], cap, maxPop(cap, stochasticity[sto], variability[var]), variability[var], "WKB_RS"))

#MTE.append(ss.mteDist(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var], "WKB_MS"))
#mteLegend.append("WKB Realspace")
#print("done!")

#MTE.append(ss.mteDist(cap, stochasticity[sto], cap, ss.maxPop(cap, stochasticity[sto], variability[var]), variability[var], "cuteAlgo"))
#mteLegend.append("QSD Algorithm")
#print("done!")

pp.multi_plot(cap, np.asarray(MTE), mteLegend, title=r"Mean Time to Extinction with $q=$"+str(stochasticity[sto])+r" and $\delta=$ "+str(variability[var]), xlab=r"capacity, $K$", ylab=r"MTE, $\tau_{MTE}$")
plt.yscale('log')
#plt.locator_params(axis='y', numticks=5)
plt.show()
"""
#---------------------Compare techniques PDF----------------------
MAX = maxPop(cap[K], stochasticity[0], variability[-1])

maximum = maxPop(cap[K], stochasticity[sto], variability[var])
population = np.linspace(1,maximum,maximum)
np.save(cwd + "/Data/population_PDF_q50d99", population)
"""
#PDF.append(pdfFP_full_normalized(population, stochasticity[sto], cap[K], variability[var]))
#pdfLegend.append("FP QSD")
#np.save(cwd + "/Data/FP_PDF_q50d99", pdfFP_full_normalized(population, stochasticity[sto], cap[K], variability[var]))

np.save(cwd + "/Data/FPGauss_PDF_q50d99", pdfFP_gaussian(population, stochasticity[sto], cap[K], variability[var]))

#PDF.append(ss.pdfFP_full_normalized(population, stochasticity[sto], cap[K], variability[var]))
#pdfLegend.append("FP WKB")

#PDF.append(distributionWKB_RS(population, cap[K], stochasticity[sto], cap[K], maximum, variability[var])/sum(distributionWKB_RS(population, cap[K], stochasticity[sto], cap[K], maximum, variability[var])))
#pdfLegend.append("WKB Realspace (normalized)")
#np.save(cwd + "/Data/wkbRS_PDF_q50d99", distributionWKB_RS(population, cap[K], stochasticity[sto], cap[K], maximum, variability[var])/sum(distributionWKB_RS(population, cap[K], stochasticity[sto], cap[K], maximum, variability[var])))


#PDF.append(ss.distributionWKB_MS(population, stochasticity[sto], cap[K], variability[var]))
#pdfLegend.append("WKB Momentumspace")

#np.save(cwd + "/Data/algoQSD_PDF_q50d99", statDistributionAlgo(cap[K], stochasticity[sto], cap[K], maximum, variability[var], 10))

pp.multi_plot(population, np.asarray(PDF), pdfLegend, title=r"Probability Distribution Function with $q=$"+str(stochasticity[sto])+r", $\delta=$ "+str(variability[var])+r" and $K=$"+str(cap[K]), xlab=r"Population", ylab=r"$Probability$")
plt.ylim(ymax=0.1)
plt.show()
"""
#-------------------varying stochasticity-------------------------
#__________________________PDF____________________________________
"""
for stochas in stochasticity_i:
    PDF.append(statDistributionAlgo(cap[K], stochas, cap[K], MAX, variability[var], 10))
    pdfLegend.append(r"$q=$"+str(stochas))
np.save(cwd + "/Data/algoQSDvsQ_PDF_d99K100", np.asarray(PDF))
np.save(cwd + "/Data/algoQSDvsQ_PDF_d99K100_legend", np.asarray(pdeLegend))
pp.multi_plot(range(1,MAX+1), np.asarray(PDF), pdfLegend, title=r"Probability Distribution Function with $\delta=$"+str(variability[var])+r" and $K=$"+str(cap[K]), xlab=r"$Population$", ylab=r"MTE, $\tau_{MTE}$")
plt.ylim(ymax=0.1)
plt.xlim(xmax=150)
plt.show()
"""
#__________________________MTE____________________________________
"""
for stochas in stochasticity_i:
    MTE.append(mteSum1D(cap[K], stochas, cap[K], maxPop(cap[K], stochas, variability), variability, "sum1d"))
    mteLegend.append(r"$q=$"+str(stochas))
    print(stochas)
np.save(cwd + "/Data/sum1dMTEvsD_K100", np.asarray(MTE))
np.save(cwd + "/Data/sum1dMTEvsD_K100_legend", np.asarray(mteLegend))
pp.multi_plot(variability, np.asarray(MTE), mteLegend, title=r"Mean Time Extinction with $K=$"+str(cap[K]), xlab=r"$\delta$", ylab=r"MTE, $\tau_{MTE}$")
plt.yscale('log')
plt.show()
"""
#-------------------varying variability-------------------------
#__________________________PDF____________________________________
"""
for delta in variability_i:
    PDF.append(statDistributionAlgo(cap[K], stochasticity[sto], cap[K], MAX, delta, 10))
    pdfLegend.append(r"$\delta=$"+str(delta))
np.save(cwd + "/Data/algoQSDvsD_PDF_q50K100", np.asarray(PDF))
np.save(cwd + "/Data/algoQSDvsD_PDF_q50K100_legend", np.asarray(pdfLegend))
pp.multi_plot(range(1,MAX+1), np.asarray(PDF), pdfLegend, title=r"Probability Distribution Function with $q=$"+str(stochasticity[sto])+r" and $K=$"+str(cap[K]), xlab=r"$\delta=$", ylab=r"MTE, $\tau_{MTE}$")
plt.ylim(ymax=0.1)
plt.xlim(xmax=150)
plt.show()
"""
#__________________________MTE____________________________________
"""
for delta in variability_i:
    MTE.append(mteSum1D(cap[K], stochasticity, cap[K], maxPop(cap[K], stochasticity, delta), delta, "sum1d"))
    mteLegend.append(r"$\delta=$"+str(delta))
    print(delta)
np.save(cwd + "/Data/sum1dMTEvsQ_K100", np.asarray(MTE))
np.save(cwd + "/Data/sum1dMTEvsQ_K100_legend", np.asarray(mteLegend))
pp.multi_plot(stochasticity, np.asarray(MTE), mteLegend, title=r"Mean Time Extinction with $K=$"+str(cap[K]), xlab=r"$q$", ylab=r"MTE, $\tau_{MTE}$")
plt.yscale('log')
plt.show()
"""

if __name__ == "__main__":
    main()
