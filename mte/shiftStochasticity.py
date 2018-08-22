import sys
import os
sys.path.insert(0, '/home/jrothschild/Research')

#import scipy as sp
import scipy.integrate as inte
#import mods.prettyplot as pp
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
    #S = np.array([deathrate(1, delta, stoch, cap)*birthrate(1, delta, stoch, cap)**(-1)])
    S = np.array([0])
    if (cap!=1):
        S = np.array([0,deathrate(1, delta, stoch, cap)*birthrate(1, delta, stoch, cap)**(-1)])

    for i in range(2,int(maxSum)+1): #Added plus 1 so that range includes maxSum
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
    return abs(float(2*cap*mp.hyp3f2(1,1,1-cap/stoch-delta*cap/(2*stoch),2,(2+delta*cap/2-2*stoch)/(stoch-2),stoch/(stoch-1))/(2+delta*cap-2*stoch)))

def mte1dsum_totaltau1(fixedPoints, stoch, cap, N, delta):
    mte = []
    for i, x in enumerate(cap):
        if(((2+delta*x/2-2*stoch)/(stoch-2))%1.0!=0.0):
            mte.append(mte1Dsum_tau1(x, stoch, x, N[i], delta))
        else:
            print("pole at K=%f"%(x))
            mte.append(0.0)
    return np.asarray(mte)

#------------small n approximation------------------ #MAB

def pdf_smalln_recursive_unnormalized(x, stoch, cap, delta=1.):
    '''
    This comes from assuming P_n>>P_{n-1} AND P_{n+1}>>P_n.
    Essentially P_{n+1} = b_n/d_{n+1} P_n.
    There should be a way to do this without that second assumption, by replacing b_n with (b_n+d_n).
    (see Meerson(?) via Ovaskainen for more details)
    This also need not be done recursively but can be written with Pochhammer symbols - it doesn't seem to compute as nicely.
    '''
    P0 = 1.E-15 #an arbitrary P0
    if x==1:
        return P0
    else:
        return (birthrate(x-1,delta,stoch,cap))*pdf_smalln_recursive_unnormalized(x-1, stoch, cap, delta)/deathrate(x,delta,stoch,cap)
#        return (birthrate(x-1,delta,stoch,cap)+deathrate(x-1,delta,stoch,cap))*pdf_smalln_recursive_unnormalized(x-1, stoch, cap, delta)/deathrate(x,delta,stoch,cap)
#        return x*(1.+delta/2.-x*stoch/cap)/(x+1)/(delta/2+(1-stoch)*(x+1)/cap)*pdf_smalln_recursive_unnormalized(x-1, stoch, cap, delta)
def pdf_smalln_recursive_normalizer(stoch, cap, delta=1.):
    normalizer = 0.
    for ex in range(1,int(cap*(1+delta)/stoch)):#ie. up to the cutoff
        normalizer+=pdf_smalln_recursive_unnormalized(ex, stoch, cap, delta)
    return normalizer
def pdf_smalln_recursive(x, stoch, cap, delta=1.):
    return pdf_smalln_recursive_unnormalized(x, stoch, cap, delta)/pdf_smalln_recursive_normalizer(stoch,cap,delta)
def pdf_smalln_recursive_list(x, stoch, cap, delta=1.):
    temp=[0. for i in x];
    for j,n in enumerate(x):
        temp[j]=pdf_smalln_recursive_unnormalized(n, stoch, cap, delta)/pdf_smalln_recursive_normalizer(stoch,cap,delta)
    return temp
def mte_smalln_recursive(stoch,cap,delta=1.):#does not like numpy.ndarray arguments
        return pdf_smalln_recursive_normalizer(stoch,cap,delta)/deathrate(1,delta,stoch,cap)
def mte_smalln_recursive_list(stoch,cap,delta=1.):
        temp=[0. for i in cap];
        for j,pacity in enumerate(cap):
            temp[j]=1./(deathrate(1,delta,stoch,pacity)*pdf_smalln_recursive(1,stoch,pacity,delta))
#            temp[j]=pdf_smalln_recursive_normalizer(stoch,pacity,delta)/deathrate(1,delta,stoch,pacity)
        return temp

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

def normalization_constant(x,stoch, cap, delta):
    return sum(pdfFP_full(x, stoch, cap, delta))
    #return inte.quad(lambda x: pdfFP_full(x,stoch,cap,delta), 0, 2.*cap)[0]

def pdfFP_full_normalized_dislikesarrays(x, stoch, cap, delta):
    return pdfFP_full(x, stoch, cap, delta)/normalization_constant(x,stoch,cap,delta)

def pdfFP_full_normalized(x, stoch, cap, delta):
    pdf = []
    for i, n in enumerate(x):
        pdf.append(pdfFP_full(n,stoch,cap,delta)/normalization_constant(stoch,cap,delta))
    return np.asarray(pdf)

def mte_from_FPpdf_dislikesarrays(fixedPoints, stoch, cap, N, delta): #MAB
    """
    Function that calculates the Mean Time for Extinction (MTE) from the Fokker-Planck quasi-stationary pdf
    n.b. the error associated with the numerical integration is dropped, but may be large
    """
    return 1./(deathrate(1, delta, stoch, cap)*pdfFP_full_normalized_dislikesarrays(1., stoch, cap, delta))

def mte_from_FPpdf(fixedPoints, stoch, cap, N, delta):
    mte = []
    for i, x in enumerate(cap):
        mte.append(mte_from_FPpdf_dislikesarrays(x, stoch, x, N[i], delta))
    return np.asarray(mte)

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
    action =  (cap/N/n/(1-stoch))*deathrate(n,delta,stoch,cap)*np.log(cap*deathrate(n,delta,stoch,cap)/n) + (cap/N/stoch/n)*birthrate(n,delta,stoch,cap)*np.log(cap*birthrate(n,delta,stoch,cap)/n) - (cap/N/n)*(deathrate(n,delta,stoch,cap)/(1-stoch) + birthrate(n,delta,stoch,cap)/stoch)
    action1 = (1/2)*np.log( deathrate(n, delta, stoch, cap)*birthrate(n, delta, stoch, cap) )

    secDerivAction = N*( (delta/2+(1-stoch)*2*n/cap)/deathrate(n, delta, stoch, cap) - (1+delta/2-stoch*2*n/cap)/birthrate(n, delta, stoch, cap) ) #ie. N*(d'/d - b'/b)

    return action, action1, secDerivAction

def distributionWKB_RS(n, stoch, cap, delta):
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
        distributionWKB_RS(population, cap[K], stochasticity[sto], cap[K], maximum, variability[var]))
    """
    fixPoints = cap
    N = maxPop(cap, stoch, delta)
    act, act1, secDerivAct = WKB_RS(n, stoch, cap, N, delta)
    actF, act1F, secDerivActF = WKB_RS(fixPoints, stoch, cap, N, delta)
    #constant term in the probability distribution found from WKB

    distribution = np.sqrt(secDerivActF/(2*np.pi*N))*np.exp(-N*(act-actF) - (act1-act1F))

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

def statDistributionAlgo(population, stoch, cap, delta, std=10):
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
    N = np.max(population)
    popDist = np.asarray([gaussian(n, cap, ) for n in range(0,N+1)])
    birthArray = np.asarray([birthrate(n, delta, stoch, cap) for n in range(0,N+1)])
    deathArray = np.asarray([deathrate(n, delta, stoch, cap) for n in range(0,N+1)])

    popDist = np.append(popDist,0)
    birthArray = np.append(birthArray,0)
    deathArray = np.append(deathArray,0)

    for i in range(0,5000):
        statDist = (birthArray[:-2]*popDist[:-2] + deathArray[2:]*popDist[2:]) / (birthArray[1:-1] + deathArray[1:-1] - popDist[2]*deathArray[2])
        #print(max(np.abs(statDist/sum(statDist))-popDist[1:-1]))
        popDist[1:-1] = np.abs(statDist/sum(statDist))

    return popDist[1:-1]

def statDistributionAlgo2(fixedPoints, stoch, cap, N, delta, std=10): #MAB
    popDist = np.asarray([gaussian(n, cap, std) for n in range(0,N+1)]) #from n=0 to n=N
    birthArray = np.asarray([birthrate(n, delta, stoch, cap) for n in range(0,N+1)])
    deathArray = np.asarray([deathrate(n, delta, stoch, cap) for n in range(0,N+1)])
    popDist = np.insert(popDist,0,0) #from n=-1 to n=N
    birthArray = np.insert(birthArray,0,0)
    deathArray = np.insert(deathArray,0,0)
    popDist = np.append(popDist,0) #from n=-1 to n=N+1
    birthArray = np.append(birthArray,0)
    deathArray = np.append(deathArray,0)
    for i in range(0,5000):
        for j in range(N):
            popDist[j+2] = (-birthArray[j]*popDist[j] + (birthArray[j+1]+deathArray[j+1])*popDist[j+1] - deathArray[2]*popDist[j+1]*popDist[2]) / (deathArray[j+2])
        popDist = popDist/sum(popDist)
    return popDist[2:-1] #since we want from n=1

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
            dist = statDistributionAlgo(maxPop(cap, stoch_i, delta), stoch_i, cap, delta, std)
            dist1 = np.append(dist1,dist[1])
    elif np.size(cap)>1:
        for i, cap_i in enumerate(cap):
            dist = statDistributionAlgo(maxPop(cap_i, stoch, delta), stoch, cap_i, delta, std)
            dist1 = np.append(dist1,dist[1])
    elif np.size(delta)>1:
        for i, delta_i in enumerate(delta):
            dist = statDistributionAlgo(maxPop(cap, stoch, delta_i), stoch, cap, delta_i, std)
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
    if technique == "WKB_RS":    #action = (1.0/N)*( n*np.log( deathrate(n, delta, stoch, cap)/birthrate(n, delta, stoch, cap)) + (cap*delta/(2*(1-stoch))) * np.log(cap*(deathrate(n, delta, stoch, cap))/n) + (cap*(1+delta/2)/stoch) * np.log(cap*birthrate(n, delta, stoch, cap))/n)
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

    maximum = maxPop(cap[K], stochasticity[sto], variability[var])
    statDistributionAlgo(maxPop(cap[K],stochasticity[sto],variability[var]), stochasticity[sto], cap[K], variability[var], std=10)
    """"
    POP = []
    for i, maxi in enumerate([int(x) for x in range(K,maximum,10)]):
        POP.append(np.linspace(1,maxi,maxi))
        PDF.append(statDistributionAlgo(maxi, stochasticity[sto], cap[K], variability[var], std=10))
        pdfLegend.append(str(maxi))

    pp.plot1D(POP[0], PDF[0], title=r"Probability Distribution Function with $q=$"+str(stochasticity[sto])+r", $\delta=$ "+str(variability[var])+r" and $K=$"+str(cap[K]), xlab=r"Population", ylab=r"Probability")
    for i, thing in enumerate(PDF):
        if i != 0:
            plt.plot(POP[i],thing)
    #plt.ylim(ymax=0.1)
    plt.yscale('log')
    plt.show()
    """
    """
    #---------------------Compare techniques PDF----------------------
    MAX = maxPop(cap[K], stochasticity[0], variability[-1])

    maximum = maxPop(cap[K], stochasticity[sto], variability[var])
    population = np.linspace(1,maximum,maximum)
    np.save("../data/population_PDF_q50d99", population)
    """

if __name__ == "__main__":
    main()
