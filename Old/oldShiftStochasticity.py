import sys
sys.path.insert(0, '/home/jrothschild/Research')

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

def mteSum1D(fixedPoints, stoch, cap, maxSum, delta, tech = "sum1d"):
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
            print(i)
            q, s = rates(stoch, cap, delta_i, np.int(maxSum[i]))
            mte.append(np.sum(q))
            if tech == "sum1d":
                for j, ratio in enumerate(s):
                    mte[i] += ratio*np.sum(q[j+1:-1])

    return np.asarray(mte)

#------------WKB Realspace------------------

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
    action = (1/N)*( n*np.log( deathrate(n, delta, stoch, cap)/birthrate(n, delta, stoch, cap)) + (cap*delta/(2*(1-stoch))) * np.log((deathrate(n, delta, stoch, cap))/n) + (cap*(1+delta/2)/stoch) * np.log(birthrate(n, delta, stoch, cap))/n)

    action1 = (1/2)*np.log( deathrate(n, delta, stoch, cap)*birthrate(n, delta, stoch, cap) )

    secDerivAction = n*N/cap*( (1-stoch)/deathrate(n, delta, stoch, cap) + stoch/birthrate(n, delta, stoch, cap) )

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

    return popDist[1:-1]

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

#------WKB MomentumSpace------

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

def genStatDistributionWKB(n, fixPoints, stoch, cap, N, delta):
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
        return 1/(deathrate(1, delta, stoch, cap)*statDistributionWKB(1, fixedPoints, stoch, cap, N, delta))
    elif technique == "cuteAlgo":
        statDist = quasiStatDist(fixedPoints, stoch, cap, N, delta)
        return 1/(deathrate(1, delta, stoch, cap)*statDist)
    else:
        print("No technique chosen for MTE")

#==============MAIN==================

def main():
    """
    Plotting the MTE as a function of the capacity given the shifting stochasticity  (q in our model) with
    birth rate = (r(1+delta/2)*n - r*stoch*n**2/cap)
    death rate = (delta/2*n + r*(1-stoch)*n**2/cap)
    """
    #=======================VARIABLES============================
    capacity = 500.0
    stochasticity = np.linspace(0.01, .99, 99)
    variability = np.linspace(0.01, 10.0, 1000)
    cap = np.linspace(1.0, capacity, num=capacity)

    stochasticity_i = np.linspace(0.01, 0.99, 10)
    variability_i = np.linspace(0.01, 9.99, 10)

    #==========================PLOTS=============================
    MTE = []
    mteLegend = []
    PDF = []
    pdfLegend = []

    #K = 300 #CHOICE OF LOW OR HIGH CAPACITY
    K = 300
    #Calculating MTE
    tech = "WKB_RS"
    comp = "delta"
    if (tech == "sum1d" or tech == "tau"):
        for i, delta in enumerate(variability):
            MTE.append(mteSum1D(cap[K], stochasticity, cap[K], maxPop(cap[K], stochasticity, delta), delta, tech))
            print(delta)
    elif (tech == "cuteAlgo" or tech == "WKB_RS"):
        if comp == "stoch":
            variable = stochasticity
            for delta in variability_i:
                print(delta)
                MTE.append(mteDist(cap[K], stochasticity, cap[K], maxPop(cap[K], stochasticity, delta), delta, tech))
                mteLegend.append(r"$\delta$ = " + str(delta))
                x_axis = r"$q$"
        elif comp == "delta":
            variable = variability
            for stochas in stochasticity_i:
                print(stochas)
                MTE.append(mteDist(cap[K], stochas, cap[K], maxPop(cap[K], stochas, variability), variability, tech))
                mteLegend.append(r"$q$ = " + str(stochas))
                x_axis = r"$\delta$"
        else:
            print("Still making capacity stuff")
            MTE = 0
    else:
        print("No technique chosen")
        MTE = 0

    #Setting a maximum and minimum to the MTE
    print(MTE)
    if np.isinf(np.log10((np.asarray(MTE)).min())) or np.isnan(np.log10((np.asarray(MTE)).min())):
        minimum = 0
    else:
        minimum = int(np.log10((np.asarray(MTE)).min()))
    if np.isinf(np.log10((np.asarray(MTE)).max())) or np.isnan(np.log10((np.asarray(MTE)).max())):
        maximum = int(np.log10(np.finfo(np.float64).max))
    else:
        maximum = int(np.log10((np.asarray(MTE)).max()))

    #Plotting
    if tech == "sum1d":
        pp.plotHeat2D(stochasticity, variability, MTE, lvls=np.logspace(minimum, maximum, maximum), title="Mean Time Extinction", xlab=r"q", ylab=r"$\delta$", zlab="$\tau_{mte}$")
        plt.colorbar(ticks=[10**minimum, 10**int((maximum-minimum)/3), 10**int((maximum-minimum)*2/3), 10**maximum])
    elif tech == "cuteAlgo" or tech == "WKB_RS":
        pp.multi_plot(variable, MTE, mteLegend, title=r"Mean Time to Extinction "+tech, xlab=x_axis, ylab=r"MTE, $\tau_{MTE}$")
        plt.locator_params(axis='y', numticks=5)

    plt.yscale('log')
    plt.show()

if __name__ == "__main__":
    main()

def plottingThings():
    #pp.plot1D(cap, mteDist(cap, stochasticity[19], cap, 1000*cap/stochasticity[19]*(1+variability[19]/2), variability[19], technique="WKB"),title="Log of the Mean Time to Extinction", xlab="capacity, K", ylab=r"log of MTE, $\tau_{MTE}$")

    #pp.plot1D(cap, mteSum1D(cap, stochasticity[19], cap, cap/stochasticity[19]*(1+variability[19]/2), variability[19]),title="Log of the Mean Time to Extinction, Sum 1D", xlab="capacity, K", ylab=r"log of MTE, $\tau_{MTE}$")

    #pp.plot1D(range(0, np.int(cap[101]/stochasticity[19]*(1+variability[19]/2))+1),statDistributionAlgo(cap[101], stochasticity[19], cap[101], np.int(cap[101]/stochasticity[19]*(1+variability[19]/2)), variability[19], 10), title="Distribution", xlab = "Population", ylab = "Probability" )

    #plt.yscale('log')
    #plt.locator_params(axis='y', numticks=5)
    #plt.show()

    #================================1D SUM======================================#
    """
    #----------------q vs delta: with capacity K----------------
    K = 300
    for delta in variability:
        MTE.append(mteSum1D(cap[K], stochasticity, cap[K], (cap[K]/stochasticity*(1+delta/2)), delta))
        print delta

    pp.plotHeat2D(stochasticity, variability, MTE, lvls=np.logspace(int(np.log10(min(min(MTE)))), int(np.log10(max(max(MTE)))), int(np.log10(max(max(MTE))))), title="Mean Time Extinction", xlab=r"q", ylab=r"$\delta$", zlab="$\tau_{mte}$")
    plt.yscale('log')
    plt.locator_params(axis='y')
    print np.log10(min(min(MTE)))
    plt.colorbar(ticks=[10**int(np.log10(min(min(MTE)))), 10**int(((np.log10(max(max(MTE))))-np.log10(min(min(MTE))))/3), 10**int(((np.log10(max(max(MTE))))-np.log10(min(min(MTE))))*2/3), 10**int(np.log10(max(max(MTE))))])
    plt.show()
    """

    #-----------capacity vs mte, Different Stochasticity--------------
    """
    for stoch in stochasticity:
        if stoch in [0.05, 0.2, 0.5, 1.0]:
            MTE.append(mteSum1D(cap, stoch, cap, cap/stoch*(1+variability[19]/2), variability[19]))
            mteLegend.append(r"$\delta$ = " + str(stoch))

    pp.multi_plot(cap, MTE, mteLegend, title=r"1D Sum Mean Time to Extinction with $q=$"+str(variability[19]), xlab=r"capacity, $K$", ylab=r"MTE, $\tau_{MTE}$")
    plt.yscale('log')
    plt.locator_params(axis='y', numticks=5)
    plt.show()
    """
    #-----------capacity vs mte, Different Variability----------------
    """
    for delta in variability:
        if delta in [0.05, 0.2, 0.5, 1.0]:
            maxPop = cap/stochasticity[19]*(1+delta/2)
            MTE.append(mteSum1D(cap, stochasticity[19], cap, maxPop, delta))
            mteLegend.append(r"$\delta$ = " + str(delta))

    pp.multi_plot(cap, MTE, mteLegend, title=r"1D Sum Mean Time to Extinction with $q=$"+str(stochasticity[19]), xlab=r"capacity, $K$", ylab=r"MTE, $\tau_{MTE}$")
    plt.yscale('log')
    plt.locator_params(axis='y', numticks=5)
    plt.show()
    """
    #-----------------------delta vs mte,cap=100,stoch=0.2-----------------
    """
    pp.plot1D(variability, mteSum1D(cap[99], stochasticity[19], cap[99], cap[99]/stochasticity*(1+variability/2), variability),title=r"1D Sum Mean Time to Extinction with $K=100$ and $\delta=$"+str(stochasticity[19]), xlab=r"$\delta$", ylab=r"MTE, $\tau_{MTE}$")
    plt.yscale('log')
    plt.locator_params(axis='y', numticks=5)
    plt.show()
    """
    #-----------------------stoch vs mte,cap=100,delta=0.2-----------------
    """
    pp.plot1D(variability, mteSum1D(cap[99], stochasticity, cap[99], cap[99]/stochasticity*(1+variability[19]/2), variability[19]),title=r"1D Sum Mean Time to Extinction with $K=100$ and $q=$"+str(variability[19]), xlab=r"$\delta$", ylab=r"MTE, $\tau_{MTE}$")
    plt.yscale('log')
    plt.locator_params(axis='y', numticks=5)
    plt.show()
    """

    #=================================QSD ALGORITHM====================================#
    #--------------pdf, Different Variability, cap=100,stoch=0.2--------------
    """
    for delta in variability:
        maxPop = np.int(cap[21]/stochasticity[19]*(1+variability[99]/2))
        maxPop = 50
        if delta in [-0.01, 0.2, 0.5, 1.0]:
            PDF.append(statDistributionAlgo(cap[21], stochasticity[19], cap[21], maxPop, delta, 5))
            pdfLegend.append(r"$\delta$ = " + str(delta))

    pp.multi_plot(range(0,maxPop+1), PDF, pdfLegend, title=r"Different Conditional Quasistationary Distributions with $q=$"+str(stochasticity[19]), xlab="Population", ylab="Probability")
    plt.show()
    """
    #--------------pdf, Different Variability, cap=100,stoch=0.2--------------
    """
    for stoch in stochasticity:
        if stoch in [0.05, 0.2, 0.6]:
            maxPop = np.int(cap[101]/stoch*(1+variability[19]/2))
            maxPop = 200
            PDF.append(statDistributionAlgo(cap[101], stoch, cap[101], maxPop, variability[19], 5))
            pdfLegend.append(r"$q$ = " + str(stoch))
            print "done"

    pp.multi_plot(range(0,maxPop+1), PDF, pdfLegend, title=r"Different Conditional Quasistationary Distributions with $\delta=$"+str(variability[19]), xlab="Population", ylab="Probability")
    plt.show()
    """

    #================================1D SUM======================================#
    """
    pp.plot1D(cap, mteDist(cap, stochasticity[19], cap, 1000*cap/stochasticity[19]*(1+variability[19]/2), variability[19], technique="cuteAlgo"),title="Log of the Mean Time to Extinction", xlab="capacity, K", ylab=r"log of MTE, $\tau_{MTE}$")
    """

    #================================1D SUM======================================#
    """
    pp.plot1D(cap, mteDist(cap, stochasticity[19], cap, 1000*cap/stochasticity[19]*(1+variability[19]/2), variability[19], technique="WKB"),title="Log of the Mean Time to Extinction", xlab="capacity, K", ylab=r"log of MTE, $\tau_{MTE}$")
    """

def goodplots():
    #================================1D SUM======================================#
    #----------------q vs delta: with capacity K----------------
    K = 300
    for delta in variability:
        MTE.append(mteSum1D(cap[K], stochasticity, cap[K], (cap[K]/stochasticity*(1+delta/2)), delta))
        print(delta)

    minimum = int(np.log10(min(min(MTE))))
    if np.isinf(np.log10(max(max(MTE)))):
        maximum = int(np.log10(np.finfo(np.float64).max))
    else:
        maximum = int(np.log10(max(max(MTE))))

    print(minimum, maximum)
    pp.plotHeat2D(stochasticity, variability, MTE, lvls=np.logspace(minimum, maximum, maximum), title="Mean Time Extinction", xlab=r"q", ylab=r"$\delta$", zlab="$\tau_{mte}$")
    plt.colorbar(ticks=[10**minimum, 10**(maximum-minimum)/3, 10**(maximum-minimum)*2/3, 10**maximum])

    #--------------
