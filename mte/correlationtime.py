"""
2018.08.16 by MattheW
created to test Roger Sidje's Krylov matrix exponentiation (from Matt Smart)
"""

import numpy as np
import scipy.sparse as spsp
from expv import expv
import matplotlib.pyplot as plt

''' yeah, I think this bit should be modular '''
def birthrate(n, delta, stoch, cap):
    return (1+delta/2)*n-stoch*n**2/cap
def deathrate(n, delta, stoch, cap):
    return delta*n/2+(1-stoch)*n**2/cap
def maxPop(cap, stoch, delta):
    maxCutoffFactor=3;
    if(stoch>1./maxCutoffFactor):#this is an artificial cutoff if the matrix is otherwise semi-infinite
        return (np.rint(cap/stoch*(1+delta))).astype(np.int64)
    else:
        return cap*maxCutoffFactor;

def calculateMean(probVec):
    nList=np.array([i+1 for i in range(len(probVec))])
    return np.dot(nList,probVec)

def calculateVariance(probVec1,probVec2):#I'm SUPER unsure how to do this
    nList=np.array([i+1 for i in range(len(probVec1))])
    return sum(nList*probVec1*nList*probVec2)

def pInit(initPop,maxPop):
    pInit=np.zeros(maxPop); pInit[initPop]=1.0;
    return pInit

def peeAtTee(t,delta,stoch,cap,initPop):
    upto=maxPop(cap,stoch,delta)
    upperdiag=np.array([birthrate(n+1,delta,stoch,cap) for n in range(upto)])#I may have mixed up upper and lower
    lowerdiag=np.array([deathrate(n+1,delta,stoch,cap) for n in range(upto)])
    maindiag=-upperdiag-lowerdiag;
    data=np.array([maindiag,upperdiag,lowerdiag])
    diags=np.array([0,+1,-1])
    M=spsp.spdiags(data,diags,upto,upto)
    peeAtTee,err,hump=expv(t,M,pInit(initPop,upto))
    del(err);del(hump);
    return peeAtTee

def autocorrelation(tmax,delta,stoch,cap,resolution=200,initPop=100):#I want initPop=cap by default but I fear I cannot
    initPop=cap #this needs to be changed!!!
    timeList=np.linspace(0,tmax,resolution)
    peeTeeList=peeAtTee(timeList,delta,stoch,cap,initPop)
    peeInitList=np.array([pInit(initPop,maxPop(cap,stoch,delta))])#this is super inefficient!!!
    kappa=calculateVariance(peeTeeList,peeInitList)-calculateMean(peeTeeList)*calculateMean(peeInitList)
    return kappa

def autocorrelationTime(tmax,delta,stoch,cap,initPop=100):
    initPop=cap #this needs to be changed!!!
    kappa=autocorrelation(tmax,delta,stoch,cap,initPop)
    lnkappa=np.log(kappa)
    m,b=fit(lnkappa,m*x+b)
    return -1./m

def main():
    """ defining constants """
    delta=0.1;
    stoch=0.9;
    cap=50;
    upto=maxPop(cap,stoch,delta);
    print(upto)
    
    """  testing the probability vector first """
    afterSomeTime=1.#?
    initialVector=pInit(cap,upto)
    finalVector=peeAtTee(afterSomeTime,delta,stoch,cap,cap)
    plt.figure()
    plt.ylabel('probability',fontsize=16)
    plt.xlabel('n-1',fontsize=16)
    plt.plot(initialVector,label="initially...")
    plt.plot(peeAtTee(afterSomeTime*.25,delta,stoch,cap,cap),label="after quarter time")
    plt.plot(peeAtTee(afterSomeTime*.50,delta,stoch,cap,cap),label="after half time")
    plt.plot(peeAtTee(afterSomeTime*.75,delta,stoch,cap,cap),label="after three quarters time")
    plt.plot(finalVector,label="after some time")
    plt.yscale('log')
    plt.legend()
    plt.show()
    #print(finalVector)
    
    del(finalVector)


if __name__ == "__main__":
    main()
