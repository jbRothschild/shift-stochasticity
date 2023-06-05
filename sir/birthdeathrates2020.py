#birth and death functions for a 2d process - specifically, SIR w/ vital
##from __future__ import division

def bn(n,m,consts):
    return consts[0]*consts[1];
def bm(n,m,consts):
    return consts[2]*m*n/consts[1];
def dn(n,m,consts):
    return n*consts[0];
def dm(n,m,consts):
    return m*consts[0];
def gg(n,m,consts):
    return bn(n,m,consts)+bm(n,m,consts)+dn(n,m,consts)+dm(n,m,consts);
def gn(n,m,consts):
    return bm(n,m,consts)+dn(n,m,consts)+dm(n,m,consts);
def gm(n,m,consts):
    return bn(n,m,consts)+dn(n,m,consts)+dm(n,m,consts);
def gnm(n,m,consts):
    return dn(n,m,consts)+dm(n,m,consts);
