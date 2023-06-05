#generates the transition matrix associated with a process (SIR w/ vital)
import birthdeathrates2020 as bdr
##from scipy.linalg import expm
from numpy import dot,ceil
import scipy.sparse as sparse

def maketrix(a,kc,maxi):
    trix = [[0 for y in range(maxi*maxi)] for x in range(maxi*maxi)]
    for i in range(maxi):
      for j in range(maxi):
        for k in range(maxi):
          for l in range(maxi):
            if(i == j and k == l and j + 1 == maxi and l+1 == maxi):
                trix[maxi*i + k][maxi*j + l]=-bdr.gnm(l+1, j + 1,a,kc);
            elif(i == j and k == l and l+1 == maxi):
                trix[maxi*i + k][maxi*j + l]=-bdr.gn(l+1, j + 1,a,kc);
            elif(i == j and k == l and j+1 == maxi):
                trix[maxi*i + k][maxi*j + l]=-bdr.gm(l+1, j + 1,a,kc);
            elif(i == j and k == l):
                trix[maxi*i + k][maxi*j + l]=-bdr.gg(l+1, j + 1,a,kc);
            elif(i == j and k-1 == l):
                trix[maxi*i + k][maxi*j + l]= bdr.bn(l+1, j + 1,a,kc);
            elif(i == j and k+1 == l):
                trix[maxi*i + k][maxi*j + l]= bdr.dn(l+1, j + 1,a,kc);
            elif(i-1 == j and k == l):
                trix[maxi*i + k][maxi*j + l]= bdr.bm(l+1, j + 1,a,kc);
            elif(i+1 == j and k == l):
                trix[maxi*i + k][maxi*j + l]= bdr.dm(l+1, j + 1,a,kc);
    return trix

def pinit(a,kc,maxi):
    n0=int(round(kc/(1.+a)));
    initcond=((n0 - 1)*maxi + n0-1)
    pinit = [0 for y in range(maxi*maxi)]; pinit[initcond]=1.;
    return pinit

##def pfinal(a,kc,maxi,time):
##    matr=maketrix(a,kc,maxi);
##    pfinal = [0 for y in range(maxi*maxi)];
##    pfinal = dot( expm([[matr[x][y]*time for x in range(maxi*maxi)] for y in range(maxi*maxi)]).T, pinit(a,kc,maxi) );
##    # I think when multiplying by the time this creates a whole new matrix
##    return pfinal

#this is most definitely NOT the best way to make it, it could be one loop or five
#NTS: check out sparse.diags
def makesparsetrix(params,maxi):
    rowvec = []; colvec = []; datvec = [];
    for i in range(maxi):
      for j in range(maxi):
        for k in range(maxi):
          for l in range(maxi):
            if(i == j and k == l and j + 1 == maxi and l+1 == maxi):
                rowvec.append(maxi*i + k);
                colvec.append(maxi*j + l);
                datvec.append(-bdr.gnm(l+1, j + 1,params));
            elif(i == j and k == l and l+1 == maxi):
                rowvec.append(maxi*i + k);
                colvec.append(maxi*j + l);
                datvec.append(-bdr.gn(l+1, j + 1,params));
            elif(i == j and k == l and j+1 == maxi):
                rowvec.append(maxi*i + k);
                colvec.append(maxi*j + l);
                datvec.append(-bdr.gm(l+1, j + 1,params));
            elif(i == j and k == l):
                rowvec.append(maxi*i + k);
                colvec.append(maxi*j + l);
                datvec.append(-bdr.gg(l+1, j + 1,params));
            elif(i == j and k-1 == l):
                rowvec.append(maxi*i + k);
                colvec.append(maxi*j + l);
                datvec.append(+bdr.bn(l+1, j + 1,params));
            elif(i == j and k+1 == l):
                rowvec.append(maxi*i + k);
                colvec.append(maxi*j + l);
                datvec.append(+bdr.dn(l+1, j + 1,params));
            elif(i-1 == j and k+1 == l):                   #I THINK this (k->k+1) is the right change to account for birth AND death
                rowvec.append(maxi*i + k);
                colvec.append(maxi*j + l);
                datvec.append(+bdr.bm(l+1, j + 1,params));
            elif(i+1 == j and k == l):
                rowvec.append(maxi*i + k);
                colvec.append(maxi*j + l);
                datvec.append(+bdr.dm(l+1, j + 1,params));
    trix = sparse.coo_matrix((datvec,(rowvec,colvec)), shape=(maxi*maxi,maxi*maxi))
    return trix

def sparsepinit(consts,maxi):
    S0=int(ceil(consts[3]*consts[1]/consts[2]));
    I0=int(ceil(consts[0]*(consts[2]-consts[3])*consts[1]/consts[2]/consts[3]));
    initcond=((I0 - 1)*maxi + S0-1)
    rowvec = [initcond]; colvec = [0]; datvec = [1.0];
    pinit = sparse.coo_matrix((datvec,(rowvec,colvec)),shape=(maxi*maxi,1));
    return pinit

def sparsevec(xIC,yIC,maxi):
    initcond=((yIC - 1)*maxi + xIC-1)
    rowvec = [initcond]; colvec = [0]; datvec = [1.0];
    pinit = sparse.coo_matrix((datvec,(rowvec,colvec)),shape=(maxi*maxi,1));
    return pinit

def exitvec(a,kc,maxi):
    rowvec = []; colvec = []; datvec = [];
    for i in range(maxi):
                rowvec.append(i);
                colvec.append(0);
                datvec.append(-bdr.dm(1,i+1,a,kc));
    pinit = sparse.coo_matrix((datvec,(rowvec,colvec)),shape=(maxi*maxi,1));
    return pinit

def negones(maxi):
    rowvec = []; colvec = []; datvec = [];
    for i in range(maxi*maxi):
                rowvec.append(i);
                colvec.append(0);
                datvec.append(-1.0);
    pinit = sparse.coo_matrix((datvec,(rowvec,colvec)),shape=(maxi*maxi,1));
    return pinit

