import maketrix2020 as mk20
from scipy.sparse.linalg import spsolve
from numpy import floor
import datetime
#OLD: K=32 takes ~3.5min [0r 7 for buffer 5?], K=64 takes ~1hour [40min?], 128 ~4 hrs

buffr = 2;
consts = [1, 200, 20, 10];
maxtrix=consts[1]*buffr;
print(datetime.datetime.now());
for kconst in range(136,137,4):
  for aa in range(0,1,2):
    #aconst = aa/10.
    homme = open("timemap-20200526.txt",'w');
    #initcondsS = if(consts[3]*consts[1]/consts[2]>1,consts[3]*consts[1]/consts[2],1.0); #n0=int(floor(kconst/(1.+aconst)));
    #initcondsI = if(consts[0]*(consts[2]-consts[3])*consts[1]/consts[2]/consts[3]>1,consts[0]*(consts[2]-consts[3])*consts[1]/consts[2]/consts[3],1.0); #initcond=((n0 - 1)*kconst*buffr + n0-1)
    
    temp3=mk20.makesparsetrix(consts,maxtrix)
    print("matrix acquired at "+str(datetime.datetime.now()))
    pickyvec = mk20.sparsepinit(consts,maxtrix);
    sumthese = spsolve(temp3.tocsc(),pickyvec.tocsc());
    tau = -sum(sumthese);
    print("with the parameters [1, 200, 20, 10] you get tau = %f" %(tau)); print();
    for elem in range(len(sumthese)):
        homme.write("%f\t"%(-sumthese[elem]));
        if((elem+1)%maxtrix==0):
            homme.write("\n");
    homme.close()
print(datetime.datetime.now());

