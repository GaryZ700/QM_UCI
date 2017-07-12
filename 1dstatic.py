import numpy as np
import scipy.linalg as scipy

lmax = 30.0
lmin = -30.0
ngrid = 1024
nwave = 3

w = 0.5
def potent(pos):
   # if (abs(pos) <= 5):
   #     return 0.0
   # return 1000000.0  # square well
    #return 150*np.sin(pos*w*2*0.5)

    return w*pos**2

diff = (lmax - lmin)/ngrid

X = np.zeros(ngrid,dtype=np.complex64)
H = np.zeros((ngrid,ngrid),dtype=np.complex64)
for pos in range(ngrid):
    X[pos] = lmin + pos*diff + 0j
    H[pos][pos] = potent(X[pos])  + 0j

for i in range(ngrid):
    for j in range(ngrid):
        if(i==j): 
            H[i,j] +=  1.0/diff**2 + 0j
            
        if(abs(i-j)==1):
            H[i,j] += - 0.5/diff**2 + 0j

E,psi = np.linalg.eigh(H)


idx = E.argsort()[::1]   
E = E[idx]
psi = psi[:,idx]

f = open("energy_eigen_val","w")
f2 = open("eigen_states","w")

for e in E:
    f.write(str(e.real) + "\n")            
for i in range(nwave):
    for pos in range(ngrid):
        f2.write(str(X[pos].real) + "    " + str(psi[pos,i].real) + "\n")
    f2.write("\n")
