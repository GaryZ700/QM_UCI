import numpy as np
import scipy.linalg as scipy

####################################################
#init variables
lmax = 30.0
lmin = -30.0
ngrid = 1024
nwave = 3
w = 0.5

diff = (lmax - lmin)/ngrid

#init position matrix and fock matrix
X = np.zeros(ngrid,dtype=np.complex64)
Fock = np.zeros((ngrid,ngrid),dtype=np.complex64)

####################################################
#important functins defined here
def potent(pos):
   #return potential

   # if (abs(pos) <= 5):
   #     return 0.0
   # return 1000000.0  # square well
    #return 150*np.sin(pos*w*2*0.5)

    return w*pos**2

####################################################
#main code goes here

for igrid in range(ngrid):
    
    X[igrid] = lmin + pos*diff + 0j
    H[igrid[igrid] = potent(X[igrid])  + 0j

for igrid in range(ngrid):
    for igrid2 in range(ngrid):
        if(igrid==igrid2): 
            H[igrid,igrid2] +=  1.0/diff**2 + 0j
            
        if(abs(igrid-igrid2)==1):
            H[igrid,igrid2] += - 0.5/diff**2 + 0j

#calculate eigen values and vectors for fock matrix
E,psi = np.linalg.eigh(fock)


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
