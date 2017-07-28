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
fock = np.zeros((ngrid,ngrid),dtype=np.complex64)

####################################################
#important functins defined here
def potent(r):
   #return potential

   # if (abs(pos) <= 5):
   #     return 0.0
   # return 1000000.0  # square well
    #return 150*np.sin(pos*w*2*0.5)

    #return w*pos**2

    

####################################################
#main code goes here

for igrid in range(ngrid):
    
    #calculate new position
    X[igrid] = lmin + igrid*diff + 0j

    #calculate diagnoals of fock matrix
    fock[igrid][igrid] = potent(X[igrid])  + 0j

for igrid in range(ngrid):
    for igrid2 in range(ngrid):
        
        #add to diagnoals of fock matrix
        if(igrid==igrid2): 
            fock[igrid,igrid2] +=  1.0/diff**2 + 0j
        
        #add to all position adjacent to diagonals
        if(abs(igrid-igrid2)==1):
            fock[igrid,igrid2] += - 0.5/diff**2 + 0j

#calculate eigen values and vectors for fock matrix
eValues,eVectors = np.linalg.eigh(fock)



#sort eigen values and eigen vectors
idx = eValues.argsort()[::1]   
eValues = eValues[idx]
eVectors = eVectors[:,idx]

#open output files
f = open("energy_eigen_val","w")
f2 = open("eigen_states","w")

#write to output files
for value in eValues:
    f.write(str(value.real) + "\n")            
for i in range(nwave):
    for igrid in range(ngrid):
        f2.write(str(X[igrid].real) + "    " + str(eVectors[igrid,i].real) + "\n")
    f2.write("\n")
