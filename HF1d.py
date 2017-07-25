#1d implementation of the hartree fock method for N numbers of particles in a box
#reference: https://www.diva-portal.org/smash/get/diva2:11442/FULLTEXT01.pdf


import math
import numpy as np

#########################################
#init main variables
#total number of grid points
ngrid = 1024

#max/min length of box along x axis
Xmax = 30
Xmin = -30

#number of particles, and basis for each particle
N = 1
basis = {
            "sigma": [0.5],
            "mu": [0],
            "MOCoeff": [1]
        }

#calculate alphas of basis gaussians
basis["alpha"] = []
for n in range(N):
    basis["alpha"].append(1 / (2 * basis["sigma"][n]**2))

#init matrices for calculation
#overlap 
S = np.zeros([N, N])
#Kinetic Energy 
T = np.zeros([N, N])
#Density Matrix
D = np.zeros([N, N])
#hamiltonian
H = np.zeros([N, N])
#Fock Matrix
F = np.zeros([N, N])

#dx for numerical integration
dx = 0.001

#maxcycles for scf iteration
maxCycle = 300

#criteria for scf convergence
convergence = 0.00001

#########################################
def primative(basis, n, x):
    #calculates gaussian values at specified position  

    return math.exp(basis["alpha"][n] * ((x - basis["mu"][n])**2))

#########################################
def primativeDPrime(basis, n, x):
    #1d laplacian operator acting on guassian function
    #equation solved on pen and paper, resulting function coded here

    return primative

#########################################
def r(igrid):
    #returns position with respect to box and igrid

    return ((Xmax - Xmin)/ngrid) * igrid 

#########################################
def KEIntegrand(basis, n, n1, x):
    #KE matrix integral equation after laplacian was applied
    
    a = basis["alpha"][n]
    a1 = basis["alpha"][n1]

    u =  basis["mu"][n]
    u1 = basis["mu"][n1]
    
    #Latex Equation
    #-2a\left(x-u\right)\mathrm{e}^{-a_1\left(x-u_1\right)^2-a\left(x-u\right)^2}
    return -2 * a1 * (x - u1) * primative(basis, n1, x) * primative(basis, n, x) 

#########################################
def KE(basis, n, n1, x):
    #create kinetic energy matrix
    
    #1d laplacian analitically solved by hand, coded below 
    a = basis["alpha"][n]
    a1 = basis["alpha"][n1]

    #numerical integration loop
    X = 0.0
    y = 0.0
    while(X < x):
        
        #y += KEIntegrand(x) * dx
        
   #     y += primative(basis, 0, X) * dx
        
        y += KEIntegrand(basis, n, n1, X) * dx
        
        print(y)

        X += dx
    
    return -0.5 * y

#########################################
def eIntegrator(basis, r, r1, N, x):
    #numerical integrator for electron integral
    
    

    X = 0.0
    y = 0.0

    while(X < x):
        y += primative(basis, N[0], r) *dx 

#########################################
def electronIntegral(basis, n, n1, N, r, r1):
    #calculates electron integral
    
    G = 0.0

    for n2 in range(N):
        for n3 in range(N):
            
            y = eIntegrator(r, r1, [n,n1,n2,n3])


            D[n2][n3] * y - 0.5 * y


#########################################
def overlap(basis, n, n1):
    #creates overlap matrix
    
    alpha = basis["alpha"][n]
    alpha1 = basis["alpha"][n1]

    p = alpha + alpha1
    q = (alpha*alpha1) / p
    Q = alpha - alpha1

    return math.exp(-q*(Q**2))

#########################################
def density(basis, n, n1, N):
    #returns density matrix value
   
    D = 0.0
   
    for iocc in range(int(math.ceil(N/2.0))):
        print("as")
        D += basis["MOCoeff"][n] * basis["MOCoeff"][n1]
    
    D *= 2.0

    return D

#########################################
def particleDensity(basis, N, D, x):
    #returns particle density at specfic point in space
    
    p = 0.0

    for n in range(N):
        for n1 in range(N):
            p+= D[n][n1] * primative(basis, n, x) * primative(basis, n, x) 

    return p

#########################################
#main codes goes here

#init main matrices
for n in range(N):
    for n1 in range(N):
    
        S[n][n1] += overlap(basis, n, n1)
        D[n][n1] += density(basis, n, n1, N)
        T[n][n1] += KE(basis, n, n1, 0)
    
                     #V not used becuase particles are not interacting with one another
        H[n][n1] = T #+ V
            
        #G[n][n1] = electronIntegral(basis, n1, n2, 0, 1) 
            
        F[n][n1] = H[n][n1] + D[n][n1]#+ G[n][n1]

#begin SCF routine
nCycle = 0

E = []
E.append(1000000)

while(nCycle < maxCycle):
   

    eValues, eVectors =  np.linalg.eig(F*np.linalg.inv(S)*D) 
   #print(eValues)
    #print(eVectors)

    E.append(np.matrix.trace(D*H) + 0.5*np.matrix.trace(D))
    
    print(E[nCycle])

    if(abs(E[nCycle] - E[nCycle-1]) <= convergence):
        print("Converged")
        break

    print(E[nCycle])

    nCycle += 1

    
