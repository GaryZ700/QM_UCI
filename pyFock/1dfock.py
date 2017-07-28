#creates fock matrix for 1d particle in a box
import math
import numpy as np

###########################################
#init main variables

#number of particles 
N = 2

#alphas and sigmas for each particles
basisGaus = {

        "mus":[0.5,0.0],
        "sigmas":[0.5,0.6],
        "alphas":[0.0,0.0]
}

#size of grid for calculation
ngrid = 1024

#max and min size of 1d box(along the only axis, X)
Xmax = 30
Xmin = -30

#grid movement rate to have box of from max to min
Rx = (Xmax - Xmin)/ngrid

#fock matrix
F = []
#overlap matrix
S = []

#Kinetic Energy Matrix
T = []

#number of cycles for SCF loop
maxCycle = 300

#difference between fock energies at which energy has converged
convergenceCriteria = 0.0000001

###########################################
#init important functions
def r(igrid):
    #returns position of particle with respect to box dimensions from grid point

    return igrid * Rx

###########################################
def basis(position, mu, sigma, alpha):
    #all basies consist of a single guassian

    gaus = math.exp((-((position - alpha)**2) * alpha))

    return gausC(sigma)*gaus

###########################################
def gausC(sigma):
    #gaussian function constant
    
    return 1/math.sqrt(2*math.pi*sigma**2) 

###########################################
def overlap(N, basisSet):
    
    S = []

    for n in range(N):
        S.append([])
        print(n)
        for n1 in range(N):
            S[n].append(gausOverlap(basisSet, n, n1))

    return np.asarray(S)

###########################################
def density(N, eVectors):
    #returns density of fock matrix, or guess matrix
    D = []
    print(N)
    for n in range(N):
	D.append([0.0 for x in range(N)])
        for n1 in range(N):
	    for n2 in range(N/2):
                D[n][n1] += eVectors[n][n2] * eVectors[n][n2]
	   
    return np.asarray(D)

###########################################
def exchange(eVectors, N):
    Exc = []

    for n in range(N):
        Exc.append([])
	for n1 in range(N):
             print(Exc)
             Exc[n] = [0.0 for x in range(N)]
             print(Exc)
	     for n2 in range(N/2):
  	         Exc[n][n1] += (eVectors[n][n2])**(8/3) 

    return np.asarray(Exc)

###########################################
def KE(N, basisSet):
    #builds kinetick energy matrix for hamiltonian
    
    T = []

    for n in range(N):
        T.append([])
        for n1 in range(N):
                
            Px = (basisSet["alphas"][n] * basisSet["mus"][n] + basisSet["alphas"][n1] * basisSet["mus"][n1]) / (basisSet["alphas"][n] + basisSet["alphas"][n1])

            t = 0.0
            t += -2 * basisSet["alphas"][n1] * gausOverlap(basisSet, n, n1)
            t += 2 * (basisSet["alphas"][n1]**2) * ((Px - basisSet["mus"][n])**2) * gausOverlap(basisSet, n, n1) 

            T[n].append(t)
    
    return np.asarray(T)

###########################################
def fockEnergy(F, T, D, N):

    E = 0.0

    for n in range(N):
	for n1 in range(N):
		E += D[n][n1] * T[n][n1] * F[n][n1]

    return E

###########################################
def gausOverlap(basisSet, n, n1):
    #uses gaussian product theory to find overlap between two particles
    #http://www.tina-vision.net/docs/memos/2003-003.pdf

        p = basisSet["alphas"][n] + basisSet["alphas"][n1]
        q = (basisSet["alphas"][n] * basisSet["alphas"][n1]) / p
        Q = basisSet["mus"][n] - basisSet["mus"][n1]
        
        return math.exp(-q*Q**2)
        
  #numerical implementation  
  #  x = -100
  #  Ymax = -float("inf")

  #  #look for max of product gaussian, which is equal to the overlap between the two guassians
  #  while(x<100):

  #      gausProduct = (basis(x, basisSet["alphas"][n],basisSet["sigmas"][n]) * basis(x, basisSet["alphas"][n1],basisSet["sigmas"][n1])) / (gausC(basisSet["sigmas"][n]) * gausC(basisSet["sigmas"][n1]))    
  #      
  #      print(gausProduct)

  #      if(gausProduct > Ymax):
  #          Ymax = gausProduct

  #      x += dx

  #  return Ymax
  #  
###########################################
def SCF(F, T, S, N, basisSet):
    
    E = []
    Exc = []
    
    eValues, eVectors = eig(S)

    #build Matrix to trasform fock matrix from mo basis
    #to orto basis
    X = eVectors * eValues**(-0.5) * np.transpose(eVectors)

    #guess starting MO coeffs
    guessVal, guessVec = eig(np.transpose(X) * T * X)

    D = density(N, guessVec)

    for n in range(N):
	Exc.append([0.0 for x in range(N)])
	

    print("####################3")

    print(D)

    #main scf loop
    nCycle = 0
    while(nCycle < maxCycle):
        
        F = T + Exc

        Fprime = np.transpose(X) * F * X

        fValues, fVectors = eig(F)

        fValues = X * fValues

        D = density(N, fVectors)
        
        Exc = exchange(fVectors, N)
	
	print(N)
        E.append(fockEnergy(F, T, D, N))

        if(E[nCycle] - E[nCycle-1] <= convergenceCriteria):
            print("SCF converged!")
            return F, E[nCycle]
        else:
            print("SCF did not converge during cycle" + str(nCycle))

    print("SCF did not converge!!!")
    return F, E[mCycle]

###########################################
def eig(matrix):
    
    #calculate eigen values and vectors for matrix
    eValues,eVectors = np.linalg.eigh(matrix)

    #sort eigen values and eigen vectors
    idx = eValues.argsort()[::1]   
    eValues = eValues[idx]
    eVectors = eVectors[:,idx]
    
    return np.asarray(eValues[idx]), np.asarray(eVectors[:,idx])

###########################################
#main code goes here

#calculate alphas for gaussians
for n in range(N):
    basisGaus["alphas"][n] = 1/(2 * basisGaus["sigmas"][n]**2)

#construct main matricies
S = overlap(N, basisGaus)
T = KE(N, basisGaus)

print(np.asarray(S))


F, E = SCF(F, T, S, N, basisGaus)

print(E)
    

