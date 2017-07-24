#creates fock matrix for 1d particle in a box
import math
import numpy as np

###########################################
#init main variables

#number of particles 
N = 3

#alphas and sigmas for each particles
basisGaus = {

        "alphas":[-2.0,3.0,1.3],
        "sigmas":[0.5,0.6,0.5]
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



gausWidthConst = 2.355
dx = 0.01

###########################################
#init important functions
def basis(position, alpha, sigma):
    #all basies consist of a single guassian

    gaus = math.exp((-((position - alpha)**2))/(2*sigma**2))

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
        for n1 in range(N):
            S[n].append(gausOverlap(basisSet, n, n1))

    return S

###########################################
def gausOverlap(basisSet, n, n1):
    #uses gaussian product theory to find overlap between two particles
    #http://www.tina-vision.net/docs/memos/2003-003.pdf
    #is numerically implemented
    
    #find width at mid height for both gaussians
#    W = basisSet["sigmas"][n] * gausWidthConst * 20
#    W1 = basisSet["sigmas"][n1] * gausWidthConst * 20
#    
#    Xr = basisSet["alphas"][n] + W 
#    X1r = basisSet["alphas"][n1] + W
#
#    Xl = basisSet["alphas"][n] - W 
#    X1l = basisSet["alphas"][n1] - W
#
#    #find which gaussian is most extends most to the right, and which to the left
#    
#    #if(basisSet["alphas"][n] == basisSet["alphas"][n1] and basisSet["sigmas"][n] == basisSet["sigmas"][n1]):
#     #   return 1
#    
#    if(Xr > X1r):
#        Xmax = Xr + 1.0
#    else: 
#        Xmax = X1r + 1.0
#
#    if(Xl < X1l):
#        x = Xl -1.0
#    else:
#        x = X1r - 1.0
#    
#    Ymax = -float("inf")
#
#    print(x)
#    print(Xmax)
#
    x = -100
    Ymax = -float("inf")

    #look for max of product gaussian, which is equal to the overlap between the two guassians
    while(x<100):

        gausProduct = (basis(x, basisSet["alphas"][n],basisSet["sigmas"][n]) * basis(x, basisSet["alphas"][n1],basisSet["sigmas"][n1])) / (gausC(basisSet["sigmas"][n]) * gausC(basisSet["sigmas"][n1]))    
        
        print(gausProduct)

        if(gausProduct > Ymax):
            Ymax = gausProduct

        x += dx

    return Ymax
    
###########################################
def SCF(F, H, S, N):
    
    E = 0.0
    
    eValues, eVectors = eig(S)

    #build Matrix to trasform fock matrix
    X = eVectors * eValues**(-0.5) * np.transpose(eVectors)
   
    print("####################3")

    print(X)

    return E

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

#construct main matricies
S = overlap(N, basisGaus)

print(np.asarray(S))

#for igrid in range(ngrid):

SCF(F, T, S, N)
    
    

