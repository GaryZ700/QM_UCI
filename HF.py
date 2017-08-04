#python implementation of the Hartree-Fock method 

import numpy as np
import math
from basis import Basis
from integrals import Integrals

Basis = Basis()
Integrals = Integrals()

##################################
#init main variables

#contains information user provided information about the  system
system = {
	#nuclear coordinates
	"R":[[0.0,0.0,0.0], [1.4,0.0,0.0]],

	#atomic number
	"Z":[1.0,1.0],
	
	#number of electron
	"N":2.0,
}

#overlap matrix
S = []

#Core hamiltonian
H = []

#Density
P = []

#Fock Matrix
F = []

#basis set 
basis = Basis.buildBasis(system["Z"], system["R"])

##################################
#main code goes here
S = Basis.buildOverlap(basis)
T = Integrals.buildKE(basis)

print(S)
print(T)
