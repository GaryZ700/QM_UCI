#python implementation of the Hartree-Fock method 

import numpy as np
import math
from basis import Basis

Basis = Basis()

##################################
#init main variables

#contains information user provided information about the  system
system = {
	
	#nuclear coordinates
	"R":[[0,0,0]],

	#atomic number
	"Z":[1],
	
	#number of electron
	"N":1,

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

print(S)
