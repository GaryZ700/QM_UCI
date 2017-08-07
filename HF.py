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
        "R":[[2.0,3.0,2.0], [1.0,2.0,0.0]],

	#atomic number
	"Z":[1,1],
	
	#number of electron
	"N":2.0,
}

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

#init main operators for first run

#nuclear nuclear repulsion
#Vnn = Integrals.buildNucNucRep(system)

#nuclear electron attraction
#Vne = Integrals.buildNucERep()

#overlap
S = Integrals.buildOverlap(basis)

#electron KE
T = Integrals.buildKE(basis)

#nuclear electron coloumb attraction 
V = Integrals.buildNuclearAttraction(basis, system)

#electron electron repulsion matrix
G = Integrals.buildElectronRepulsion(basis)

print(G)
