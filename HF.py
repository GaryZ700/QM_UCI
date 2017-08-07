#python implementation of the Hartree-Fock method 

import numpy as np
import math
from basis import Basis
from integrals import Integrals
from scf import SCF

Basis = Basis()
Integrals = Integrals()
SCF = SCF()

##################################
#init main variables

#contains information user provided information about the  system
system = {
	#nuclear coordinates
	"R":[[0.0,0.0,0.0], [0.0,0.0,2.0]],

	#atomic number
	"Z":[1,1],
	
	#number of electron
	"N":2.0,
}

#basis set 
basis = Basis.buildBasis(system["Z"], system["R"])

#scf energy loop threshold
Ediff = 1.0 * (10.0**(-6))

#maximum number of scf cycles to run
maxCycle = 10

##################################
#main code goes here

#init main operators for first run

#overlap
S = Integrals.buildOverlap(basis)

#electron KE
T = Integrals.buildKE(basis)

#nuclear electron coloumb attraction 
Vext = Integrals.buildNuclearAttraction(basis, system)

#electron electron repulsion matrix
G = Integrals.buildElectronRepulsion(basis)

print(G)
print("---------")
print(Vext)

print("---------")
print(T)

print("---------")
print(S)
#number of basis functions
basisNumber = len(basis["alphas"])

#init core hamiltonian 
#KE operator plus Nuclear attraction operator
HCore = T + Vext

#init variables for scf loop

#AO to MO transformation operator
X = scf.getTransform(S)

#init system energy list
#index refers to cycle number
E = [-float("inf")]



#main scf loop
#for cycle in range(maxCycle):
#	
#	
#
#	#break statement
#	#if convergence has occured
#	if( abs(E[cycle] - E[cycle-1]) < Ediff ):
#		print("SCF convergence occured in " + str(cycle) + " cycles. \n")
#		break
#
##if convergence did not occur
#if(cycle = maxCycle-1):
#	print("Energy convergence did not occur after " + str(maxCycle) + " cycles. \n")

