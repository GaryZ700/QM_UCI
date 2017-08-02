#python implementation of the Hartree-Fock Method

import numpy as np
import math

############################
#init main variables

#nucleus information
nucli = {
	
	"position":[[0,0,0]]
	"atomicNumbers":[[1]]

}

#atomic orbitals information
AO = {

       "position":Nucli["position"]
       "alphas":[1]
}

############################
#init important functions

