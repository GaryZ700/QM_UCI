#contains functions for constructing hartree basis functions
import constants

###############################################
class basis:
    
    def constructBasis(atomicNumbers, coord):
    #build basis set for each atom

    basis = []

    for atom in range(len(atomicNumbers)):
	basis.append([])
	for alpha in const.alpha[atom]:
            basis["gaussian"][atom].append(sOrbital(coord[atomicNumbers[atom]]),alpha)

###############################################
   def sOrbital(coord, alpha):
   #return s orbital dictionary 
      orbital = {
          
          "coord":coord
	  "alpha":alpha
          "amplitude":((2*alpha)/math.pi)**(0.75)
}      

     return orbital

###############################################
#ef overlap():
#construct ovelap matrix 
