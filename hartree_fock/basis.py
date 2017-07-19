#contains functions for constructing hartree basis functions
import constants as consts
import math

###############################################
class basis:
    
    def constructBasis(self, atomicNumbers, coord):
        #build basis set for each atom

        basis = []

        for atom in range(len(atomicNumbers)):
            basis.append([])
            for alpha in consts.alpha[atom]:
                basis[atom].append(self.sOrbital(coord[atom],alpha))

        return basis

###############################################
    def sOrbital(self, coord, alpha):
   #return s orbital dictionary 
	orbital = {
          
            "coord":coord,
	    "alpha":alpha,
            "amplitude":((2*alpha)/math.pi)**(0.75)
}      

	return orbital

###############################################
    def overlap(self,basis):
    #constructs overlap function
    

###############################################
    def guasProduct(alpha1, alpha2, x1, x2):
    #use guassian product theorem to return product of two guassians
    
    p = alpha1 + alpha2
    q = (alpha1*alpha2)/p
    
    P = ((alpha1** x1) + (alpha2**x2))/p
    Q = x1 - x2
 
    return math.exp(-q*(Q**2))
###############################################
   def gausOverlap(gaus1,gaus2):
       
        E = []  
      
	for dimension in range(3):
	     E.append(guasProduct()) 
