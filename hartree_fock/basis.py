#contains functions for constructing hartree basis functions
import constants as consts
import math

###############################################
class basis:
    
    def constructBasis(self, atomicNumbers, coord):
        #build basis set for each atom

        basis = []

        basisCounter = -1

        for atom in range(len(atomicNumbers)):
            if(atomicNumbers[atom] == 1):
                #if hydrogen atom
                
                basis.append([])
                basisCounter+=1

                for alpha in consts.alpha[1]:
                    basis[basisCounter].append(self.sOrbital(coord[atom], alpha))
                
                basisCounter += 1
                basis.append([])
                basis[basisCounter].append(self.sOrbital(coord[atom], 0.1612778))
        
        print(basis)

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
    def gausOverlap(self, gaus1, gaus2):
       
        overLap = []
    
        for dimension in range(3):
            overLap.append(self.gausProduct(gaus1["alpha"], gaus2["alpha"], gaus1["coord"][dimension], gaus2["coord"][dimension]))
    

        Overlap = overLap[0]*overLap[1]*overLap[2]*(((math.pi/(gaus1["alpha"]+gaus2["alpha"]))**0.5)**3) * gaus1["amplitude"]*gaus2["amplitude"] 
        
        print("###################3")
        print(Overlap)
        return float(Overlap)


###############################################
    def overlap(self, basis):
    #constructs overlap function
    
        S = []

        for basis1 in range(len(basis)):
            S.append([])
            for basis2 in range(len(basis)):
                S[basis1].append(0.0)
                
                for gaus1 in range(len(basis[basis1])):
                    print(gaus1)
                    for gaus2 in range(len(basis[basis2])):
                     #   print(gaus2)
                      #  print("!!!!!!!!!!!!!!!!!")
                       # print((basis1==basis2))

                        print("&&&&&&&&&&&&&&&&")
                        print(basis[basis2][gaus2])

                        Overlap = self.gausOverlap(basis[basis1][gaus1],basis[basis2][gaus2])
                        normalization = 1 #basis[basis1][gaus1]["alpha"] * basis[basis2][gaus2]["alpha"]

                        S[basis1][basis2] = S[basis1][basis2] + (Overlap*normalization)

        return S

###############################################
    def gausProduct(self, alpha1, alpha2, x1, x2):
    #use guassian product theorem to return product of two guassians
    
        p = alpha1 + alpha2
        q = (alpha1*alpha2)/p
        
        P = ((alpha1** x1) + (alpha2**x2))/p
        Q = x1 - x2
     
        return math.exp(-q*(Q**2))
                                                                                                    
