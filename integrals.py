#integral calculator for python implementation of the Hartree-Fock Method
#referece material: "Modern Quantom Chemistry Introduction to Advanced Electronic Structure Theory" by 

import numpy as np
import math
from scipy import special as advMath

#################################
class Integrals:

#################################
	def constants(self, basis, b1, b2, p1, p2):
	#init constants used to evaluate integrals
	        
		#init integral semi-constants dictionary 
		C = {}

		#get important guassian data
		C["a1"] = basis["alphas"][b1][p1]
		C["a2"] = basis["alphas"][b2][p2]

                #center of guassians
		C["r1"] = basis["r"][b1]
		C["r2"] = basis["r"][b2]

                #normalization of each guassians
                C["N1"] = basis["N"][b1][p1]
                C["N2"] = basis["N"][b2][p2]

                #guassian contraction coeffs
		C["c1"] = basis["coeffs"][b1][p1]
		C["c2"] = basis["coeffs"][b2][p2]
		
		#calculate overlap values
                #https://youtu.be/I27vnSDtyiI?t=3m58s
                #https://youtu.be/I27vnSDtyiI?t=7m36s
		C["p"] = C["a1"] + C["a2"]
                C["P"] = [ (C["a1"] * C["r1"][dim] + C["a2"] * C["r2"][dim] ) / C["p"] for dim in range(3) ]
                C["m"] = C["a1"]*C["a2"]
	
                C["q"] = C["m"]/C["p"]
		C["Q"] = [ (C["r1"][dim] - C["r2"][dim])**2 for dim in range(3) ]
	
                C["overlap"] = 1.0

                for dim in range(3):
                    C["overlap"] *= math.exp(-C["q"]*C["Q"][dim])#[ math.exp(-C["q"]*C["Q"][dim]) for dim in range(3) ] 

                #calculate 3D analytical integral
                C["integrand"] = math.sqrt( math.pi / (C["p"]) ) ** 3.0
               
                #calculate product of contraction coeffs and normalization constants
                C["c12"] = C["c1"] * C["c2"]
                C["N12"] = C["N1"] * C["N2"]

		return C

#################################
        def boys(self, n, x):
            #n = order, x = position
            #implementation of the boy's function
            #https://youtu.be/N7A_o0TL_ho?t=1m56s 
            
            #if gamma function goes to inf
            if( x <= 0.0):
                
                n = (2*n) + 1.0
        
                return 1.0 / n
        
            else:
                
                n = n + 0.5
                xP = x ** n
                
                return (advMath.gamma(n) * advMath.gammainc(x, n)) / (2.0*xP)
        
#################################
	def buildKE(self, basis):
		#builds KE matrix T

		#number that represents the number of basis functions being used for the system
		basisNumber = len(basis["alphas"])

#lpha*g1.x0 + g2.alpha*g2.x0)/p;
			
		#init empty overlap matrix
		T  = np.zeros([basisNumber, basisNumber])
		test  = [ [x for x in range(basisNumber)] for y in range(basisNumber) ]

		#iterate over atom basis twice
		for b1 in range(basisNumber):
			for b2 in range(basisNumber):
				test[b1][b2] = [len(basis["alphas"][b1]), len(basis["alphas"][b2])]
				#iterate over primatives used in basis twice
				for p1 in range(len(basis["alphas"][b1])):
					for p2 in range(len(basis["alphas"][b2])):
						
						#get integral constants
						C = self.constants(basis, b1, b2, p1, p2)
                                                

						#https://youtu.be/RHkWFlIhNHo?t=2m4s
                                                #https://youtu.be/W6zfFHE5zIE?t=9m27s
                                                term1 = 3.0 * C["overlap"] * C["a2"] * C["c12"]
                                                
                                                d = [ ((C["P"][dim] - C["r2"][dim]) ** 2) + (1.0/(2.0*C["p"])) for dim in range(3) ]
					        a2 = -2.0 * (C["a2"]**2)             

                                                term2 = [ a2 * d[dim] * C["overlap"] * C["c12"] for dim in range(3) ]
                        
                                                T[b1][b2] += term1 
                                                for dim in range(3):
                                                    T[b1][b2] += term2[dim]

						    
		print(test)                                				
		return T
	
#################################
	def buildOverlap(self, basis):
		#builds overlap matrix S

		#number that represents the number of basis functions being used for the system
		basisNumber = len(basis["alphas"])

		#init empty overlap matrix
		S = np.zeros([basisNumber, basisNumber])
		
		#iterate over atom basis twice
		for b1 in range(basisNumber):
			for b2 in range(basisNumber):
				
				print(b1)
				print(b2)
				print("-------------")
				
				#iterate over primatives used in basis twice
				
				print(basis["coeffs"])
				
				for p1 in range(len(basis["alphas"][b1])):
					for p2 in range(len(basis["alphas"][b2])):
						
						
						#get integral constants
						C = self.constants(basis, b1, b2, p1, p2)							
						#3D overlap times integral portion times normalization constants
                                                S[b1][b2] += C["overlap"] * C["integrand"] * C["N12"] * C["c12"] 

		return S				

#################################
        def buildNuclearAttraction(self, basis, system):
                #builds coulomb matrix for couloumb force between each nucli and electron
                #https://youtu.be/N7A_o0TL_ho?t=3m9s
                
                #number that represents the number of basis functions being used for the system
                basisNumber = len(basis["alphas"])
        
                #init empty overlap matrix
                V = np.zeros([basisNumber, basisNumber])
        
                #loop over all nucli in system
                for nucli in system["Z"]:
        
                    #loop over all orbital baisis functions twice
                    for b1 in range(basisNumber):
                        for b2 in range(basisNumber):
                            
                            #loop over all primative guassians in basis twice
                            for p1 in range(len(basis["alphas"][b1])):
                                for p2 in range(len(basis["alphas"][b2])):
        
                                    #init integral constants
                                    C = self.constants(basis, b1, b2, p1, p2)
                                   
                                    RPA2 = [ ((system["Z"][x] - C["P"][x]) ** 2) for x in range(len(system["Z"])) ] 
                                    RPA2 = sum(RPA2) * C["p"]
        
                                    #calculate coulomb nuclear attraction
                                    V[b1][b2] += -2.0 * nucli * self.boys(0, RPA2) * C["N12"] * C["c12"]   
        
                return V

#################################
        def buildElectronRepulsion(self, basis):
            #builds electron repulsion matrix
            #https://youtu.be/vlNxTTF1kK4?t=1m25s
        
            #number of basis functions used in total
            basisNumber = len(basis["alphas"])
        
            #init electron repulsion matrix
            G = np.zeros( [basisNumber, basisNumber, basisNumber, basisNumber] )
        
            #loop over basis functions twice
            for b1 in range(basisNumber):
                for b2 in range(basisNumber):
                    
                    #loop over primative guassians twice
                    for p1 in range(len(basis["alphas"][b1])):
                        for p2 in range(len(basis["alphas"][b2])):
        
                            #init integral contants 
                            C = self.constants(basis, b1, b2, p1, p2)
        
                            #amplitude of two basis 
                            A12 = C["overlap"] * C["c12"] * C["N12"]
        
                            #loop over basis functions twice more
                            for b3 in range(basisNumber):
                                for b4 in range(basisNumber):
        
                                    #loop over primative guassians twice more
                                    for p3 in range(len(basis["alphas"][b3])):
                                        print("BBBBBBBBBBb")
                                        print(basisNumber)
                                        for p4 in range(len(basis["alphas"][b4])):
                                            
                                            #init 2nd integral constants
                                            C2 = self.constants(basis, b3, b4, p3, p4)
        
                                            #amplitude of 2nd level basis
                                            A34 = C2["overlap"] * C2["c12"] * C2["N12"]
                                            
                                            p = C["p"] + C2["p"]
                                            m = C["p"] * C2["p"]
        
                                            R = 0.0
                                            for dim in range(3):
                                                R += (C["P"][dim] - C["P"][dim]) ** 2
        
                                            alpha =  m / p
                                                                        
                                            const = (2.0 * ( math.pi ** 2.5 )) / ( m * math.sqrt(p))  
        
                                            G[b1][b2][b3][b4] += A12 * A34 * self.boys(0, alpha * R) * const  
            return G
                                
