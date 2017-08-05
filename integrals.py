#integral calculator for python implementation of the Hartree-Fock Method
#referece material refers to textbook "Modern Quantom Chemistry Introduction to Advanced Electronic Structure Theory" by 

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
	
                C["overlap"] = 1

                for dim in range(3):
                    C["overlap"] *= math.exp(-C["q"]*C["Q"][dim])#[ math.exp(-C["q"]*C["Q"][dim]) for dim in range(3) ] 

                #calculate 3D analytical integral
                C["integrand"] = math.sqrt( math.pi / (C["p"]) ) ** 3
               
                print("-----------")
                print(C["overlap"])

                #calculate product of contraction coeffs and normalization constants
                C["c12"] = C["c1"] * C["c2"]
                C["N12"] = C["N1"] * C["N2"]

                print("GGGGGGGGGGGGG")
                print(C["overlap"])

		return C

#################################
	def buildKE(self, basis):
		#builds KE matrix T

		#number that represents the number of basis functions being used for the system
		basisNumber = len(basis["alphas"])

		#init empty overlap matrix
		T  = np.zeros([basisNumber, basisNumber])
		
		#iterate over atom basis twice
		for b1 in range(basisNumber):
			for b2 in range(basisNumber):
				
				#iterate over primatives used in basis twice
				for p1 in range(len(basis["alphas"][b1])):
					for p2 in range(len(basis["alphas"][b2])):
						
						#get integral constants
						C = self.constants(basis, b1, b2, p1, p2)						
		                                				
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
		#builds overlap matrix S
		#number that represents the number of basis functions being used for the system
		basisNumber = len(basis["alphas"])

		#init empty overlap matrix
		V = np.zeros([basisNumber, basisNumber])
		
		#iterate over atom basis twice
		for b1 in range(basisNumber):
			for b2 in range(basisNumber):
				
				#iterate over primatives used in basis twice	
				for atom in range(len(system["Z"])):
					for p1 in range(len(basis["alphas"][b1])):
						for p2 in range(len(basis["alphas"][b2])):
	
							#get integral constants
							C = self.constants(basis, b1, b2, p1, p2)						
							t = C["q"]*C["Q"]
							
							P = [ ((C["a1"]*C["mu1"][x] + C["a2"]*C["mu2"][x])/C["p"]) for x in range(3)]
							
																			     
							 
							RPA2 = sum((np.asarray(system["R"][atom]) - np.asarray(P))) * C["p"]
			
							term1 = (2*math.pi)/C["p"] 
							if(RPA2 == 0):
								boys = 0
							else:
								print(RPA2)
								boys = (advMath.gamma(0.5) * advMath.gammainc(0, 0.5)) / (2 * (RPA2**(0.5)))
								print("GGGGGGGGG")
								print(advMath.gammainc(RPA2,0.5))	
					
							V[b1][b2] = -system["Z"][atom] * term1 * boys * C["c12"]											
	
			
		return V
