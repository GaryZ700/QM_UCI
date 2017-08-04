#integral calculator for python implementation of the Hartree-Fock Method

import numpy as np
import math

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

		C["mu1"] = basis["mus"][b1][0]
		C["mu2"] = basis["mus"][b2][0]

		C["c1"] = basis["coeffs"][b1][p1]
		C["c2"] = basis["coeffs"][b2][p2]
		
		#calculate integral values
		C["p"] = C["a1"] + C["a2"]
		C["P"] = C["a1"]*C["a2"]
		C["Pp"] = ( C["a1"]*C["mu1"] + C["a2"]*C["mu2"] ) / C["p"]
	
		C["q"] = C["P"]/C["p"]
		C["Q"] = abs(C["mu1"] - C["mu2"])**2

		C["c12"] = C["c1"] * C["c2"]
		
		C["overlap"] = math.exp(-C["q"]*C["Q"])		

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
						
						#Calculate integral
						c1 = 3 - (2*C["q"]*C["Q"])
						c2 = (math.pi / C["p"]) ** (3/2) 		
						T[b1][b2] += C["c12"] * C["q"] * c1 * c2  * C["overlap"]
						#T[b1][b2] += 3 * C["a2"] * C["overlap"] * C["c12"]
						#T[b1][b2] += -2 * ( ((C["Pp"] - C["mu2"])**2) + (1/(2*C["p"])) ) * C["overlap"] * C["c12"]
					#	term1 = C["a2"] * C["overlap"]
					#	term2 =  -2 * C["a2"]**2 * C["overlap"]
					#	term3 = -0.5 * C["overlap"]
						
					#	T[b1][b2] += term1 + term2 + term3
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
						S[b1][b2] += C["c12"] * C["overlap"]

		
		return S				

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
						S[b1][b2] += C["c12"] * C["overlap"]

		
		return S				
