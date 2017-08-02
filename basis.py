#basis functions for python implementation of Hartree-Fock Method

import numpy as np
import math

#################################
class Basis:
	
	#guassian constants used for each atom, Atomic number - 1 = index coresspoding to that atom
	constant = {
		
		#guassian alpha values
		"alphas":[[0.168856, 0.623913, 3.42525]],
		
		#contraction coeffs for each guassian
		"coeffs":[[0.444635, 0.535328, 0.154329]]

}
#################################
	def buildBasis(self, Z, R):
		#builds basis based on atomic number and coordinates
		
		#init alphas and contraction coeffs holders
		A = []
		CC = []
		
		#init alphas and contraction coeffs
		a = self.constant["alphas"]
		cc = self.constant["coeffs"]		
		
		print(Z[0])		
	
		#for each atom 
		for atom in range(len(Z)):		
			#for each guassian used to represent the atom
			for gaussian in range(len(a)):
				#append all gaussian data to lists
				A[atom].append(a[Z[atom]-1])
				CC[atom].append(cc[Z[atom]-1])
		
		return { "alphas":A, "coeffs": CC, "mus":R }
				
################################
	def overlap(self, basis, b1, b2, p1, p2):
		#implements guassian product theory to calculate Hartree-Fock integrals
		
		#get important guassian data
		a1 = basis["alphas"][b1][p1]
		a2 = basis["alphas"][b2][p2]

		mu1 = basis["mus"][b1][0]
		mu2 = basis["mus"][b2][0]
		
		#new guassian alpha
		A = a + b
		
		#find constant for new guassian
		#divide equation for constant into three parts
	#	a12 = a1 * a2 
	#	c1 = (2 * a12) / (( p * math.pi ) ** ( 3/4 ) ) 
	#	c2 = math.exp( ( -a12 / (a12 * (abs(mu1-mu2)**2)))		
	#	C = a1*a2

		#center of new guassian
	#	Mu = (a1*mu1 + a2*mu2) / A
		
		q = (a1*a2)/A
		Q = abs(mu1 - mu2)**2

		return math.exp(-qQ)
		
################################
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
						S[b1][b2] += basis["coeffs"][b1][p1]*basis["coeffs"][b2][p2]  * self.overlap(basis, b1, b2, p1, p2 )

		
		return S				
		


