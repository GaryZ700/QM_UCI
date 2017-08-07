#basis functions for python implementation of Hartree-Fock Method

import numpy as np
import math

#################################
class Basis:
	
	#guassian constants used for each atom, Atomic number - 1 = index coresspoding to that atom
	constant = {
	
                #Hydrogen uses STO-3G from Basis Set Exchange
                #https://bse.pnl.gov/bse/portal

		#guassian alpha values
	#	"alphas":[[3.42525091, 0.62391373, 0.16885540]],
		
		#contraction coeffs for each guassian
	#	"coeffs":[[0.15432897, 0.53532814, 0.44463454]]

		"alphas": [ [ [18.7311370, 2.8253937, 0.6401217], [0.1612778]  ] ], 
		"coeffs":[ [ [0.03349460, 0.23472695, 0.81375733], [1.0]  ] ]

}
#################################
	def buildBasis(self, Z, R):
		#builds basis based on atomic number and coordinates
		
		#init alphas, contraction coeffs, normalization constants, and guassian center holders
    		A = []
                CC = []
                N = []
		r = []

		#init alphas and contraction coeffs
		a = self.constant["alphas"]
                print(a)
		cc = self.constant["coeffs"]			
	    	
		#for each atom 
		for atom in range(len(Z)):	

			#for each basis used to represent the atom
			for basis in range(len(a[Z[atom]-1])):
					
				#for each guassian used to represent the atom
				A.append(a[Z[atom]-1][basis])
				CC.append(cc[Z[atom]-1][basis])
				r.append(R[atom])

		
               	#https://youtu.be/RHkWFlIhNHo?t=2m4s
		for primative in range(len(A)):
			N.append( [ ((2*alpha)/math.pi) ** (0.75) for alpha in A[primative] ] )
		
		#list of guassian data, alphas,
                #contraction coeffs
                #r: positions
                #N: normalization value
                #to access data, pass string of requsted value, and then internal atom number
                return { "alphas":A, "coeffs": CC, "r":r, "N":N }
				
################################
