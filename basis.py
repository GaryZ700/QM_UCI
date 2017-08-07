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
		"alphas":[[3.42525091, 0.62391373, 0.16885540]],
		
		#contraction coeffs for each guassian
		"coeffs":[[0.15432897, 0.53532814, 0.44463454]]

}
#################################
	def buildBasis(self, Z, R):
		#builds basis based on atomic number and coordinates
		
		#init alphas, contraction coeffs, and normalization constants holders
    		A = []
                CC = []
                N = []

		#init alphas and contraction coeffs
		a = self.constant["alphas"]
                print(a)
		cc = self.constant["coeffs"]			
	    	
		#for each atom 
		for atom in range(len(Z)):		
			#for each guassian used to represent the atom
			A.append(a[Z[atom]-1])
			CC.append(cc[Z[atom]-1])

                        #https://youtu.be/RHkWFlIhNHo?t=2m4s
                        N.append( [ ((2*alpha)/math.pi) ** (0.75) for alpha in A[atom] ] )

                        print("+++++++++++++")
                        print(N)


		#list of guassian data, alphas,
                #contraction coeffs
                #r: positions
                #N: normalization value
                #to access data, pass string of requsted value, and then internal atom number
                return { "alphas":A, "coeffs": CC, "r":R, "N":N }
				
################################
