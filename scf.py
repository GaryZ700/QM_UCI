#Self Consistant Field Method, Python Implementation
#operator names with a p mean that the operator is in MO basis, all other operator are in AO basis


import numpy as np

##################################
class SCF:
	
	#dictionary to store SCF data
	SCFData = {
		
		#electron repulsion operator
		"G":[],
		
		#density matrix
		"D":[],

		#fock matrix in AO basis
		"F":[],
		
		#fock matrix in MO basis
		"Fp":[],

		#MO coeffs
		"C":[],
		
		#transformation matrix for AO to MO basis
		"X":[]

}

##################################
	def sort(eig):
		#sorts eigen values and vectors from numpy.linalg.eig output

		#separate eigen values and vectors
		eVal = eig[0]
		eVec = eig[0]

		#get sorting index 
		#sorts from largest to smallest 
		index = eVal.argsort[::-1]
		
		return eVal[idx], eVec[:, idx]

##################################
	def getTransform(S):
		#creates transform operator to move from AO to MO basis
		
		#get electron density matrix
		eVal, D = np.linalg.eig(S)
		
		D = D**(-0.5)
		
		return eVal  * D * eVal.getH()

##################################
	def scf(SCFData):
		#single scf iteration

		#get eigen values and vectors of fock matrix
		#and sort the values and vectors
	        eVal, eVec = self.sort(np.linalg.eig(SCFData["F"]))	

		SCFData["Fp"] = SCFData["X"].getH() * SCFData["F"] * SCFData["X"]
		
		return SCFData
		
##################################
