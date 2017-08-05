#basis implementation of a hamiltonian constructor for a system of M nucli and N electrons

import numpy as np
import math

#############################################################
#init important variables
#nucli coordinates
R = [[0,0,0]]
#electron coordinates
r = [[1,0,0]]

#############################################################
#init important functions
def dist(x1, x2):
    #returns distance between two pointts in 3D space
 
    d = 0.0	

    for dim in range(3):
	d+= math.abs((x1[dim] - x2[dim])**2)

    return math.sqrt(d)

#############################################################
def gaussian(alpha, r):
    #1s orbital basis guassian

    return 0
#############################################################
def overlap():
    #returns overlap matrix of basis guassians
    
    S = []    
    
    return S

#############################################################
#main code goes here
