#hartree fock master control file

import basis

###############################################
#important functions go here
###############################################
def ui():
#gets atom charges and coordinates from user
    atomicNumbers = []
    coord = []
    charges = []
    option = 0
    
    while(option != 2):
        print("Choose an option. \n")
        print("1. Add Atom")
        print("2. Exit")   
    
        try:
            option = int(raw_input())  
        except:
            print("Invalid input! \n \n")
            continue
        
        if(option == 1):
            data = addAtom()
            
            atomicNumbers.append(data[0])
            coord.append(data[1])
            charges.append(data[2])
       
    return atomicNumbers, coord, charges

###############################################
def addAtom():
#allows user to interactively add atoms
    
    coord = []    

    print("\n \n")
    print("Enter atomic number of atom:")
    
    atomicNumber = int(raw_input())
    
    print("Enter atom charge:")
    
    charge = int(raw_input())

    print("Enter coordinates of atom in x y z format:")
    
    testCoord = raw_input()

    for pos in testCoord:
	if(pos != " "):
            coord.append(pos)
      
    return atomicNumber, coord, charge
     
###############################################3
#main code goes here

#get user input
atomicNumbers, coord, charges = ui()

basis = constructBasis(atomicNumbers, coord)



print(atomicNumbers)
print(coord)
print(charges)
print(basis)
#nuclearRepulsion = calcNuclearRepulsion(atomicNumbers, coord)
#basis = buildBasis(atomicNumbers, coord)
#S = buildOverlap(bais)
#T = buildKinetic(basis)
#nuclearAttraction = calcNuclearAtraction(atomicNumbers, coord, basis)
#electronicRepulsion = calcElectronicRepulsion()
