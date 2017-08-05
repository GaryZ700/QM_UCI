import math

def basis(position, alpha, sigma):
    #all basies consist of a single guassian


    gaus = math.exp((-((position - alpha)**2))/(2*sigma**2))

    return c(sigma)*gaus

def c(sigma):
    return 1/math.sqrt(2*math.pi*sigma**2) 


gaus1 = [4.0,1.0]
gaus2 = [1.0,2.0]

dx = 0.0001

x = -20

y1, y2, yo =0,0,0 

#
#while(x < 20):
#    
#    yo += dx*((basis(x,gaus1[0],gaus1[1])*basis(x,gaus1[0],gaus1[1]))/(c(gaus1[1])*c(gaus1[1])))
#    y2 += dx*basis(x,gaus1[0],gaus1[1])
#    y1 += dx * basis(x,gaus1[0],gaus1[1]) 
#
#    #y+= 2*x*dx
#    print(x)
#        
#    x+=dx
#
#print(yo)
#print(y1)
#print(y2)
#
#z = yo/(y1+y2)
#print("))))))))))") 
#print(z)

ym = -float("inf")
while(x<20):

    test = (basis(x,gaus1[0],gaus1[1])*basis(x,gaus1[0],gaus1[1]))/(c(gaus1[1])*c(gaus1[1]))
    if(test > ym):
        ym = test

    print(x)
    x+=dx

print(ym)
