import time
import numpy as np 
np.random.seed(10)



x = np.loadtxt("data/potentialDataVR_100.dat")

distance = x[:,0]
potential=x[:,1]
efield = x[:,2]
counter=1

# boolarray = (distance>0.5 & distance<0.5)
print(np.amin(distance),np.amax(distance))
test = distance[(distance>0.5) & (distance<1.5)]
print(np.amin(test),np.amax(test))