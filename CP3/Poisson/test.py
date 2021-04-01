import numpy as np 

array = np.loadtxt("data/SOR/SOR_100_0.0001.dat")

np.savetxt("data/SOR/SOR_100_0.0001.dat",np.transpose(array),fmt='%.4f')
