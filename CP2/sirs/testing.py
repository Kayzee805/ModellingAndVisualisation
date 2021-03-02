import numpy as np


a = np.zeros(10)
b =np.zeros(10)
c=np.zeros(10)
x = np.array((a,b,c))
np.savetxt("test.dat",np.transpose(x))