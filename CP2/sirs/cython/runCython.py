import test
import logging
import threading
import time
import numpy as np
size=50
sweeps=10000
pS=0.8
pI=0.1
pR=0.01
#100    0.8  0.1    0.001
#test.main(size,sweeps,pS,pI,pR,runAnim=False,genData=False,task5DO=True)



name = input("Enter file name")

array2 = test.task5(size,sweeps,name)
array2 = np.array(array2)
np.savetxt("data/"+name+".dat",np.transpose(array2),fmt='%.7f')