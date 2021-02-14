import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import colors
import matplotlib.patches as mpatches
from scipy.stats import sem
import time
import random
import sys
matplotlib.use('TKAgg')


array1 = np.transpose(np.loadtxt("data/one.dat"))
array2 = np.transpose(np.loadtxt("data/two.dat"))
array3 =np.transpose( np.loadtxt("data/three.dat"))
array4 = np.transpose(np.loadtxt("data/four.dat"))
array5 = np.transpose(np.loadtxt("data/four.dat"))



finalArray=[]
errors=[]
pImmune = np.linspace(0,1,101)
print(f"Lengths= {len(array1)},{len(array2)},{len(array3)},{len(array4)},{len(array5)}")
for i in range(len(array1)):
    temp = [array1[i],array2[i],array3[i],array4[i],array5[i]]
    finalArray.append(np.mean(temp))
    errors.append(sem(temp))

plt.scatter(pImmune,finalArray,s=10,color='k')
plt.plot(pImmune,finalArray,color='b')
plt.errorbar(pImmune,finalArray,yerr=errors,ecolor='r')

plt.xlabel("immune")
plt.ylabel("Infection")
plt.title("immune vs infection")
plt.savefig("figures/immune.png")
plt.show()