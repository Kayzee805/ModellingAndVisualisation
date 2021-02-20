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

array = np.loadtxt("data/task3ProcessedData.dat")
p1s=array[:,0]
p3s=array[:,1]
infected = array[:,2]
variance = array[:,3]

print(f"Lengh.  p1s={p1s.size}  p2s={p3s.size} infected={len(infected)}")
p1s = p1s.reshape((21,21))
p3s=p3s.reshape((21,21))
infected = infected.reshape((21,21))
variance = variance.reshape((21,21))

plt.figure()
#  CS=plt.contour(p1s,p3s,infected)
CS=plt.contourf(p1s,p3s,infected,cmap='magma')
#plt.clabel(CS,fontsize=8,colors='k')
cbar=plt.colorbar(CS)
plt.xticks(np.linspace(0,1,11))
plt.yticks(np.linspace(0,1,11))

plt.xlabel("P1")
plt.ylabel("P3")
plt.title("Average infected of p1-p3 plane")
plt.savefig("figures/infectedContourv2.png")
plt.show()

#plt.cla()

plt.figure()
#  CS=plt.contour(p1s,p3s,infected)
CS=plt.contourf(p1s,p3s,variance,cmap='magma')
#plt.clabel(CS,fontsize=8,colors='k')
cbar=plt.colorbar(CS)
plt.xticks(np.linspace(0,1,11))
plt.yticks(np.linspace(0,1,11))

plt.xlabel("P1")
plt.ylabel("P3")
plt.title("Scaled variance of p1-p3 plane")
plt.savefig("figures/varianceContourv2.png")
plt.show()

print("stop here")

# array1 = np.transpose(np.loadtxt("data/100/1v2.dat"))
# array2 = np.transpose(np.loadtxt("data/100/2v2.dat"))
# array3 =np.transpose( np.loadtxt("data/100/3v2.dat"))
# array4 = np.transpose(np.loadtxt("data/100/4v2.dat"))
# array5 = np.transpose(np.loadtxt("data/100/5v2.dat"))

array1 = np.transpose(np.loadtxt("data/1v2.dat"))
array2 = np.transpose(np.loadtxt("data/2v2.dat"))
array3 =np.transpose( np.loadtxt("data/3v2.dat"))
array4 = np.transpose(np.loadtxt("data/4v2.dat"))
array5 = np.transpose(np.loadtxt("data/5v2.dat"))



finalArray=[]
errors=[]
pImmune = np.linspace(0,1,101)
allArray= np.zeros((101,5))
print(f"Lengths= {len(array1)},{len(array2)},{len(array3)},{len(array4)},{len(array5)}")
for i in range(len(array1)):
    temp = [array1[i],array2[i],array3[i],array4[i],array5[i]]
    allArray[i,0]=array1[i]
    allArray[i,1]=array2[i]
    allArray[i,2]=array3[i]
    allArray[i,3]=array4[i]
    allArray[i,4]=array5[i]
    finalArray.append(np.mean(temp))
    errors.append(sem(temp))

combined=np.array((pImmune,finalArray,errors))
np.savetxt('data/task5ProcessedData.dat',np.transpose(combined),fmt='%.6f')
np.savetxt('data/task5RawData.dat',np.transpose(combined),fmt='%.6f')

plt.scatter(pImmune,finalArray,s=10,color='k')
plt.plot(pImmune,finalArray,color='b')
plt.errorbar(pImmune,finalArray,yerr=errors,ecolor='r')

plt.xlabel("immune")
plt.ylabel("Infection")
plt.title("immune vs infection fo 100x100")
plt.savefig("figures/immune50v2.png")
plt.show()