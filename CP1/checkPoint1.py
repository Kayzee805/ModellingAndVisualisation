import matplotlib
matplotlib.use('TKAgg')

import sys
import math
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

J=1.0
nstep=10000

#input

if(len(sys.argv) != 3):
    print("Usage python ising.animation.py N T")
    sys.exit()



lx=int(sys.argv[1]) 
ly=lx 
kT=float(sys.argv[2]) 

spin=np.zeros((lx,ly),dtype=float)

#initialise spins randomly

for i in range(lx):
    for j in range(ly):
        r=random.random()
        if(r<0.5): spin[i,j]=-1
        if(r>=0.5): spin[i,j]=1

fig = plt.figure()
im=plt.imshow(spin, animated=True)
#print(spin)

#update loop here - for Glauber dynamics
def delta_e(lattice, i,j):
    si = lattice[i,j]
    top = lattice[i,(j-1)%ly]
    bottom = lattice[i,(j+1)%ly]
    left = lattice[(i-1)%lx,j]
    right = lattice[(i+1)%lx,j]
    e = 2*J*si*(top+bottom+left+right)
    return e

def calculateProbability(energyDiff):
    kb=1
    T=1 #for now
    p = np.exp((-energyDiff)/(kb*kT))
    return p

for n in range(nstep):
    for i in range(lx):
        for j in range(ly):

#select spin randomly
            itrial=np.random.randint(0,lx)
            jtrial=np.random.randint(0,ly)
            spin_new=-spin[itrial,jtrial]

           

#compute delta E eg via function (account for periodic BC)
            #change it for now then check if I want to accept it or not
            #spin[itrial,jtrial] =spin_new
            energyD =delta_e(spin,itrial,jtrial) 
           # print(f"Change in E = {energyD} at iteration = {n}")
            if(energyD>0):
                prob = calculateProbability(energyD)
                randomNumber = np.random.rand(1)
                if(prob>randomNumber[0]):
                    spin[itrial,jtrial] = spin_new
                 #   print("FLIPPING")
            else:
                spin[itrial,jtrial]=spin_new
#perform metropolis test
                
#occasionally plot or update measurements, eg every 10 sweeps
    if(n%10==0): 
#       update measurements
        energy=0.0
        for i in range(lx):
            for j in range(ly):
                iup=i+1
                if(i==lx-1):iup=0
                jup=j+1
                if(j==ly-1):jup=0
                energy += -J*spin[i,j]*(spin[iup,j]+spin[i,jup])
        print(f"Iteration = {n}  Energy = {energy}")
#       dump output
        
        f=open('spins.dat','w')
        for i in range(lx):
            for j in range(ly):
                f.write('%d %d %lf\n'%(i,j,spin[i,j]))
              #  print('%d %d %lf\n'%(i,j,spin[i,j]))

        f.close()
#       show animation
        plt.cla()
        im=plt.imshow(spin, animated=True,vmin=-1,vmax=1)
        plt.draw()
        plt.pause(0.0001)
   # print(f"Iteration = {n}")

plt.show()