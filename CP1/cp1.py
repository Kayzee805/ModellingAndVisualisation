import matplotlib
matplotlib.use('TKAgg')

import sys
import math
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from Lattice import Lattice
from Glauber import Glauber
from Kawasaki import Kawasaki
import time

def generateGlauberData(n,T,nSteps):
    TemperatureValues = 21
    #not sure how many values of T I need
    tempList = np.linspace(1,T,TemperatureValues)
    averageMagnetisation =np.zeros(TemperatureValues)
    susceptibility=np.zeros(TemperatureValues)
    specificHeat = np.zeros(TemperatureValues)
    averageEnergy =np.zeros(TemperatureValues)
  
    specificHeatError = np.zeros(TemperatureValues)
    glauber = Glauber(n,0)

    #setting spin/lattice to all up spins
    glauber.allOnes()
    for t in range(TemperatureValues):
        glauber.setTemperature(tempList[t])
        magnetisation = []
        energy = []
        start = time.time()
        for n in range(nSteps):
            #now need to flip?
            start = time.time()
            glauber.update()
            print(f"Time to update: {time.time()-start}")
            if(n%10==0):
                if(n>=100):
                #    print("hello")
                    #this is the wait for equilibration?
                    energy.append(glauber.totalEnergy())
                    #magetisation needs to be the abs value
                    magnetisation.append(glauber.totalMagnetisation())
        #use np.var for varaince
       # print(f"size = {len(magnetisation)}")
        averageMagnetisation[t] = np.mean(magnetisation)
        averageEnergy[t] = np.mean(energy)
        print(f"t ={t+1}  Time taken ={time.time()-start}")
        susceptibility[t] = glauber.calculuateSusceptibility(np.var(magnetisation))
        specificHeat[t] = glauber.calculateHeatCapacity(np.var(energy))
        specificHeatError[t] = glauber.jacknife(energy,specificHeat[t])
    #now save each of the arrays as a text file
    #with temp as the initial column
    np.savetxt('data/Glauber_AverageMagnetisation.dat',(tempList,averageMagnetisation))
    np.savetxt('data/Glauber_AverageTotalEnergy.dat',(tempList,averageEnergy))
    np.savetxt('data/Glauber_SpecificHeat.dat',(tempList,specificHeat))
    np.savetxt('data/Glauber_Susceptibility.dat',(tempList,susceptibility))
    print(f"Finished generating data")

    plt.plot(tempList,specificHeat)
    plt.xlabel("Temp")
    plt.ylabel("Specific heat")
    plt.title("Temp vs C")
    plt.show()

def generateKawasakiData(n,T,nSteps):
    TemperatureValues = 21
    #not sure how many values of T I need
    tempList = np.linspace(1,T,TemperatureValues)
    specificHeat = np.zeros(TemperatureValues)
    averageEnergy =np.zeros(TemperatureValues)
    specificHeatError = np.zeros(TemperatureValues)
    kawasaki = Kawasaki(n,0)
    kawasaki.halfHalf()
    print(f"Magnetisation = {kawasaki.totalMagnetisation()}")
    for t in range(TemperatureValues):
        kawasaki.setTemperature(tempList[t])
        energy=[]

        for n in range(nSteps):
            start=time.time()

            kawasaki.update()#
            print(f"Time taken for update: {time.time()-start}")
            if(n%10==0 and n>=100):
                energy.append(kawasaki.totalEnergy())

        print(f"Time taken for T={t}  is {time.time()-start}s")
        print(f"Energy size = {len(energy)}")
        averageEnergy[t] = np.mean(energy)
        specificHeat[t] = kawasaki.calculateHeatCapacity(np.var(energy))
        specificHeatError[t] = kawasaki.jacknife(energy,specificHeat[t])
    np.savetxt('data/Kawasaki_AverageTotalEnergy.dat',(tempList,averageEnergy))
    np.savetxt('data/Kawasaki_SpecificHeat.dat',(tempList,specificHeat))
    np.savetxt('data/Kawasaki_SpecificHeatWithError.dat',(tempList,specificHeat,specificHeatError))

    plt.plot(tempList,averageEnergy)
    plt.xlabel("Temp")
    plt.ylabel("Specific heat")
    plt.title("Temp vs C")
    plt.show()



def plotGraphs():
    fileNames = np.loadtxt('data/plotNames.dat',dtype='str')
    #I can use this method if all files have the same amount of columns
    #cannot use this for the error ones :( )
    for name in fileNames:
        fileName = name[0]
        xLabel = name[1]
        yLabel = name[2]
        array = np.loadtxt("data/"+fileName)
        #print(array.shape)
        array=np.transpose(array)
        #print(array.shape)
        xVals = array[:,0]
        yVals = array[:,1]
        plt.cla()
        plt.plot(xVals,yVals)
        plt.xlabel(xLabel)
        plt.ylabel(yLabel)
        title = xLabel+" against "+yLabel
        plt.title(title)
        plt.scatter(xVals,yVals,s=8)
        plt.savefig("figures/"+fileName[:len(fileName)-4]+"_Figure.png")
        plt.show()
      #  print(fileName,xVals[0],yVals[0],len(xVals),len(yVals))




def animate(n,T,nSteps):
    
    while True:
        method = input("0 for Glauber and 1 for kawasaki")
        try:
            method = int(method)
        except ValueError:
            print("Must be an Integer")
            continue
        if method==0 or method==1:
            break
        else:
            print("Please enter 0 for Glauber and 1 for Kawasaki")


    if(method==0):
        #loop and update glauber
        d = Glauber(n,T)
    else:
        #loop and update kawasadki
        d = Kawasaki(n,T)

    #now do the updating and animating here
    for n in range(nSteps):
        d.update()
        if(n%10==0): 
#       update measurements
            energy = d.totalEnergy()
            print(f"Iteration = {n}  Energy = {energy} Magnetisation = {d.totalMagnetisation()}")
            f=open('spins.dat','w')
            for i in range(lx):
                for j in range(ly):
                    f.write('%d %d %lf\n'%(i,j,d.spin[i,j]))
                #  print('%d %d %lf\n'%(i,j,spin[i,j]))

        f.close()
#       show animation
        plt.cla()
        im=plt.imshow(d.spin, animated=True,vmin=-1,vmax=1)
        plt.draw()
        plt.pause(0.0001)
   # print(f"Iteration = {n}")

    plt.show()





if __name__ == "__main__":
    if(len(sys.argv) != 3):
        print("Usage python ising.animation.py N T")
        sys.exit()
    t1=time.time()
    print("Started")
    lx = int(sys.argv[1])
    ly = lx
    kT=float(sys.argv[2]) 
    nSteps= 10000
    #animate(lx,kT,nSteps)
    #generateGlauberData(lx,kT,nSteps)
    generateKawasakiData(lx,kT,nSteps)
   # plotGraphs()
    t2=time.time()
    print(f"Time taken for everything= {t2-t1}")

    

