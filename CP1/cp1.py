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
            glauber.update()
            if(n%500==0):
                print(n)
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
        susceptibility[t] = glauber.calculuateSusceptibility(np.var(magnetisation))
        specificHeat[t] = glauber.calculateHeatCapacity(np.var(energy))
        specificHeatError[t] = glauber.jacknife(energy,specificHeat[t])
        print(f"t ={t+1}  Time taken ={time.time()-start}   Specific heat = {specificHeat[t]}")

    #now save each of the arrays as a text file
    #with temp as the initial column
    np.savetxt('data/Glauber_AverageMagnetisation.dat',(tempList,averageMagnetisation))
    np.savetxt('data/Glauber_AverageTotalEnergy.dat',(tempList,averageEnergy))
    np.savetxt('data/Glauber_SpecificHeat.dat',(tempList,specificHeat))
    np.savetxt('data/Glauber_Susceptibility.dat',(tempList,susceptibility))
    np.savetxt('data/Glauber_SpecificHeatWithError.dat',(tempList,specificHeat,specificHeatError))

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
        start = time.time()
        for n in range(nSteps):

            kawasaki.update()
            if(n%500==0):
                print(n)
            if(n%10==0 and n>=100):
                energy.append(kawasaki.totalEnergy())

        print(f"Time taken for T={t+1}  is {time.time()-start}s   {nSteps}")
        averageEnergy[t] = np.mean(energy)
        specificHeat[t] = kawasaki.calculateHeatCapacity(np.var(energy))
        specificHeatError[t] = kawasaki.jacknife(energy,specificHeat[t])
        print(f"Specific heat at T={tempList[t]} = {specificHeat[t]}")
    np.savetxt('data/Kawasaki_AverageTotalEnergy.dat',(tempList,averageEnergy))
    
    np.savetxt('data/Kawasaki_SpecificHeat.dat',(tempList,specificHeat))
    np.savetxt('data/Kawasaki_SpecificHeatWithError.dat',(tempList,specificHeat,specificHeatError))

    plt.plot(tempList,averageEnergy)
    plt.xlabel("Temp")
    plt.ylabel("Energy")
    plt.title("Temp vs C")
    plt.show()

    plt.plot(tempList,specificHeat)
    plt.errorbar(tempList,specificHeat,linewidth=0.9,yerr=specificHeatError)
    plt.xlabel("Temp")
    plt.ylabel("Specific heat")
    plt.title("Temp vs C")
    plt.show()



def plotGraphs():
    array = np.loadtxt("data/allData.dat")  
    array = np.transpose(array)
    temp = array[:,0]
    gMagnetisation = array[:,1]
    gEnergy = array[:,2]
    gHeatCapacity = array[:,3]
    gHeatError = array[:,4]
    gSusceptibility = array[:,5]
    kEnergy = array[:,6]
    kHeatCapacity = array[:,7]
    kHeatError = array[:,8]

    #Temperature and Glauber Magnetisation
    plt.plot(temp,gMagnetisation)
    plt.scatter(temp,gMagnetisation,c='r',marker='x')
    plt.xlabel("Temperature")
    plt.ylabel("Magnetisation")
    plt.title("Temperature against Glauber Magnetisation")
    plt.savefig("figures/Glauber_AverageMagnetisation_Figure.png")
    plt.show()

    #Temperature and Glauber Total Energy
    plt.plot(temp,gEnergy)
    plt.scatter(temp,gEnergy,c='r',marker='x')
    plt.xlabel("Temperature")
    plt.ylabel("Total Energy")
    plt.title("Temperature against Glauber Energy")
    plt.savefig("figures/Glauber_AverageTotalEnergy_Figure.png")
    plt.show()

    #Temperature and Glauber Heat Capacity
    plt.plot(temp,gHeatCapacity,linewidth= 0.3,color='b')
    plt.scatter(temp,gHeatCapacity,c='r',marker='x')
    plt.errorbar(temp,gHeatCapacity,yerr=gHeatError,ecolor='k')
    plt.xlabel("Temperature")
    plt.ylabel("Specific Heat Capacity")
    plt.title("Temperature against Glauber Specific Heat Capacity with Error")
    plt.savefig("figures/Glauber_SpecificHeat_Figure.png")
    plt.show()

    #Temperature and Glauber Susceptibility
    plt.plot(temp,gSusceptibility)
    plt.scatter(temp,gSusceptibility,c='r',marker='x')
    plt.xlabel("Temperature")
    plt.ylabel("Susceptibility")
    plt.title("Temperature against Glauber Susceptibility")
    plt.savefig("figures/Glauber_Susceptibility_Figure.png")
    plt.show()

    #Temperature and Kawasaki energy
    plt.plot(temp,kEnergy)
    plt.scatter(temp,kEnergy,c='r',marker='x')
    plt.xlabel("Temperature")
    plt.ylabel("Total Energy")
    plt.title("Temperature against Kawasaki Energy")
    plt.savefig("figures/Kawasaki_AverageTotalEnergy_Figure.png")
    plt.show()


    #Temperature and Kawasaki Heat Capacity
    plt.plot(temp,kHeatCapacity,linewidth=0.3,color='b')
    plt.scatter(temp,kHeatCapacity,c='r',marker='x')
    plt.errorbar(temp,kHeatCapacity,yerr=kHeatError,ecolor='k')
    plt.xlabel("Temperature")
    plt.ylabel("Specific Heat Capacity")
    plt.title("Temperature against Kawasaki Specific Heat Capacity with Error")
    plt.savefig("figures/Kawasaki_SpecificHeat_Figure.png")
    plt.show()

    print("Done plotting")

def plotGraphs2():
    fileNames = np.loadtxt('data/plotNames.dat',dtype='str')
    #I can use this method if all files have the same amount of columns
    #cannot use this for the error ones :( )#
    '''
    0 Temperature
    1 GlauberMagnetisation
    2 GlauberEnergy
    3 Glauber HeatCapacity
    4 Galuber Heatcap Error
    5 Glauber Susceptibility
    6 Kawasaki energy
    7 kawasaki heatCapacity
    8 Kawasaki heat cap error
    '''
    array = np.zeros((7))
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
    print(f"Nsteps = {nSteps}")
    dynamics = int(input("0=GenerateGlauber\n1=Generatekawasaki\n2=animate"))
    if(dynamics==0):
        generateGlauberData(lx,kT,nSteps)
    elif(dynamics==1):
        generateKawasakiData(lx,kT,nSteps)
    elif(dynamics==2):
        animate(lx,kT,nSteps)
    else:
        print("Just going to plot")
    #plotGraphs2()
    plotGraphs()
    t2=time.time()
    print(f"Time taken for everything= {t2-t1}")

    

