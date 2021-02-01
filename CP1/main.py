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
import matplotlib
matplotlib.use('TKAgg')
#plt.rcParams['figure.figsize'] = (9, 9)

def generateGlauberData(systemSize,Temperature,nSteps):
    '''
    Runs the simulation for glauber dynamic rule 
    then saves the generated data in a file called
    'glauberData.dat' in a folder called data
    Parameters:
    -----------
    systemSize :  Integer
                  Size of the system
    T:  Float
        Temperature of the system, here we generate data so T is max temp
        it will start from T=1 to T=T
    
    nSteps: Integer
            Sweep size per temperature
    '''
    TemperatureValues = 21
    #Creating arrays to store data. Multiple arrays instead of 1 big array

    tempList = np.linspace(1,Temperature,TemperatureValues)
    averageMagnetisation =np.zeros(TemperatureValues)
    susceptibility=np.zeros(TemperatureValues)
    specificHeat = np.zeros(TemperatureValues)
    averageEnergy =np.zeros(TemperatureValues)
    specificHeatError = np.zeros(TemperatureValues)
    
    #create a glauber object, which contains the spin lattice 2d array
    glauber = Glauber(systemSize,1)

    #setting spin/lattice to all up spins
    glauber.allOnes()

    #looping over all temperatures
    for t in range(TemperatureValues):
        #change the temperature of the glauber object
        glauber.setTemperature(tempList[t])

        #reset arrays to empty arrays
        magnetisation = []
        energy = []
        start = time.time()
        
        #loop for the sweep
        for n in range(nSteps):
            
            #for each sweep, update the lattice
            glauber.update()

            #uncomment below to keep track of how far into the sweep we are
            # if(n%1000==0):
            #     print(n)

            #only taking in every 10th measurement after the equilibration wait of 100
            if(n%10==0 and n>=100):
                
                #append the measurements to their respective arrays
                energy.append(glauber.totalEnergy())
                magnetisation.append(glauber.totalMagnetisation())

        #for each temperature after the sweeps, do the calculations and index it to the right array
        averageMagnetisation[t] = np.mean(magnetisation)
        averageEnergy[t] = np.mean(energy)
        susceptibility[t] = glauber.calculuateSusceptibility(glauber.calculateVariance(magnetisation))
        specificHeat[t] = glauber.calculateHeatCapacity(glauber.calculateVariance(energy))
        
        #call jacknife method to calculate the error 
        specificHeatError[t] = glauber.jacknife(energy)
        print(f"t ={tempList[t]}  Time taken ={time.time()-start}   Specific heat = {specificHeat[t]}")


    #after simulation done for all arrays, save it in a big array to write it to a file
    combinedArray = np.array((tempList,averageMagnetisation,averageEnergy,specificHeat,specificHeatError,susceptibility))
    #write a transposed array to make it easier to read.
    np.savetxt('data/glauberData.dat',np.transpose(combinedArray),fmt='%.7f')
    print(f"Finished generating data")

def generateKawasakiData(systemSize,Temperature,nSteps):
    '''
    Runs the simulation for kawasaki dynamic rule 
    then saves the generated data in a file called
    'kawasakiData.dat' in a folder called data
    Parameters:
    -----------
    systemSize :  Integer
                  Size of the system
    Temperature:  Float
        Temperature of the system, here we generate data so T is max temp
        it will start from T=1 to T=T
    
    nSteps: Integer
            Sweep size per temperature
    '''

    TemperatureValues = 21
    #Creating arrays to store data. Multiple arrays instead of 1 big array

    tempList = np.linspace(1,Temperature,TemperatureValues)
    specificHeat = np.zeros(TemperatureValues)
    averageEnergy =np.zeros(TemperatureValues)
    specificHeatError = np.zeros(TemperatureValues)

    #create a Kawasaki object, which contains the spin lattice 2d array
    kawasaki = Kawasaki(systemSize,Temperature)
    #initialise the kawasaki spins to half up and half downs so magnetisation=0
    kawasaki.halfHalf()

    print(f"Magnetisation = {kawasaki.totalMagnetisation()}")

    #looping over all temperatures
    for t in range(TemperatureValues):
        #change the temperature of the kawasaki object

        kawasaki.setTemperature(tempList[t])

        #reset the energy array for each emperature
        energy=[]
        #loop for nSteps sweeps
        start =time.time()
        for n in range(nSteps):
            #for each sweep, update the spin lattice
            kawasaki.update()

            #uncomment below to keep track of how far into the sweep we are
            # if(n%1000==0):
            #     print(n)

            #take every 10th measurement after the equilibriation wait 
            if(n%10==0 and n>100):
                #add each measurement to the energy array
                energy.append(kawasaki.totalEnergy())
         

        #For each temperature calculate the speacific heat and average energy
        averageEnergy[t] = np.mean(energy)
        specificHeat[t] = kawasaki.calculateHeatCapacity(kawasaki.calculateVariance(energy))
        #call jacknife method to calculate the error 
        specificHeatError[t] = kawasaki.jacknife(energy)
        print(f"Time taken for T={tempList[t]}  is {time.time()-start}s   specificHeat = {specificHeat[t]}  error={specificHeatError[t]}")
    
    #after simulation done for all arrays, save it in a big array to write it to a file
    combinedArray = np.array((tempList,averageEnergy,specificHeat,specificHeatError))

    #write a transposed array to make it easier to read
    np.savetxt('data/kawasakiData.dat',np.transpose((combinedArray)),fmt='%.7f')

def plotGraphs():
    '''
    Plots all the necessary plots then save the figures in a folder
    called figures
    '''
    glauberData = np.loadtxt("data/glauberData.dat")  
    temp = glauberData[:,0] 
    gMagnetisation = glauberData[:,1]
    gEnergy = glauberData[:,2]
    gHeatCapacity = glauberData[:,3]
    gHeatError = glauberData[:,4]
    gSusceptibility = glauberData[:,5]

    kawasakiData= np.loadtxt("data/kawasakiData.dat")
    kEnergy = kawasakiData[:,1]
    kHeatCapacity = kawasakiData[:,2]
    kHeatError = kawasakiData[:,3]

    #Temperature and Glauber Magnetisation
    plt.plot(temp,gMagnetisation)
    plt.scatter(temp,gMagnetisation,c='r',marker='x')
    plt.xlabel("Temperature")
    plt.ylabel("Average Absolute Magnetisation")
    plt.title("Temperature against Glauber Magnetisation")
    plt.savefig("figures/Glauber_AverageMagnetisation_Figure.png")
    plt.show()

    #Temperature and Glauber Total Energy
    plt.plot(temp,gEnergy)
    plt.scatter(temp,gEnergy,c='r',marker='x')
    plt.xlabel("Temperature")
    plt.ylabel("Average Total Energy")
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
    plt.ylabel("Average Total Energy")
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



def animate(n,T,nSteps):
    '''
    Parameters:
    -----------
    systemSize: Integer
                dimensions of the system.
    Temperature: float
                 Temperature of the system
    nSteps: Integer
            Number of sweeps for the animation
    '''

    #takes in user input on which dynamic rule to animate
    while True:
        method = input("0 for Glauber animation\n1 for kawasaki animation")
        try:
            method = int(method)
        except ValueError:
            print("Must be an Integer")
            continue
        if method==0 or method==1:
            break
        else:
            print("Please enter 0 for Glauber and 1 for Kawasaki or Ctrl+c to exit")


    #initialising which lattice sub class to use according to user input
    if(method==0):
        d = Glauber(n,T)
    else:
        d = Kawasaki(n,T)

    #initalise the lattice with random spins
    d.randomInit()
    #Running the animation for nSteps sweeps
    for n in range(nSteps):
        #updating the lattice for each sweep
        d.update()

        if(n%10==0): 
        #update measurements every 10th sweep
            energy = d.totalEnergy()
            #print the total magnetisation and the energy of the system at nth sweep
            print(f"Iteration = {n}  Energy = {energy} Magnetisation = {d.totalMagnetisation()}")
            
            #writing the files and using the code provided in the pseudocode
            f=open('spins.dat','w')
            for i in range(lx):
                for j in range(ly):
                    f.write('%d %d %lf\n'%(i,j,d.spin[i,j]))
            

        f.close()
        #show animation for each sweep
        plt.cla()
        im=plt.imshow(d.spin, animated=True,vmin=-1,vmax=1)
        plt.draw()
        plt.pause(0.0001)
    plt.show()





if __name__ == "__main__":

    #Running the system
    if(len(sys.argv) != 3):
        print("Usage python main.py N T")
        sys.exit()
    
    t1=time.time()
    print("Started")
    #storing teh respective arguments to respective variables
    lx = int(sys.argv[1])
    ly = lx
    kT=float(sys.argv[2]) 
    nSteps= 10000
    print(f"Nsteps = {nSteps}")
    
    dynamics = int(input("0 to Generate Glauber data\n1 to Generate kawasaki data\n2 to show animation\n3 to Plot and save figures\nAnything else to exit"))

    if(dynamics==0):
        generateGlauberData(lx,kT,nSteps)
        #plotGraphs()
    elif(dynamics==1):
        generateKawasakiData(lx,kT,nSteps)
        plotGraphs()
    elif(dynamics==2):
        animate(lx,kT,nSteps)
    elif(dynamics==3):
        plotGraphs()
    else:
        print("Exiting")
    t2=time.time()
    print(f"Time taken for everything= {t2-t1}")
    sys.exit()

    

# Time taken for T=1  is 14.966191530227661s   400
# Specific heat at T=1.0 = 0.07939555555555561