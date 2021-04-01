import numpy as np 
from cahnHilliard import cahnHilliard
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from tqdm import tqdm
import sys,time
def animate(size,sweeps,phi0,dx,dt,noise):
    '''
    Animation for the cahn hilliard update equation. Animating every 100th
    frame.
    '''
    #initialise a cahnHilliard object
    
    model = cahnHilliard(size,phi0,dx,dt,noise)
    colour='magma'
    plt.figure()
    im=plt.imshow(model.lattice,animated=True,cmap=colour)
    plt.colorbar()
    plt.clim(-1,1)
    plt.draw()
    plt.pause(0.001)
    
    allSweeps = np.linspace(0,sweeps,sweeps+1,dtype=int)
    for i in (allSweeps):
        #call the update method
        model.update()

        #then plot every 100th frame to the animation
        if(i%100==0):
            plt.cla()
            im=plt.imshow(model.lattice,animated=True,cmap=colour)
            plt.draw()
            plt.pause(0.0001)

            #print every 1000th free energy
            currentEnergy = model.calculateFreeEnergy()
            if(i%1000==0):
                print(f"Free energy ={currentEnergy} at {i}")


def generateFreeEnergy(size,sweeps,phi0,dx,dt,noise):
    '''
    Generate free energy of a system and save the free energy vs iteration data.
    '''

    #intiliase a cahnHilliard object
    model = cahnHilliard(size,phi0,dx,dt,noise)
    freeEnergy=[]    
    times=[]
    allSweeps = np.linspace(0,sweeps,sweeps+1,dtype=int)

    #then update it allSweeps times
    for i in tqdm(allSweeps):
        #call the update method of the object
        model.update()

        #calculate every Free energy and save the time/iteration number and the free energy to a list
        if(i%1==0):
            currentEnergy = model.calculateFreeEnergy()
            freeEnergy.append(currentEnergy)
          #  print(f"Free energy ={currentEnergy} at {i}")
            times.append(i)
        
    #Save the data to a data file
    combined = np.array((freeEnergy,times))
    np.savetxt(f"data/freeEnergy{phi0}_{size}.dat",np.transpose(combined),fmt='%.4f')


    # plt.plot(times,freeEnergy,color='b')
    # plt.scatter(times,freeEnergy,color='k',s=1)

    # plt.xlabel("Time (iterations)")
    # plt.ylabel("Free energy")
    # plt.title(f"Plot of freeEnergy vs time. for phi0={phi0} n={size}")
    # plt.savefig(f"freeEnergy{phi0}_{size}.png")
    # plt.show()

def plotGraphs(size,phi0):
    array = np.loadtxt(f"data/freeEnergy{phi0}_{size}.dat")
    iterations = array[:,1]
    freeEnergy = array[:,0]

   # plt.plot(iterations,freeEnergy,color='b')
    plt.scatter(iterations[::10],freeEnergy[::10],color='k',s=1)
    plt.xlabel("Time (iterations)")
    plt.ylabel("Free energy")
    plt.title(f"Plot of free energy vs time for phi0={phi0} and size={size}")
    plt.savefig(f"figures/freeEnergy{phi0}_{size}.png")
    plt.show()

if __name__=="__main__":
    if(len(sys.argv)!=3):
        print("Usage python main.py N phi0")
        exit(0)
    

    size = int(sys.argv[1])
    phi0 = float(sys.argv[2])

    #these values are preset but can be changed here. 
    sweeps = 100000
    dx=1
    dt=1
    noise = 0.1
    task = int(input("Input 0: Animate\nInput 1: to generate Free energy data\nInput 2: to Plot free energy Data\nDefault = 2\n\nInput= "))
    t1=time.time()
    if(task==0):
        animate(size,sweeps,phi0,dx,dt,noise=noise)
    elif(task==1):
        generateFreeEnergy(size,sweeps,phi0,dx,dt,noise)
        plotGraphs(size,phi0)

    else:
        plotGraphs(size,phi0)
    t2=time.time()
    print(f"Time taken to run {t2-t1}s")