import numpy as np 
from cahnHilliard import cahnHilliard
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from tqdm import tqdm

def animate(size,sweeps,phi0,dx,dt,noise):
    model = cahnHilliard(size,phi0,dx,dt,noise)
    plt.figure()
    im=plt.imshow(model.lattice,animated=True)
    plt.colorbar
    plt.draw()
    plt.pause(0.001)
    plt.colorbar()
    freeEnergy=[]    
    times=[]
    allSweeps = np.linspace(0,sweeps,sweeps+1,dtype=int)
    for i in (allSweeps):
        model.update()
        if(i%10==0):
            plt.cla()
            im=plt.imshow(model.lattice,animated=True)
            plt.draw()
            plt.pause(0.00001)
            currentEnergy = model.calculateFreeEnergy()
            freeEnergy.append(currentEnergy)
            print(f"Free energy ={currentEnergy} at {i}")
            times.append(i)
        
        
    combined = np.array((freeEnergy,times))
    np.savetxt("freeEnergy0.5.dat",np.transpose(combined),fmt='%.4f')
    plt.show()
    plt.plot(times,freeEnergy,color='b')
    plt.xlabel("Time (iterations)")
    plt.ylabel("Free energy")
    plt.title("Plot of freeEnergy vs time. for phi0=0")
    plt.savefig("freeEnergy0.5.png")
    plt.show()


def generateFreeEnergy(size,sweeps,phi0,dx,dt,noise):
    model = cahnHilliard(size,phi0,dx,dt,noise)
    freeEnergy=[]    
    times=[]
    allSweeps = np.linspace(0,sweeps,sweeps+1,dtype=int)
    for i in tqdm(allSweeps):
        model.update()
        if(i%1==0):
            currentEnergy = model.calculateFreeEnergy()
            freeEnergy.append(currentEnergy)
          #  print(f"Free energy ={currentEnergy} at {i}")
            times.append(i)
        
        
    combined = np.array((freeEnergy,times))
    np.savetxt(f"freeEnergy{phi0}.dat",np.transpose(combined),fmt='%.4f')
    plt.plot(times,freeEnergy,color='b')
    plt.scatter(times,freeEnergy,color='k',s=1)

    plt.xlabel("Time (iterations)")
    plt.ylabel("Free energy")
    plt.title(f"Plot of freeEnergy vs time. for phi0={phi0}")
    plt.savefig(f"freeEnergy{phi0}.png")
    plt.show()

def plotFreeEnergy():
    array = np.loadtxt("freeEnergy.dat")
    energy = array[:,0]
    times = array[:,1]
    plt.plot(times,energy)
    plt.xlabel("Time (iterations)")
    plt.ylabel("Free energy")
    plt.title("Plot of freeEnergy vs time. for phi0=0")
    plt.savefig("freeEnergy.png")
    plt.show()
if __name__=="__main__":
    size = 50
    sweeps = 100000
    dx=1
    dt=1
    noise = 0.1
    phi0=0.5
   # animate(size,sweeps,phi0,dx,dt,noise)
   # generateFreeEnergy(size,sweeps,0,dx,dt,noise)
    generateFreeEnergy(size,sweeps,0,dx,dt,noise)

   # plotFreeEnergy()
