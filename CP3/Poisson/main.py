import numpy as np 
from Poisson import poisson
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#epsilon = 0.01 and 0.001

def monopole(n,method,epsilon):
    model = poisson(n,method,epsilon)
    model.setBoundaries()
    model.setPointCharge()

    if(method=="gauss"):
        model.gaussSeidelUpdate()
    else:
        model.jacobiUpdate()
    potential = model.getPotential()
    im=plt.imshow(potential)
    plt.colorbar(im)
    plt.savefig(f"monopole_{method}.png")
    plt.show()
    np.savetxt(f"monopole_{method}.dat",potential)
    model.plotEField()

def wire(n,method,epsilon):
    model = poisson(n,method,epsilon)
    model.setChargedWire()
    model.setBoundaries()

    if(method=="gauss"):
        model.gaussSeidelUpdate()
    else:
        model.jacobiUpdate()
    potential = model.getPotential()
    im=plt.imshow(potential)
    plt.colorbar(im)
    plt.savefig(f"wire_{method}.png")
    plt.show()
    np.savetxt(f"wire_{method}.dat",potential)


def plotPotential(method):
    potential =np.loadtxt(f"monopole_{method}.dat")
    plt.imshow(potential,cmap="magma")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.title("Potential plot for charged wire")
    plt.colorbar()
    plt.savefig(f"monopole_{method}.png")
    plt.show()





def plotPotentialWire(method):
    potential =np.loadtxt(f"wire_{method}.dat")
    plt.imshow(potential,cmap="magma")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.title("Potential plot for charged wire")
    plt.colorbar()
    plt.savefig(f"wire_{method}.png")
    plt.show()



if __name__=='__main__':
    print("starting")
    n=50
    method="jacobi"
    epsilon=0.001
    monopole(n,method,epsilon)
    #wire(n,method,epsilon)
    #plotPotential(method)
    #plotPotentialWire(method)
   # plotPotential()
   # plotPotentialWire()
    print("finished")
