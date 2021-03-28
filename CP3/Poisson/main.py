import numpy as np 
from Poisson import poisson
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
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
    model.getMagneticFiled()

def SOR(n,epsilon,minW,maxW):

    model = poisson(n,epsilon=epsilon,minW=minW,maxW=maxW)
    model.setPointCharge()

    model.overRelaxationUpdate_all()

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

def genSor(n,epsilon):
    model = poisson(n,epsilon=epsilon)
    model.generate_SOR()


if __name__=='__main__':
    print("starting")
    t1=time.time()
    n=100
    method="jacobi"
    epsilon=0.001
    monopole(n,method,epsilon)
    #wire(n,method,epsilon)
    #SOR(50,epsilon,1,2)
    #plotPotential(method)
    #plotPotentialWire(method)
   # plotPotential()
   # plotPotentialWire()
   # genSor(n,epsilon)
    print(f"finished at {time.time()-t1}")
