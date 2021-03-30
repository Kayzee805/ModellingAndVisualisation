from Poisson import poisson
import numpy as np
import matplotlib.pyplot as plt
import time

def monopole(n,method,epsilon):
    model = poisson(n,epsilon)
    model.setPointCharge3D()
    model.setBoundaries3D(model.lattice)
    if(method=="gauss"):
        model.gaussSeidel_CheckerBoard()
        #model.gaussSeidelUpdate_roll()
        model.generateElectricData()
    else:
        model.jacobiUpdate()
        model.generateElectricData()

    potential = model.getPotential()
    im=plt.imshow(potential)
    plt.colorbar(im)
    plt.title(f"Potential plot for monopole size={n}")
    plt.savefig(f"figures/ElectricField/monopole_{n}.png")
    plt.show()
# np.savetxt(f"monopole_{method}.dat",potential)
    #model.plotEField()

def chargedWire(n,method,epsilon):
    model=poisson(n,epsilon)
    model.setChargedWire3D()
    model.setBoundaries3D(model.lattice)
    if(method=="gauss"):
        print("starting gauss")
        model.gaussSeidel_CheckerBoard()
        #model.gaussSeidelUpdate_roll()
       # model.gaussSeidelUpdateOriginal()
        model.generateMagneticData()

    else:
        model.jacobiUpdate()
        model.generateMagneticData()

    potential = model.getPotential()
    im=plt.imshow(potential,cmap="magma")
    plt.colorbar(im)
  # // plt.savefig(f"monopole_{method}.png")
    plt.title(f"Potential plot for charged wire size={n}")
    plt.savefig(f"figures/magneticField/potentialChargedWire_{n}.png")
    plt.show()

def sor(n,epsilon):
    model = poisson(n,epsilon)
    model.generate_SOR()
    #model.generate_SOR3D()
    #model.generate_SOR_Point()
  #  model.generate_SOR_Point3D()

def getFit(n):
    print("Seperate fit\n\n")
    array = np.loadtxt(f"data/potentialDataVR_{n}.dat",dtype=float)
    distance = np.log2(array[:,0])
    potential = np.log2(array[:,1])
    efield = np.log2(array[:,2])
    print(f"distance = {len(distance)}  potetial={len(potential)} efield={len(efield)}")
    newDistance = distance[(distance>0.5) & (distance<1.6)]
    newPotential = potential[(distance>0.5) & (distance<1.6)]
    xfit,xin = np.polyfit(newDistance,newPotential,1)
    print(f"Potential Fit = xfit={xfit}  and xin={xin}")
    print(f"Using potential {newPotential[0]} {array[0,2]}")

    newDistance = distance[(distance>0.5) & (distance<2.1)]
    newElectric = efield[(distance>0.5) & (distance<2.1)]
    xfit,xin = np.polyfit(newDistance,newElectric,1)
    print(f"ElectricField  Fit = xfit={xfit}  and xin={xin}")

    print("Seperate fit\n\n")
    array = np.loadtxt(f"data/potentialDataVRMagnetic_{n}.dat",dtype=float)
    distance = np.log2(array[:,0])
    potential = (array[:,1])
    efield = np.log2(array[:,2])
 
    print(f"distance = {len(distance)}  potetial={len(potential)} magnetic field={len(efield)}")
    newDistance = distance[(distance>0.5) & (distance<1.5)]
    newPotential = potential[(distance>0.5) & (distance<1.5)]
    xfit,xin = np.polyfit(newDistance,newPotential,1)

    print(f"Potential Fit Magnetic = xfit={xfit}  and xin={xin}")

    newDistance = distance[(distance>1) & (distance<2.5)]
    newElectric = efield[(distance>1) & (distance<2.5)]
    xfit,xin = np.polyfit(newDistance,newElectric,1)
    print(f"magnetic  Fit = xfit={xfit}  and xin={xin}")

def main():
    n =100
    epsilon=0.001
    method="gauss"
  #  monopole(n,method,epsilon)
   # chargedWire(n,method,epsilon)
   # sor(n,epsilon)
    getFit(n)

t1=time.time()
main()
t2=time.time()
print(f"Time taken ={t2-t1}")