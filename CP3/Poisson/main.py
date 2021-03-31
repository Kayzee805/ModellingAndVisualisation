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
    plt.savefig(f"figures/ElectricField/size{n}/monopole_{n}.png")
    plt.show()
# np.savetxt(f"monopole_{method}.dat",potential)
    #model.plotEField()

def chargedWire(n,method,epsilon):
    model=poisson(n,epsilon)
    model.setChargedWire3D()
    model.setBoundaries3D(model.lattice)
    if(method=="gauss"):
        print("starting gauss")
       # model.gaussSeidel_CheckerBoard()
        model.gaussSeidelUpdate_roll()
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
    plt.savefig(f"figures/magneticField/size{n}/potentialChargedWire_{n}.png")
    plt.show()

def sor(n,epsilon):
    model = poisson(n,epsilon)
    model.generate_SOR()
    #model.generate_SOR3D()
    #model.generate_SOR_Point()
  #  model.generate_SOR_Point3D()



def getElectricFit(n,plotGraphs):
    if(n==50):
        distanceMax = 1.7
        distanceMin = 0.7
        electricMax = 2.5
        electricMin=1.5
    else:
        distanceMax = 1.7
        distanceMin = 0.7
        electricMax = 2.5
        electricMin=1.5
    array = np.loadtxt(f"data/electricField/potentialDataVR_{n}.dat")
    distance = np.log(array[:,0])
    potential =np.log(array[:,1])
    field = np.log(array[:,2])

    copyDistance = np.copy(distance)
    distance = distance[(copyDistance>0)]
    potential = potential[(copyDistance>0)]
    field = field[(copyDistance>0)]
    #distance vs potential
    newDistance = distance[(distance>distanceMin) &(distance<distanceMax)]
    newPotential = potential[(distance>distanceMin) &(distance<distanceMax)]
    xfit,xin = np.polyfit(newDistance,newPotential,1)

    print(f"monopole potential Fit = {xfit} for n={n}")
    if(plotGraphs):
        plt.scatter(distance,(potential),s=5)
        plt.plot(distance,xfit*distance+xin,color='r',linewidth=0.4)
        plt.xlabel("log_(distance)")
        plt.ylabel("log_(potential)")
        plt.title(f"Log plot of distance v potential for monopole n={n}")
        plt.savefig(f"figures/electricField/size{n}/distanceVsPotential_{n}.png")
        plt.show()

    #electricFit
    newDistance = distance[(distance>distanceMin) &(distance<distanceMax)]
    newField = field[(distance>distanceMin) &(distance<distanceMax)]
    xfit,xin = np.polyfit(newDistance,newField,1)
    print(f"electric fit = {xfit} for n={n}")
    if(plotGraphs):
        plt.scatter(distance,(field),s=5)
        plt.plot(distance,xfit*distance+xin,color='r',linewidth=0.4)
        plt.xlabel("log_(distance)")
        plt.ylabel("log_(electricField)")
        plt.title(f"Log plot of distance v electricField for monopole n={n}")
        plt.savefig(f"figures/electricField/size{n}/distanceVsElectricField_{n}.png")
        plt.show()

def getMagnetic(n,plotGraphs):
    if(n==50):
        distanceMax = 1.7
        distanceMin = 0.7
        magneticMax = 2.5
        magneticMin=1.5
    else:
        distanceMax = 1.7
        distanceMin = 0.7
        magneticMax = 2.4
        magneticMin=1.5
    array = np.loadtxt(f"data/magneticField/potentialDataVR_{n}.dat")
    distance = np.log(array[:,0])
    potential =(array[:,1])
    field = np.log(array[:,2])
    copyDistance = np.copy(distance)

    distance = distance[(copyDistance>0)]
    potential = potential[(copyDistance>0)]
    field = field[(copyDistance>0)]
    for x in distance:
        if(x==0):
            print("TESTINGSGDSGF")
    #distance vs potential
    newDistance = distance[(distance>distanceMin) &(distance<distanceMax) &(distance>0)]
    newPotential = potential[(distance>distanceMin) &(distance<distanceMax) & (distance>0)]
    xfit,xin = np.polyfit(newDistance,newPotential,1)

    print(f"Charged wire potential Fit = {xfit} for n={n}")
    if(plotGraphs):
        plt.scatter(distance,(potential),s=5)
        plt.plot(distance,xfit*distance+xin,color='r',linewidth=0.4)
        plt.xlabel("log_(distance)")
        plt.ylabel("Potential")
        plt.title(f"Log linear plot of distance v potential for chargedWire n={n}")
        plt.savefig(f"figures/magneticField/size{n}/distanceVsPotentialWire_{n}.png")
        plt.show()
    #distance vs field
    newDistance = distance[(distance>magneticMin) &(distance<magneticMax)]
    newField = field[(distance>magneticMin) &(distance<magneticMax)]
    xfit,xin = np.polyfit(newDistance,newField,1)

    print(f"magnetic field Fit = {xfit} for n={n}")

    if(plotGraphs):
        plt.scatter(distance,(field),s=5)
        plt.plot(distance,xfit*distance+xin,color='r',linewidth=0.4)
        plt.xlabel("log_(distance)")
        plt.ylabel("log_(magneticField)")
        plt.title(f"Log plot of distance v magneticField for chargedWire n={n}")
        plt.savefig(f"figures/magneticField/size{n}/distanceVsMagneticWire_{n}.png")
        plt.show()







def main():
    n =100
    epsilon=0.001
    method="gauss"
   # monopole(n,method,epsilon)
  #  chargedWire(n,method,epsilon)
   # sor(n,epsilon)
   # getFit(n)
    plotGraphs=True
    getElectricFit(n,plotGraphs)
    getMagnetic(n,plotGraphs)
t1=time.time()
main()
t2=time.time()
print(f"Time taken ={t2-t1}")

-2.025924129840355


'''
electric field
50x50= potential = [0.5,2.4]
100x100 potential = [0.7.1.2]

50x50 field = [1.5,2.5]
100x100 field=[1.5,2.5]

'''