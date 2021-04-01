from Poisson import poisson
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
def monopole(n,method,epsilon,checkerBoard=True):
    '''
    Calls the update method to the model and calls the method to
    plot and save data for a monopole charge in 3D.
    Parameters:
    -----------
    n:      Type Integer
            size of the model
    method: String
            either jacobi update or gauss-seidel update, if invalid jacobi is used
    epsilon:Type float
            precision for model to converge
    '''
    model = poisson(n,epsilon)
    model.setPointCharge3D()
    model.setBoundaries3D(model.lattice)

    #carry out updates for a point charge in 3D, monopole charge
    if(method=="gauss"):
        if(checkerBoard):
            model.gaussSeidel_CheckerBoard()
        else:
            model.gaussSeidelUpdate_roll()
    
    else:
        model.jacobiUpdate()

    #generate and calculate the required data files
    model.generateElectricData()

    #plot the potential of the monopole
    potential = model.getPotential()
    im=plt.imshow(potential)
    plt.colorbar(im)
    plt.title(f"Potential plot for monopole size={n}")
    plt.savefig(f"figures/ElectricField/size{n}/monopole_{n}.png")
    plt.show()

def chargedWire(n,method,epsilon,checkerBoard=True):
    '''
    Calls the update method to the model and calls the method to
    plot and save data for a charged wire aligned in the z axis. 
    Parameters:
    -----------
    n:      Type Integer
            size of the model
    method: String
            either jacobi update or gauss-seidel update, if invalid jacobi is used
    epsilon:Type float
            precision for model to converge
    '''
    model=poisson(n,epsilon)
    model.setChargedWire3D()
    model.setBoundaries3D(model.lattice)

    #update the model using specified methods
    if(method=="gauss"):
        if(checkerBoard):
            model.gaussSeidel_CheckerBoard()
        else:
            model.gaussSeidelUpdate_roll()

    else:
        model.jacobiUpdate()
    
    #generate the required data files
    model.generateMagneticData()

    #plot the potential of cut of charged wire
    potential = model.getPotential()
    im=plt.imshow(potential,cmap="magma")
    plt.colorbar(im)
    plt.title(f"Potential plot for charged wire size={n}")
    plt.savefig(f"figures/magneticField/size{n}/potentialChargedWire_{n}.png")
    plt.show()

def sor(n,epsilon,checkerBoard=True):
    '''
    Calls the update method to the model and calls the method to
    plot and save data for a cut of charged wire algined in z axis.
    Parameters:
    -----------
    n:      Type Integer
            size of the model
    epsilon:Type float
            precision for model to converge
    '''
    model = poisson(n,epsilon)
    #update the model for the sepcified method.
    if(checkerBoard):
        model.generate_SOR_CheckerBoard()
    else:
        model.generate_SOR()


def getElectricFit(n,plotGraphs):
    '''
    Plots and calculates the fit of potential and electric field against distance
    for a monopole charge
    '''
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

    print(f"\nmonopole potential Fit = {xfit} for n={n}")
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
    print(f"\nelectric fit = {xfit} for n={n}")
    if(plotGraphs):
        plt.scatter(distance,(field),s=5)
        plt.plot(distance,xfit*distance+xin,color='r',linewidth=0.4)
        plt.xlabel("log_(distance)")
        plt.ylabel("log_(electricField)")
        plt.title(f"Log plot of distance v electricField for monopole n={n}")
        plt.savefig(f"figures/electricField/size{n}/distanceVsElectricField_{n}.png")
        plt.show()

def getMagnetic(n,plotGraphs):
    '''
    Plots and calculates the fit of potential and magnetic field against distance
    for a cut of charged wire algined in the z axis.
    '''
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

    #distance vs potential
    newDistance = distance[(distance>distanceMin) &(distance<distanceMax) &(distance>0)]
    newPotential = potential[(distance>distanceMin) &(distance<distanceMax) & (distance>0)]
    xfit,xin = np.polyfit(newDistance,newPotential,1)

    print(f"\nCharged wire potential Fit = {xfit} for n={n}")
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

    print(f"\nmagnetic field Fit = {xfit} for n={n}\n")

    if(plotGraphs):
        plt.scatter(distance,(field),s=5)
        plt.plot(distance,xfit*distance+xin,color='r',linewidth=0.4)
        plt.xlabel("log_(distance)")
        plt.ylabel("log_(magneticField)")
        plt.title(f"Log plot of distance v magneticField for chargedWire n={n}")
        plt.savefig(f"figures/magneticField/size{n}/distanceVsMagneticWire_{n}.png")
        plt.show()

def main():
    if(len(sys.argv)!=3):
        print("usage python main.py N epsilon")
        exit(0)
    size = int(sys.argv[1])
    epsilon = float(sys.argv[2])

    task = int(input(("Input 0:to generate data\nInput 1: to plot exisiting data\nDefault: 1\nInput: ")))

    #default value of checkerboard is true, so will use checkerboard method
    checkerBoard=False
    if(task==0):
        taskNumber = int(input("Input 0: generate monopole data\n1: generate charged wire data\n2: generate SOR data\nDefault=0\nInput: "))
        method = int(input("Input 0 to use jacobi method\nInput 1 to use Gauss-seidel method\nDefault=1\nInput: "))
        if(method==1):
            useMethod = "gauss"
        else:
            useBoard="jacobi"
        
        
        if(taskNumber==1):
            chargedWire(size,method,epsilon,checkerBoard)
        elif(taskNumber==2):
            sor(size,epsilon,checkerBoard)
        else:
            monopole(size,useMethod,epsilon,checkerBoard)
    
    else:
        plotGraphs=False
        getElectricFit(size,plotGraphs)
        getMagnetic(size,plotGraphs)



t1=time.time()
main()
t2=time.time()
print(f"Time taken ={t2-t1}")



