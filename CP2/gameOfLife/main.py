import numpy as np
from gameOfLife import gameOfLife
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
def animate(size,sweeps,initialisation):
    '''
    Plots the game of life animation for the given input.
    Parameters:
    -----------
    size: Type Integer
          Size of the system
    sweeps: Type Integer
            Number of sweeps for the animation. 
    initialisation: Type string
                    Type of initialisation, if given wrong initialisation, random is used.
    '''

    model = gameOfLife(size,initialisation)

    #calculating the number of live cells 
    model.currentActiveSite=model.activeSites()
    im=plt.imshow(model.lattice,animated=True,vmin=-1,vmax=1,cmap='magma')
    plt.draw()
    plt.pause(0.001)

    
    for i in range(sweeps):
        #updating the model for each sweep
        model.update()
        #updating the number of live cells for each update
        #tracking number of live cells for testing purposes
        newActiveSite = model.activeSites()
        if(newActiveSite==model.currentActiveSite):
            model.activeSitesCounter+=1
        else:
            model.activeSitesCounter=1
            model.currentActiveSite=newActiveSite
        plt.cla()
        im=plt.imshow(model.lattice,animated=True,vmin=-1,vmax=1,cmap='magma')
        plt.draw()
        plt.pause(0.001)
      #  print(f"iteration {i} active state= {model.activeSites()}")
    plt.show()

def generateHistogram(size,sweeps,initialisation):
    '''
    Generates and saves a dat file for the time taken to reach the absorbing state
    for 1000 simulations. Time is measured in number of sweeps here.
    Parameters:
    -----------
    size: Type Integer
          Size of the system
    sweeps: Type Integer
            Number of sweeps for the animation. 
    initialisation: Type string
                    Type of initialisation, if given wrong initialisation, random is used.
    
    '''
    absorbingState = []

    for i in range(100):
        model = gameOfLife(size,initialisation)
        model.currentActiveSite=model.activeSites()
        print(f"Iteration = {i}")
        for j in range(sweeps):
            model.update()
            newActiveSite = model.activeSites()
            if(newActiveSite==model.currentActiveSite):
                model.activeSitesCounter+=1
            else:
                model.activeSitesCounter=1
                model.currentActiveSite= newActiveSite
            if(model.activeSitesCounter==10):
                #break and add to absorbing state
                absorbingState.append(j-9)
                break
   #print(f"Absorbing state = {absorbingState}")
    np.savetxt("data/histogram1000.dat",np.transpose(absorbingState),fmt='%.4f')
      
def generateCom(size,sweeps,initialisation):
    '''
    Calculates tthe centre of mass for the x and y component for upto 300sweeps.
    centre of mass is only calculated for glider that does not crossses the boundary
    of the lattice.
    The centre of mass is saved in a dat file.
    Parameters:
    -----------
    size: Type Integer
          Size of the system
    sweeps: Type Integer
            Number of sweeps for the animation. 
    initialisation: Type string
                    Type of initialisation, if given wrong initialisation, random is used.
    
    '''
    xCom =[]
    yCom=[]
    t=[]
    model = gameOfLife(size,initialisation)
    print(f"Starting to generate centre of mass.")
    for i in range(500):
        model.update()
       # print(len(np.where(model.lattice==1)[0]))
        com = model.centreOfMass()
        #com[0]= x, com[1] = y
        if(com[0]!=-1 and com[1]!=-1):
            #add com
            xCom.append(com[0])
            yCom.append(com[1])
            t.append(i)
        #model.update()

    arrayCombined = np.array((t,xCom,yCom))
   # print(arrayCombined)
    np.savetxt("data/centreOfMass.dat",np.transpose(arrayCombined),fmt='%.4f')
#
def getVelocity(allArray):
    t=allArray[:,0]
    x=allArray[:,1]
    y=allArray[:,2]
    newX=[]
    newY=[]
    newT = []
    for i in range(len(t)):
        if(i>150 and i<250):
            newX.append(x[i])
            newY.append(y[i])
            newT.append(i)
    
    newX=np.asarray(newX)
    newY=np.asarray(newY)
    newT=np.asarray(newT)

    xfit,xin = np.polyfit(newT,newX,1)
    yfit,fin = np.polyfit(newT,newY,1)

    print(f"Xvel = {xfit}\nyVel = {yfit}")

    vel = np.sqrt(xfit**2+yfit**2)
    print(f"Velocity = {vel}")
    return vel


def plotAll():

    #data 1
    data1 = np.loadtxt("data/histogram1000.dat")
    plt.hist(data1,bins=25,density=False,facecolor='none',edgecolor='k')
    plt.title("Histogram for the time taken to reach absorbing state")
    plt.xlabel("Time step till absorbing state")
    plt.ylabel("Frequency")
    plt.savefig("figures/Histogram.png")
    plt.show()

    
    centreOfMass = np.loadtxt("data/centreOfMass.dat")
    t=centreOfMass[:,0]
    xCom=centreOfMass[:,1]
    yCom=centreOfMass[:,2]

    plt.scatter(t[::10],xCom[::10],s=10)
    plt.title("Scatter for xCom")
    plt.ylabel("centre of mass X")
    plt.xlabel("Time (Sweeps)")
    plt.savefig("figures/xCom.png")
    plt.show()

    plt.scatter(t[::10],yCom[::10],s=10)
    plt.title("Scatter for yCom")
    plt.ylabel("centre of mass Y")
    plt.xlabel("Time (Sweeps)")
    plt.savefig("figures/yCom.png")
    plt.show()

    plt.scatter(t[::10],xCom[::10],s=10,label="x")
    plt.scatter(t[::10],yCom[::10],s=10,label="y")
    plt.title("Scatter of both com")
    plt.xlabel("Time (Sweeps)")
    plt.ylabel("Centre of mass")
    plt.legend()
    plt.savefig("figures/bothCOM.png")
    plt.show()


if __name__=="__main__":
    size = 50
    sweeps = 10000
    initialisation="glider"
    t1=time.time()
   # animate(size,sweeps,initialisation)
   # generateHistogram(size,sweeps,initialisation)
   # generateCom(size,sweeps,initialisation)
    centreOfMass = np.loadtxt("data/centreOfMass.dat")
    vel=getVelocity(centreOfMass)
    # t2=time.time()
    # print(f"Time taken for all = {t2-t1}s")

    plotAll()