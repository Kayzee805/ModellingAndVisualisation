from test import gameOfLife
import numpy as np
from gameOfLife import gameOfLife
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time

def animate(size,sweeps,initialisation):
    #for now I only have random initialisation

    #need a switch case here for glider and blinker
    model = gameOfLife(size,initialisation)

    model.currentActiveSite=model.activeSites()
    im=plt.imshow(model.lattice,animated=True,vmin=-1,vmax=1)
    plt.draw()
    plt.pause(0.1)

    
    for i in range(sweeps):
        model.update()
        newActiveSite = model.activeSites()
        if(newActiveSite==model.currentActiveSite):
            model.activeSitesCounter+=1
        else:
            model.activeSitesCounter=1
            model.currentActiveSite=newActiveSite
       # model.centreOfMass()
        plt.cla()
        im=plt.imshow(model.lattice,animated=True,vmin=-1,vmax=1)
        plt.draw()
        plt.pause(0.1)
        print(model.activeSites())
  

    plt.show()

def generateHistogram(size,sweeps,initialisation):

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
    print(f"Absorbing state = {absorbingState}")
    np.savetxt("data/normalHistogram2.dat",np.transpose(absorbingState),fmt='%.4f')
    #plt.hist(absorbingState)
   #plt.show()
   # np.savetxt('data/Normalhistogram3.dat',np.transpose(absorbingState),fmt='%.4f')

def generateCom(size,sweeps,initialisation):
    xCom =[]
    yCom=[]
    t=[]
    model = gameOfLife(size,initialisation)
    for i in range(300):
        if(i%100==0):print(f"Iteration = {i}")
        com = model.centreOfMass()
        #com[0]= x, com[1] = y
        if(com[0]!=-1 and com[1]!=-1):
            #add com
            #i can probs keep it how it is and just slice an array later?
            xCom.append(com[0])
            yCom.append(com[1])
            t.append(i)
        model.update()

    arrayCombined = np.array((t,xCom,yCom))
    np.savetxt("data/centreOfMass2.dat",np.transpose(arrayCombined),fmt='%.4f')
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
    data1 = np.loadtxt("data/Cythonhistogram2.dat")
    plt.hist(data1,bins=25)
    plt.title("Cython histogram 1")
    plt.xlabel("Time step till absorbing state")
    plt.ylabel("Frequency")
    plt.savefig("figures/Histogram2.png")
    plt.show()

    
    centreOfMass = np.loadtxt("data/centreOfMass2.dat")
    t=centreOfMass[:,0]
    xCom=centreOfMass[:,1]
    yCom=centreOfMass[:,2]

    plt.scatter(t,xCom,s=3)
    plt.title("Scatter for xCom")
    plt.ylabel("centre of mass X")
    plt.xlabel("Time (Sweeps)")
    plt.savefig("figures/xCom2.png")
    plt.show()

    plt.scatter(t,yCom,s=3)
    plt.title("Scatter for yCom")
    plt.ylabel("centre of mass Y")
    plt.xlabel("Time (Sweeps)")
    plt.savefig("figures/yCom2.png")
    plt.show()

    plt.scatter(t,xCom,s=3,label="x")
    plt.scatter(t,yCom,s=3,label="y")
    plt.title("Scatter of both com")
    plt.xlabel("Time (Sweeps)")
    plt.ylabel("Centre of mass")
    plt.savefig("figures/bothCOM2.png")
    plt.show()


if __name__=="__main__":
    size = 50
    sweeps = 10000
    initialisation="random"
    t1=time.time()
   # animate(size,sweeps,initialisation)
    generateHistogram(size,sweeps,initialisation)
    #generateCom(size,sweeps,initialisation)

    t2=time.time()
    print(f"Time taken for all = {t2-t1}s")

    #plotAll()