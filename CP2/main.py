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
    #plt.hist(absorbingState)
   #plt.show()
    np.savetxt('data/Normalhistogram3.dat',np.transpose(absorbingState),fmt='%.4f')



def plotAll():

    #data 1
    data1 = np.loadtxt("data/Cythonhistogram.dat")
    plt.hist(data1,bins=25)
    plt.title("Cython histogram 1")
    plt.ylabel("Time step till absorbing state")
    plt.xlabel("Frequency")
    plt.savefig("figures/CythonHistogram.png")
    plt.show()
    #data 2
    data2 = np.loadtxt("data/Cythonhistogram2.dat")
    plt.hist(data2,bins=25)
    plt.title("Cython histogram 2")
    plt.ylabel("Time step till absorbing state")
    plt.xlabel("Frequency")
    plt.savefig("figures/CythonHistogram2.png")
    plt.show()
    #data 3
    data3 = np.loadtxt("data/Cythonhistogram3.dat")
    plt.hist(data3,bins=25)
    plt.title("Cython histogram 3")
    plt.ylabel("Time step till absorbing state")
    plt.xlabel("Frequency")
    plt.savefig("figures/CythonHistogram3.png")
    plt.show()
    #data 4
    data4 = np.loadtxt("data/Normalhistogram1.dat")
    plt.hist(data4,bins=25)
    plt.title("Normal histogram 1")
    plt.ylabel("Time step till absorbing state")
    plt.xlabel("Frequency")
    plt.savefig("figures/normalHistogram1.png")
    plt.show()

    #data 5
    data5 = np.loadtxt("data/Normalhistogram2.dat")
    plt.hist(data5,bins=25)
    plt.title("Normal histogram 2")
    plt.ylabel("Time step till absorbing state")
    plt.xlabel("Frequency")
    plt.savefig("figures/normalHistogram2.png")
    plt.show()

    #data 4
    data6 = np.loadtxt("data/Normalhistogram3.dat")
    plt.hist(data6,bins=25)
    plt.title("Normal histogram 3")
    plt.ylabel("Time step till absorbing state")
    plt.xlabel("Frequency")
    plt.savefig("figures/normalHistogram3.png")
    plt.show()
        #data 4
    data7 = np.loadtxt("data/cythonhistogram4.dat")
    plt.hist(data7,bins=40)
    plt.title("200 dat ")
    plt.ylabel("Time step till absorbing state")
    plt.xlabel("Frequency")
    plt.savefig("figures/200data.png")
    plt.show()


if __name__=="__main__":
#     size = 50
#     sweeps = 10000
#     initialisation="random"
#     t1=time.time()
#    # animate(size,sweeps,initialisation)
#     generateHistogram(size,sweeps,"random")
#     t2=time.time()
#     print(f"Time taken for all = {t2-t1}s")
    data7 = np.loadtxt("data/cythonhistogram5.dat")
    plt.hist(data7,bins=40)
    plt.title("200 dat ")
    plt.ylabel("Time step till absorbing state")
    plt.xlabel("Frequency")
    #plt.savefig("figures/200data.png")
    plt.show()

    #plotAll()