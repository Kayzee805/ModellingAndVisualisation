from test import gameOfLife
import numpy as np
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

    for i in range(200):
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
    plt.hist(absorbingState)
    plt.title("Histogram of time step taken to reach stable state")
    plt.ylabel("Freqeuncy")
    plt.xlabel("Time step")
    plt.savefig("figures/200V2Data.png")
    plt.show()
    np.savetxt('data/Cythonhistogram5.dat',np.transpose(absorbingState),fmt='%.4f')


if __name__=="__main__":
    size = 50
    sweeps = 10000
    initialisation="random"
    t1=time.time()
   # animate(size,sweeps,initialisation)
    generateHistogram(size,sweeps,"random")
    t2=time.time()
    print(f"Time taken for all = {t2-t1}s")

