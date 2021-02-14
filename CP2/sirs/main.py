from sirs import sirs
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import colors
import matplotlib.patches as mpatches
from scipy.stats import sem
import time
import random
import sys
matplotlib.use('TKAgg')

def animate(size,sweeps,pS,pI,pR):
    '''
    pink = sus
    black= infected
    white = rec
    '''
    #for now I only have random initialisation

    #need a switch case here for glider and blinker
    model = sirs(size,pS,pI,pR)
    cMap='magma'
    im=plt.imshow(model.lattice,cmap=cMap,animated=True,vmin=-1,vmax=1)
    plt.draw()
    plt.pause(0.0001)
    for i in range(sweeps):
        model.update()
        plt.cla()
        im=plt.imshow(model.lattice,cmap=cMap,animated=True,vmin=-1,vmax=1)
        plt.draw()
        plt.pause(0.0001)
        if(i%1000==0):
            print(i)
        if(model.infected==0):
            print(f"Infected = 0 at {i}")

    plt.show()



def task3(size):
    p2=0.5
    #basically need an array of p1 and p3. np.linspace 
    #then for each combination? find the number of infected after 1000 sweeps?
    #use 100 sweeps as equilibration time, so only take data after 100 sweeps
    #measure every sweep#


    #if Infection ==0. stop
    '''
    p1. p3. averageInfected averageofSq variance
    '''
    p1s= np.linspace(0,1,21)
    p3s= np.linspace(0,1,21)
    #maybe I only save average infected and varaince
    allArray = np.zeros((21*21,4))
    N=size*size
    t1=time.time()
    counter =0
    for i in range(10,21):
        #print(i)
        start=time.time()

        for j in range(5,21):

            model = sirs(size,p1s[i],p2,p3s[j])
            model.setImmune(0.1)
            infectedRate=[]
            for n in range(1000):
                model.update()
                if(n>100):
                    infected = model.infected
                   # temp = model.meanInfected()
                    infectedRate.append(infected)
                    if(infected==0):
                       # print(f"i={i} j={j}")
                        break
            if(infected==0):
                averageInfected=0
                variance=0
            else:
                infectedRate=np.asarray(infectedRate)
                variance = np.var((infectedRate))/N
                averageInfected = np.mean(infectedRate)/N
            #my indexing is wrong
            allArray[counter]=[p1s[i],p3s[j],averageInfected,variance]
            counter+=1
        print(f"Finished {i} in time {time.time()-start}s")


    print(f"time taken = {time.time()-t1}")
    np.savetxt("data/task32.dat",allArray,fmt='%.7f')
            
def task4(size):
    '''
    so I plot x axis = p1s
    y axmis = varance?
    '''
    t1=time.time()
    p3=0.5
    p2=0.5
    p1s=np.linspace(0.2,0.5,31)
  
    varianceArray=[]
    varianceError = []
    N=size*size
    print("starting")
    print(p1s)
    for i in range(len(p1s)):
        model=sirs(size,p1s[i],p2,p3)
        infectedRate=[]
        for n in range(10000):
            model.update()
            if(n>=100):
                infected=model.infected
                infectedRate.append(infected)
                if(infected==0):
                    print(f"Breaking at n={n}")
                    break
        
        if(infected==0):
            variance=0
            vError=0
        else:
            #infectedRate=np.array(infectedRate)
            variance = np.var(infectedRate)/N
            vError = model.jacknife(infectedRate)
        print(f"{p1s[i]} variance={variance} error = {vError}")

        varianceArray.append(variance)
        varianceError.append(vError)
    plt.scatter(p1s,varianceArray,s=5)
    plt.plot(p1s,varianceArray)
    plt.errorbar(p1s,varianceArray,yerr=varianceError)
    plt.show()
    combined = np.array((p1s,varianceArray,varianceError))
    np.savetxt("data/task4Variance.dat",np.transpose(combined),fmt='%.6f')
    t2=time.time()
    print(f"Time taken = {t2-t1}s")


def task5(size,sweeps):
    '''
    1000-10,000
    y-axis is avg infected
    x axis is avg immune
    '''
    p1=p2=p3=0.5
    pImmune = np.linspace(0,1,3) 
    precision = len(pImmune)
    N=size*size
    infectionArray = np.zeros((precision,5))
    
    for s in range(5):
        print(s)
        t1=time.time()

        for i in range(len(pImmune)):
            model=sirs(size,p1,p2,p3)
            model.isImmune(pImmune[i])
            infectedRate=[]

            for n in range(sweeps):
                model.update()

                if(n>=100):
                    infected=model.infected
                    infectedRate.append(infected)
                    if(infected==0):
                        break
            infectionArray[i,s]=np.mean(infectedRate)/N 
        print(f"Time taken {s} == {time.time()-t1}s")
    finalArray = []
    errors = []
    for i in range(precision):
        finalArray.append(np.mean(infectionArray[i]))
        errors.append(sem(infectionArray[i]))
    
    combined = np.array((pImmune,finalArray,errors))
    np.savetxt('data/TestTASK5.dat',np.transpose(combined),fmt='%.6f')

    plt.scatter(pImmune,finalArray,s=5)
    plt.plot(pImmune,finalArray)
    plt.xlabel("Immune")
    plt.ylabel("Infection")
    plt.title("immune vs infection")
    plt.savefig("figures/immune.png")
    plt.show()

    
def pltContour():
    array = np.loadtxt("data/task32.dat")
    p1s=array[:,0]
    p3s=array[:,1]
    infected = array[:,2]
    variance = array[:,3]

    print(f"Lengh.  p1s={p1s.size}  p2s={p3s.size} infected={len(infected)}")
    p1s = p1s.reshape((21,21))
    p3s=p3s.reshape((21,21))
    infected = infected.reshape((21,21))
    variance = variance.reshape((21,21))

    plt.figure()
  #  CS=plt.contour(p1s,p3s,infected)
    CS=plt.contourf(p1s,p3s,infected,cmap='magma')
    #plt.clabel(CS,fontsize=8,colors='k')
    cbar=plt.colorbar(CS)
    plt.xticks(np.linspace(0,1,11))
    plt.yticks(np.linspace(0,1,11))

    plt.xlabel("P1")
    plt.ylabel("P3")
    plt.title("Average infected of p1-p3 plane")
    plt.savefig("figures/infectedContour.png")
    plt.show()

    #plt.cla()

    plt.figure()
  #  CS=plt.contour(p1s,p3s,infected)
    CS=plt.contourf(p1s,p3s,variance,cmap='magma')
    #plt.clabel(CS,fontsize=8,colors='k')
    cbar=plt.colorbar(CS)
    plt.xticks(np.linspace(0,1,11))
    plt.yticks(np.linspace(0,1,11))

    plt.xlabel("P1")
    plt.ylabel("P3")
    plt.title("Scaled variance of p1-p3 plane")
    plt.savefig("figures/varianceContour.png")
    plt.show()


    array2 = np.loadtxt("data/task4Variance.dat")
    xAxis = array2[:,0]
    yAxis=array2[:,1]
    yError=array2[:,2]
    plt.scatter(xAxis,yAxis,s=5,color='k')
    plt.plot(xAxis,yAxis)
    plt.errorbar(xAxis,yAxis,yerr=yError,ecolor='k')

    plt.xlabel("P1")
    plt.ylabel("Variance")
    plt.title("Variance vs P1 with p2=p3=0.5")
    plt.savefig("figures/variance_p1.png")
    plt.show()
if __name__=="__main__":
    if(len(sys.argv)!=5):
        print("usage python main.py N p1 p2 p3")
        exit(0)
    size = int(sys.argv[1])
    pS = float(sys.argv[2])
    pI = float(sys.argv[3])
    pR = float(sys.argv[4])
    sweeps = 10000

    #animate(size,sweeps,pS,pI,pR)
    #task3(size)
    #ask4(size)
    task5(size,sweeps)
  #  pltContour()
#