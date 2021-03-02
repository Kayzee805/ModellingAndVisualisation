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
import seaborn as sns
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
    #basically need an array of p1 and p3. np.linspace 
    #then for each combination? find the number of infected after 1000 sweeps?
    #use 100 sweeps as equilibration time, so only take data after 100 sweeps
    #measure every sweep#


    #if Infection ==0. stop
    '''
    p1. p3. averageInfected averageofSq variance
    '''
    p2=0.5
    p1s= np.linspace(0,1,21)
    p3s= np.linspace(0,1,21)
    #maybe I only save average infected and varaince
    allArray = np.zeros((21*21,4))
    N=size*size
    t1=time.time()
    counter =0
    for i in range(21):
        #print(i)
        start=time.time()

        for j in range(21):

            model = sirs(size,p1s[i],p2,p3s[j])
            infectedRate=[]
            for n in range(1000):
                model.update()
                if(n>=100):
                    infected = model.infected
                    infectedRate.append(infected)
                    if(infected==0):
                        break
            if(infected==0):
                averageInfected=0
                variance=0
            else:
                infectedRate=np.asarray(infectedRate)
                variance = np.var((infectedRate))/N
                averageInfected = np.mean(infectedRate)/N
            allArray[counter]=[p1s[i],p3s[j],averageInfected,variance]
            counter+=1
        print(f"Finished {i} in time {time.time()-start}s")


    print(f"time taken = {time.time()-t1}")
    #np.savetxt("data/Task3ProcessedData.dat",allArray,fmt='%.7f')
    np.savetxt("data/Task3ProcessedData100x100.dat",allArray,fmt='%.7f')

    array = np.loadtxt("data/Task3ProcessedData100x100.dat")
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

    CS=plt.contourf(p1s,p3s,variance,cmap='magma')
    #plt.clabel(CS,fontsize=8,colors='k')
    cbar=plt.colorbar(CS)
    plt.xticks(np.linspace(0,1,11))
    plt.yticks(np.linspace(0,1,11))

    plt.xlabel("P1")
    plt.ylabel("P3")
    plt.title("Scaled variance of p1-p3 plane")
    plt.savefig("figures/Task3_ScaledVariance100x100.png")
    plt.show()


def task4(size,sweeps):
    '''
    so I plot x axis = p1s
    y axmis = varance?
    '''
    t1=time.time()
    p3=0.5
    p2=0.5
    p1s=np.linspace(0.2,0.5,31)
    N=size*size
    allArray = np.zeros((len(p1s),1))
    
    for s in range(5):
        print(f"Starting {s}")
        for i in range(len(p1s)):
            start = time.time()

            model = sirs(size,p1s[i],p2,p3)
            infectedRate=[]
            for n in range(sweeps):
                model.update()
                if(n>=100):
                    infected=model.infected
                    infectedRate.append(infected)
                    if(infected==0):
                        break
            if(infected==0):
                variance=0
            else:
                variance=np.var(infectedRate)/N
            
            allArray[i,s]=variance
            print(f"Time for {s} {i} at time = {time.time()-start}")
    
    np.savetxt(f"data/Task4_RawDataCombined.dat",np.transpose(allArray),fmt='%.6f')
    finalArray =[]
    errors = []
    for i in range(len(p1s)):
       finalArray.append(np.mean(allArray[i]))
       errors.append(sem(allArray[i]))
    combined = np.array((p1s,finalArray,errors))
    np.savetxt("data/Task4_ProcessedData.dat",np.transpose(combined),fmt='%.6f')
    print("DONEEE")   


def task5(size,sweeps):
    '''
    1000-10,000
    y-axis is avg infected
    x axis is avg immune
    '''
    p1=p2=p3=0.5
    pImmune = np.linspace(0,1,101) 
    precision = len(pImmune)
    N=size*size
    times=5
    infectionArray = np.zeros((precision,times))
    for s in range(times):
        print(s)
        t1=time.time()

        for i in range(len(pImmune)):
            print(pImmune[i])
            model=sirs(size,p1,p2,p3,isImmune=True,immuneProbability=pImmune[i])
            infectedRate=[]
            for n in range(sweeps):
                model.update()
                if(n>=100):
                    infected=model.infected
                    infectedRate.append(infected)
                    if(infected==0):
                        break
            if(infected==0):
                infectionArray[i,s]=0
            else:
                infectionArray[i,s]=np.mean(infectedRate)/N 
        print(f"Time taken {s} == {time.time()-t1}s")
    finalArray = []
    errors = []
    for i in range(precision):
        finalArray.append(np.mean(infectionArray[i]))
        errors.append(sem(infectionArray[i]))
    
    combined = np.array((pImmune,finalArray,errors))
    np.savetxt('data/Task5_RawDataCombined.dat',(infectionArray),fmt='%.6f')
    np.savetxt('data/Task5_ProcessedData.dat',np.transpose(combined),fmt='%.6f')
def task7(size,sweeps):
    '''
    1000-10,000
    y-axis is avg infected
    x axis is avg immune
    '''
    p1=0.8
    p2=0.1
    p3=.02
    pImmune = np.linspace(0,1,101) 
    precision = len(pImmune)
    N=size*size
    times=5
    infectionArray = np.zeros((precision,times))
    for s in range(times):
        print(s)
        t1=time.time()
        for i in range(len(pImmune)):
            model=sirs(size,p1,p2,p3,isImmune=True,immuneProbability=pImmune[i])
            infectedRate=[]
            for n in range(sweeps):
                model.update()

                if(n>=100):
                    infected=model.infected
                    infectedRate.append(infected)
                    if(infected==0):
                        break
            if(infected==0):
                infectionArray[i,s] =0
            else:
                infectionArray[i,s]=np.mean(infectedRate)/N 
        print(f"Time taken {s} == {time.time()-t1}s")

    finalArray = []
    errors = []
    for i in range(precision):
        finalArray.append(np.mean(infectionArray[i]))
        errors.append(sem(infectionArray[i]))
    
    combined = np.array((pImmune,finalArray,errors))
    np.savetxt('data/Task7_RawDataCombined.dat',(infectionArray),fmt='%.6ft')
    np.savetxt('data/Task7_ProcessedData.dat',np.transpose(combined),fmt='%.6f')

    # plt.scatter(pImmune,finalArray,s=5)
    # plt.plot(pImmune,finalArray)
    # plt.xlabel("Immune")
    # plt.ylabel("Infection")
    # plt.title("immune vs infection")
    # plt.savefig("figures/immune.png")
    # plt.show()


def plotAll():
    array = np.loadtxt("data/task3ProcessedData.dat")
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
    plt.savefig("figures/Task3_AverageInfected.png")
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
    plt.savefig("figures/Task3_ScaledVariance.png")
    plt.show()

    variance = np.asarray(variance)
    #print(f"Max variance = {np.amax(variance)}")
    xTicks = np.round(np.linspace(0,1,11),2)
    plt.imshow(np.transpose(variance),cmap='magma',extent=[0,1,1,0])
    plt.xlabel("P1")
    plt.ylabel("P3")
    plt.xticks(xTicks)
    plt.yticks(xTicks)
    plt.title("Scaled variance of p1-p3 plane")
    plt.colorbar(label="Scaled variance")
    plt.gca().invert_yaxis()
    plt.savefig("figures/TASK3_ScaledVariance_imshow.png")
    plt.show()
     
    xTicks = np.round(np.linspace(0,1,11),2)
    plt.imshow(np.transpose(infected),cmap='magma',extent=[0,1,1,0])
    plt.xlabel("P1")
    plt.ylabel("P3")
    plt.xticks(xTicks)
    plt.yticks(xTicks)
    plt.title("Average infected of p1-p3 plane")
    plt.colorbar(label="Ratio of average infected")
    plt.gca().invert_yaxis()

    plt.savefig("figures/TASK3_AverageInfected_imshow.png")
    plt.show()




    array2 = np.loadtxt("data/Task4_ProcessedData.dat")
    xAxis = array2[:,0]
    yAxis=array2[:,1]
    yError=array2[:,2]
    plt.scatter(xAxis,yAxis,s=5,color='r')
    plt.plot(xAxis,yAxis,color='b')
    plt.errorbar(xAxis,yAxis,yerr=yError,ecolor='k')

    plt.xlabel("P1")
    plt.ylabel("Variance")
    plt.title("Variance vs P1 with p2=p3=0.5")
    plt.savefig("figures/Task4_ScaledVariance_cutP3.png")
    plt.show()


    #TASK 5
    task5Array = np.loadtxt("data/Task5_ProcessedData.dat")
    immuneProbability = task5Array[:,0]
    task5Variance = task5Array[:,1]
    task5Error = task5Array[:,2]

    plt.scatter(immuneProbability,task5Variance,s=10,color='k')
    plt.plot(immuneProbability,task5Variance,color='b')
    plt.errorbar(immuneProbability,task5Variance,yerr=task5Error,ecolor='r')
    plt.xlabel("Immune probability")
    plt.ylabel("Average infected")
    plt.title("Immune vs Infected for 50x50 p1=p2=p3=0.5")
    plt.savefig("figures/Task5_50x50.png")
    plt.show()

        #TASK 7
    task7Array = np.loadtxt("data/Task5Part2_ProcessedData.dat")
    immuneProbability7 = task7Array[:,0]
    task7Variance = task7Array[:,1]
    task7Error = task7Array[:,2]

    plt.scatter(immuneProbability7,task7Variance,s=10,color='k')
    plt.plot(immuneProbability7,task7Variance,color='b')
    plt.errorbar(immuneProbability7,task7Variance,yerr=task7Error,ecolor='r')
    plt.xlabel("Immune probability")
    plt.ylabel("Average infected")
    plt.title("Immune vs Infected for 100x100 p1=0.8,p2=0.1,p3=0.02")
    plt.savefig("figures/Task5_100x100.png")
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

   # animate(size,sweeps,pS,pI,pR)
    task3(size)
   # task4(size,sweeps)
   # task5(size,sweeps)
    #plotAll()
  ##  pltContour()
#