import numpy as np
import matplotlib.pyplot as plt 
import random as r
import matplotlib.animation as animation
from tqdm import tqdm
from scipy.stats import sem
import math

'''
cell either infected or healthy
if random cell is healthy, nothing happens
if random cell is infected
    It has a probability 1-p to become healthy
    and probability p to infect a random nn
        If chosen neighbour is infected already, nothing happens
    
'''



#1 =healthy
#0=unhealthy
class Model(object):

    def __init__(self,n,p,initialProbability=0.5,taskf=False):
        self.n=n
        self.p = p
        self.initialProbability =initialProbability
        self.lattice = np.zeros((self.n,self.n))
        self.infectedCounter=0
        self.initialiseLattice()

        if(taskf):
            self.lattice=np.ones((self.n,self.n))
            self.singleActive()



    def initialiseLattice(self):
        counter=0
        for i in range(self.n):
            for j in range(self.n):
                randNumber = r.random()
                if (randNumber>self.initialProbability):
                    self.lattice[i,j]=1
                    counter+=1
        
        self.infectedCounter= (self.n*self.n)-counter

    def singleActive(self):
        i = r.randint(0,self.n-1)
        j = r.randint(0,self.n-1)
        self.lattice[i,j]=0
        self.infectedCounter=1

    def updateRandomNeighbour(self,i,j):
        #basically have some system to choose one neighbours out of the 4 neighbours
        #1,2,3,4 = left,right,up,down
        neighbour = r.randint(1,4)

        if(neighbour==1):
            #check for left
            if(self.lattice[i,(j-1)%self.n]==1):
                self.lattice[i,(j-1)%self.n]=0
                self.infectedCounter+=1
        elif(neighbour==2):
            #check for right
            if(self.lattice[i,(j+1)%self.n]==1):
                self.lattice[i,(j+1)%self.n]=0
                self.infectedCounter+=1

        elif(neighbour==3):
            #check for up
            if(self.lattice[(i-1)%self.n,j]==1):
                self.lattice[(i-1)%self.n,j]=0
                self.infectedCounter+=1

        elif(neighbour==4):
            #check for down
            if(self.lattice[(i+1)%self.n,j]==1):
                self.lattice[(i+1)%self.n,j]=0
                self.infectedCounter+=1

        else:
            print(f"It should not be here. Neighbour {neighbour} does not exist")
            exit(0)

    def update(self):
       # print("starting update")
        for i in range(self.n):
            for j in range(self.n):
                randomI = r.randint(0,self.n-1)
                randomJ = r.randint(0,self.n-1)
                if(self.lattice[randomI,randomJ]==1):
                    #is healthy
                    continue
                else:
                    if(r.random()>self.p):
                        #so higher than 1-p, so healthy
                        self.lattice[randomI,randomJ]=1 
                        self.infectedCounter-=1

                    else:
                        self.updateRandomNeighbour(randomI,randomJ)


def animate(n,p):
    model = Model(n,p)
    sweeps = 10000
    colour= 'magma'
    plt.figure()
    im = plt.imshow(model.lattice,animated=True,cmap=colour,vmin=-1,vmax=1)
    plt.colorbar(im)
    plt.draw()
    plt.pause(0.001)
    for i in range(sweeps):
        model.update()
        if i%10==0:
            plt.clf()
            im=plt.imshow(model.lattice,animated=True,cmap=colour,vmin=-1,vmax=1)
            plt.draw()
            plt.pause(0.0001)
    
    plt.show()


def taskB(n):
    p1 = 0.6
    p2=0.7

    model1 = Model(n,p1)
    model2 = Model(n,p2)

    p1Infected = []
    p2Infected=[]
    sweeps = 2000
    steps = np.linspace(0,sweeps,sweeps+1)

    N=n*n
    for i in tqdm(steps):
        model1.update()
        model2.update()
        
        p1Infected.append(model1.infectedCounter)
        p2Infected.append(model2.infectedCounter)


    plt.plot(steps,p1Infected,label="p=0.6",color='r')
    plt.plot(steps,p2Infected,label="p=0.7",color='b')
    plt.xlabel("Time(sweeps)")
    plt.ylabel("Number of infected sites")
    plt.legend()
    plt.savefig(f"figures/taskB.png")
    plt.show()   


    print(len(steps),len(p1Infected),len(p2Infected))
    np.savetxt("data/taskB.dat",np.transpose(np.array((steps,p1Infected,p2Infected))),fmt='%.0f') 


def taskC(n):
    allP = np.arange(0.55,0.7,0.005)
    N = n*n
    sweeps = 1000

    averageInfected=[]
    for i in range(len(allP)):
        model = Model(n,allP[i])
        for j in range(sweeps):
            model.update()

            if(model.infectedCounter==0):
                print(allP[i])
                break 
        average = model.infectedCounter/N
        averageInfected.append(average)
        print(f"Finished p={allP[i]} with average = {average}")


    plt.scatter(allP,averageInfected,s=10,marker='+',color='r')
    plt.plot(allP,averageInfected,color='b')
    plt.xlabel("p values")
    plt.ylabel("Average Infected")
    plt.title("Plot of average infected for p=0.55->0.7")
    plt.savefig("figures/taskC.png")
    plt.show()
    np.savetxt("data/taskC.dat",np.transpose(np.array((allP,averageInfected))),fmt='%.4f')



def taskD(n):
    allP = np.arange(0.55,0.7,0.005)
    N = n*n
    sweeps = 800

    allArray = np.zeros((len(allP),5))

    for x in range(5):
        for i in range(len(allP)):
            model = Model(n,allP[i])

            total=[]

            for j in range(sweeps):
                model.update()
                total.append(model.infectedCounter)
                if(model.infectedCounter==0):
                    break
            if(model.infectedCounter==0):
                variance =0
            else:
                variance = np.var(total)/N
            
            allArray[i,x]=variance
            print(x,allP[i])

    averageVar = []
    varError = []

    for i in range(len(allP)):
        averageVar.append(np.average(allArray[i]))
        varError.append(sem(allArray[i]))
    

    plt.scatter(allP,averageVar,s=5,color='k')
    plt.plot(allP,averageVar,color='b')
    plt.errorbar(allP,averageVar,yerr=varError,ecolor='r')
    plt.xlabel("Probability")
    plt.ylabel("Variance")
    plt.title("Plot of varaince over p=0.55->0.7")
    plt.savefig(f"figures/taskD.png")
    plt.show()

    np.savetxt("data/taskD.dat",np.transpose(np.array((allP,averageVar,varError))),fmt='%.5f')


def taskF(n,p):
    '''
    Plot survival probability as a function of time?

    x axis = time
    y axis = average survivial probability over many runs
    '''
    sweeps = 300
    attempts = 100
    steps =np.linspace(1,sweeps,sweeps)
    print(steps)
    allArray = np.zeros((sweeps,attempts))

    for i in range(attempts):
        model = Model(n,p,taskf=True)
        print(i)
        for j in range(sweeps):
            allArray[j,i]=model.infectedCounter
            model.update()

    
    allArray = allArray/(n*n)
    averages = []
    for i in range(len(allArray)):
        averages.append(np.mean(allArray[i]))

    plt.plot(steps,averages)
    plt.xlabel("time(sweeps)")
    plt.ylabel("Survival probability")
    plt.title(f"Survival rate for p = {p}")
    plt.savefig(f"figures/taskf_{p}.png")
    plt.show()
    np.savetxt(f"data/taskf_{p}.dat",np.transpose(np.array((steps,averages))),fmt='%.7f')

    logSteps = np.log10(steps)
    logAverages = np.log10(averages)

    mask = np.isinf(logAverages)
    logAverages[mask]=0
    plt.plot(logSteps,logAverages)
    plt.xlabel("log_10(time(Sweeps))")
    plt.ylabel("log_10(Survival probability)")
    plt.title(f"Log log plot of Surivavl rate for p={p}")
    plt.savefig(f"figures/taskf_loglog_{p}.png")
    plt.show()
    np.savetxt(f"data/taskf_loglog_{p}.dat",np.transpose(np.array((logSteps,logAverages))),fmt='%.5f')


    plt.plot(steps,averages)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("time(sweeps)")
    plt.ylabel("Survival probability")
    plt.title(f"Survival rate for p = {p}")
    plt.savefig(f"figures/test{p}.png")
    plt.show()
if __name__=="__main__":
    n=50
    p=0.6
    #animate(n,p)
    #taskB(n)
   # taskC(n)
  #  taskD(n)
    taskF(n,0.6)


