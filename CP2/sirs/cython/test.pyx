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

class sirs(object):
    def __init__(self,int size,double pS,double pI,double pR):
        '''
        pS=p1. probability of S->I
        pI=p2. probability of I->R
        pR=p3, probability of R->S
        '''
        self.size=size
        self.p1=pS
        self.p2=pI
        self.p3=pR
        self.lattice = np.zeros((size,size))
        self.infected=0
        self.setRandom()
        self.immune=0

    def setRandom(self):
        cdef int counter
        cdef double r
        counter= 0
        for i in range(self.size):
            for j in range(self.size):
                r= random.random()
                if(r<0.33):
                    self.lattice[i,j]=0   #S
                elif(r<0.66):
                    self.lattice[i,j]=-1 #I
                    counter+=1
                else:#
                    self.lattice[i,j]=1  #R

        self.infected=counter

    def setImmune(self,double pImmune):
        cdef double pRest, r
        cdef int counter,immCounter
        pRest = (1-pImmune)/3
        counter=0
        immCounter=0
        for i in range(self.size):
            for j in range(self.size):
                r=random.random()
                if(r<pImmune):
                    #2 is immune
                    self.lattice[i,j]=2
                    immCounter+=1
                 #   print("hello")
                if(r<(pRest+pImmune)):
                    self.lattice[i,j]=-1
                    counter+=1
                if(r<(pRest*2+pImmune)):
                    self.lattice[i,j]=1
                else:
                    self.lattice[i,j]=0
        self.infected=counter
        self.immune=immCounter

    def isImmune(self,double pI):
        cdef int counter
        cdef double r
        counter=0
        for i in range(self.size):
            for j in range(self.size):
                r = random.random()
                if(r<=pI):
                    self.lattice[i,j]=2
                    counter+=1
        self.immune=counter


    def nearestNeighboursInfected(self,int i,int j):
        #Here I check if the nearest neighbour is infected or not
        #top bot left and right
        #returning true means is infected
        #returning false means not infected, so Sus or Rec
        cdef int size,top,bot,left,right
        size = self.size
        top=self.lattice[i,(j-1)%size]
        bot= self.lattice[i,(j+1)%size]
        left = self.lattice[(i-1)%size,j]
        right = self.lattice[(i+1)%size,j]
        if(top==-1 or bot ==-1 or left==-1 or right==-1):
            return True
        return False

    def update(self):
        cdef int iTrial,jTrial,value, counter
        counter=0
        for i in range(self.size):
            for j in range(self.size):
                iTrial = random.randint(0,self.size-1)
                jTrial = random.randint(0,self.size-1)
                #if 0, check for infected states and carry out test
                #if 1, generate random number then check <= p2
                #if 2, generate random number then check <=p3
                r = random.random()
                value = self.lattice[iTrial,jTrial]
                if(value==0):
                    #is sus, so check for nn
                    if(self.nearestNeighboursInfected(iTrial,jTrial)):
                        #nn is infected
                        if(r<=self.p1):
                            self.lattice[iTrial,jTrial]=-1
                            counter+=1
                elif(value==-1):
                    #is infected
                    if(r<=self.p2):
                        self.lattice[iTrial,jTrial]=1
                    else:
                        counter+=1
                elif(value==1):
                    #in recovery
                    if(r<=self.p3):
                        self.lattice[iTrial,jTrial]=0
                else:
                    continue
                    #its immune
        self.infected = counter



    def jacknifeError(self,data):
        N=self.size*self.size
       # resamples = jackknife_resampling(data)
        length=len(data)
        resamples = np.zeros((len(data),len(data)-1))


        for i in range(length):
            resamples[i] = np.delete(data,i)


        x_resamples = np.zeros((len(data)))
        for i in range(len(resamples)):
            x_resamples[i] = np.var(resamples[i])/N
        
        v = np.var(data)/N
        result=0

        for i in range(len(x_resamples)):
            result+=(x_resamples[i]-v)**2


        return np.sqrt(result)
       
def animate(int size,int sweeps,double pS,double pI,double pR):
    '''
    pink = sus
    black= infected
    white = rec
    '''
    #for now I only have random initialisation

    #need a switch case here for glider and blinker
    model = sirs(size,pS,pI,pR)
    print(f"INfected = {model.infected}")
    cMap='magma'
    im=plt.imshow(model.lattice,cmap=cMap,animated=True,vmin=-1,vmax=1)
    plt.draw()
    plt.pause(0.0001)
    uniqueVals = [-1,0,1]
    colors = [im.cmap(im.norm(value)) for value in uniqueVals]
    patches = [ mpatches.Patch(color=colors[i], label="Level {l}".format(l=uniqueVals[i]) ) for i in range(len(uniqueVals)) ]
    plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
    for i in range(sweeps):
        model.update()
        plt.cla()
        im=plt.imshow(model.lattice,cmap=cMap,animated=True,vmin=-1,vmax=1)
        plt.draw()
        plt.pause(0.0001)

    plt.show()



def task3(int size):
    cdef int N,counter, infected
    cdef double variance, averageInfected,p2
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
    print("Starting now")
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
                   # temp = model.meanInfected()
                    infectedRate.append(infected)
                    if(infected==0):
                       # print(f"i={i} j={j}")
                        break

            if(infected==0):
                variance=0
                averageInfected=0
            else:
                infectedRate=np.asarray(infectedRate)
                variance = np.var((infectedRate))/N
                averageInfected = np.mean(infectedRate)/N

            #my indexing is wrong
            allArray[counter]=[p1s[i],p3s[j],averageInfected,variance]
            counter+=1
        print(f"Finished {i} in time {time.time()-start}s")


    print(f"time taken = {time.time()-t1}")
    np.savetxt("data/task3ProcessedData.dat",allArray,fmt='%.7f')
            

def calculateVariance(data):
    meanSq = np.mean(np.square(data))
    squareMean = np.square(np.mean(data))
    print(f"Mean square ={meanSq}  squareMean ={squareMean} ")
    return (np.mean(np.square(data))-(np.square(np.mean(data))))
           
def task4(int size):
    '''
    so I plot x axis = p1s
    y axmis = varance?
    '''
    cdef int infected,i,name
    cdef double t1,p3,p2,variance,vError
    t1=time.time()
    p3=0.5
    p2=0.5
    p1s=np.linspace(0.2,0.5,31)
    size=50
    varianceArray=[]
    varianceError = []
    N=size*size
    print("starting")
    sweeps=10000
    wait =300
    for i in range(len(p1s)):
        model=sirs(50,p1s[i],0.5,0.5)
        infectedRate=np.zeros(sweeps-wait)
        for n in range(sweeps):
            model.update()
            if(n>=wait):
                infected=model.infected
                infectedRate[n-wait]=infected
                if(infected==0):
                # print(f"Breaking at n=d{n}")
                    break
   
            #infectedRate=np.array(infectedRate)
       # autoVar = np.var(infectedRa)
        variance = calculateVariance(infectedRate)/N
        vError = model.jacknifeError(infectedRate)
        print(f"{p1s[i]} variance={variance}  error = {vError} mean = {np.mean(infectedRate)}")
        varianceArray.append(variance)
        varianceError.append(vError)

    combined = np.array((p1s,varianceArray,varianceError))
    np.savetxt("data/task4ProcessedData.dat",np.transpose(combined),fmt='%.6f')

    plt.scatter(p1s,varianceArray,s=5)
    plt.plot(p1s,varianceArray)
    plt.errorbar(p1s,varianceArray,yerr=varianceError)
    plt.show()
    t2=time.time()
    print(f"Time taken = {t2-t1}s")

def task5(int size,int sweeps, name):
    '''
    1000-10,000
    y-axis is avg infected
    x axis is avg immune
    '''
    cdef double p1,p2,p3,t1
    cdef int precision,N,infected,s,i,n
    p1=p2=p3=0.5
    pImmune = np.linspace(0,1,101) 
    precision = len(pImmune)
    N=size*size
    times = 1
    infectionArray = np.zeros((precision,times))
    print(f"Starting {name}")
    for s in range(times):
        t1=time.time()

        for i in range(len(pImmune)):
            model=sirs(size,p1,p2,p3)
            model.isImmune(pImmune[i])
            infectedRate=np.zeros(sweeps-100)
            t3=time.time()
            for n in range(sweeps): 
                model.update()

                if(n>=100):
                    infected=model.infected
                    infectedRate[n-100]=infected
                    if(infected==0):
                        break
            infectionArray[i,s]=np.mean(infectedRate)/N 
            print(f"{name} at i {i}Time taken {s} == {time.time()-t3}s")
        print(f"{name}Time taken {s} == {time.time()-t1}s")
    
    return infectionArray

def task7(int size,int sweeps, name):
    '''
    1000-10,000
    y-axis is avg infected
    x axis is avg immune
    '''
    cdef double p1,p2,p3,t1
    cdef int precision,N,infected,s,i,n
    p1=0.8
    p2=0.1
    p3=0.02
    pImmune = np.linspace(0,1,101) 
    precision = len(pImmune)
    N=size*size
    times = 1
    infectionArray = np.zeros((precision,times))
    print(f"Starting {name}")
    for s in range(times):
        t1=time.time()

        for i in range(len(pImmune)):
            model=sirs(size,p1,p2,p3)
            model.isImmune(pImmune[i])
            infectedRate=np.zeros(sweeps-100)
            t3=time.time()
            for n in range(sweeps): 
                model.update()

                if(n>=100):
                    infected=model.infected
                    infectedRate[n-100]=infected
                    if(infected==0):
                        break
            infectionArray[i,s]=np.mean(infectedRate)/N 
            print(f"{name} at i {i}Time taken {s} == {time.time()-t3}s")
        print(f"{name}Time taken {s} == {time.time()-t1}s")
    
    return infectionArray
    # finalArray = []
    # errors = []
    # for i in range(precision):
    #     finalArray.append(np.mean(infectionArray[i]))
    #     errors.append(sem(infectionArray[i]))
    
    # combined = np.array((pImmune,finalArray,errors))
    # np.savetxt('data/'+name+'.dat',np.transpose(combined),fmt='%.6f')

    # plt.scatter(pImmune,finalArray,s=5)
    # plt.plot(pImmune,finalArray)
    # plt.xlabel("Immune")
    # plt.ylabel("Infection")
    # plt.title("immune vs infection")
    # plt.savefig("figures/immune.png")
    # plt.show()

    
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

def main(int size,int sweeps,double pS,double pI,double pR,runAnim=False,genData=False,task5DO=False):

    if(runAnim):
        animate(size,sweeps,pS,pI,pR)
    if(genData):
        task3(size)
    if(task5DO):
        task5(size,sweeps,"hello")
    pltContour()