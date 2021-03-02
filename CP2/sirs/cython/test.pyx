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
    def __init__(self,int size,double pS,double pI,double pR, isImmune = False,double immuneProbability=0):
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
        self.immune=0
        if(isImmune):
            self.setImmune(immuneProbability)
        else:
            self.setRandom()

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

    def setImmune(self,double immuneProbability):
        cdef double pRest, r
        cdef int infectedCounter,immuneCounter
        pRest = (1-immuneProbability)/3
        infectedCounter=0
        immuneCounter=0
        for i in range(self.size):
            for j in range(self.size):
                r=random.random()
                if(r<immuneProbability):
                    #2 is immune
                    self.lattice[i,j]=2
                    immuneCounter+=1
                 #   print("hello")
                elif(r<(pRest+immuneProbability)):
                    self.lattice[i,j]=-1
                    infectedCounter+=1
                elif(r<(pRest*2+immuneProbability)):
                    self.lattice[i,j]=1
                else:
                    self.lattice[i,j]=0
        self.infected=infectedCounter
        self.immune=immuneCounter


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
        cdef int iTrial,jTrial,value, infectedCounter
        infectedCounter=self.infected
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
                            infectedCounter+=1
                elif(value==-1):
                    #is infected
                    if(r<=self.p2):
                        self.lattice[iTrial,jTrial]=1
                        infectedCounter-=1
                elif(value==1):
                    #in recovery
                    if(r<=self.p3):
                        self.lattice[iTrial,jTrial]=0
                else:
                    continue
                    #its immune
        self.infected = infectedCounter



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



def task3(int size):
    cdef int N,counter, infected
    cdef double variance, averageInfected,p2

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
        print(i)
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
    np.savetxt("data/Task3ProcessedData.dat",allArray,fmt='%.7f')
            

def calculateVariance(data):
    meanSq = np.mean(np.square(data))
    squareMean = np.square(np.mean(data))
    print(f"Mean square ={meanSq}  squareMean ={squareMean} ")
    return (np.mean(np.square(data))-(np.square(np.mean(data))))
           
def task4(int size,int sweeps):
    '''
    so I plot x axis = p1s
    y axmis = varance?
    '''
    fileName = input("Enter teh name of the file")
    print(fileName)
    cdef int infected,i,name,N
    cdef double t1,p3,p2,variance,vError,start
    t1=time.time()
    p3=0.5
    p2=0.5
    p1s=np.linspace(0.2,0.5,31)
    N=size*size
    allArray = np.zeros((len(p1s),1))
    
    for s in range(1):
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
    
    np.savetxt(f"data/Task4_RawData{fileName}.dat",np.transpose(allArray),fmt='%.6f')
    #finalArray =[]
    #errors = []
    #for i in range(len(p1s)):
      #  finalArray.append(np.mean(allArray[i]))
     #   errors.append(sem(allArray[i]))
    #combined = np.array((p1s,finalArray,errors))
    #np.savetxt("data/Task4_ProcessedData.dat",np.transpose(combined),fmt='%.6f')
    print("DONEEE")    

def task5(int size,int sweeps):
    '''
    1000-10,000
    y-axis is avg infected
    x axis is avg immune
    '''
    cdef double p1,p2,p3,t1
    cdef int precision,N,infected,s,i,n,times
    times=1
    name = input("Enter file name")
    p1=p2=p3=0.5
    pImmune = np.linspace(0,1,101) 
    precision = len(pImmune)
    N=size*size
    infectionArray = np.zeros((precision,times))
    for s in range(times):
        for i in range(len(pImmune)):
            print(pImmune[i])
            t1=time.time()

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
            print(f"Time taken {i} == {time.time()-t1}s")
    
    np.savetxt(f"data/task5/Task5_RawData{name}.dat",infectionArray,fmt='%.6f')
    #finalArray = []
  #  errors = []
   # for i in range(precision):
       # finalArray.append(np.mean(infectionArray[i]))
       # errors.append(sem(infectionArray[i]))
    
    #combined = np.array((pImmune,finalArray,errors))

def task7(int size,int sweeps):
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
    size=100
    times=1
    name = input("Enter file name")
    pImmune = np.linspace(0,1,101) 
    precision = len(pImmune)
    N=size*size
    infectionArray = np.zeros((precision,times))
    for s in range(times):
        for i in range(len(pImmune)):
            print(pImmune[i])
            t1=time.time()

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
            print(f"Time taken {i} == {time.time()-t1}s")
    
    np.savetxt(f"data/task7/Task7_RawData{name}.dat",infectionArray,fmt='%.6f')
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