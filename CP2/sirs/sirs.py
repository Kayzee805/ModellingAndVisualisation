import numpy as np
import random
from astropy.stats import jackknife_resampling
from astropy.stats import jackknife_stats

class sirs(object):
    def __init__(self,size,pS,pI,pR):
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
        self.setRandom()
        

    def setRandom(self):
        counter=0
        for i in range(self.size):
            for j in range(self.size):
                r= random.random()
                if(r<1/3):self.lattice[i,j]=0   #S
                elif(r<2/3):
                    self.lattice[i,j]=-1 #I
                    counter+=1
                else:self.lattice[i,j]=1  #R
        self.infected=counter
    
    def setImmune(self,pImmune):
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



    def nearestNeighboursInfected(self,i,j):
        #Here I check if the nearest neighbour is infected or not
        #top bot left and right
        #returning true means is infected
        #returning false means not infected, so Sus or Rec
        size = self.size
        top=self.lattice[i,(j-1)%size]
        bot= self.lattice[i,(j+1)%size]
        left = self.lattice[(i-1)%size,j]
        right = self.lattice[(i+1)%size,j]
        if(top==-1 or bot ==-1 or left==-1 or right==-1):
            return True
        return False
    def isImmune(self,pI):
        counter=0
        for i in range(self.size):
            for j in range(self.size):
                r = random.random()
                if(r<=pI):
                    self.lattice[i,j]=2
                    counter+=1
        self.immune=counter

    def update(self):
        counter=0

        for i in range(self.size):
            for j in range(self.size):
                if(self.lattice[i,j]==-1):
                    counter+=1
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
                elif(value==-1):
                    #is infected
                    if(r<=self.p2):
                        self.lattice[iTrial,jTrial]=1
                elif(value==1):
                    #in recovery
                    if(r<=self.p3):
                        self.lattice[iTrial,jTrial]=0
                else:
                #    print("Immune so skip?")
                    continue
        self.infected = counter
    
    def jacknife(self,data):
        N=self.size*self.size
        length=len(data)
        newData = np.zeros(length)
        resamples=np.empty([length,length-1])
        for i in range(length):
            resamples[i] = np.delete(data,i)
            newData[i] = np.var(resamples[i])/N
        result=0.0

        originalC = np.var(data)/N   
        for i in newData:
            result+= (i-originalC)**2
        return np.sqrt(result)

    def jackknife2(self,data, n):
        resamples = jackknife_resampling(data)
        x_resamples = [np.var(resamples[i])/(n*n) for i in range (len(resamples))]
        v = np.var(data)/(n*n)
        return np.sqrt(np.sum([(x_resamples[i] - v) * (x_resamples[i]-v) for i in range (len(x_resamples))]))
'''
onnly 4 nearest neighbours. top,bot,left and right

1. S-> I if, one of the nn is infected. after 1 update
2.  I->R if, after x updates? ??
3. R->S if, after y updates. 

Parrallel Deterministics:
-------------------------
All  agents are updated at the same time, similar to gol
so cant have a cell change and affect near by cells

Sequenctial deterministic version:
----------------------------------
Agents change turn wise turn, which means that changing one agent can change the 
outcome of the next agent

Randoms equenctial algorithm:
------------------------------

1. Pick a random cell in the lattice with uniform probability
2. Rules:
    a. S->I, with proability p1 if one nn is infected
    b. I->R, with probability p2.
    c. R->S, with probability p3

What are the possible steady states/behaviours?


System takes 4 args. e.g. main.py N p1 p2 p3

1. Absorbing state:
------------------
    if all cells are S, I=0, nothing happens
    <fraction of infected cells> = <I>/N. <I> = num of infected cells?

2. Dynamic equilibrium:
-----------------------
Basically repeats?
<I>/N != 0.  
but constant over time/ will just fluctuate a tiny bit

p1=p2=p3 does not hold


3. Travelling waves of infection spreading:
------------------------------------------

similar to glider? because of oscillator.
need n.n to spread infection
<I>/N !=0 but not constant over time. 
so waves. 



'''