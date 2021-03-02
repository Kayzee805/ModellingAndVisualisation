import numpy as np
import random
from astropy.stats import jackknife_resampling
from astropy.stats import jackknife_stats

class sirs(object):

    '''
    A sirs object class that contains initialisation and update methods for the 
    sirs model
    '''
    def __init__(self,size,pS,pI,pR,isImmune=False,immuneProbability=0):
        '''
        pS=p1. probability of S->I
        pI=p2. probability of I->R
        pR=p3, probability of R->S
        '''
        '''
        Intialises the lattice for the sirs model and some variables here
        some args has default assignments to it as it wont be used for all
        methods. i.e. immune probability and isImmune is only used for 
        models with immunity in the lattice
        '''
        self.size=size
        self.p1=pS
        self.p2=pI
        self.p3=pR
        self.lattice = np.zeros((size,size))
        self.infected=0
        self.immune=0
        if(isImmune):
            #set lattice with immune fraction
            self.setImmune(immuneProbability)
        else:
            self.setRandom()
            
        

    def setRandom(self):
        '''
        Initialises the lattice for random values of 
        Susceptibility, Infection and Recovered
        It also updates the number of infected in the model
        '''
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
    
    def setImmune(self,immuneProbability):
        '''
        Initialises the lattice model with a fraction of it having a probability
        of being immune. The infected and immune counter is also udpated
        '''
        pRest = (1-immuneProbability)/3
        infectedCounter=0
        immuneCounter=0
      #  print(f"SETTTING IMMUNE = {immuneProbability}")
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
       # print(f"INFECTED = {infectedCounter} immune = {immuneCounter}")



    def nearestNeighboursInfected(self,i,j):
        '''
        Returns true if any nearest neighbours of a point in the lattice is infected.
        Top,Bottom, Left and Right of the point are considered as nearest neighbours.
        Return False if no nearest neighbour is infected.
        '''
        size = self.size
        top=self.lattice[i,(j-1)%size]
        bot= self.lattice[i,(j+1)%size]
        left = self.lattice[(i-1)%size,j]
        right = self.lattice[(i+1)%size,j]
        if(top==-1 or bot ==-1 or left==-1 or right==-1):
            return True
        return False


    def isImmune(self,pI):
        '''
        not being used rn
        '''
        counter=0
        for i in range(self.size):
            for j in range(self.size):
                r = random.random()
                if(r<=pI):
                    self.lattice[i,j]=2
                    counter+=1
        self.immune=counter

    def update(self):
       # counter=0
        infectedCounter =self.infected
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
                        #    counter+=1
                            infectedCounter+=1
                elif(value==-1):
                    #is infected
                    if(r<=self.p2):
                        self.lattice[iTrial,jTrial]=1
                        infectedCounter-=1
                      #  counter-=1
                elif(value==1):
                    #in recovery
                    if(r<=self.p3):
                        self.lattice[iTrial,jTrial]=0
                        
                else:
                    continue


        self.infected = infectedCounter
        #print(f"update ={infectedCounter} infectedCounter = {self.infected}")


        

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