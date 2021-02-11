import numpy as np
import random
class gameOfLife(object):
    #random.seed(12345)

    def __init__(self,size, initialisation):
        self.lattice = np.zeros((size,size))
        self.size = size
        if(initialisation=="blinker"):
            self.setBlinker()
            print("set blinker here")

        elif(initialisation=="glider"):
            self.setGlider()
            print("Set glider here")
        else:
            self.setRandom()
            print("set random ehre")
        self.nextLattice = np.zeros((size,size))
       # random.seed(10012)
        self.activeSitesCounter=1
        self.currentActiveSite=0
 
    def newState(self,i,j,nn):
        '''
        Live cell: 2< or >3 dies
                    only live cell with 2 or 3 live nn lives
        dead cell: only if nn= 3, cell goes to live
        nn==nearest neighbours
        '''
        cell = self.lattice[i,j]
        if(cell==0):
            #dead
            if(nn==3):
                return 1
        elif(cell==1):
            #alive
            if(nn==2 or nn==3):
                return 1
        else:
            print("Cell should be either 0 or 1")
            exit(0)
        return 0
        
    def activeSites(self):
        return np.sum(self.lattice)

    def update(self):
        '''
        I can check the time to move arrays so lattice = new lattice
        and just creating a whole new array and compare the time
        
        '''
        for i in range(self.size):
            for j in range(self.size):
                nn = self.nearestNeighbours(i,j)
              #  if(nn>0):
                  #  print(i,j)
                self.nextLattice [  i,j]=self.newState(i,j,nn)
        self.lattice = self.nextLattice.copy()

    def setLattice(self,newLattice):
        self.lattice = newLattice
      #  print(f"size = {np.sum(self.lattice)}")


    def setBlinker(self):
        '''
        need to figure out how to do this
        '''
        size =int( self.size/2)
        # self.lattice[size,size]=1
        # self.lattice[size,size-1]=1
        # self.lattice[size,size+1]=1
        for i in range(size,size+3):
            for j in range(size,size+7):
                self.lattice[i,j]=1
    
    def setGlider(self):
        size = int(self.size/2)
        # self.lattice[size+1,size]=1
        # for i in range(size,size+3):
        #     self.lattice[i,size+2]=1
        # self.lattice[size+2,size+1]=1
        for i in range(size,size+3):
            self.lattice[size,i]=1
        self.lattice[size-2,size+1]=1
        self.lattice[size-1,size+2]=1
    
    def setRandom(self):
        for i in range(self.size):
            for j in range(self.size):
                r= random.random()
                if(r<0.5):self.lattice[i,j]=1
                else:self.lattice[i,j]=0

    def nearestNeighbours(self,i,j):
        #I guess I return the sum of the nearest neighbours 
        size = self.size
        top=self.lattice[i,(j-1)%size]
        topRight=self.lattice[(i+1)%size,(j-1)%size]
        topLeft = self.lattice[(i-1)%size,(j-1)%size]
        bot= self.lattice[i,(j+1)%size]
        botRight = self.lattice[(i+1)%size,(j+1)%size]
        botLeft = self.lattice[(i-1)%size,(j+1)%size]
        left = self.lattice[(i-1)%size,j]
        right = self.lattice[(i+1)%size,j]

        return top+topRight+topLeft+bot+botRight+botLeft+left+right


    def centreOfMass(self):
     #   print("carry out method for it")
        #should return all alive states?
        #this returns two arrays, one for all x and one for all y
        #[[x1,x2,x4,x5],[y1,y2,y3,y4,y5]]
        allAlive = np.where(self.lattice==1)
        xLength = int(len(self.lattice[0])/2)
        yLength = int(len(self.lattice[1])/2)
        activeStates =len(allAlive[0])
        #or can use the function that I made..

        #returning negative values because 
        #x boundary check
        if(abs(np.amax(allAlive[0])-np.amin(allAlive[0]))>xLength):
            return [-1,-1]
        
        #y boundary check
        if(abs(np.amax(allAlive[1])-np.amin(allAlive[1]))>yLength):
            return [-1,-1]
        
        #can calculate the com now, for x and y
        #because the numerator is just the sum of the indexs in the lattice?

        '''
        Only question is what if indexing starts from 1 in the equation?
        but here it starts from 0
        so Do i need to add 1? or add active states to compensate for that?
        '''
        xCom = np.sum(allAlive[0])/activeStates
        yCom = np.sum(allAlive[1])/activeStates

        return [xCom,yCom]

'''
when is a system stable?
so if number or active states over the last 10 iteration == const?

if number of active states doesnt change and stays the same
for "10" steps.
it is an absorbant state

if I get glider, which changes the number of states
we can remove that from the result?

'''