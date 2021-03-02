import numpy as np
import random
class gameOfLife(object):
    #random.seed(12345)

    def __init__(self,size, initialisation):
        '''
        Initialisation of a game of life object. A lattice can be initialisaed randomly,
        for a glider or a blinker.#
        Parameters:
        -------------
        Size: Type integer
              Size of the system/lattice, system will be of Size x Size
        initialisation: Type Integer
                        Initialisation method for the system, can be 2=glider,1= blinker or 0=random
                        Default:0 random
        '''

        #Lattice of an object, of size (size*size)
        self.lattice = np.zeros((size,size))
        self.size = size

        #Initialising the lattice over different states. Default: random
        if(initialisation==1):
            self.setBlinker()
            print("Initialise blinker lattice")

        elif(initialisation==2):
            self.setGlider()
            print("Initialise glider lattice")
        else:
            self.setRandom()
            print("Initialise random lattice")
        
        #Since system is fully determinsitic, we go index by index, we have another lattice
        #which stores the values of the next lattice, after update has been done.
        #this saves us from creating a new size by size array every single update
        self.nextLattice = np.zeros((size,size))

        #variables to keep track of the state of the system.        
        self.activeSitesCounter=1
        self.currentActiveSite=0
 
    def newState(self,i,j,nn):
        '''
        Determines if a state is alive or dead by following the rules of game of life.
        if cell is alive, it will only survive for next turn/state if it has 2 or 3 neighbours.
        If cell is dead, it will come back alive if it has 3 neighbours.
        Neighbours are cells surrounding a cell at index i,j.
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
        '''
        Returns:
        --------
        The number of live cells in the lattice
        As live cells have a value of 1, and dead 0,
        we just return the sum of the lattice.
        '''
        return np.sum(self.lattice)

    def update(self):
        '''
        Updates the lattice in a fully deterministic way.
        Uses nextLattice to store the new values of each state, then assigns 
        the original lattice to this new nextLattice.
        '''
        for i in range(self.size):
            for j in range(self.size):
                #calculates the amout of neighbours.
                nn = self.nearestNeighbours(i,j)
                self.nextLattice[i,j]=self.newState(i,j,nn)
        
        #setting the nextLattice to lattice. 
        self.lattice = self.nextLattice.copy()


    def setBlinker(self):
        '''
        Sets the lattice into a small system so the final state of the system is a 
        blinker/oscillator. 
        '''
        size =int( self.size/2)
        for i in range(size,size+3):
            for j in range(size,size+7):
                self.lattice[i,j]=1
    
    def setGlider(self):
        '''
        Initialising the lattice so there is one glider that glides across the system
        '''
        size = int(self.size/2)
        for i in range(size,size+3):
            self.lattice[size,i]=1
        self.lattice[size-2,size+1]=1
        self.lattice[size-1,size+2]=1
    
    def setRandom(self):
        '''
        Randomly initialising the lattice with live and dead cells. Each have a probability
        of 0.5.
        '''
        for i in range(self.size):
            for j in range(self.size):
                r= random.random()
                if(r<0.5):self.lattice[i,j]=1
                else:self.lattice[i,j]=0

    def nearestNeighbours(self,i,j):
        '''
        Parameters:
        -----------
        i,j: Type Integer
             ith and jth index of a lattice
        Returns:
        --------
        The sum of the neighbours. Neighbours consists of diagonal neighbours as well.
        So topleft,topRight,botLeft and botRight and the normal neighbours, top,bot,left and right
        Boundary conditions have been taken core of by taking the modulus, %. 
        '''
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
        '''
        Returns the centre of mass for both x and y components of a glider.
        The method mentioned in the lectures has been used to compute the centre of mass.
        '''

        #this returns two arrays, one for all x and one for all y
        #[[x1,x2,x4,x5],[y1,y2,y3,y4,y5]]
        #the shape of the array would be (2,numberOfActiveStates), where 2 is for x and y
        allAlive = np.where(self.lattice==1)

        xLength = int(len(self.lattice[0])/2)
        yLength = int(len(self.lattice[1])/2)
        #print(f"Length = {xLength}")
        #total amount of active states.
        activeStates =len(allAlive[0])

 
        #x boundary check, return [-1,-1] if glider crossing through boundaries
        if(abs(np.amax(allAlive[0])-np.amin(allAlive[0]))>xLength):
            return [-1,-1]
        
        #y boundary check, return [-1,-1] if glider crossing through boundaries
        if(abs(np.amax(allAlive[1])-np.amin(allAlive[1]))>yLength):
            return [-1,-1]
        
        #can calculate the com now, for x and y

        #method to calculate the centre of mass for x and y components.
        xCom = np.sum(allAlive[0])/activeStates
        yCom = np.sum(allAlive[1])/activeStates

        return [xCom,yCom]