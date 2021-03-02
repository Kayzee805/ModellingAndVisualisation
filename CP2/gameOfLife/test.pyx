import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import random
class gameOfLife(object):
    #random.seed(12345)
    def __init__(self,size, initialisation):
        self.lattice = np.zeros((size,size))
        self.size = size
        if(initialisation=="blinker"):
            self.setBlinker()
        elif(initialisation=="glider"):
            self.setGlider()
        else:
            self.setRandom()
        self.nextLattice = np.zeros((size,size))
       # random.seed(10012)
        self.activeSitesCounter=1
        self.currentActiveSite=0
 
    def newState(self,int i,int j,int nn):
        '''
        Live cell: 2< or >3 dies
                    only live cell with 2 or 3 live nn lives
        dead cell: only if nn= 3, cell goes to live
        nn==nearest neighbours
        '''
        cdef int cell = self.lattice[i,j]
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
        cdef int nn
        for i in range(self.size):
            for j in range(self.size):
                nn = self.nearestNeighbours(i,j)
              #  if(nn>0):
                  #  print(i,j)
                self.nextLattice [i,j]=self.newState(i,j,nn)
        self.lattice = self.nextLattice.copy()

    def setLattice(self,newLattice):
        self.lattice = newLattice
      #  print(f"size = {np.sum(self.lattice)}")


    def setBlinker(self):
        print("set blinker here")
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
        print("self lattice to glider")
        size = int(self.size/2)
        self.lattice[size+1,size]=1
        for i in range(size,size+3):
            self.lattice[i,size+2]=1
        self.lattice[size+2,size+1]=1

    
    def setRandom(self):
        for i in range(self.size):
            for j in range(self.size):
                r= random.random()
                if(r<0.5):self.lattice[i,j]=1
                else:self.lattice[i,j]=0

    def nearestNeighbours(self,int i,int j):
        #I guess I return the sum of the nearest neighbours 
        cdef int size,top,topRight,topLeft,bot,botRight,botLeft,left,right
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
        cdef int xLength = int(len(self.lattice[0])/2)
        cdef int yLength = int(len(self.lattice[1])/2)
        cdef int activeStates =len(allAlive[0])
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
        cdef float xCom = np.sum(allAlive[0])/activeStates
        cdef float yCom = np.sum(allAlive[1])/activeStates

        return [xCom,yCom]


def animate(int size,int sweeps,str initialisation):
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
       # model.centreOfMass()
        plt.cla()
        im=plt.imshow(model.lattice,animated=True,vmin=-1,vmax=1)
        plt.draw()
        plt.pause(0.1)
        print(model.activeSites())
  

    plt.show()

def generateHistogram(int size,int sweeps,str initialisation):
    cdef double t1
    cdef int newActiveSite, i,j
    absorbingState = []
    newActiveSite=0

    for i in range(1000):
        t1=time.time()
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
        print(f"Time taken for {i} == {time.time()-t1}s")       
    #print(f"Absorbing state = {absorbingState}")
    np.savetxt("data/histogramCython.dat",np.transpose(absorbingState),fmt='%.4f')
    plt.hist(absorbingState,bins=100)
    plt.show()

def generateCom(int size,int sweeps,str initialisation):
    xCom =[]
    yCom=[]
    t=[]
    model = gameOfLife(size,initialisation)
    for i in range(300):
        if(i%100==0):print(f"Iteration = {i}")
        com = model.centreOfMass()
        #com[0]= x, com[1] = y
        if(com[0]!=-1 and com[1]!=-1):
            #add com
            #i can probs keep it how it is and just slice an array later?
            xCom.append(com[0])
            yCom.append(com[1])
            t.append(i)
        model.update()

    arrayCombined = np.array((t,xCom,yCom))
    np.savetxt("data/centreOfMassCYTHON.dat",np.transpose(arrayCombined),fmt='%.4f')
#
def getVelocity(allArray):
    t=allArray[:,0]
    x=allArray[:,1]
    y=allArray[:,2]
    newX=[]
    newY=[]
    newT = []
    for i in range(len(t)):
        if(i>150 and i<250):
            newX.append(x[i])
            newY.append(y[i])
            newT.append(i)
    
    newX=np.asarray(newX)
    newY=np.asarray(newY)
    newT=np.asarray(newT)

    xfit,xin = np.polyfit(newT,newX,1)
    yfit,fin = np.polyfit(newT,newY,1)

    print(f"Xvel = {xfit}\nyVel = {yfit}")

    vel = np.sqrt(xfit**2+yfit**2)
    print(f"Velocity = {vel}")
    return vel


def plotAll():

    #data 1
    data1 = np.loadtxt("data/histogramCython.dat")
    plt.hist(data1,bins=25)
    plt.title("Cython histogram 1")
    plt.xlabel("Time step till absorbing state")
    plt.ylabel("Frequency")
    plt.savefig("figures/HistogramCython.png")
    plt.show()

    
    centreOfMass = np.loadtxt("data/centreOfMass3Cython.dat")
    t=centreOfMass[:,0]
    xCom=centreOfMass[:,1]
    yCom=centreOfMass[:,2]

    plt.scatter(t[::10],xCom[::10],s=8)
    plt.title("Scatter for xCom")
    plt.ylabel("centre of mass X")
    plt.xlabel("Time (Sweeps)")
    plt.savefig("figures/xComCython.png")
    plt.show()

    plt.scatter(t[::10],yCom[::10],s=8)
    plt.title("Scatter for yCom")
    plt.ylabel("centre of mass Y")
    plt.xlabel("Time (Sweeps)")
    plt.savefig("figures/yComCython.png")
    plt.show()

    plt.scatter(t[::10],xCom[::10],s=8,label="x")
    plt.scatter(t[::10],yCom[::10],s=8,label="y")
    plt.title("Scatter of both com")
    plt.xlabel("Time (Sweeps)")
    plt.ylabel("Centre of mass")
    plt.legend()
    plt.savefig("figures/bothCOMCython.png")
    plt.show()



'''
when is a system stable?
so if number or active states over the last 10 iteration == const?

if number of active states doesnt change and stays the same
for "10" steps.
it is an absorbant state

if I get glider, which changes the number of states
we can remove that from the result?

'''