import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from IPython.core.display import clear_output
from tqdm import tqdm
from scipy.stats import sem


'''
phi=1 for r<R. r = distance to centre and R is a given Radius
phi = 0 for r>=R

'''
class Model(object):

    def __init__(self,n,R,dt,D=1,alpha=1,dx=1,twoD=True,noise=False):
        self.n=n
        self.R = R
        self.dt = dt
        self.D = D
        self.alpha = alpha
        self.dx = dx 
        self.mid = int(self.n/2)
        if twoD:
            self.lattice = np.zeros((n,n))
            self.initialiseLattice()
        
        if noise:
            self.lattice=np.ones((n,n))
            noiseArray = np.random.uniform(-0.1,0.1,(n,n))
            self.lattice+= noiseArray
        
    
    def calculateDistance(self,i,j):
        return np.sqrt((i-self.mid)**2+(j-self.mid)**2)
    def initialiseLattice(self):
        for i in range(self.n):
            for j in range(self.n):
                distance = self.calculateDistance(i,j)
                if(distance<self.R):
                    self.lattice[i,j]=1
                #else its already 0
    
    def boundaryConditions(self):
        self.lattice[0]=1
        self.lattice[-1] = self.lattice[-2]
    def initialise1DLattice(self):
        self.lattice = np.zeros((self.n))
        x0 =int( self.n/10)

        for i in range(self.n):
            if i<= x0:
                self.lattice[i]=1
        
        self.boundaryConditions()
        

    def initialiseTaskD(self,k):
        self.lattice = np.zeros((self.n))
        self.lattice[0]=1
        for i in range(1,self.n):
            self.lattice[i] = np.exp(-k*i)

    def updateSingle(self):
        newLattice = np.zeros((self.n))
        for i in range(1,self.n-1):
            grad = self.D*(self.lattice[(i+1)%self.n]+self.lattice[(i-1)%self.n]-4*self.lattice[i])
            rhs = self.alpha*self.lattice[i]*(1-self.lattice[i])
            newLattice[i] = self.lattice[i]+(self.dt/(self.dx**2))*(grad+rhs)
        
        self.lattice = newLattice.copy()
        self.boundaryConditions()
        # grad = self.D*(np.roll(self.lattice,1,0)+np.roll(self.lattice,-1,0)-4*self.lattice)
        # rhs = self.alpha*self.lattice*(1-self.lattice)

        # self.lattice = self.lattice+ ((self.dt)/(self.dx**2))*(grad+rhs)
        # self.boundaryConditions()
        

    def calculateNearestNeighbours(self,i,j):
        allNeighbours = self.lattice[(i+1)%self.n,j]+self.lattice[(i-1)%self.n,j]+self.lattice[i,(j+1)%self.n]+self.lattice[i,(j-1)%self.n]-4*self.lattice[i,j]
        return allNeighbours
    def update(self):
        newArray = np.zeros((self.n,self.n))
        for i in range(self.n):
            for j in range(self.n):
                 newArray[i,j] = self.lattice[i,j]+((self.dt)/(self.dx)**2)*((self.D*self.calculateNearestNeighbours(i,j))+(self.alpha*self.lattice[i,j]*(1-self.lattice[i,j])))
        
        self.lattice= np.copy(newArray)
    def updateRoll(self):
        grad = self.D*(np.roll(self.lattice,1,0)+np.roll(self.lattice,-1,0)+np.roll(self.lattice,1,1)+np.roll(self.lattice,-1,1)-4*self.lattice)
        
        self.lattice = self.lattice+ (self.dt/(self.dx)**2)*(grad+((self.alpha*self.lattice)*(1-self.lattice)))

    def updateCahnHilliard(self,M,a,k):
        mu = a*self.lattice*(self.lattice-1)*(self.lattice-2)-(k/(self.dx**2))*(np.roll(self.lattice,1,0)+np.roll(self.lattice,-1,0)+np.roll(self.lattice,1,1)+np.roll(self.lattice,-1,1)-4*self.lattice)

        allCombined = (M/(self.dx**2))*(np.roll(mu,1,0)+np.roll(mu,-1,0)+np.roll(mu,1,1)+np.roll(mu,-1,1)-4*mu)+self.alpha*self.lattice*(1-self.lattice)

        self.lattice = self.lattice+ (self.dt)*(allCombined)
    
        
def taskE():
    n=30
    M=a=k=0.1
    alpha = 0.0003
    R=10
    dt=0.1
    model = Model(n,R,dt,alpha=alpha,noise=True,twoD=False)

    convergence=1
    before = np.sum(model.lattice)
    steps = np.linspace(1,500000,500000)
    for i in tqdm(steps):
        model.updateCahnHilliard(M,a,k)
  
    
    plt.figure()
    colour='magma'
    im = plt.imshow(model.lattice,cmap=colour)
    plt.colorbar(im)
    plt.savefig("figures/taskE.png")
    plt.show()

    np.savetxt("data/taskE.dat",model.lattice,fmt='%.4f')

def taskD():
    allK = np.arange(0.1,1,0.1)
    N=1000
    R=10
    dt=0.01
    velocities = []

    for i in range(len(allK)):
        print(i)
        steps = []
        integral=[]
        model = Model(N,R,dt,twoD=False)
        model.initialiseTaskD(allK[i])
        for j in range(500):
            model.updateSingle()
            steps.append(j)
            integral.append(np.sum(model.lattice))
        xfit,xin = np.polyfit(steps,integral,1)
        velocities.append(xfit)
    
    plt.scatter(allK,velocities,s=20,marker='+',color='r')
    plt.plot(allK,velocities)
    plt.xlabel("K")
    plt.ylabel("Velocity")
    plt.savefig("figures/taskD.png")
    plt.show()


def taskC():
    N=1000
    R=10
    dt= 0.01
    model = Model(N,R,dt,twoD=False)
    model.initialise1DLattice()

    steps = []
    integral = []
    for i in range(500):
        model.updateSingle()
        integral.append(np.sum(model.lattice)*model.dt)
        steps.append(i)
        if(i%100==0):
            print(i)
    
    plt.plot(steps,integral)
    plt.xlabel("steps")
    plt.ylabel("integral")
    plt.savefig(f"figures/taskC.png")
    plt.show()
    xfit,xin = np.polyfit(steps[0:100],integral[0:100],1)
    print(xfit,xin)

 


def animate():
    n=50
    R=10
    dt=0.01
    model = Model(n,R,dt)
    colour='magma'
    plt.figure()
    im = plt.imshow(model.lattice,animated=True,cmap=colour)
    plt.colorbar(im)
    plt.draw()
    plt.pause(0.0001)
    for i in range(10000):
        model.updateRoll()
        if(i%10==0):
            plt.clf()
            im = plt.imshow(model.lattice,animated=True,cmap=colour)
            plt.colorbar(im)
            plt.draw()
            plt.pause(0.0001)
            print(f"i={i}  sum={np.sum(model.lattice)}")
        
            if(i==0 or i==250 or i==1000 or i==2000):
                plt.savefig(f"figures/taskB_{i}.png")
    
    plt.show()


def main():
   # animate()
   #taskC()
   #taskD()
   taskE()

main()


# def animate(n,D1,D2,R,F,dt):
#     model= Model(n,D1,D2,R,F,dt)
#     colour='magma'
#     plt.figure()
#     im = plt.imshow(model.V,animated=True,cmap=colour)
#     plt.colorbar(im)
#     plt.draw()
#     plt.pause(0.0001)

#     print("starting animation")
#     convergence = 1
#     error =0.0001
#     before = np.sum(model.U)
#     i = 0
#     while convergence>error:
#         model.update()
#         if i%100==0:
#             plt.clf()
#             im=plt.imshow(model.V,animated=True,cmap=colour)
#             plt.colorbar(im)
#             plt.draw()
#             plt.pause(0.0001)
#             print(f"coutner={i} and convergence = {convergence}")
        
#         i+=1
#         after = np.sum(model.U)
#         convergence=abs(after-before)
#         before=after

#     plt.show()
'''

i=0  sum=305.0
i=10  sum=305.7048739175764
i=20  sum=307.3457922703377
i=30  sum=309.6467859534556
i=40  sum=312.4683125085087
i=50  sum=315.73063736224486
i=60  sum=319.3837043965842
i=70  sum=323.3937419422798
'''