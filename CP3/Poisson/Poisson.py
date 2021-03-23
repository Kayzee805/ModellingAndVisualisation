import numpy as np 
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation
class poisson(object):
    def __init__(self,n,method,epsilon=1,dx=1,sweeps=100000):
        super().__init__()
        self.n=n
        self.epsilon = epsilon
        self.dx=dx
        self.dt = (dx**2)/6
        self.rho = np.zeros((n,n,n))
        self.sweeps=sweeps
        self.lattice=np.zeros((n,n,n))
        self.method=method
        self.nextLattice=np.zeros((n,n,n))
    

    def setBoundaries(self):
        #set boundaries to 0?
        self.lattice[0,:,:]= 0
        self.lattice[:,0,:]= 0
        self.lattice[:,:,0]= 0
        self.lattice[-1,:,:]= 0
        self.lattice[:,-1,:]= 0
        self.lattice[:,:,-1]= 0

    
    def setPointCharge(self):
        midPoint = int(self.n/2)
        self.rho[midPoint,midPoint,midPoint]=1

    def setChargedWire(self):
        midPoint = int(self.n/2)
        for z in range(1,self.n-1):
            self.rho[midPoint,midPoint,z]=1
        

    def setBoundariesForArray(self,array):
        array[0,:,:]= 0
        array[:,0,:]= 0
        array[:,:,0]= 0
        array[-1,:,:]= 0
        array[:,-1,:]= 0
        array[:,:,-1]= 0
    def jacobiUpdate(self):
        #does update over n sweeps

        convergence=self.epsilon+1
        counter=0
        print(f"Starting jacobi update for {self.method}")
        mid = int(self.n/2)
        print(f"At mid = {self.rho[mid,mid,mid]} sum ={np.sum(self.lattice)}")
            
        while(convergence>=self.epsilon):
            t1=time.time()
            self.setBoundaries()
     
            sumBefore = np.sum((self.lattice))

            self.nextLattice= (self.rho+(
            np.roll(self.lattice,1,axis=0)+np.roll(self.lattice,-1,axis=0)
            +np.roll(self.lattice,1,axis=1)+np.roll(self.lattice,-1,axis=1)+
            np.roll(self.lattice,1,axis=2)+np.roll(self.lattice,-1,axis=2)))/6

            self.setBoundariesForArray(self.nextLattice)

            sumAfter = np.sum((self.nextLattice))
            convergence = np.abs(sumBefore-sumAfter)
            if(counter%100==0):
                print(f"counter={counter}  convergence={convergence}  sumBefore={sumBefore}  sumAfter={sumAfter} ")
            self.lattice = np.copy(self.nextLattice)
            counter+=1
        
        print("done updating")

        # while(convergence>=self.epsilon):
        #     #allNeighbours = np.roll(self.lattice,1,axis=0)+
        #     t1=time.time()
        #     for i in range(1,self.n-1):
        #         for j in range(1,self.n-1):
        #             for k in range(1,self.n-1):
        #                 #do update here
        #                 nn = self.calculateNearestNeighbours(i,j,k)
        #                 rho = self.rho[i,j,k]
        #                 self.nextLattice[i,j,k] = (nn+rho)/6
            
        #     convergence = abs(np.sum(self.lattice)-np.sum(self.nextLattice))
         
        #     if(counter%1==0):
        #          print(f"counter={counter}  convergence={convergence}  sumBefore={np.sum(self.lattice)}  sumAfter={np.sum(self.nextLattice)}")
        #     self.lattice= self.nextLattice.copy()
        #     counter+=1

    def gaussSeidelUpdate(self):
        #does update for one gauss?
        print()
        convergence =self.epsilon+1
        #do it in a while loop
        counter=0
        mid = int(self.n/2)
        while(convergence>=self.epsilon):
            t1=time.time()
            sumBefore = np.sum(self.lattice)
            for i in range(1,self.n-1):
                for j in range(1,self.n-1):
                    for k in range(1,self.n-1):
                        nn = self.calculateNearestNeighbours(i,j,k)
                        rho = self.rho[i,j,k]*self.dx**2
                        self.lattice[i,j,k] = (nn+rho)/6
            sumAfter =np.sum(self.lattice)
            convergence=np.abs(sumAfter-sumBefore)
            counter+=1
            if(counter%10==0):
                print(f"Counter={counter}  convergence={convergence}  time={time.time()-t1}   mid={self.lattice[mid,mid,mid]}   midAfter={self.nextLattice[mid,mid,mid]} midAfter={self.rho[mid,mid,mid]}")


            #do break checks here?


    def overRelaxationUpdate(self):
        #update for over relaxation
        for s in range(self.sweeps):
            print()

    def calculateNearestNeighbours(self,i,j,k):
        #find the 6 neighbours, so +- for all i,j,k
        n=self.n
        return self.lattice[(i+1)%n,j,k]+self.lattice[(i-1)%n,j,k]+self.lattice[i,(j-1)%n,k]+self.lattice[i,(j+1)%n,k]+self.lattice[i,j,(k-1)%n]+self.lattice[i,j,(k+1)%n]

    def getEx(self,i,j,k):
        gradient = self.lattice[(i+1)%self.n,j,k]- self.lattice[(i-1)%self.n,j,k]
        return -gradient/2

    def getEy(self,i,j,k):
        gradient = self.lattice[i,(j+1)%self.n,k]- self.lattice[i,(j-1)%self.n,k]
        return -gradient/2

    def getEz(self,i,j,k):
        gradient = self.lattice[i,j,(k+1)%self.n]- self.lattice[i,j,(k-1)%self.n]
        return -gradient/2

    def getAllEfield(self):
        n=self.n
        gradient = -(np.roll(self.lattice,1,axis=0)+np.roll(self.lattice,-1,axis=0)+np.roll(self.lattice,1,axis=1)+np.roll(self.lattice,-1,axis=1)+np.roll(self.lattice,1,axis=2)+np.roll(self.lattice,-1,axis=2))/2


    
    def getPotential(self):
        return self.lattice[:,:,int(self.n/2)]

    def plotEField(self):
        xGradient = np.zeros((self.n,self.n))
        yGradient = np.zeros((self.n,self.n))
        mid =int(self.n/2)
        for i in range(self.n):
            for j in range(self.n):
                xGradient[i,j] = self.getEx(i,j,mid)
                yGradient[i,j] = self.getEy(i,j,mid)

        #now normalise it?    
        xNorm = np.zeros((self.n,self.n))
        yNorm = np.zeros((self.n,self.n))
        for i in range(self.n):
            for j in range(self.n):
                norm = np.sqrt((xGradient[i,j]**2)+(yGradient[i,j]**2))
                xNorm[i,j] =xGradient[i,j]/norm
                yNorm[i,j]=yGradient[i,j]/norm
        
        ranges = np.arange(0,self.n,1)

        x,y = np.meshgrid(ranges,ranges)
        plt.quiver(y,x,xNorm,yNorm,linewidth=0.5)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("Electric field for one charge")
        plt.savefig("eFiledForCharge.png")
        plt.show()

    def normaliseArray(self,x,y):
        for i in range(self.n):
            for j in range(self.n):
                norm = np.sqrt((x[i,j]**2)+(y[i,j]**2)) 
                x[i,j]/=norm
                y[i,j]/=norm   
        return x,y
                
                
    def getMagneticFiled(self):
        print()

'''
Needs phi or the normal derivative of the phi?

phi: Dirichet
derivative: Neumann


To solve Delta sq. phi=-rho/e
convert into a suitable initial value probelmn: so look at dphi/dt = lap ^2 phi + p/e. Initial value problem 
can be solved through timestep? 

when we get the steady state, dphi/dt = 0
so solve the eq that we want. (the rhs) eq is solved in steady state

Since rhs = 0, we can use the version with - or +.  +=well behaved. - = will blow up? so wont go to steady state



Jacobi method = Time step the equation, with as large dt as possible. Larger the dt, quicker we get to the steady state
jacobi has slow convergence. takes N^2 steps, each step will take N^d? for new potential in each lattice

find error????


Algorithms:

    1. Gauss-seidel algorithm. 

        which is the same except that the same array of values is used for the
        previous estimate (on the rhs of the previous equation) and the next estimate (on the lhs). That is, the array is updated in-place, as opposed
        to creating a new array with values based on the old one. You should
        convince yourself that this algorithm converges to the same solution.


Jacobian

need phiOld and phiNew, need to do np.copy 

in gauss, only need one array



    2. Successive over-relaxation (SOR)

         LOOK AT ODE rather than pde
         dx/dt = g
         x(n+1)-x(n) dt g(x(n))
        x(n+1) = f(x(n))
        plotting this gives


'''