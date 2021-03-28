import numpy as np 
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation
class poisson(object):
    def __init__(self,n,method='gauss',epsilon=1,dx=1,sweeps=100000,minW=0,maxW=0):
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
        self.minW=minW
        self.maxW=maxW
    

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
            if(counter%1==0):
                print(f"counter={counter}  convergence={convergence}  sumBefore={sumBefore}  sumAfter={sumAfter} ")
            self.lattice = np.copy(self.nextLattice)
            counter+=1
        
        print(f"Finished jacobi at {convergence} at step={counter}")
        
        print("done updating")

    def gaussSeidelUpdate(self):
        #does update for one gauss?
        print()
        convergence =self.epsilon+1
        #do it in a while loop
        counter=0
        mid = int(self.n/2)
        print("staritng gasuss seidel")
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
            if(counter%1==0):
                print(f"time = {time.time()-t1}s Counter={counter}  convergence={convergence}  time={time.time()-t1}   mid={self.lattice[mid,mid,mid]}   midAfter={self.nextLattice[mid,mid,mid]} midAfter={self.rho[mid,mid,mid]}")


    def gaussSeidelUpdateV2(self):
        #so generate all the +1??
        convergene=self.epsilon+1
        counter=0


    def overRelaxationUpdate(self,omega):
        #update for over relaxation
        convergence = self.epsilon+1
        counter=0

        while(convergence>=self.epsilon):
            sumBefore = np.sum(self.lattice)

            for i in range(1,self.n-1):
                for j in range(1,self.n-1):
                    for k in range(1,self.n-1):
                        #omega between 0 and 2
                        #lattice[n+1] = (1-w)lattice[n] +w*lattice[n+1]
                        before=self.lattice[i,j,k]
                        nn = self.calculateNearestNeighbours(i,j,k)
                        rho = self.rho[i,j,k]*self.dx**2
                        after= (nn+rho)/6
                        self.lattice[i,j,k] = (((1-omega)*before)+(omega*after))

            sumAfter = np.sum(self.lattice)
            convergence = np.abs(sumAfter-sumBefore)
            if(counter%10==0):
                print(f"counter={counter}  convergence={convergence} before={sumBefore}  after={sumAfter}")
            counter+=1

        print(f"DONEEEE convergence= {convergence}  epsilon={self.epsilon}")
        return counter

    def overRelaxationUpdate_all(self):
        amount=11
        w = np.linspace(self.minW,self.maxW,amount)
        print(f"Omega = {w}")
        stableSweep=[]
        for i in range(amount):
            print(f"Starting iteration={i}")
            self.lattice=np.zeros((self.n,self.n,self.n))
            sweeps = self.overRelaxationUpdate(w[i])
            stableSweep.append(sweeps)
            print(f"i={i} sweeps={sweeps}  w={w[i]}")
        
        np.savetxt("SOR.dat",np.array((w,stableSweep)))

        plt.scatter(w,stableSweep,s=10,color='k')
        plt.plot(w,stableSweep,color='b')
        plt.title(f"SOR plot for Omega={self.minW}->{self.maxW}")
        plt.xlabel("omega")
        plt.ylabel("sweeps till stable")
        plt.savefig("sor.png")
        plt.show()

    
    def set2D_point(self):
        self.lattice=np.zeros((self.n,self.n))
        self.rho = np.zeros((self.n,self.n))
        mid = int(self.n/2)
        self.rho[mid,mid]=1
    
    def neighbours_2D(self,i,j):
        n=self.n
        return 1/4*(self.lattice[(i+1)%n,j]+self.lattice[(i-1)%n,j]+self.lattice[i,(j+1)%n]+self.lattice[i,(j-1)%n]+self.rho[i,j])
    def generate_SOR(self):
        w= np.arange(1,2,0.01)
        stableSweeps=[]
        print(w)
        self.n=100
        for x in w:
            self.set2D_point()
            convergence=self.epsilon+1
            counter=0
            print(f"Starting for w={x}")
            while convergence>=self.epsilon:
                convergence=0
                for i in range(1,self.n-1):
                    for j in range(1,self.n-1):
                        before = self.lattice[i,j]
                        after = self.neighbours_2D(i,j)
                        self.lattice[i,j] = (1-x)*before+x*after
                        convergence+=abs(self.lattice[i,j]-before)

                if(counter%1000==0):
                    print(f"counter={counter} convergence={convergence}")
                counter+=1
            print(f"Finished w={x} in {counter}sweeps")
            stableSweeps.append(counter)
        

        np.savetxt("SOR_DATA.dat",np.array((w,stableSweeps)))
        plt.scatter(w,stableSweeps,s=10,color='k')
        plt.plot(w,stableSweeps,color='b')
        plt.title(f"SOR plot for Omega={self.minW}->{self.maxW}")
        plt.xlabel("omega")
        plt.ylabel("sweeps till stable")
        plt.savefig("sor.png")
        plt.show()
                
                


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
        for i in range(1,self.n-1):
            for j in range(1,self.n-1):
                norm = np.sqrt((xGradient[i,j]**2)+(yGradient[i,j]**2))
               # print(norm)
                xNorm[i,j] =xGradient[i,j]/norm
                yNorm[i,j]=yGradient[i,j]/norm
        
        ranges = np.arange(0,self.n,1)

        x,y = np.meshgrid(ranges,ranges)
        plt.quiver(y,x,xNorm,yNorm,linewidth=0.5)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("Electric field for one charge")
        plt.savefig(f"eFiledForCharge_{self.method}.png")
        plt.show()

    def normaliseArray(self,x,y):
        for i in range(1,self.n-1):
            for j in range(1,self.n-1):
                norm = np.sqrt((x[i,j]**2)+(y[i,j]**2)) 
                x[i,j]/=norm
                y[i,j]/=norm   
        return x,y
                
                
    def getMagneticFiled(self):
        potential = self.getPotential()

        bx= np.zeros((self.n,self.n))
        by= np.zeros((self.n,self.n))
        for i in range(1,self.n-1):
            for j in range(1,self.n-1):
                #now do the gradient thing here
                bx[i,j] = (potential[i,(j+1)%self.n]-potential[i,(j-1)%self.n])/(2)
                by[i,j] = -(potential[(i+1)%self.n,j]-potential[(i-1)%self.n,j])/(2)
        
        bxNorm,byNorm = self.normaliseArray(bx,by)

        ranges = np.arange(0,self.n,1)

        x,y = np.meshgrid(ranges,ranges)
        plt.quiver(y,x,bxNorm,byNorm,linewidth=0.5)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("magnetic field for line of charge")
        plt.savefig(f"BFiledForWire_{self.method}.png")
        plt.show()

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