import numpy as np 

class poisson(object):
    def __init__(self,n,sweeps,method,epsilon=1,dx=1):
        super().__init__()
        self.n=n
        self.epsilon = epsilon
        self.dx=dx
        self.dt = (dx**2)/6
        self.rho = np.zeros((n,n,n))
        self.sweeps=sweeps
        self.lattice=np.zeros((n,n,n))
        self.method=method
        if(method=="jacobi"):
            self.nextLattice=np.zeros((n,n,n))
    

    def setBoundaries(self):
        #set boundaries to 0?
        self.lattice[0,:,:]=0
        self.lattice[self.n-1,:,:]=0
        self.lattice[:0,:]=0
        self.lattice[:self.n-1,:]=0
        self.lattice[::0,]=0
        self.lattice[::self.n-1,]=0

    
    def setPointCharge(self):
        midPoint = int(self.n/2)
        self.rho[midPoint,midPoint,midPoint]=1

    def setChargedWire(self):
        midPoint = int(self.n/2)
        for z in range(1,self.n-1):
            self.rho[midPoint,midPoint,z]=1
        

    def jacobiUpdate(self):
        #does update over n sweeps
        for s in range(self.sweeps):

            for i in range(1,self.n-1):
                for j in range(1,self.n-1):
                    for k in range(1,self.n-1):
                        #do update here
                        nn = self.calculateNearestNeighbours(i,j,k)
                        rho = self.rho[i,j,k]
                        self.nextLattice[i,j,k] = (nn+rho)/6
            self.lattice= self.nextLattice.copy()

            #do break checks here


    def gaussSeidelUpdate(self):
        #does update for one gauss?
        print()
        for s in range(self.sweeps):

            for i in range(1,self.n-1):
                for j in range(1,self.n-1):
                    for k in range(1,self.n-1):
                        nn = self.calculateNearestNeighbours(i,j,k)
                        rho = self.rho[i,j,k]*self.dx**2
                        self.lattice[i,j,k] = (nn+rho)/6


            #do break checks here?


    def overRelaxationUpdate(self):
        #update for over relaxation
        for s in range(self.sweeps):
            print()

    def calculateNearestNeighbours(self,i,j,k):
        #find the 6 neighbours, so +- for all i,j,k
        n=self.n
        return self.lattice[(i+1)%n,j,k]+self.lattice[(i-1)%n,j,k]+self.lattice[i,(j-1)%n,k]+self.lattice[i,(j+1)%n,k]+self.lattice[i,j,(k-1)%n]+self.lattice[i,j,(k+1)%n]


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