import numpy as np 
import time
import matplotlib.pyplot as plt

class poisson(object):
    def __init__(self,n,epsilon=1,dx=1,minW=1,maxW=2):
        self.n=n
        self.epsilon=epsilon
        self.dx=dx
        self.rho = np.zeros((n,n,n))
        self.lattice=np.zeros((n,n,n))
        self.nextLattice=np.zeros((n,n,n))
        self.minW=minW
        self.maxW=maxW

    


    def setBoundaries3D(self,Array):
        #set boundaries to 0?
        Array[0,:,:]= 0
        Array[:,0,:]= 0
        Array[:,:,0]= 0
        Array[-1,:,:]= 0
        Array[:,-1,:]= 0
        Array[:,:,-1]= 0
    
    def setBoundaries2D(self,Array):
        Array[0,:]=0
        Array[-1,:]=0
        Array[:,0]=0
        Array[:,-1]=0 
    
    def setPointCharge2D(self):
        mid=int(self.n/2)
        self.rho[mid,mid]=1
    
    def setPointCharge3D(self):
        mid=int(self.n/2)
        self.rho[mid,mid,mid]=1
    
    def setChargedWire3D(self):
        mid=int(self.n/2)
        for z in range(1,self.n-1):
            self.rho[mid,mid,z]=1
    
    def calculateNearestNeighbours3D(self,i,j,k):
        #find the 6 neighbours, so +- for all i,j,k
        n=self.n
        return self.lattice[(i+1)%n,j,k]+self.lattice[(i-1)%n,j,k]+self.lattice[i,(j-1)%n,k]+self.lattice[i,(j+1)%n,k]+self.lattice[i,j,(k-1)%n]+self.lattice[i,j,(k+1)%n]

    def dependentNeighbours(self,i,j,k):
        n=self.n
        return self.lattice[(i-1)%n,j,k]+self.lattice[i,(j-1)%n,k]+self.lattice[i,j,(k-1)%n]

    def jacobiUpdate(self):
        print("Starting jacobi")
        convergence = self.epsilon+1
        counter=0
        sumBefore = np.sum(self.lattice)
        while(convergence>=self.epsilon):
            t1=time.time()
            self.setBoundaries3D(self.lattice)
            
            self.nextLattice= (self.rho+(
            np.roll(self.lattice,1,axis=0)+np.roll(self.lattice,-1,axis=0)
            +np.roll(self.lattice,1,axis=1)+np.roll(self.lattice,-1,axis=1)+
            np.roll(self.lattice,1,axis=2)+np.roll(self.lattice,-1,axis=2)))/6

            self.setBoundaries3D(self.nextLattice)
            sumAfter=np.sum(self.nextLattice)
            convergence = abs(sumAfter-sumBefore)
            self.lattice = np.copy(self.nextLattice)
            sumBefore=sumAfter
            counter+=1

            if(counter%1==0):
                print(f"Counter={counter} convergence={convergence} time={time.time()-t1}s")

        print("finished jacobi")

    
    def gaussSeidelUpdateOriginal(self):
        print("Starting gauss seidel")
        convergence = self.epsilon+1
        counter=0
        sumBefore = np.sum(self.lattice)
        while(convergence>=self.epsilon):
            t1=time.time()
            for i in range(1,self.n-1):
                for j in range(1,self.n-1):
                    for k in range(1,self.n-1):
                        nn = self.calculateNearestNeighbours3D(i,j,k)
                        self.lattice[i,j,k]=(self.rho[i,j,k]+nn)/6
            sumAfter=np.sum(self.lattice)
            convergence=abs(sumAfter-sumBefore)
            counter+=1
            sumBefore=sumAfter

            if(counter%1==0):
                print(f"Counter={counter} convergence={convergence} time={time.time()-t1}s")

        print("Finished gauss seidel")

    
    def gaussSeidelUpdate_roll(self):
        convergence = self.epsilon+1
        counter=0
        sumBefore = np.sum(self.lattice)
        while(convergence>=self.epsilon):
            t1=time.time()
            independentNeighbours = np.roll(self.lattice,-1,axis=0)+np.roll(self.lattice,-1,axis=1)+np.roll(self.lattice,-1,axis=2)+self.rho

            for i in range(1,self.n-1):
                for j in range(1,self.n-1):
                    for k in range(1,self.n-1):
                        self.lattice[i,j,k] = (independentNeighbours[i,j,k]+self.dependentNeighbours[i,j,k])/6
            sumAfter=np.sum(self.lattice)
            convergence = abs(sumAfter-sumBefore)
            counter+=1
            sumBefore=sumAfter
            
            if(counter%1==0):
                print(f"Counter={counter} convergence={convergence} time={time.time()-t1}s")