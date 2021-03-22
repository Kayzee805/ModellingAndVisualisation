import numpy as np
class cahnHilliard(object):

    def __init__(self,size,phi0,dx,dt,noise=0):
        self.n=size
        self.phi0=phi0
        self.dx = dx
        self.dt = dt
        self.a = 0.1
        self.M = 0.1
        self.k =0.1
        self.lattice = np.full((size,size),phi0)
        if(noise!=0):
            noiseArray=np.random.uniform(-0.1,0.1,(size,size))
            self.lattice = np.add(self.lattice,noiseArray)
        
        self.kdx = self.k/(self.dx**2)
        self.constants = (self.M*self.dt)/(self.dx**2)

    def update(self):
        #mu = np.zeros((50,50))
        mu = -self.a*self.lattice+self.a*self.lattice**3-self.kdx*(np.roll(self.lattice,1,axis=0)+np.roll(self.lattice,-1,axis=0)+np.roll(self.lattice,1,axis=1)+np.roll(self.lattice,-1,axis=1)-(4*self.lattice))
        self.lattice += self.constants*(np.roll(mu,1,axis=0)+np.roll(mu,-1,axis=0)+np.roll(mu,1,axis=1)+np.roll(mu,-1,axis=1)-(4*mu))
    
    def calculateFreeEnergy(self):
        freeEnergyDensity = -(self.a/2)*self.lattice**2+(self.a/4)*self.lattice**4+(self.k/(8*self.dx*self.dx))*(np.roll(self.lattice,1,axis=0)-np.roll(self.lattice,-1,axis=0)+np.roll(self.lattice,1,axis=1)-np.roll(self.lattice,-1,axis=1))**2
        freeEnergy = np.sum(freeEnergyDensity)
        return freeEnergy
