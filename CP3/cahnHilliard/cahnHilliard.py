import numpy as np
#a cahnHilliard object class that will initialise the lattice and all the needed variables
class cahnHilliard(object):

    def __init__(self,size,phi0,dx,dt,noise=0,a=0.1,M=0.1,k=0.1):
        '''
        Constructor for a cahnHilliard object. Will initialise the lattice and all teh variables
        and will add the noise, if noise arguement in not zero.
        '''
        self.n=size
        self.phi0=phi0
        self.dx = dx
        self.dt = dt
        self.a = a
        self.M = M
        self.k =k
        #initialising the lattice with phi0
        self.lattice = np.full((size,size),phi0)

        #adding noise to the lattice
        if(noise!=0):
            noiseArray=np.random.uniform(-0.1,0.1,(size,size))
            self.lattice = np.add(self.lattice,noiseArray)
        #precomputing some cosntant calculations
        self.kdx = self.k/(self.dx**2)
        self.constants = (self.M*self.dt)/(self.dx**2)

    def update(self):
        '''
        cahnHilliard update method using np.roll
        First we precompute the cehcmical potential array
        Then use that precomputed chemical potential array alongside the numerical method of cahnHilliard method
        to update the lattice/model.
        
        '''
        mu = -self.a*self.lattice+self.a*self.lattice**3-self.kdx*(np.roll(self.lattice,1,axis=0)
            +np.roll(self.lattice,-1,axis=0)+np.roll(self.lattice,1,axis=1)+np.roll(self.lattice,-1,axis=1)
            -(4*self.lattice))


        self.lattice += self.constants*(np.roll(mu,1,axis=0)+np.roll(mu,-1,axis=0)+np.roll(mu,1,axis=1)+np.roll(mu,-1,axis=1)-(4*mu))
    
    def calculateFreeEnergy(self):
        '''
        Use np.roll to calculate the free energy Density and a discretised form of 'f' to calculate the free energy density
        then take the sum of it to calculate the freeEnergy of the lattice/model
        '''
        freeEnergyDensity = -(self.a/2)*self.lattice**2+(self.a/4)*self.lattice**4+(self.k/(8*self.dx*self.dx))*(np.roll(self.lattice,1,axis=0)-np.roll(self.lattice,-1,axis=0)+np.roll(self.lattice,1,axis=1)-np.roll(self.lattice,-1,axis=1))**2
        freeEnergy = np.sum(freeEnergyDensity)
        return freeEnergy
