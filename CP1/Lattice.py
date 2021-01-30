import numpy as np
import math
import random

class Lattice(object):
   # random.seed(10)

    def __init__(self,n,T,J=1.0, kb=1.0):
        '''
        Sets up the lattice model given some arguements. For this check point
        we are assuming a square lattice. However, can be changed by adding more
        arguments to the function.
        Parameters:
        ------------
        n: Type integer
           width or the height of the lattice.

        T: Type float/Integer
           Temperature of the system.
        J: Type float/Integer
           Default value is 1.0 unless stated otherwise
        
        kb:Type float
           Boltzmanns constant
           Default value is 1.0 unless stated otherwise
        '''
        #initialising object variables
        self.n = n      
        self.lx = n
        self.ly = n
        self.T = T      
        self.J = J   
        self.kb = kb
        #Initialising the spin lattice of the class object randoml

        self.spin =[] #spin states of the lattice

    def randomInit(self):
        spin = np.zeros((self.lx,self.ly),dtype=float)
        for i in range(self.lx):
            for j in range(self.ly):
                r = random.random()
                if(r<0.5):spin[i,j]=-1
                else:spin[i,j]=1
        self.spin = spin
    def setTemperature(self,newTemp):
        '''
        Setter for the temperature. Sets the temperature of the
        current object to newTemp.
        Parameters:
        ------------
        newTemp: Type float
                 Temperature of the system
        '''
        self.T = newTemp
        
    def deltaE_Glauber(self,i,j):
        '''
        Calculates the change in energy of a lattice after a flip.
        Made for glauber dynamic rule.
        
        Parameters:
        ------------
        spin: Type 2darray
              Each element in the array representing spin up or down of the system.
        
        i,j:    Type Integer
                indicies for an element in a 2d array.
        Returns:
        --------
        energyDiff: Type float
                    energy difference before and after flip.
        '''

        si = self.spin[i,j]
        
        #Finding the nearest neighbours of an element. 
        # '%self.lx' takes care of the periodic boundary conditions
        top =  self.spin[i,(j-1)%self.ly]
        bottom =  self.spin[i,(j+1)%self.ly]
        left =  self.spin[(i-1)%self.lx ,j]
        right =  self.spin[(i+1)%self.lx ,j]

        #Method shown in the lecture to avoid unnecessary computations
        energyDiff = 2*self.J*si*(top+bottom+left+right)
        return energyDiff
        
    def energy_Kawasaki(self,i,j):
        '''
        Calculates the energy of spin by using its nearest neighbours.
        Made for kawasaki dynamic rule.
        
        Parameters:
        ------------
        spin: Type 2darray
              Each element in the array representing spin up or down of the system.
        
        i,j:    Type Integer
                indicies for an element in a 2d array.
        
        Returns:
        --------
        energy:  Type float
                 Energy of a spin
        '''
        si =  self.spin[i,j]
        
        #nearast neighbours of spin[i,j] 
        #follows the same method as deltaE_Glauber
        top =  self.spin[i,(j-1)%self.ly]
        bottom =  self.spin[i,(j+1)%self.ly]
        left =  self.spin[(i-1)%self.lx ,j]
        right =  self.spin[(i+1)%self.lx ,j]
        
        #This method returns the energy of a spin, not the change in energy
        energy = -self.J*si*(top+bottom+left+right)
        return energy    
    
    def calculateProbability(self,energyDifference):
        '''
        Calculates the probability of a spin swap or flip.
        Parameters:
        ------------
        energyDifference: Type float
                          The energydifference before and after flip or swap
        
        Returns:
        --------
        Type float
        The probability of a spin swap or flip
        '''
        return np.exp((-energyDifference)/(self.kb*self.T))
    
    def totalMagnetisation(self):
        '''
        Calculates the absolute total magnetisation of a system.
        
        Returns:
        --------
        Type Integer
        The absolute magnetisation of the system.
        '''
        return np.abs(np.sum(self.spin))
    
    def calculuateSusceptibility(self,MagnetisationVariance):
        '''
        Calculates the Susceptibility of a system
        Parameters:
        ------------
        MagnetisationVariace:  Type float
                               The variance of magnetisation for a temperature. (Usually 10,000 sweeps)
        Returns:
        --------
        Type float
        Susceptibility of a system
        '''
        return (1/(self.lx*self.ly*self.kb*self.T))*MagnetisationVariance

    def calculateHeatCapacity(self,EnergyVariance):
        '''
        Calculates the specific heat capacity of a system
        Parameters:
        ------------
        EnergyVariance:  Type float
                               The variance of average total energy for a temperature. (Usually 10,000 sweeps)
        Returns:
        --------
        Type float
        specific heat capacity of a system
        '''
        return (1/(self.lx*self.ly*self.kb*(self.T**2)))*EnergyVariance

    def calculateVariance(self,data):
        '''
        Calculates the variance of an array.
        Parameters:
        ------------
        data:  Type float/integer array
               The array to calculate the variance of.
        Returns:
        --------
        Type float
        Variance of the input data
        '''
        return (np.mean(np.square(data))-(np.square(np.mean(data))))

    def totalEnergy(self):
        '''
        Calculates the total energy of the system.
        Returns:
        --------
        energy: Type float
                the total energy of the system
        '''
        energy = 0.0
        for i in range(self.lx):
            for j in range(self.ly):
                iup=i+1
                if(i==self.lx-1):iup=0
                jup=j+1
                if(j==self.ly-1):jup=0
                energy+= -self.J*self.spin[i,j]*(self.spin[iup,j]+self.spin[i,jup])


                # #using only the bottom and right nearest neighbour as mentioned in the lectures
                # bottom = self.spin[i,(j+1)%self.ly]
                # right = self.spin[(i+1)%self.lx,j]
                # #appending each energy to the total energy.
                # energy+= -self.J*self.spin[i,j]*(bottom+right)
        return energy


    def jacknife(self,energy):
        '''
        Calculates the error with the jacknife resampling method
        Parameters:
        ------------
        energy:  Type float/integer array
                 The array to calculate the error for.
        Returns:
        --------
        Type float/integer 
        Returns an error containing error for each element in the input array. 

        '''
        length = len(energy)

        #create a new array of the same size
        newC = np.zeros(length)
        #create a 2d array of the same length as energy but with one less element in each array. 
        resamples = np.empty([length,length-1])

        for i in range(length):
            #delete the ith element from the array 
            resamples[i]=np.delete(energy,i)
            #then calculate the new heat capacity without the ith element
            newC[i] = self.calculateHeatCapacity(self.calculateVariance(resamples[i]))
        

        result = 0.0
        #calculate the error
        meanNewC = np.mean(newC)
        for i in newC:
            result+= (i-meanNewC)**2

        #return the square root of the error
        return np.sqrt(result)


