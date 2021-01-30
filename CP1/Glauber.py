
from Lattice import Lattice
import numpy as np
import random
class Glauber(Lattice):
    '''
    Glauber is a sub/child class of Lattice. It chas access to all methods in 
    lattice object and also, some additional ones as defined in the subclass.
    It has access to all of the variables as well, such as spin lattice and temperatures.
    It takes  the same arugements as lattice as well.
    '''
    def update(self):
        '''
        We randomly select a spin/point, then check the energyDifference if it was flipped
        if energyDifference <=0, we accept the flip else, we compare it with the 
        probability on it being accepted. 

        The following method will be called for each sweep, similar to the pseudocode provided
        '''
        for i in range(self.lx):
            for j in range(self.ly):

                #select a random spin in the lattice
                iTrial = random.randint(0,self.lx-1)
                jTrial = random.randint(0,self.ly-1)
                spin_new = -self.spin[iTrial,jTrial]

                #calculate the energy difference if flipped
                energyDiff = self.deltaE_Glauber(iTrial,jTrial)

                
                if(energyDiff>0):
                    #energy>0 therefore need to check the probability to see if accepted
                    probability = self.calculateProbability(energyDiff)
                    randomNumber = random.random()
                    if(probability>randomNumber):
                        #we accept the flip
                        self.spin[iTrial,jTrial] = spin_new
                    #else not accepted so dont do anything
                else:
                    #energy<=0 therefore is accepted
                    self.spin[iTrial,jTrial] = spin_new   

    def allOnes(self):
        '''
        Initialises the spin lattice to be all up spins.
        '''
        self.spin = np.ones((self.n,self.n))                 
                
