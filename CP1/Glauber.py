
from Lattice import Lattice
import numpy as np
import random
class Glauber(Lattice):
    #using Algorithm as object so I dont have to initialise spin array
    #can use Algorithm for both glauber and kawasaki classes 
    #this way I only have to initialise spin array once
    #glauber will take same args as algorithm
    #np.random.seed(10)

    def update(self):
        for i in range(self.lx):
            for j in range(self.ly):
                iTrial = random.randint(0,self.lx-1)
                jTrial = random.randint(0,self.ly-1)
                spin_new = -self.spin[iTrial,jTrial]

                #not doing the new spin thing here because 
                #the diff in energy becomes 2*J*si
                #so -si is not needed
                energyDiff = self.deltaE(self.spin,iTrial,jTrial)
                #if energyDiff<= 0, can carryout swap else no
                if(energyDiff>0):
                    #energy>0 therefore need to check the probability to see if accepted
                    probability = self.calculateProbability(energyDiff)
                    randomNumber = random.random()
                    if(probability>randomNumber):
                        self.spin[iTrial,jTrial] = spin_new
                else:
                    #energy<=0 therefore is accepted
                    self.spin[iTrial,jTrial] = spin_new   

    def allOnes(self):
        self.spin = np.ones((self.n,self.n))                 
                
