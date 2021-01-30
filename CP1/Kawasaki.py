
from Lattice import Lattice
import numpy as np
import time
import random
class Kawasaki(Lattice):

    '''
    Subclass of lattice and has methods that only kawasaki dynamic rule uses
    '''
    def swap(self,i,j,x,y):
        '''
        similar to being given coords of two points and we swap it
        however, since we know the elements are either +1 or -1, we 
        just multiply by -1
        '''
        #update the spin lattice with the points swapped
        self.spin[i,j] *= -1
        self.spin[x,y] *=-1
    
    def halfHalf(self):
        '''
        Initialises the first half of the spin lattice +1
        and the other half -1
        '''
        self.spin = np.ones((self.lx,self.ly))
        halfY = int(self.lx/2)


        for i in range(halfY,self.lx):
            for j in range(self.lx):
                self.spin[i,j] = -1

    def update(self):
        '''
        We randomly select two different points in the lattice,
        then compare the energy before swap and energy after swap
        then if energydiff is less than or equal to 0, we accept the swap
        else try for a probability 
        '''
        #method is carried out for each sweep.
        for i in range(self.lx):

            for j in range(self.ly):
                #we repeat until we get two distinct points in the lattice
                iOne = random.randint(0,self.lx-1)
                jOne = random.randint(0,self.ly-1)
                while True:
                    #we only look for one different point
                    iTwo = random.randint(0,self.lx-1)
                    jTwo = random.randint(0,self.ly-1)
                    #we only accept if spins are different and are point at two different points
                    #dont need to check if pointing at the same index because that would mean spin sum!=0
                    if((self.spin[iOne,jOne]+self.spin[iTwo,jTwo]==0)):
                        break
                    #else repeat
  
                #energy total before the swap
                energyOneBefore = self.energy_Kawasaki(iOne,jOne)
                energyTwoBefore = self.energy_Kawasaki(iTwo,jTwo)
                energyBefore = energyOneBefore+energyTwoBefore  

                #carry out swap
                self.swap(iOne,jOne,iTwo,jTwo)

                #check if nearest neighbours
                jDiff = np.abs(jOne-jTwo)%(self.ly-2)
                iDiff = np.abs(iOne-iTwo)%(self.ly-2)

                if(((iOne==iTwo) and jDiff==1) or ((iDiff==1) and (jOne==jTwo))):
                    #if nearestNeighbours, then we calculate the energy after swap
                    energyOneAfter = self.energy_Kawasaki(iOne,jOne)
                    energyTwoAfter = self.energy_Kawasaki(iTwo,jTwo)
                    energyAfter = energyOneAfter+energyTwoAfter

                    #then calculate energy difference
                    energyDiff = energyAfter-energyBefore 
                else:
                    #if not nearest neighbours, we simiplify the calculations by
                    #following a similar method to glauber
                    energyDiff = -2* energyBefore


                #no need to check for <= 0 because its already flipped
                if(energyDiff>0):
                    #carry out probability check
                    probability = self.calculateProbability(energyDiff)
                    randomNumber = random.random()
                    #it is already flipped, so we unflip if prob<randomNumber
                    #else: keep it flipped, so do nothing
                    if(probability<=randomNumber):
                        self.swap(iOne,jOne,iTwo,jTwo)
           
