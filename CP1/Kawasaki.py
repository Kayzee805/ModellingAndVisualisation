
from Lattice import Lattice
import numpy as np
import time
import random
class Kawasaki(Lattice):
    '''
    Below is glauber, change it to update of kawasaki
    '''
    def flip(self,i,j,x,y):
        self.spin[i,j] *= -1
        self.spin[x,y] *=-1
    
    def halfHalf(self):
        '''
        Make half of spin +1
        and the other half -1
        '''
        self.spin = np.ones((self.lx,self.ly))
        halfY = int(self.lx/2)


        for i in range(halfY,self.lx):
            for j in range(self.lx):
                self.spin[i,j] = -1

    def update(self):
        for i in range(self.lx):
            for j in range(self.ly):

                iOne = random.randint(0,self.lx-1)
                jOne = random.randint(0,self.ly-1)
                iTwo = random.randint(0,self.lx-1)
                jTwo = random.randint(0,self.ly-1)
                if((self.spin[iOne,jOne]+self.spin[iTwo,jTwo]!=0) or (iOne ==iTwo and jOne==jTwo)):
                    continue

  
                energyOneBefore = self.deltaEKawasaki(self.spin,iOne,jOne)
                energyTwoBefore = self.deltaEKawasaki(self.spin,iTwo,jTwo)
                energyBefore = energyOneBefore+energyTwoBefore  #+correction?
                self.flip(iOne,jOne,iTwo,jTwo)
                energyOneAfter = self.deltaEKawasaki(self.spin,iOne,jOne)
                energyTwoAfter = self.deltaEKawasaki(self.spin,iTwo,jTwo)
                energyAfter = energyOneAfter+energyTwoAfter

                energyDiff=energyAfter-energyBefore
                if((energyDiff)>0):
                    probability = self.calculateProbability(energyDiff)
                    randomNumber = random.random()
                    #it is already flipped, so we unflip if prob<randomNumber
                    #else: keep it flipped
                    if(probability<=randomNumber):
                        self.flip(iOne,jOne,iTwo,jTwo)  

                # energyDiff = energyAfter-energyBefore 
                # jDiff = np.abs(jOne-jTwo)%(self.ly-2)
                # iDiff = np.abs(iOne-iTwo)%(self.ly-2)
                # if(((iOne==iTwo) and jDiff==1) or ((iDiff==1) and (jOne==jTwo))):
                #    # print(iOne,iTwo,jOne,jTwo)
                #     energyOneAfter = self.deltaEKawasaki(self.spin,iOne,jOne)
                #     energyTwoAfter = self.deltaEKawasaki(self.spin,iTwo,jTwo)
                #     energyAfter = energyOneAfter+energyTwoAfter

                #     energyDiff = energyAfter-energyBefore 
                # else:
                #     energyDiff = -2* energyBefore

                # #no need to check for <= 0 because its already flipped
                # if(energyDiff>0):
                #     #need to do probability check
                #     probability = self.calculateProbability(energyDiff)
                #     randomNumber = random.random()
                #     #it is already flipped, so we unflip if prob<randomNumber
                #     #else: keep it flipped
                #     if(probability<=randomNumber):
                #         self.flip(iOne,jOne,iTwo,jTwo)

'''
energyDiff == energy before
need to flip and calc energy again
then take the diff 
thats energy diff
# '''
                # #check if neighbours
                # iDiff = np.abs(iOne-iTwo)%self.n
                # jDiff = np.abs(jOne-jTwo)%self.n   #need to mod self.n due to boundary

                # if(((iOne==iTwo) and (jDiff==1)) or ((iDiff==1) and (jOne==jTwo))):
                #     energyDiff -= -4  #whatever the correction is
                #     #think its -4 because we encourage swaps of nn
                #     #and each nn contributes to -2 in energy
                #     #as si = +-1 and corresponding nn = -+1

              
                

