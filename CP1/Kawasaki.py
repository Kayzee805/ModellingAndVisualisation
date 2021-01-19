
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
        for i in range(self.lx):
            for j in range(halfY):
                self.spin[i,j] = -1

    def update(self):
        for i in range(self.lx):
            for j in range(self.ly):

                iOne = np.random.randint(0,self.lx)
                jOne = np.random.randint(0,self.ly)
                '''
                Not sure if i need to do it until i find unique points
                or can I just go to next?
                '''
                iTwo = np.random.randint(0,self.lx)
                jTwo = np.random.randint(0,self.ly)
                if((self.spin[iOne,jOne]+self.spin[iTwo,jTwo]!=0) or (iOne ==iTwo and jOne==jTwo)):
                    continue

                # while True:
                #     #getting two different states

                #     iTwo = np.random.randint(0,self.lx)
                #     jTwo = np.random.randint(0,self.ly)
                    
                #     #making sure its not the same state
                #     if(iOne != iTwo and jOne!=jTwo):

                #         if((lhs+self.spin[iTwo,jTwo]) ==0):
                #             #making sure its not the same spin
                #             #valid random i and js
                #             break
                #print(f"Taken {count} tries")
                energyOneBefore = self.deltaEKawasaki(self.spin,iOne,jOne)
                energyTwoBefore = self.deltaEKawasaki(self.spin,iTwo,jTwo)
                energyBefore = energyOneBefore+energyTwoBefore  #+correction?
                
                self.flip(iOne,jOne,iTwo,jTwo)

                energyOneAfter = self.deltaEKawasaki(self.spin,iOne,jOne)
                energyTwoAfter = self.deltaEKawasaki(self.spin,iTwo,jTwo)
                energyAfter = energyOneAfter+energyTwoAfter

                energyDiff = energyAfter-energyBefore 
                #no need to check for <= 0 because its already flipped
                if(energyDiff>0):
                    #need to do probability check
                    probability = self.calculateProbability(energyDiff)
                    randomNumber = random.random()
                    #it is already flipped, so we unflip if prob<randomNumber
                    #else: keep it flipped
                    if(probability<=randomNumber):
                        self.flip(iOne,jOne,iTwo,jTwo)


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

              
                

