import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from IPython.core.display import clear_output
from tqdm import tqdm
from scipy.stats import sem
import random as r
'''
description:

select a random cell, and flip the sign
move is accepted or rejected acootding to the metropolis test

metropolis  p = exp(-deltaE/(kb*T))

if deltaE<0 then move is accepted
else:
    accepted if random number<p


average and the variance of the magnetisation

average and the variance of the stagered magnetisation

average of the energy of the system? 
'''

class Model(object):

    def __init__(self,n,h,T=1,kb=1,J=-1,equality=0.5):
        self.n = n
        self.h = h
        self.T=T
        self.kb = kb
        self.J=J
        self.equality = equality
        self.initialiseLattice()
        self.initialiseSngArray()

    def initialiseSngArray(self):
        self.staggered = np.ones((self.n,self.n))
        for i in range(self.n):
            for j in range(self.n):
                sgn = (-1)**((i+1)+(j+1))
                self.staggered[i,j] = sgn

    def initialiseLattice(self):
        self.lattice = np.ones((self.n,self.n))
        for i in range(self.n):
            for j in range(self.n):
                randNumber = r.random()
                if(randNumber<self.equality):
                    self.lattice[i,j]=-1
    
    def calculateDeltaE(self,i,j):
        top = self.lattice[(i-1)%self.n,j]
        bot = self.lattice[(i+1)%self.n,j]
        right = self.lattice[i,(j+1)%self.n]
        left = self.lattice[i,(j-1)%self.n]

        hBefore = -self.h*self.lattice[i,j]
        hAfter = -self.h*self.lattice[i,j]*-1
        hTotal = hAfter-hBefore
        return 2*self.J*self.lattice[i,j]*(top+bot+right+left)+hTotal
    
    def calculateProbability(self,energyDiff):
        return np.exp(-(energyDiff)/(self.kb*self.T))
    
    def calculateVariance(self,array):
        squareAverage = np.average(np.square(array))
        averageSquare = np.square(np.average(array))
        return squareAverage-averageSquare
    
    def calculateMagnetisation(self):
        return np.sum(self.lattice)
    
    def calcualteStaggeredMagnetisation(self):


        sgnArray = self.staggered*self.lattice
        test2= np.sum(sgnArray)
        return test2

    def calculateTotalEnergy(self):
        total=0
        for i in range(self.n):
            for j in range(self.n):
                bot = self.lattice[(i+1)%self.n,j]
                right = self.lattice[i,(j+1)%self.n]
                total+= -self.J*self.lattice[i,j] - self.h*self.lattice[i,j]
        return total

    
    def calculateField(self,P,tau,h0,x,y,n):
        x+=1
        y+=1
        h = h0*np.cos((2*np.pi*x)/P)*np.cos((2*np.pi*y)/P)*np.sin((2*np.pi*n)/tau)
        return h


    def updateD(self,P,tau,h0,n):
        
        for i in range(self.n):
            for j in range(self.n):
                jRandom = r.randint(0,self.n-1)
                iRandom = r.randint(0,self.n-1)
                newSpin =  -self.lattice[iRandom,jRandom]
                self.h = self.calculateField(P,tau,h0,iRandom,jRandom,n)
                energyDiff = self.calculateDeltaE(iRandom,jRandom)

                if(energyDiff>0):
                    #do metropolis check
                    probability = self.calculateProbability(energyDiff)
                    randNumber = r.random()
                    if(probability>randNumber):
                        self.lattice[iRandom,jRandom] = newSpin
                else:
                    self.lattice[iRandom,jRandom]=newSpin
    
    def update(self):
        
        for i in range(self.n):
            for j in range(self.n):
                jRandom = r.randint(0,self.n-1)
                iRandom = r.randint(0,self.n-1)
                newSpin =  -self.lattice[iRandom,jRandom]

                energyDiff = self.calculateDeltaE(iRandom,jRandom)

                if(energyDiff>0):
                    #do metropolis check
                    probability = self.calculateProbability(energyDiff)
                    randNumber = r.random()
                    if(probability>randNumber):
                        self.lattice[iRandom,jRandom] = newSpin
                else:
                    self.lattice[iRandom,jRandom]=newSpin
    

def generateDataE():
    n=50
    h0=10
    tau=10000
    P=25
    sweeps=10000
    model = Model(n,1)
    steps = np.linspace(1,sweeps,sweeps)
    i=0
    staggered = []
    fieldStrength=[]
    allSteps= []
    print(f"Task e for p={P}")
    for z in tqdm(steps):
        model.updateD(P,tau,h0,i)
        if i%1==0:
            staggered.append(model.calcualteStaggeredMagnetisation())
            f = h0*np.sin((2*np.pi*i)/tau)
            fieldStrength.append(f)
            allSteps.append(i)
        i+=1

    plt.plot(allSteps,fieldStrength)
    plt.xlabel("time(sweeps)")
    plt.ylabel("Field strength")
    plt.title("Field strength over time")
    plt.savefig(f"figures/taskE_fieldStrength.png")
    plt.show()

    plt.plot(allSteps,staggered)
    plt.xlabel("time(sweeps)")
    plt.ylabel("staggered magnetisation")
    plt.title(f"staggered magnetisation over time for P={P}")
    plt.savefig(f"figures/taskE_staggeredUpdateV2_P{P}.png")
    plt.show()


    plt.plot(fieldStrength,staggered)
   # plt.scatter(fieldStrength,staggered,s=1)
    plt.xlabel("time(sweeps)")
    plt.ylabel("staggered magnetisation")
    plt.title(f"staggered magnetisation over fieldStrength for P={P}")
    plt.savefig(f"figures/taskE_staggeredVfield_{P}.png")
    plt.show()
def generateDataD():
    P=25
    tau = 10000
    h0=10
    n=50
    sweeps=10000
    
    model = Model(n,1)
    steps = np.linspace(1,sweeps,sweeps)
    i =0
    allChecked=0
    magnetisation=[]
    allSteps=[]
    for z in tqdm(steps):
        model.updateD(P,tau,h0,i)
        '''
        when sin(2pi/tau) == -1,0,1
        -1 when 2*n/tau = 1/3?
        0 when 2*n/tau = 0 or just 1
        1 when 2*n/tau = 1/2
        '''
        if(i%1==0 and i>100):
            magnetisation.append(model.calculateMagnetisation())
            allSteps.append(i)

        # if((2*i)/tau)==0:
        #     plt.figure()
        #     im = plt.imshow(model.lattice,cmap='magma')
        #     plt.colorbar(im)
        #     plt.title("For sin(2*pi*n/tau)==0")
        #     plt.savefig(f"figures/taskD_0.png")
        #     plt.show(block=False)
        #     plt.pause(2)
        #     plt.close()
        # elif((2*i)/tau)==1/2:
        #     plt.figure()
        #     im = plt.imshow(model.lattice,cmap='magma')
        #     plt.colorbar(im)
        #     plt.title("For sin(2*pi*n/tau)==1")
        #     plt.savefig(f"figures/taskD_1.png")
        #     plt.show(block=False)
        #     plt.pause(2)
        #     plt.close()
        
        # elif((2*i)/tau)== 3/2:
        #     plt.figure()
        #     im = plt.imshow(model.lattice,cmap='magma')
        #     plt.colorbar(im)
        #     plt.title("For sin(2*pi*n/tau)==-1")
        #     plt.savefig(f"figures/taskD_-1.png")
        #     plt.show(block=False)
        #     plt.pause(2)
        #     plt.close()
        i+=1

    
    plt.plot(allSteps,magnetisation)
    plt.xlabel("time(sweeps)")
    plt.ylabel("magnetisation")
    plt.title("Task D spin pattern as a function of time")
    plt.savefig("figures/taskD_spinPatter.png")
    plt.show()
def generateDataC():
    allH = np.arange(0,10.5,0.5)
    magnetisationVaraince = []
    staggeredVariance = []
    totalEnergy=[]
    n=50
    sweeps=5000
    counter=0
    averageMagnetisation=[]
    averageStaggered=[]
    for z in tqdm(allH):
        h = allH[counter]
        model = Model(n,h)

        magnetisation=[]
        staggeredMagnetisation=[]
        energy=[]
        print(f"current h = {h}")
        for i in  range(sweeps):
            model.update()
            if(i%10==0 and i>100):
                magnetisation.append(model.calculateMagnetisation())
                staggeredMagnetisation.append(model.calcualteStaggeredMagnetisation())
                energy.append(model.calculateTotalEnergy())

        print("\n\n")
        totalEnergy.append(np.average(energy))
        staggeredVariance.append(model.calculateVariance(staggeredMagnetisation))
        magnetisationVaraince.append(model.calculateVariance(magnetisation))
        averageMagnetisation.append(np.average(magnetisation))
        averageStaggered.append(np.average(staggeredMagnetisation))
        counter+=1
    
    plt.scatter(allH,magnetisationVaraince,s=10,color='k')
    plt.plot(allH,magnetisationVaraince)
    plt.xlabel("h")
    plt.ylabel("magnetisation variance")
    plt.title("h vs magnetistaiton variance")
    plt.savefig(f"figures/magnetisationVariance.png")
    plt.show()

    plt.scatter(allH,staggeredVariance,s=10,color='k')
    plt.plot(allH,staggeredVariance)
    plt.xlabel("h")
    plt.ylabel("Staggered magnetisation variance")
    plt.title("h vs Staggeredmagnetistaiton variance")
    plt.savefig(f"figures/StaggeredmagnetisationVariance.png")
    plt.show()

    plt.scatter(allH,averageMagnetisation,s=10,color='k')
    plt.plot(allH,averageMagnetisation)
    plt.xlabel("h")
    plt.ylabel("average Magnetisation ")
    plt.title("h vs averageMagnetisation ")
    plt.savefig(f"figures/averageMagnetisation.png")
    plt.show()

    plt.scatter(allH,averageStaggered,s=10,color='k')
    plt.plot(allH,averageStaggered)
    plt.xlabel("h")
    plt.ylabel("average Staggered ")
    plt.title("h vs averageStaggered ")
    plt.savefig(f"figures/averageStaggered.png")
    plt.show()

    plt.scatter(allH,totalEnergy,s=10,color='k')
    plt.plot(allH,totalEnergy)
    plt.xlabel("h")
    plt.ylabel("average total energy ")
    plt.title("h vs totalEnergy ")
    plt.savefig(f"figures/totalEnergy.png")
    plt.show()

    np.savetxt(f"data/taskC_magnetisation.dat",np.transpose(np.array((allH,averageMagnetisation,magnetisationVaraince))),fmt='%.5f')
    np.savetxt(f"data/taskC_staggered.dat",np.transpose(np.array((allH,averageStaggered,staggeredVariance))),fmt='%.5f')
    np.savetxt(f"data/taskC_totalEnergy.dat",np.transpose(np.array((allH,totalEnergy))),fmt='%.5f')




def main():
    #generateDataC()
    #generateDataD()
    generateDataE()
main()