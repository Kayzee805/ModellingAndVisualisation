import numpy as np
import math
import random
from astropy.stats import jackknife_resampling, bootstrap
class Lattice(object):
    def __init__(self,n,T):
        self.n = n      #lattice dimension. e.g. 50= 50x50
        self.lx = n
        self.ly = n
        self.T = T      #temperature
        self.J = 1.0    
        self.kb = 1.0
        #I can probably do the fol  lowing in one line with numpy
        #but keep it like this Until i have rest of the stuff working
        spin = np.zeros((self.lx,self.ly),dtype=float)
        for i in range(self.lx):
            for j in range(self.ly):
                r = random.random()
                if(r<0.5):spin[i,j]=-1
                else:spin[i,j]=1

        #print(spin)
        self.spin =spin #spin states of the lattice

    def setTemperature(self,newTemp):
        self.T = newTemp
        
    def deltaE(self,spin,i,j):
        #this only works for a squre lattice
        #as I use slef.n for both lx and ly
        #can easily change it by introducing self.lx and self.ly
        si = spin[i,j]
        top = spin[i,(j-1)%self.ly]
        bottom = spin[i,(j+1)%self.ly]
        left = spin[(i-1)%self.lx ,j]
        right = spin[(i+1)%self.lx ,j]
        e = 2*self.J*si*(top+bottom+left+right)
        return e
        
    def deltaEKawasaki(self,spin,i,j):
        #this only works for a squre lattice
        #as I use slef.n for both lx and ly
        #can easily change it by introducing self.lx and self.ly
        si = spin[i,j]
        top = spin[i,(j-1)%self.ly]
        bottom = spin[i,(j+1)%self.ly]
        left = spin[(i-1)%self.lx ,j]
        right = spin[(i+1)%self.lx ,j]
        e = -self.J*si*(top+bottom+left+right)
        return e    
    
    def calculateProbability(self,energyDifference):
        return np.exp((-energyDifference)/(self.kb*self.T))
    
    def totalMagnetisation(self):
        #apparently they want teh abs value??
        return np.abs(np.sum(self.spin))
    
    def calculuateSusceptibility(self,MagnetisationVariance):
        #here N is the number of spins? 
        return (1/(self.lx*self.ly*self.kb*self.T))*MagnetisationVariance

    def calculateHeatCapacity(self,EnergyVariance):
        
        x= (1/(self.lx*self.ly*self.kb*self.T*self.T))*EnergyVariance
       # print(f"Temperature {self.T}  and c = {x}")
        return x


    def totalEnergy(self):
        '''
        Two methods.
        One that I followed from the lectures where i add top and right
        but in his code, he used bottom and right?? 
        '''
        energy = 0.0
        for i in range(self.lx):
            for j in range(self.ly):
               # lecture
                iup=i+1
                if(i==self.lx-1):iup=0
                jup=j+1
                if(j==self.ly-1):jup=0
                energy+= -self.J*self.spin[i,j]*(self.spin[iup,j]+self.spin[i,jup])

                #mine
                # bottom = self.spin[i,(j+1)%self.ly]
                # right = self.spin[(i+1)%self.lx,j]
                # energy+= -self.J*self.spin[i,j]*(bottom+right)
        #can prob reduce computation by multiplying by -self.J once at the end?
        return energy

    def jacknife(self,energy):
        length = len(energy)
        newC = np.zeros(length)
        resamples = np.empty([length,length-1])
        for i in range(length):
            resamples[i]=np.delete(energy,i)
            newC[i] = self.calculateHeatCapacity(np.var(resamples[i]))
        result = 0.0
        meanNewC = np.mean(newC)
        for i in newC:
            result+= (i-meanNewC)**2

        print(f"Error = {np.sqrt(result)} at temp = {self.T}")
        return np.sqrt(result)




    # def jacknife(self,energy,c):
    #     #c here is the original heatcapacity?
    #     #print("TESTINGG")
    #     resamples = jackknife_resampling(np.asarray(energy))
    #     newC = np.zeros(len(resamples))
    #     for i in range(len(resamples)):
    #         newC[i] = self.calculateHeatCapacity(np.var(resamples[i]))

    #     result = 0
    #     for i in range(len(energy)):
    #         result += (newC[i]-c)**2

    #     print(np.sqrt(result))
    #     return np.sqrt(result)

