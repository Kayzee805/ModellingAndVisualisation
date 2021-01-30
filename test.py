import numpy as np 
import random
import time


def deltaE_Glauber(spin,i,j):
    si=spin[i,j]
    length=len(spin)
    top = spin[i,(j-1)%length]
    bottom = spin[i,(j+1)%length]
    left = spin[(i-1)%length,j]
    right=spin[(i+1)%length,j]
    return 2*si*(top+bottom+left+right)

def energy_Kawasaki(spin,i,j):
    si=spin[i,j]
    length=len(spin)
    top = spin[i,(j-1)%length]
    bottom = spin[i,(j+1)%length]
    left = spin[(i-1)%length,j]
    right=spin[(i+1)%length,j]
    return -si*(top+bottom+left+right)

def calculateProbability(energyDiff,T):
    return np.exp((-energyDiff)/T)

def totalMagnetisation(spin):
    return np.abs(np.sum(spin))
def calculateHeatCapacity(energyVariance,T,size):
    return (1/(size*size*T**2))*energyVariance

def calculuateSusceptibility(MagnetisationVariance,T,size):
    return (1/(size*size*T))*MagnetisationVariance
def calculateVariance(data):
    return (np.mean(np.square(data))-(np.square(np.mean(data))))

def totalEnergy(spin):
    energy = 0.0
    length=len(spin)
    for i in range(length):
        for j in range(length):
            iup=i+1
            if(i==length-1):iup=0
            jup=j+1
            if(j==length-1):jup=0
            energy+= -1*spin[i,j]*(spin[iup,j]+spin[i,jup])

    return energy


    
def halfHalf():
    '''
    Initialises the first half of the spin lattice +1
    and the other half -1
    '''
    spin = np.ones((50,50))
    halfY = int(50/2)


    for i in range(halfY,50):
        for j in range(50):
            spin[i,j] = -1
    return spin

def swap(spin,i,j,a,b):
    spin[i,j] *=-1
    spin[a,b]*=-1

def updateKawasaki(spin):
    length=len(spin)
    for i in range():

        for j in range(length):
            #we repeat until we get two distinct points in the lattice
            iOne = random.randint(0,length-1)
            jOne = random.randint(0,length-1)
            while True:
                #we only look for one different point
                iTwo = random.randint(0,length-1)
                jTwo = random.randint(0,length-1)
                #we only accept if spins are different and are point at two different points
                #dont need to check if pointing at the same index because that would mean spin sum!=0
                if((spin[iOne,jOne]+spin[iTwo,jTwo]==0)):
                    break
                #else repeat

            #energy total before the swap
            energyOneBefore = energy_Kawasaki(spin,iOne,jOne)
            energyTwoBefore = energy_Kawasaki(spin,iTwo,jTwo)
            energyBefore = energyOneBefore+energyTwoBefore  

            #carry out swap
            swap(spin,iOne,jOne,iTwo,jTwo)

            #check if nearest neighbours
            jDiff = np.abs(jOne-jTwo)%(length-2)
            iDiff = np.abs(iOne-iTwo)%(length-2)

            if(((iOne==iTwo) and jDiff==1) or ((iDiff==1) and (jOne==jTwo))):
                #if nearestNeighbours, then we calculate the energy after swap
                energyOneAfter = energy_Kawasaki(spin,iOne,jOne)
                energyTwoAfter = energy_Kawasaki(spin,iTwo,jTwo)
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
                probability = calculateProbability(energyDiff)
                randomNumber = random.random()
                #it is already flipped, so we unflip if prob<randomNumber
                #else: keep it flipped, so do nothing
                if(probability<=randomNumber):
                    swap(spin,iOne,jOne,iTwo,jTwo)
        
def generateKawasaki(systemSize,Temperature,nSteps):
    TemperatureValues = 21
    tempList = np.linsapce(1,Temperature,TemperatureValues)
    specificHeat = np.zeros(TemperatureValues)
    averageEnergy =np.zeros(TemperatureValues)    
    spin = halfHalf()
    print(f"Magnetisation = {totalMagnetisation(spin)}")

    for t in range(TemperatureValues):
        currentTemp = tempList[t]
        energy=[]
        start = time.time()
        for n in range(nSteps):
            updateKawasaki(spin)
            if(n%500==0):
                print(n)
            if(n%10==0 and n>100):
                energy.append(totalEnergy(spin))
        
        averageEnergy[t] =np.mean(energy)
        specificHeat[t] = calculateHeatCapacity(calculateVariance(energy),currentTemp,50)
        print(f"Timet aken for T={currentTemp} is {time.time()-start} specificHeat = {specificHeat[t]}")

    combinedArray = np.array((tempList,averageEnergy,specificHeat))

    #write a transposed array to make it easier to read
    np.savetxt('data/kawasakiData.dat',np.transpose((combinedArray)),fmt='%.7f')   





