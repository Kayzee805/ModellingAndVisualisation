import numpy as np
import time
import random
from astropy.stats import bootstrap
data= [1,2,3,4,5,6,7,1,25,12,5,1,2,5]
data = np.asarray(data)
def variance(data):
    mean = np.mean(data)
    sq = np.mean(data**2)
    meansq = mean ** 2
    v = sq - meansq
    return v#



def bootstrap2(data):
    newC = np.zeros(len(data))
    for i in range(len(data)):
        newC[i] = np.var(np.random.choice(data,len(data)))
    
    lhs = np.mean(np.square(newC))
    rhs = np.square(np.mean(newC))

    sigma = np.sqrt(lhs-rhs)
    return sigma

error = bootstrap2(data)
print(f"Error manual = {error}\n\n")




def bootstrapAstropy(data):
    return bootstrap(data,len(data),bootfunc=np.mean)
error2 = bootstrapAstropy(data)
print(f"Error auto = {error2}")

# fileNames = np.loadtxt('plotNames.dat',dtype='str')
# for x in fileNames:
# #     print((x))

# np.savetxt('data/test.dat',ran)



# test = [1,2,3,4,5,6,1]
# from astropy.stats import jackknife_resampling
# resample = jackknife_resampling(np.asarray(test))
# print(resample)
# print(np.var(test))


# test=[0.00222,0.021111,0.124125,0.1251251]
# dummy = [float('%.3f'%elem) for elem in test]

# print(dummy)
'''
Glauber
    Plot of average abs mag against T
    Average total energy at each T?
    specific heat plot 
    Susceptibility 
    Need to have data files for each? or one but know which stuff to plot


Kawasaki
    Average energy
    Specific heat
    need data files for each plot

Extra
    account for error bars with bootsrap or jacknife for specific heat

Need plot for average abs Mag agaisnt T Glauber

Susceptibility against T Glauber


'''