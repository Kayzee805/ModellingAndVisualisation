import numpy as np
data= [1,2,3,4,5,6,7]
data = np.asarray(data)
def variance(data):
    mean = np.mean(data)
    sq = np.mean(data**2)
    meansq = mean ** 2
    v = sq - meansq
    return v#

print(f"Manual = {variance(data)}")
print(f"np = {np.var(data)}")

test = np.zeros(5)
ran=[0,1,2,3,4,5]
print(np.mean(ran))

# fileNames = np.loadtxt('plotNames.dat',dtype='str')
# for x in fileNames:
# #     print((x))

# np.savetxt('data/test.dat',ran)



# test = [1,2,3,4,5,6,1]
# from astropy.stats import jackknife_resampling
# resample = jackknife_resampling(np.asarray(test))
# print(resample)
# print(np.var(test))

test = np.linspace(1,3,21)
test1 = np.linspace(1,3,21)
test2 = np.linspace(1,3,21)

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