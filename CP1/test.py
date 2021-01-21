import numpy as np
import time
import random
from astropy.stats import bootstrap
random.seed(10)
np.random.seed(10)


test = np.ones((50,50))
mid = int(50/2)
for i in range(mid,50):
    for j in range(0,50):
        test[i][j]=-1

print(np.sum(test))

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