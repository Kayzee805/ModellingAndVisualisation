import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from scipy.stats import sem
import random as r
import sys

class Model(object):

    def __init__(self,n):
        self.n=n

    


    def update(self):

        for i in range(self.n):
            for j in range(self.n):
                print("update")
    



def animate(n,sweeps):

    model = Model(n)
    plt.figure()
    colour='magma'
    im=plt.imshow(model.lattice,animated=True,cmap=colour)
    plt.colorbar(im)
    plt.draw()
    plt.pause(0.0001)

    for i in range(sweeps):
        model.update()
        if i%10==0:
            plt.clf()
            im=plt.imshow(model.lattice,animated=True,cmap=colour)
            plt.colorbar(im)
            plt.draw()
            plt.pause(0.0001)
   
    plt.show()    

def main():
    if(len(sys.argv)!=3):
        print("Usage python main.py N")
        sys.exit()
    
    n=int(sys.argv[1])
    sweeps = 10000

   # animate(n,sweeps)
