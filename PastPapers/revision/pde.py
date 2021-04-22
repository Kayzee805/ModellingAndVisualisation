import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from IPython.core.display import clear_output
from tqdm import tqdm
class PDE(object):

    def __init__(self,n,noise,dt=0.2,dx=1,rho=0.5,sigma=10,k=0.01,v0=0.01):
        self.n = n
        self.noise = noise 
        self.sigma = sigma 
        self.k = k
        self.D = 1 
        self.dx = dx
        self.dt = dt
        self.lattice = np.full((n,n),rho)
        noiseArray=np.random.uniform(-noise,noise,(n,n))
        self.lattice+= noiseArray
        self.p = np.zeros((n,n))
        self.initialiseSourceTerm()
        self.v0 = v0

        self.velocity = np.zeros((n,n))
        self.initialiseVelocity()

    def initialiseVelocity(self):
        for i in range(self.n):
            for j in range(self.n):
                self.velocity[i,j] = -self.v0*np.sin((np.pi*2*j)/self.n)

    def initialiseSourceTerm(self):
        for i in range(self.n):
            for j in range(self.n):
                self.p[i,j] = self.calculateSourceTerm(i,j)
    
    def calculateSourceTerm(self,x,y):
        mid = int(self.n/2)
        distance = np.sqrt((x-mid)**2+(y-mid)**2)
        return np.exp(-((distance**2)/(self.sigma**2)))

    def update(self):
        grad =  (self.D)*(np.roll(self.lattice,1,axis=0)+np.roll(self.lattice,-1,axis=0)+np.roll(self.lattice,1,axis=1)+np.roll(self.lattice,-1,axis=1)-4*self.lattice)
        self.lattice += (grad + self.p - self.k*self.lattice)*self.dt
        

    def velocityUpdate(self):
        grad =  (self.D)*(np.roll(self.lattice,1,axis=0)+np.roll(self.lattice,-1,axis=0)+np.roll(self.lattice,1,axis=1)+np.roll(self.lattice,-1,axis=1)-4*self.lattice)
        vel = self.velocity*((1/4)*(np.roll(self.lattice,-1,axis=0)-np.roll(self.lattice,1,axis=0)))

        self.lattice += (grad + self.p - self.k*self.lattice)*self.dt -vel*self.dt




def animate(n,noise,dt,dx=1,rho=0.5,sigma=10,k=0.01):
    model = PDE(n,noise)
    colour = 'magma'
    plt.figure()
    im=plt.imshow(model.lattice,animated=True,cmap=colour)
    plt.colorbar(im)
    plt.draw()
    plt.pause(0.001)
    for i in range(10000):
        model.update()
        if i%100==0:
            plt.clf()
            im=plt.imshow(model.lattice,animated=True,cmap=colour)
            c=plt.colorbar(im)
            plt.draw()
            c.update_normal(im)
            clear_output(wait=True)
            plt.pause(0.0001)
            print(i)
    
    plt.show()

def generateData3(n,noise,dt,dx=1,rho=0.5,sigma=10,k=0.01):
    model = PDE(n,noise,dt,dx,rho,sigma,k)
    sweeps=10000
    steps=[]
    averageValues =[]
    for i in range(sweeps):
        model.update()
        averageValues.append(np.mean(model.lattice))
        steps.append(i)

    plt.plot(steps,averageValues)
    plt.xlabel("time")
    plt.ylabel("average value")
    plt.savefig(f"figures/averageValues.png")
    plt.show()
    np.savetxt(f"data/averageValues.dat",np.transpose(np.array((steps,averageValues))),fmt='%.8f')


def generateData4(n,noise,dt,dx=1,rho=0.5,sigma=10,k=0.01):
    model = PDE(n,noise,dt,dx,rho,sigma,k)
    sweeps=3000

    for i in range(sweeps):
        model.update()
    
    mid=int(n/2)
    values=[]
    distance=[]
    for i in range(n):
        for j in range(n):  
            values.append(model.lattice[i,j])
            tempDistance = np.sqrt((mid-i)**2+(mid-j)**2)
            distance.append(tempDistance)

    distance = np.log(distance)
    values=np.log(values)
    plt.scatter(distance,values,s=5)
    plt.xlabel("distance to centre")
    plt.ylabel("value")
    plt.title("plot of value vs distance")
    plt.savefig("figures/distVsValue.png")
    plt.show()

    np.savetxt(f"data/valueVsDistance.dat",np.transpose(np.array((distance,values))),fmt='%.8f')
    #no fit it?

    minimum=1.8
    maximum=3
    newDistance = distance[(distance>minimum) & (distance<maximum)]
    newValues = values[(distance>minimum)&(distance<maximum)]
    xfit,xin = np.polyfit(newDistance,newValues,1)
    print(xfit,xin)
    plt.scatter(distance,values,s=5)
    plt.plot(distance,distance*xfit+xin,color='r',linewidth=0.4,label="fit")
    plt.xlabel("log(r)")
    plt.ylabel("log(psi")
    plt.title("log(psi) vs log(r)")
    plt.legend()
    plt.savefig(f"data/valueVsDistanceFIT.png")
    plt.show()

def generateData6(n,noise,dt,dx=1,rho=0.5,sigma=10,k=0.01,v0=0.01):
    model = PDE(n,noise,dt,dx,rho,sigma,k,v0=v0)
    sweeps=10000
    steps=[]
    averageValues =[]
    for i in range(sweeps):
        model.velocityUpdate()

    colour = 'magma'
    im=plt.imshow(model.lattice,cmap=colour)
    plt.colorbar(im)
    plt.title(f"Contour for v0={v0}")
    plt.xlim([0,50])
    plt.ylim([0,50])
    plt.savefig(f"figures/contour_{v0}.png")
    plt.show()

def main():
    n = 50
    noise=0.1
    rho = 0.5
    dt=0.1
 #  animate(n,noise,dt)
 #   generateData3(n,noise,dt)
    generateData4(n,noise,dt)
   # generateData6(n,noise,dt,v0=0.5)
main()
