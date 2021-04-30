import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import random as r
import sys
from tqdm import tqdm

'''
Description:

phi= particle concentration
m= local magnetisation

phi = phi0+ noise of -0.01 to 0.01

m = 0+- noise -0.01 to 0.01

'''

class Model(object):

    def __init__(self,n,phi0,X,dt,dx=1,a=0.1,c=0.1,k=0.1,M=0.1,D=1,noise=0.01):
        self.n=n
        self.phi0 = phi0
        self.noise = noise
        self.X=X
        self.a=a
        self.c=c
        self.k=k
        self.M=M
        self.D= D
        self.dx=dx
        self.dt = dt
        self.initialiseM()
        self.initialisePhi()

    
    def initialisePhi(self):
        noiseArray = np.random.uniform(-self.noise,self.noise,(self.n,self.n))
        self.lattice = np.full((self.n,self.n),self.phi0)
        self.lattice= np.add(self.lattice,noiseArray)
    
    def initialiseM(self):
        noiseArray = np.random.uniform(-self.noise,self.noise,(self.n,self.n))
        self.m = np.zeros((self.n,self.n))
        self.m+=noiseArray

    
    def update(self):
        
        mu = -(self.a*self.lattice)+ (self.a*(self.lattice**3))-((self.X/2)*(self.m**2)) -(self.k/(self.dx**2))*(np.roll(self.lattice,1,0)+np.roll(self.lattice,-1,0)+np.roll(self.lattice,1,1)+np.roll(self.lattice,-1,1)-4*self.lattice)

        phiMu = (np.roll(mu,1,0)+np.roll(mu,-1,0)+np.roll(mu,1,1)+np.roll(mu,-1,1)-4*mu)/(self.dx**2)
        self.lattice += (self.dt)*(self.M*phiMu)


        mRHS = (self.c-(self.X*self.lattice))*self.m+ (self.c*(self.m**3))

        mLap = (np.roll(self.m,1,0)+np.roll(self.m,-1,0)+np.roll(self.m,1,1)+np.roll(self.m,-1,1)-4*self.m)/(self.dx**2)

        self.m += self.dt*((self.D*mLap)-mRHS)

    
    def updateD(self,alpha,phiBar):
        
        mu = -(self.a*self.lattice)+ (self.a*(self.lattice**3))-((self.X/2)*(self.m**2)) -(self.k/(self.dx**2))*(np.roll(self.lattice,1,0)+np.roll(self.lattice,-1,0)+np.roll(self.lattice,1,1)+np.roll(self.lattice,-1,1)-4*self.lattice)

        phiMu = (np.roll(mu,1,0)+np.roll(mu,-1,0)+np.roll(mu,1,1)+np.roll(mu,-1,1)-4*mu)/(self.dx**2)

        self.lattice += (self.dt)*(self.M*phiMu)-self.dt*(alpha*(self.lattice-phiBar))

        mRHS = (self.c-(self.X*self.lattice))*self.m+ (self.c*(self.m**3))

        mLap = (np.roll(self.m,1,0)+np.roll(self.m,-1,0)+np.roll(self.m,1,1)+np.roll(self.m,-1,1)-4*self.m)/(self.dx**2)

        self.m += self.dt*((self.D*mLap)-mRHS)

def animate(n,sweeps,phi0,X,dt):
    print(f"Animating for phi0={phi0} and x={X}")
    model = Model(n,phi0,X,dt)
    plt.figure()
    colour='magma'
    im=plt.imshow(model.m,animated=True,cmap=colour)
    plt.colorbar(im)
    plt.draw()
    plt.pause(0.0001)
    i=0
    steps= np.linspace(1,sweeps,sweeps)
    for z in tqdm(steps):
        model.update()
        if i%100==0:
            plt.clf()
            im=plt.imshow(model.m,animated=True,cmap=colour)
            plt.colorbar(im)
            plt.title(f"taskB phi0={phi0} and X={X}")
            plt.draw()
            plt.pause(0.0001)
            print(i)
        i+=1
    plt.show()

def generateTaskB(n,sweeps,phi0,X,dt):
    print(f"generating taskB for phi0={phi0} and x={X}")
    model = Model(n,phi0,X,dt)
    colour='magma'
    steps = np.linspace(1,sweeps,sweeps)
    for i in tqdm(steps):
        model.update()
     
    
    plt.figure()
    im = plt.imshow(model.lattice,cmap=colour)
    plt.title(f"taskB phi for phi0={phi0} and X={X} for sweeps={sweeps}")
    plt.colorbar(im)
    plt.draw()
    plt.savefig(f"figures/taskB/sweeps{sweeps}/phi/taskB_phi_{phi0}_{X}_{sweeps}.png")
    plt.show() 

    np.savetxt(f"data/taskB/taskB_phi_{phi0}_{X}_{sweeps}.dat",np.array(model.lattice),fmt='%.6f')

    plt.figure()
    im = plt.imshow(model.m,cmap=colour)
    plt.colorbar(im)
    plt.title(f"taskB m for  phi0={phi0} and X={X} for sweeps={sweeps}")
    plt.draw()
    plt.savefig(f"figures/taskB/sweeps{sweeps}/m/taskB_m_{phi0}_{X}_{sweeps}.png")
    plt.show() 

    np.savetxt(f"data/taskB/taskB_m_{phi0}_{X}_{sweeps}.dat",np.array(model.m),fmt='%.6f')

def generateTaskC(n,sweeps,phi0,X,dt):
    print(f"Generating task C for phi0={phi0} and X={X}")
    model = Model(n,phi0,X,dt)
    phiAverage = []
    mAverage = []
    steps = []
    i=0
    sweepsArray=np.linspace(1,sweeps,sweeps)
    for z in tqdm(sweepsArray):
        model.update()
        if(i%100==0):
            phiAverage.append(np.mean(model.lattice))
            mAverage.append(np.mean(model.m))
            steps.append(i)
        i+=1

    #plt.scatter(steps,phiAverage,s=2,color='k')
    plt.plot(steps,phiAverage)
    plt.xlabel("time(sweeps)")
    plt.ylabel("phi Average")
    plt.title(f"Average phi over time phi={phi0} X={X}")
    plt.savefig(f"figures/taskC/2taskC_phi_{phi0}_{X}.png")
    plt.show()

   # plt.scatter(steps,phiAverage,s=2,color='k')
    plt.plot(steps,mAverage)
    plt.xlabel("time(sweeps)")
    plt.ylabel("m Average")
    plt.title(f"Average m over time, phi={phi0} X={X}")
    plt.savefig(f"figures/taskC/2taskC_m_{phi0}_{X}.png")
    plt.show()

    np.savetxt(f"data/taskC/2averagePhi_{phi0}_{X}.dat",np.transpose(np.array((steps,phiAverage))),fmt='%.6f')
    np.savetxt(f"data/taskC/2averageM_{phi0}_{X}.dat",np.transpose(np.array((steps,mAverage))),fmt='%.6f')

def generateTaskD(n,sweeps,phi0,X,dt,alpha,phiBar):
    print(f"Generating task d for alpha={alpha}")
    model = Model(n,phi0,X,dt)
    steps=np.linspace(1,sweeps,sweeps)
    for i in tqdm(steps):
        model.updateD(alpha,phiBar)

        
    colour='magma'
    plt.figure()
    im = plt.imshow(model.lattice,cmap=colour)
    plt.colorbar(im)
    plt.title(f"taskD phi for alpha={alpha}")
    plt.draw()
    plt.savefig(f"figures/taskD/phi/taskD_phi_{alpha}.png")
    plt.show() 

    np.savetxt(f"data/taskD/taskD_phi_{alpha}.dat",np.array(model.lattice),fmt='%.6f')

    plt.figure()
    im = plt.imshow(model.m,cmap=colour)
    plt.colorbar(im)
    plt.title(f"taskD m for alpha={alpha}")
    plt.draw()
    plt.savefig(f"figures/taskD/m/taskD_m_{alpha}.png")
    plt.show() 

    np.savetxt(f"data/taskD/taskD_m_{alpha}.dat",np.array(model.m),fmt='%.6f')


def calculateVar(array):
    array = np.array(array)
    meanSquare = np.mean(np.square(array))
    squareMean = np.square(np.mean(array))

    return meanSquare-squareMean
def generateTaskE(n,sweeps,phi0,X,dt,phiBar):
    allAlpha = np.arange(0.0005,0.005+0.0005,0.0005)
    averageM = []
    varianceM=[]
    z=0
    steps = np.linspace(1,sweeps,sweeps)
    for j in tqdm((allAlpha)):
        model= Model(n,phi0,X,dt)
       # print(allAlpha[z])
        for i in tqdm(steps):
            model.updateD(allAlpha[z],phiBar)
        averageM.append(np.mean(model.m))   
        varianceM.append(calculateVar(model.m))
        z+=1

    plt.plot(allAlpha,averageM)
    plt.xlabel("alpha")
    plt.ylabel("Average m")
    plt.title(f"average m vs alpha plot taskE")
    plt.savefig("figures/taskE/averageM.png")
    plt.show()

    plt.plot(allAlpha,varianceM)
    plt.xlabel("alpha")
    plt.ylabel("variance m")
    plt.title(f"variance m vs alpha plot taskE")
    plt.savefig("figures/taskE/varianceM.png")
    plt.show()
    np.savetxt(f"data/taskE/averageM.dat",np.transpose(np.array((allAlpha,averageM))))
    np.savetxt(f"data/taskE/varianceM.dat",np.transpose(np.array((allAlpha,varianceM))))


def main():
    if(len(sys.argv)!=5):
        print("Usage python main.py N phi0 chi dt")
        sys.exit()
    
    n=int(sys.argv[1])
    phi0 = float(sys.argv[2])
    X = float(sys.argv[3])
    dt = float(sys.argv[4])
    sweeps = 100000
    phiBar=phi0

    task = int(input("Which task to carry out\n0: animate m\n1: generateTaskb\n2: generate task c\n3: generate task D\n4: generate taskE\ndefault: animate"))

    if(task==1):
        generateTaskB(n,sweeps,phi0,X,dt)
    elif(task==2):
        generateTaskC(n,sweeps,phi0,X,dt)
    elif(task==3):
        alpha = float(input("enter alpha value: "))
        generateTaskD(n,sweeps,phi0,X,dt,alpha,phiBar)
    elif(task==4):
        generateTaskE(50,sweeps,0.5,0.3,0.1,0.5)
    else:
        animate(n,sweeps,phi0,X,dt)

main()