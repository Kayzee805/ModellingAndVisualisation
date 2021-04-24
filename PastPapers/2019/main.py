import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from IPython.core.display import clear_output
from tqdm import tqdm
from scipy.stats import sem

class Model(object):

    def __init__(self,n,D1,D2,R,F,dt,dx=1,k=0.06,noise=0.01):
        self.n=n
        self.D1=D1
        self.D2 = D2
        self.R = R
        self.dt = dt
        self.dx = dx
        self.k=k
        self.noise = noise
        self.mid = int(self.n/2)
        self.U = self.initialiseLattice(1,0.5)
        self.V = self.initialiseLattice (0.01,0.25)
        self.F = F
 
    def initialiseLattice(self,greater,less):
        temporaryArray = np.zeros((self.n,self.n),dtype=np.float64)
        noiseArray=np.random.uniform(-self.noise,self.noise,(self.n,self.n))
        for i in range(self.n):
            for j in range(self.n):
                distance = self.distanceCalculator(i,j)
                if(distance>self.R):
                    temporaryArray[i,j]=greater
                else:
                    temporaryArray[i,j]=less
        temporaryArray+= noiseArray
        return np.copy(temporaryArray)

    def distanceCalculator(self,i,j):
        return np.sqrt((i-self.mid)**2+(j-self.mid)**2)
    
    
    def update(self):
        #left,right,up,down - 4 then multiply att by dt

        grad = (np.roll(self.U,1,axis=0)+np.roll(self.U,-1,axis=0)+np.roll(self.U,1,axis=1)+np.roll(self.U,-1,axis=1)-4*self.U)/(self.dx**2)
        
        newU= (self.dt)*((self.D1*grad)-(self.U*(np.square(self.V)))+(self.F*(1-self.U)))

        gradV =(np.roll(self.V,1,axis=0)+np.roll(self.V,-1,axis=0)+np.roll(self.V,1,axis=1)+np.roll(self.V,-1,axis=1)-(4*self.V))/(self.dx**2)

        newV = (self.dt)*((self.D2*gradV)+(self.U*(np.square(self.V)))-((self.F+self.k)*self.V))
        self.U+=newU
        self.V+=newV





def animate(n,D1,D2,R,F,dt):
    model= Model(n,D1,D2,R,F,dt)
    colour='magma'
    plt.figure()
    im = plt.imshow(model.V,animated=True,cmap=colour)
    plt.colorbar(im)
    plt.draw()
    plt.pause(0.0001)

    print("starting animation")
    convergence = 1
    error =0.0001
    before = np.sum(model.U)
    steps = np.linspace(1,50000,50000)
    for i in tqdm(steps):
        model.update()
        if i%100==0:
            plt.clf()
            im=plt.imshow(model.U,animated=True,cmap=colour)
            plt.colorbar(im)
            plt.draw()
            plt.pause(0.0001)
   
    plt.show()

    
def taskB(n,d1,d2,r,f,dt):
    model = Model(n=n,D1=d1,D2=d2,R=r,F=f,dt=dt)
    colour='magma'
    sweeps=50000
    steps = np.linspace(1,sweeps,sweeps)
    
    convergence = 1
    error = 0.001

    counter=0
    for i in tqdm(steps):
        model.update()
    plt.figure()
    im = plt.imshow(model.U,cmap=colour)
    plt.colorbar(im)
    plt.title(f"Snap shot of U at time ={sweeps} for F={f}")
    plt.savefig(f"figures/taskB/taskB_U_{f}.png")
    plt.show()

    plt.figure()
    im = plt.imshow(model.V,cmap=colour)
    plt.colorbar(im)
    plt.title(f"Snap shot of V at time ={sweeps} for F={f}")
    plt.savefig(f"figures/taskB/taskB_V_{f}.png")
    plt.show()

    array = np.array((model.U))
    array2 = np.array((model.V))
    np.savetxt(f"data/taskB/U_{f}.dat",array)
    np.savetxt(f"data/taskB/V_{f}.dat",array2)



def calculateVariance(array):
    averageSquare = np.average(np.square(array))
    squareAverage = np.square(np.average(array))
    return averageSquare-squareAverage

def taskC(n,d1,d2,r,dt):
    n=50
    allF = np.arange(0.02,0.055,0.005)
    allVariance = []
    sweeps=100000
    N=n*n
    dt=0.01
    colour='magma'
    counter=0
    steps = np.linspace(1,sweeps,sweeps)
    for i in tqdm(allF):
        model=Model(n,d1,d2,r,allF[counter],dt)

        
        arrayConvergence=[]
        for j in tqdm(steps):
            model.update()
            arrayConvergence.append(calculateVariance(model.U))
     
        plt.figure()
        im = plt.imshow(model.U,cmap=colour)
        plt.colorbar(im)
        plt.title(f"Snap shot of U at time ={sweeps} for F={round(allF[counter],3)}")
        plt.savefig(f"figures/taskC/taskC_{round(allF[counter],3)}.png")
        plt.show(block=False)
        plt.pause(3)
        plt.close()
        # for j in range(sweeps):
        #     model.update()
        
        allVariance.append(np.mean(arrayConvergence))
        counter+=1
    plt.plot(allF,allVariance)
    plt.title("Variance over F")
    plt.xlabel("F")
    plt.ylabel("Variance")
    plt.savefig(f"figures/taskC.png")
    plt.show()

    np.savetxt(f"data/taskC.dat",np.transpose(np.array((allF,allVariance))),fmt='%.6f')
        

def taskD(F,variancePlot = False):
    n=50
    dx=2
    k=0.05
    d1=0.2
    d2=0.1
    dt=0.01
    
    model = Model(n,d1,d2,R=20,F=F,dt=dt,dx=dx,k=k)
    sweeps=100000
    steps = []
    variance = []
    steps = np.linspace(0,sweeps-1,sweeps)
    for i in tqdm(steps):
        model.update()
        variance.append(calculateVariance(model.U))

    
    if(variancePlot):
        plt.plot(steps,variance)
        plt.title(f"Variance plot for f={F} over time")
        plt.xlabel("sweeps")
        plt.ylabel("Variance")
        plt.savefig(f"figures/taskD/taskD_{F}.png")
        plt.show()
    
    colour='magma'
    plt.figure()
    im = plt.imshow(model.U,cmap=colour)
    plt.colorbar(im)
    plt.title(f"Snap shot of  U for F={F}")
    plt.savefig(f"figures/taskD/taskD_U_{F}.png")
    plt.show()

    colour='magma'
    plt.figure()
    im = plt.imshow(model.V,cmap=colour)
    plt.colorbar(im)
    plt.title(f"Snap shot of  V for F={F}")
    plt.savefig(f"figures/taskD/taskD_v_{F}.png")
    plt.show()


def taskE():
    n=50
    dx=2
    k=0.05
    d1=0.2
    d2=0.1
    dt=0.01
    R=20
    allF = np.arange(0.005,0.035,0.005)

    counter=0
    sweeps=100000
    variance = np.zeros((len(allF),5))
    steps = np.linspace(1,sweeps-1,sweeps)
    errors = []
    average=[]
    
    for z in range(5):
        counter=0
        for i in tqdm(allF):
            model=Model(n=n,D1=d1,D2=d2,R=R,F=allF[counter],dt=dt,dx=dx,k=k)

            for j in tqdm(steps):
                model.update()
            
            variance[counter,z]=calculateVariance(model.U)
            counter+=1

    for i in range(len(allF)):
        errors.append(sem(variance[i]))
        average.append(np.mean(variance[i]))

    plt.scatter(allF,average,s=20,marker='+',color='k')
    plt.plot(allF,average,color='b')
    plt.errorbar(allF,average,yerr=errors,ecolor='r')
    plt.xlabel("F")
    plt.ylabel("variance")
    plt.title("variance vs F with standard errors")
    plt.savefig("figures/taskE.png")
    plt.show()
    np.savetxt(f"data/taskE.dat",np.transpose(np.array((allF,average,errors))))

def main():
    n=100
    R=20
    D1=0.2
    D2=0.1
    F=0.049
    dt=0.1
 #   animate(n,D1,D2,R,F,dt)
  #  taskB(n,D1,D2,R,F,dt)
    #taskC(n,D1,D2,R,dt)
    #taskD(0.01,variancePlot=True)
    taskE()
main()




