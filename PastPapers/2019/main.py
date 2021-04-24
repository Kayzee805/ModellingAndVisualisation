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
        temporaryArray = np.zeros((self.n,self.n))
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

        grad = np.roll(self.U,1,axis=0)+np.roll(self.U,-1,axis=0)+np.roll(self.U,1,axis=1)+np.roll(self.U,-1,axis=1)-4*self.U
        
        newU= (self.dt)*((self.D1*grad)-(self.U*(np.square(self.V)))+(self.F*(1-self.U)))

        gradV =np.roll(self.V,1,axis=0)+np.roll(self.V,-1,axis=0)+np.roll(self.V,1,axis=1)+np.roll(self.V,-1,axis=1)-(4*self.V)

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
    i = 0
    while convergence>error:
        model.update()
        if i%100==0:
            plt.clf()
            im=plt.imshow(model.V,animated=True,cmap=colour)
            plt.colorbar(im)
            plt.draw()
            plt.pause(0.0001)
            print(f"coutner={i} and convergence = {convergence}")
        
        i+=1
        after = np.sum(model.U)
        convergence=abs(after-before)
        before=after

    plt.show()

    
def taskB(n,d1,d2,r,f,dt):
    model = Model(n,d1,d2,r,f,dt)
    colour='magma'
    sweeps=100000
    steps = np.linspace(1,sweeps,sweeps)
    
    convergence = 1
    error = 0.001

    counter=0
    while convergence>=error:
        before = np.sum(model.U)
        model.update()
        after=np.sum(model.U)
        convergence = abs(after-before)
        counter+=1

        if(counter%100==0):
            print(f"Counter={counter}  and convergence = {convergence}")
    
    print(counter,convergence)

    # for i in tqdm(steps):
    #     model.update()

    
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
    counter=0
    N=n*n
    dt=0.0001
    for i in range(len(allF)):
        model=Model(n,d1,d2,r,allF[counter],dt)
        convergence = 1
        tolerance = 0.001
        before = np.sum(model.U)
        counter2=0
        arrayConvergence=[]
        for j in range(50000):
            model.update()
            after = np.sum(model.U)
            convergence=abs(after-before)
            before=after
            counter2+=1
            if(counter2%10000==0):
                print(f"Counter={counter2} convergence={convergence}")

            arrayConvergence.append(calculateVariance(model.U))
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
    dt=0.0001
    
    model = Model(n,d1,d2,R=20,F=F,dt=dt,dx=dx,k=k)
    sweeps=100000
    steps = []
    variance = []
    error = 0.001
    convergence = 1
    before = np.sum(model.U)
    counter=0
    while convergence>error:
        model.update()
        after=np.sum(model.U)
        convergence=abs(after-before)
        before=after
        variance.append(calculateVariance(model.U))
        steps.append(counter)
        if(counter%1000==0):
            print(f"Counter={counter} convergence={convergence}")
        
        counter+=1

    # for i in tqdm(steps):
    #     model.update()
    #     variance.append(calculateVariance(model.U))
    
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
    dt=0.05
    R=20
    allF = np.arange(0.005,0.03,0.005)

    allVariance = []    

    
    for i in range(len(allF)):
        model = Model(n=n,D1=d1,D2=d2,R=R,F=allF[i],dt=dt,dx=dx,k=k)
        error = 1
        convergence = 10
        before= np.sum(model.U)
        variance = []
        print(allF[i])
        counter=0
        for j in range(50000):
            model.update()
            variance.append(calculateVariance(model.U))
            after = np.sum(model.U)
            convergence= abs(after-before)
            after=before
            if(counter%10000==0):
                print(f"allF={allF[i]}  Counter={counter} convergence={convergence}")
            counter+=1
        allVariance.append(calculateVariance(model.U))
    

    plt.scatter(allF,allVariance,s=5,color='r')
    plt.plot(allF,allVariance)
    plt.xlabel("F")
    plt.ylabel("Variance")
    plt.title("Variance vs F plot for task e")
    plt.savefig("figures/taskE.png")
    plt.show()


    
   

def main():
    n=100
    R=20
    D1=0.2
    D2=0.1
    F=0.005
    dt=0.03
  #  animate(n,D1,D2,R,F,dt)
   # taskB(n,D1,D2,R,F,dt)
    taskC(n,D1,D2,R,dt)
   # taskD(0.03,variancePlot=True)
   # taskE()
main()




