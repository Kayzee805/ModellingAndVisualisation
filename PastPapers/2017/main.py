import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from IPython.core.display import clear_output
from tqdm import tqdm
from scipy.stats import sem
import random as r

'''
Each cell is either green,red or blue. 

cell i schosen randomly. randomly select one of its nn. call this cell j

cell i converts cell j with a probability of p(Si,Sj)

P(R,G)=P(G,B)=P(B,R) = p1
P(G,R) = P(B,G)= P(R,B)=p2

if Si==Sj do nothing

R=-1
B=0
G=1
'''
class Model(object):

    def __init__(self,n,p1,p2):
        self.n = n
        self.p1 = p1
        self.p2 = p2
        self.initialiseEqual()

        # self.redCounter=0
        # self.blueCounter=0
        # self.greenCounter=0

    def initialiseEqual(self):
        self.lattice = np.ones((self.n,self.n))
        self.redCounter=0
        self.blueCounter=0
        self.greenCounter=0
        for i in range(self.n):
            for j in range(self.n):
                randNumber = r.random()
                if(randNumber<(1/3)):
                    #its red
                    self.lattice[i,j]=-1
                    self.redCounter+=1
                elif(randNumber<(2/3)):
                    #its blue
                    self.lattice[i,j]=0
                    self.blueCounter+=1
                else:
                    #its already Green
                    self.greenCounter+=1
                    continue

    def update(self):

        for i in range(self.n):
            for j in range(self.n):
                iOne = r.randint(0,self.n-1)                   
                jOne = r.randint(0,self.n-1)
                
                self.updateNeighbour(iOne,jOne)

    def updateNeighbour(self,i,j):
        #left,right,top,bot = 1,2,3,4
        randInt = r.randint(1,4)
        neighbourindex = []
        if(randInt==1):
            neighbour = self.lattice[i,(j-1)%self.n]
            neighbourindex.append(i)
            neighbourindex.append((j-1)%self.n)
        elif(randInt==2):
            neighbour = self.lattice[i,(j+1)%self.n]
            neighbourindex.append(i)
            neighbourindex.append((j+1)%self.n)
        elif(randInt==3):
            neighbour = self.lattice[(i-1)%self.n,j]
            neighbourindex.append((i-1)%self.n)
            neighbourindex.append(j)
        else:
            neighbour = self.lattice[(i+1)%self.n,j]
            neighbourindex.append((i+1)%self.n)
            neighbourindex.append(j)
        
        if(self.lattice[i,j]==neighbour):
            ##both the same so ignore
            return
        
        if(self.lattice[i,j]==-1):
            #red
            #check for blue and green
            if(neighbour==0 and r.random()<self.p2):
                self.lattice[neighbourindex[0],neighbourindex[1]]=-1
                self.redCounter+=1
                self.blueCounter-=1
            else:
                if(r.random()<self.p1):
                    self.lattice[neighbourindex[0],neighbourindex[1]]=-1
                    self.redCounter+=1
                    self.greenCounter-=1


        elif(self.lattice[i,j]==0):
            #blue
            if(neighbour==-1 and r.random()<self.p1):
                self.lattice[neighbourindex[0],neighbourindex[1]]=0
                self.blueCounter+=1
                self.redCounter-=1

            else:
                #is green
                if(r.random()<self.p2):
                    self.lattice[neighbourindex[0],neighbourindex[1]]=0
                    self.blueCounter+=1
                    self.greenCounter-=1

        else:
            if(neighbour==-1 and r.random()<self.p2):
                self.lattice[neighbourindex[0],neighbourindex[1]]=1
                self.redCounter-=1
                self.greenCounter+=1
            else:
                #is blue
                if(r.random()<self.p1):
                    self.lattice[neighbourindex[0],neighbourindex[1]]=1
                    self.greenCounter+=1
                    self.blueCounter-=1



def animate():
    n=50
    p1=1
    p2=1
    model = Model(n,p1,p2)

    plt.figure()
    colour='magma'
    im=plt.imshow(model.lattice,animated=True,cmap=colour)
    plt.colorbar(im)
    plt.draw()
    plt.pause(0.0001)

    for i in range(50000):
        model.update()
        if i%10==0:
            plt.clf()
            im=plt.imshow(model.lattice,animated=True,cmap=colour)
            plt.colorbar(im)
            plt.draw()
            plt.pause(0.0001)
   
    plt.show()    


def taskB():
    n=50
    p1=1
    p2=1
    model=Model(n,p1,p2)
    fr = []
    fb=[]
    fg=[]
    sweeps=5000
    steps= np.linspace(0,sweeps-1,sweeps)
    N=n*n
    for i in tqdm(steps):
        model.update()
        fr.append(model.redCounter/(N))
        fb.append(model.blueCounter/N)
        fg.append(model.greenCounter/N)
    

    plt.plot(steps,fr,color='r',linestyle=':',label='Fr')
    plt.plot(steps,fb,color='b',linestyle='dashdot',label='Fb')
    plt.plot(steps,fg,color='g',linestyle='--',label='Fg')
    plt.xlabel("Steps")
    plt.ylabel("Fraction")
    plt.title("Fr,Fb,Fg against time")
    plt.legend()
    plt.savefig("figures/taskB.png")
    plt.show()
    np.savetxt(f"data/taskB.dat",np.transpose(np.array(steps,fr,fb,fg)),fmt='%.5f')

def taskC():
    '''
    set p1=1 and vary p2
    '''
    tries = 5
    allP = np.arange(0.5,1.05,0.05)
    allArray = np.zeros((len(allP),tries))
    p1=1
    n=50 

    steps = np.linspace(0,tries-1,tries)
    sweeps=5000
    stepCounter=0
    for i in tqdm(steps):
        
        pCounter=0
        for j in range(len(allP)):
            model= Model(n,p1,allP[j])

            counter=0
            blue = model.blueCounter
            red = model.redCounter
            green = model.greenCounter
            absorbing = 0
            for z in range(sweeps):
                model.update()
                if(model.redCounter==red and model.blueCounter==blue and model.greenCounter==green):
                    counter+=1
                    if(counter==10):
                        absorbing=1
                        break
                else:
                    counter=0
                red = model.redCounter
                blue=model.blueCounter
                green = model.greenCounter

            
            allArray[j,stepCounter]=absorbing
        stepCounter+=1
    
    probabilities = []
    for i in range(len(allP)):
        probabilities.append(np.average(allArray[i]))
    plt.plot(allP,probabilities)
    plt.xlabel("P2")
    plt.ylabel("Probability of absorbing state")
    plt.title("Absoribing state probability for varying P2")
    plt.savefig("figures/taskC.png")
    plt.show()

    np.savetxt(f"data/taskC.dat",np.transpose(np.array((allP,probabilities))),fmt='%.5f')

            
            


def main():
   # animate()
   #taskB()
    taskC()

t1=time.time()
main()
t2=time.time()
print(f"Time taken = {t2-t1}s")
