import numpy as np
from gameOfLife import gameOfLife
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time,sys
def animate(size,sweeps,initialisation):
    '''
    Plots the game of life animation for the given input.
    Parameters:
    -----------
    size: Type Integer
          Size of the system
    sweeps: Type Integer
            Number of sweeps for the animation. 
    initialisation: Type string
                    Type of initialisation, if given wrong initialisation, random is used.
    '''

    model = gameOfLife(size,initialisation)

    #calculating the number of live cells 
    model.currentActiveSite=model.activeSites()
    im=plt.imshow(model.lattice,animated=True,vmin=-1,vmax=1,cmap='magma')
    plt.draw()
    plt.pause(0.1)


    for i in range(sweeps):
        #updating the model for each sweep
        model.update()
        #tracking number of live cells for testing purposes
        newActiveSite = model.activeSites()
        if(newActiveSite==model.currentActiveSite):
            model.activeSitesCounter+=1
        else:
            model.activeSitesCounter=1
            model.currentActiveSite=newActiveSite
        plt.cla()
        im=plt.imshow(model.lattice,animated=True,vmin=-1,vmax=1,cmap='magma')
        plt.xlabel("X axis")
        plt.ylabel("Y axis")
        plt.title("Animation for game of life model")
        plt.draw()
        plt.pause(0.1)
      #  print(f"iteration {i} active state= {model.activeSites()}")
    plt.show()

def generateHistogram(size,sweeps,initialisation):
    '''
    Generates and saves a dat file for the time taken to reach the absorbing state
    for 1000 simulations. Time is measured in number of sweeps here.
    Parameters:
    -----------
    size: Type Integer
          Size of the system
    sweeps: Type Integer
            Number of sweeps for the animation. 
    initialisation: Type string
                    Type of initialisation, if given wrong initialisation, random is used.
    
    '''
    absorbingState = []

    for i in range(100):
        model = gameOfLife(size,initialisation)
        model.currentActiveSite=model.activeSites()
        print(f"Iteration = {i}")
        for j in range(sweeps):
            model.update()
            newActiveSite = model.activeSites()
            if(newActiveSite==model.currentActiveSite):
                model.activeSitesCounter+=1
            else:
                model.activeSitesCounter=1
                model.currentActiveSite= newActiveSite
            if(model.activeSitesCounter==10):
                #break and add to absorbing state
                absorbingState.append(j-9)
                break
   #print(f"Absorbing state = {absorbingState}")
    np.savetxt("data/histogram1000.dat",np.transpose(absorbingState),fmt='%.4f')
      
def generateCom(size,sweeps,initialisation):
    '''
    Calculates tthe centre of mass for the x and y component for upto 300sweeps.
    centre of mass is only calculated for glider that does not crossses the boundary
    of the lattice.
    The centre of mass is saved in a dat file.
    Parameters:
    -----------
    size: Type Integer
          Size of the system
    sweeps: Type Integer
            Number of sweeps for the animation. 
    initialisation: Type string
                    Type of initialisation, if given wrong initialisation, random is used.
    
    '''
    xCom =[]
    yCom=[]
    t=[]
    model = gameOfLife(size,initialisation)
    print(f"Starting to generate centre of mass.")
    for i in range(500):
        model.update()
       # print(len(np.where(model.lattice==1)[0]))
        com = model.centreOfMass()
        #com[0]= x, com[1] = y
        if(com[0]!=-1 and com[1]!=-1):
            #add com
            xCom.append(com[0])
            yCom.append(com[1])
            t.append(i)
        #model.update()

    arrayCombined = np.array((t,xCom,yCom))
   # print(arrayCombined)
    arrayCombined= np.transpose(arrayCombined)
    np.savetxt("data/centreOfMass.dat",(arrayCombined),fmt='%.4f')

    velocities = getVelocity(arrayCombined)
    print(f"X velocity ={velocities[0]}\nY velocity = {velocities[1]}\nTotal velocity = {velocities[2]}")
    np.savetxt("data/velocity.dat",np.transpose(velocities),fmt='%.4f')
    print("Done")


def getVelocity(allArray):
    '''
    Using numpy poly fit to calculate the velocity of the glider 
    Parameters:
    -----------
    allArray: Type ndArray
               MultiD array that contains the centre of mass for x and y components and the
               corresponding time.
    '''
    t=allArray[:,0]
    x=allArray[:,1]
    y=allArray[:,2]
    newX=[]
    newY=[]
    newT = []
    for i in range(len(t)):
        #adding only the elements between 150 and 250, to calculate the velocity
        if(i>150 and i<250):
            newX.append(x[i])
            newY.append(y[i])
            newT.append(i)
    
    newX=np.asarray(newX)
    newY=np.asarray(newY)
    newT=np.asarray(newT)

    xfit,xin = np.polyfit(newT,newX,1)
    yfit,yin = np.polyfit(newT,newY,1)

    #print(f"X velocity = {xfit}\nY velocity = {yfit}")

    vel = np.sqrt(xfit**2+yfit**2)
   # print(f"Total Velocity = {vel}")
    return [xfit,yfit,vel]


def plotAll():

    #data 1
    data1 = np.loadtxt("data/histogram1000.dat")
    plt.hist(data1,bins=25,density=False,facecolor='none',edgecolor='k')
    plt.title("Histogram for the time taken to reach absorbing state")
    plt.xlabel("Time step till absorbing state")
    plt.ylabel("Frequency")
    plt.savefig("figures/Histogram.png")
    plt.show()

    
    centreOfMass = np.loadtxt("data/centreOfMass.dat")
    t=centreOfMass[:,0]
    xCom=centreOfMass[:,1]
    yCom=centreOfMass[:,2]

    plt.scatter(t[::10],xCom[::10],s=10,color='k',marker='x')
    plt.title("Centre of mass for the x component, plotting every 10th data")
    plt.ylabel("centre of mass for the X component")
    plt.xlabel("Time (Sweeps)")
    plt.savefig("figures/xCom.png")
    plt.show()

    plt.scatter(t[::10],yCom[::10],s=10,color='r')
    plt.title("Centre of mass for the y component, plotting every 10th data")
    plt.ylabel("centre of mass for the y component")
    plt.xlabel("Time (Sweeps)")
    plt.savefig("figures/yCom.png")
    plt.show()

    plt.scatter(t[::10],xCom[::10],s=10,color='k',marker='x',label="X component")
    plt.scatter(t[::10],yCom[::10],s=10,color='r',label="y component")
    plt.title("Centre of mass for the x and y component, plotting every 10th data")
    plt.xlabel("Time (Sweeps)")
    plt.ylabel("Centre of mass")
    plt.legend()
    plt.savefig("figures/bothCOM.png")
    plt.show()


if __name__=="__main__":
    if(len(sys.argv)!=3):
        print("usage python main.py N sweeps")
    size = int(sys.argv[1])
    sweeps = int(sys.argv[2])
    while True:
        initialisation = int(input("0 for Random initialisation\n1 for Blinker\n2 for Glider"))
        if(initialisation==0 or initialisation==1 or initialisation==2):
            break
        else:
            print("Please enter a valid input")
            
    
    toDo = int(input("\n\n1 to animate\n2 to generate histogram data\n3 generate centre of mass data and generate velocity data: Recommended for Glider initialisation"
    +"\n4 Plot graphs"))

    if(toDo==1):
        animate(size,sweeps,initialisation)
    elif(toDo==2):
        generateHistogram(size,sweeps,initialisation)
    elif(toDo==3):
        generateCom(size,sweeps,initialisation)
        #centreOfMassData = np.loadtxt("data/centreOfMass.dat")
       # getVelocity(centreOfMassData)
    elif(toDo==4):
        plotAll()
    else:
        print("Exiting as invalid input")
