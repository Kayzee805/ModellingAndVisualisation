from sirs import sirs
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import colors
import matplotlib.patches as mpatches
from scipy.stats import sem
import time
import random
import sys
matplotlib.use('TKAgg')
def animate(size,sweeps,pS,pI,pR):
    '''
    pink = sus
    black= infected
    white = rec
    Animation for the sirs model for the given input.
    '''
    model = sirs(size,pS,pI,pR)
    cMap='magma'
    im=plt.imshow(model.lattice,cmap=cMap,animated=True,vmin=-1,vmax=1)
    plt.draw()
    plt.pause(0.0001)
    for i in range(sweeps):
        model.update()
        plt.cla()
        im=plt.imshow(model.lattice,cmap=cMap,animated=True,vmin=-1,vmax=1)
        plt.draw()
        plt.pause(0.0001)
        if(i%1000==0):
            print(i)
        # if(model.infected==0):
        #     print(f"Infected = 0 at {i}")

    plt.show()



def task3(size,sweeps):
    '''
    Runs the simulation for p2=0.5 over p1-p3 plane and saves the average infected
    and the scaled variance in a dat file.
    Parameters:
    -----------
    size:  Type Integer
           Size of the system
    sweeps: Type Integer
            Number of sweeps for the simulation.
    Returns:
    --------
    Saves a dat file of the average infected and scaled variance over p1-p3 plane.
    '''
    print("Starting Task3")
    p2=0.5
    p1s= np.linspace(0,1,21)
    p3s= np.linspace(0,1,21)
    #hard setting sweep to 1,000 for this task
    sweeps=1000
    #Saving the data at each point in allArray
    allArray = np.zeros((21*21,4))
    N=size*size
    t1=time.time()
    #keeping track of the index
    counter =0
    for i in range(21):
        print(i)
        start=time.time()

        for j in range(21):

            #create a sirs object which initialises a lattice with the given arguments 
            model = sirs(size,p1s[i],p2,p3s[j])

            #empty array that will contain the number of infected over sweeps
            infectedTotal=[]

            #iterate for sweeps
            for n in range(sweeps):
                #update the lattice for each sweep
                model.update()
                
                #wait 100 sweeps, equilibirum wait
                if(n>=100):
                    #then add the number of infected in a system to the infectedTotal array
                    infected = model.infected
                    infectedTotal.append(infected)

                    #stop sweeps if system reaches absorbing state
                    if(infected==0):
                        #has reached absorbing state so break out of sweeps loop
                        break
            if(infected==0):
                #is absorbing state, so set <I> and variance to 0
                averageInfected=0
                variance=0
            else:
                #not absorbing state so we calculate the average infected and scaled variance
                infectedTotal=np.asarray(infectedTotal)
                variance = np.var((infectedTotal))/N
                averageInfected = np.mean(infectedTotal)/N
            
            #storing the data into a huge array, which we will save as a dat file
            allArray[counter]=[p1s[i],p3s[j],averageInfected,variance]
            #incrementing counter for the array index
            counter+=1
        print(f"Finished {i} in time {time.time()-start}s")

    print(f"time taken = {time.time()-t1}")
    np.savetxt("data/Task3ProcessedData.dat",allArray,fmt='%.7f')
    #np.savetxt("data/Task3ProcessedData100x100.dat",allArray,fmt='%.7f')

    #test plot for a 100x100 system
#     array = np.loadtxt("data/Task3ProcessedData100x100.dat")
#     p1s=array[:,0]
#     p3s=array[:,1]
#     infected = array[:,2]
#     variance = array[:,3]

#     print(f"Lengh.  p1s={p1s.size}  p2s={p3s.size} infected={len(infected)}")
#     p1s = p1s.reshape((21,21))
#     p3s=p3s.reshape((21,21))
#     infected = infected.reshape((21,21))
#     variance = variance.reshape((21,21))
#     plt.figure()
#   #  CS=plt.contour(p1s,p3s,infected)

#     CS=plt.contourf(p1s,p3s,variance,cmap='magma')
#     #plt.clabel(CS,fontsize=8,colors='k')
#     cbar=plt.colorbar(CS)
#     plt.xticks(np.linspace(0,1,11))
#     plt.yticks(np.linspace(0,1,11))

#     plt.xlabel("P1")
#     plt.ylabel("P3")
#     plt.title("Scaled variance of p1-p3 plane")
#     plt.savefig("figures/Task3_ScaledVariance100x100.png")
#     plt.show()


def task4(size,sweeps):
    '''
    Runs the simulation for p2=p3=0.5 over varying p1 and saves the scaled variance in a dat file.
    Parameters:
    -----------
    size:  Type Integer
           Size of the system
    sweeps: Type Integer
            Number of sweeps for the simulation.
    Returns:
    --------
    Saves a dat file of the scaled variance over a cut of p2,p3.
    '''
    print("Starting Task4")
    t1=time.time()
    p3=0.5
    p2=0.5
    p1s=np.linspace(0.2,0.5,31)
    N=size*size
    times=5
    allArray = np.zeros((len(p1s),times))
    
    for s in range(times):
        print(f"Starting {s}")
        for i in range(len(p1s)):
            start = time.time()
            #create a sirs object which initialises a lattice with the given arguments 
            model = sirs(size,p1s[i],p2,p3)

            #empty array that will contain the number of infected over sweeps
            infectedTotal=[]
            for n in range(sweeps):
                #update system for each sweep
                model.update()

                #wait 100 sweeps, equilibirum wait
                if(n>=100):
                    #then add the number of infected in a system to the infectedTotal array
                    infected=model.infected
                    infectedTotal.append(infected)
                                        
                    #stop sweeps if system reaches absorbing state
                    if(infected==0):
                        #has reached absorbing state so break out of sweeps loop
                        break
            if(infected==0):
                #is absorbing state so set variance to 0
                variance=0
            else:
                #is not absorbing state so caluclate the scaled variance
                variance=np.var(infectedTotal)/N
            
            #storing the variance in the array for index i,s
            allArray[i,s]=variance
            print(f"Time for {s} {i} at time = {time.time()-start}")
    
    #saving the array as a dat file
    np.savetxt(f"data/Task4_RawDataCombined.dat",np.transpose(allArray),fmt='%.6f')
    
    #calculating the standard error and the average scaled varaince over 5 simulations
    averageInfected =[]
    standardErrors = []
    for i in range(len(p1s)):
       averageInfected.append(np.mean(allArray[i]))
       standardErrors.append(sem(allArray[i]))
    combined = np.array((p1s,averageInfected,standardErrors))
    np.savetxt("data/Task4_ProcessedData.dat",np.transpose(combined),fmt='%.6f')
    print("DONEEE")   


def task5(size,sweeps):
    '''
    Runs the simulation for p1=p2=p3=0.5 for varying fraction of immunity for a 50x50 system.
    -----------
    size:  Type Integer
           Size of the system
    sweeps: Type Integer
            Number of sweeps for the simulation.
    Returns:
    --------
    Saves a dat file of the average infected at a coressponding fraction of immunity.
    '''
    print("Starting Task5 part one")

    #parameter intialisation
    p1=p2=p3=0.5
    immuneProbabilities = np.linspace(0,1,101) 
    precision = len(immuneProbabilities)
    N=size*size
    simulations=5

    #creating a 101x5 array to store data
    infectionArray = np.zeros((precision,simulations))

    #run the program for multiple simulations
    for s in range(simulations):
        print(s)
        t1=time.time()
        #looping for immune probability 0->1 at an interval of 0.01
        for i in range(len(immuneProbabilities)):
            print(f"Starting immune fraction {immuneProbabilities[i]}")
            #initialising a sirs object with new immune fraction.
            model=sirs(size,p1,p2,p3,isImmune=True,immuneProbability=immuneProbabilities[i])
            infectedTotal=[]
            for n in range(sweeps):
               # t2=time.time()
               #update system for each sweep
                model.update()
                #take measurements after the equilibirum wait
                if(n>=100):
                    infected=model.infected      #calling for the number of infected in the system
                    infectedTotal.append(infected)
                    if(infected==0):
                        #is absorbing state so break from sweeps loop
                        break
            
            #sets average infection to 0 if system in absoribing state
            if(infected==0):
                infectionArray[i,s]=0
            #else take the averag scaled infection
            else:
                infectionArray[i,s]=np.mean(infectedTotal)/N 
            
        print(f"Time taken {s}== {time.time()-t1}s")

    #calculating the standard error and the average scaled varaince over 5 simulations
    finalArray = []
    errors = []
    for i in range(precision):
        finalArray.append(np.mean(infectionArray[i]))
        errors.append(sem(infectionArray[i]))
    
    combined = np.array((immuneProbabilities,finalArray,errors))

    #storing the data as a dat file
    np.savetxt('data/Task5_RawDataCombined.dat',(infectionArray),fmt='%.6f')
    np.savetxt('data/Task5_ProcessedData.dat',np.transpose(combined),fmt='%.6f')
def task5_Part2(size,sweeps):
    '''
    Runs the simulation for p1=0.8, p2=0.1,p3=0.02 for varying fraction of immunity for a 100x100 system.
    -----------
    size:  Type Integer
           Size of the system
    sweeps: Type Integer
            Number of sweeps for the simulation.
    Returns:
    --------
    Saves a dat file of the average infected at a coressponding fraction of immunity.
    '''
    print("Starting Task5 part 2")
    #initialising variables to fit the task
    p1=0.8
    p2=0.1
    p3=.02
    immuneProbabilities = np.linspace(0,1,101) 
    precision = len(immuneProbabilities)
    N=size*size
    simulations=5

    #creating a 101x5 array to store data
    infectionArray = np.zeros((precision,simulations))

    #run the program for multiple simulations
    for s in range(simulations):
        print(s)
        t1=time.time()
        #looping for immune probability 0->1 at an interval of 0.01
        for i in range(len(immuneProbabilities)):
            print(f"Starting immune fraction {immuneProbabilities[i]}")
            #initialising a sirs object with new immune fraction.
            model=sirs(size,p1,p2,p3,isImmune=True,immuneProbability=immuneProbabilities[i])
            infectedTotal=[]
            for n in range(sweeps):
                #update system for each sweep
                model.update()

                #take measurements after the equilibirum wait
                if(n>=100):
                    infected=model.infected
                    infectedTotal.append(infected)
                    if(infected==0):
                        #is absorbing state so break from sweeps loop
                        break

            #sets average infection to 0 if system in absoribing state
            if(infected==0):
                infectionArray[i,s] =0
            #else take the averag scaled infection
            else:
                infectionArray[i,s]=np.mean(infectedTotal)/N 
        print(f"Time taken {s} == {time.time()-t1}s")

    #calculating the standard error and mean scaled variance over 5 simulations
    finalArray = []
    errors = []
    #using scipy.stats.sem to calculate standard error
    for i in range(precision):
        finalArray.append(np.mean(infectionArray[i]))
        errors.append(sem(infectionArray[i]))
    
    #saving arrays as a dat file for plotting
    combined = np.array((immuneProbabilities,finalArray,errors))
    np.savetxt('data/Task7_RawDataCombined.dat',(infectionArray),fmt='%.6ft')
    np.savetxt('data/Task7_ProcessedData.dat',np.transpose(combined),fmt='%.6f')



def plotAll():
    #Plotting the heatmaps/contours for task 3
    array = np.loadtxt("data/Task3ProcessedData.dat")
    p1s=array[:,0]
    p3s=array[:,1]
    infected = array[:,2]
    variance = array[:,3]

    print(f"Lengh.  p1s={p1s.size}  p2s={p3s.size} infected={len(infected)}")
    p1s = p1s.reshape((21,21))
    p3s=p3s.reshape((21,21))
    infected = infected.reshape((21,21))
    variance = variance.reshape((21,21))

    #contour of average infected
    plt.figure()
    CS=plt.contourf(p1s,p3s,infected,cmap='magma')
    cbar=plt.colorbar(CS)
    plt.xticks(np.linspace(0,1,11))
    plt.yticks(np.linspace(0,1,11))

    plt.xlabel("P1")
    plt.ylabel("P3")
    plt.title("Average infected of p1-p3 plane")
    plt.savefig("figures/Task3_AverageInfected.png")
    plt.show()


    #contour of scaled variance
    plt.figure()
    CS=plt.contourf(p1s,p3s,variance,cmap='magma')
    cbar=plt.colorbar(CS)
    plt.xticks(np.linspace(0,1,11))
    plt.yticks(np.linspace(0,1,11))

    plt.xlabel("P1")
    plt.ylabel("P3")
    plt.title("Scaled variance of p1-p3 plane")
    plt.savefig("figures/Task3_ScaledVariance.png")
    plt.show()

    #heatmap of scaled variance
    variance = np.asarray(variance)
    #print(f"Max variance = {np.amax(variance)}")
    xTicks = np.round(np.linspace(0,1,11),2)
    plt.imshow(np.transpose(variance),cmap='magma',extent=[0,1,1,0])
    plt.xlabel("P1")
    plt.ylabel("P3")
    plt.xticks(xTicks)
    plt.yticks(xTicks)
    plt.title("Scaled variance of p1-p3 plane")
    plt.colorbar(label="Scaled variance")
    plt.gca().invert_yaxis()
    plt.savefig("figures/TASK3_ScaledVariance_imshow.png")
    plt.show()
     
     #heatmap of average infected
    xTicks = np.round(np.linspace(0,1,11),2)
    plt.imshow(np.transpose(infected),cmap='magma',extent=[0,1,1,0])
    plt.xlabel("P1")
    plt.ylabel("P3")
    plt.xticks(xTicks)
    plt.yticks(xTicks)
    plt.title("Average infected of p1-p3 plane")
    plt.colorbar(label="Ratio of average infected")
    plt.gca().invert_yaxis()

    plt.savefig("figures/TASK3_AverageInfected_imshow.png")
    plt.show()



    #TASK4
    array2 = np.loadtxt("data/Task4_ProcessedData.dat")
    xAxis = array2[:,0]
    yAxis=array2[:,1]
    yError=array2[:,2]
    plt.scatter(xAxis,yAxis,s=5,color='r')
    plt.plot(xAxis,yAxis,color='b')
    plt.errorbar(xAxis,yAxis,yerr=yError,ecolor='k')

    #Variance plot for the cut with p2=p3=0.5
    plt.xlabel("P1")
    plt.ylabel("Variance")
    plt.title("Variance vs P1 with p2=p3=0.5")
    plt.savefig("figures/Task4_ScaledVariance_cutP3.png")
    plt.show()


    #TASK 5
    #plotting the immunefraction vs infected fraction for 50x50
    task5Array = np.loadtxt("data/Task5_ProcessedData.dat")
    immuneProbability = task5Array[:,0]
    task5Infected = task5Array[:,1]
    task5Error = task5Array[:,2]

    plt.scatter(immuneProbability,task5Infected,s=10,color='k')
    plt.plot(immuneProbability,task5Infected,color='b')
    plt.errorbar(immuneProbability,task5Infected,yerr=task5Error,ecolor='r')
    plt.xlabel("Immune probability")
    plt.ylabel("Average infected")
    plt.title("Immune vs Infected for 50x50 p1=p2=p3=0.5")
    plt.savefig("figures/Task5_50x50.png")
    plt.show()

    #TASK 7 part2
    #plotting the immunefraction vs infected fraction for 100x100
    task5Part2Array = np.loadtxt("data/Task5Part2_ProcessedData.dat")
    immuneProbability7 = task5Part2Array[:,0]
    task7Infected = task5Part2Array[:,1]
    task7Error = task5Part2Array[:,2]

    plt.scatter(immuneProbability7,task7Infected,s=10,color='k')
    plt.plot(immuneProbability7,task7Infected,color='b')
    plt.errorbar(immuneProbability7,task7Infected,yerr=task7Error,ecolor='r')
    plt.xlabel("Immune probability")
    plt.ylabel("Average infected")
    plt.title("Immune vs Infected for 100x100 p1=0.8,p2=0.1,p3=0.02")
    plt.savefig("figures/Task5_100x100.png")
    plt.show()


if __name__=="__main__":
    if(len(sys.argv)!=5):
        print("usage python main.py N p1 p2 p3")
        exit(0)
    size = int(sys.argv[1])
    pS = float(sys.argv[2])
    pI = float(sys.argv[3])
    pR = float(sys.argv[4])
    #sweeps will be 10,000 for all tasks but 1 for which it has been changed in the method
    sweeps = 10000

    toDo = int(input("\n0 to animate\n1 to plot graphs\n2 to generate data for Task 3"+
    "\n3 to generate data for Task 4\n4 to generate data for Task5"))
    
    if(toDo==0):
        animate(size,sweeps,pS,pI,pR)
    elif(toDo==1):
        plotAll()
    elif(toDo==2):
        task3(size,sweeps)
    elif(toDo==3):
        task4(size,sweeps)
    elif(toDo==4):
        part = int(input("\n\n0 for part one, 50x50\n1 for part two, 100x100"))
        if(part==0):
            task5(size,sweeps)
        elif(part==1):
            task5_Part2(size,sweeps)
        else:
            print("Please try again for a valid number for task 5.")
            exit(0)
    else:
        print("Please enter a valud task to carry out")
        exit(0)
    
    print("Finished running task")
