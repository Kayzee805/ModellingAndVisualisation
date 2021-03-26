import numpy as np
import time
np.random.seed(10)

n=50
x=np.zeros((n,n,n))
noise = np.random.uniform(-1,1,(n,n,n))
mid = int(n/2)
rho=np.zeros((n,n,n))
rho[mid,mid,mid]=1
# x=np.add(x,noise)

# print(x)

epsilon = 0.001
def setBoundaries3D(Array):
    #set boundaries to 0?
    Array[0,:,:]= 0
    Array[:,0,:]= 0
    Array[:,:,0]= 0
    Array[-1,:,:]= 0
    Array[:,-1,:]= 0
    Array[:,:,-1]= 0
def calculateNearestNeighbours(array,i,j,k):
    a= array[(i+1)%n,j,k]+array[(i-1)%n,j,k]+array[i,(j-1)%n,k]+array[i,(j+1)%n,k]+array[i,j,(k-1)%n]+array[i,j,(k+1)%n]
    return a

def calculatePrevious(array,i,j,k):
    return array[(i-1)%n,j,k]+array[i,(j-1)%n,k]+array[i,j,(k-1)%n]

def normalGaussSeidel(x):
    convergence = epsilon+1
    counter=0
    convergenceArray=[]
    print("Original")
    while(convergence>=epsilon):
        sumBefore = np.sum(x)
        t1=time.time()
        for i in range(1,n-1):
            for j in range(1,n-1):
                for k in range(1,n-1):
                    nn=calculateNearestNeighbours(x,i,j,k)
                    dummy = rho[i,j,k]
                    x[i,j,k] = (nn+dummy)/6
        sumAfter = np.sum(x)
        convergence = np.abs(sumAfter-sumBefore)
        counter+=1
        convergenceArray.append(convergence)
      #  print(f"Counter={counter} convergence={convergence}  timetaken={time.time()-t1}")
        print(f"finished original at counter={counter} and convergence= {convergence}")
    return convergenceArray

def myUpdate(x):
    convergence=epsilon+1
    counter=0
    convergenceArray=[]
    print("My version")
    while(convergence>=epsilon):
        t1=time.time()
        sumBefore = np.sum(x)
        independent = np.roll(x,-1,axis=0)+np.roll(x,-1,axis=1)+np.roll(x,-1,axis=2)+rho
        for i in range(1,n-1):
            for j in range(1,n-1):
                for k in range(1,n-1):
                    x[i,j,k] = (independent[i,j,k]+calculatePrevious(x,i,j,k))/6
        
        sumAfter = np.sum(x)
        convergence=np.abs(sumAfter-sumBefore)
        counter+=1
        convergenceArray.append(convergence)
       # print(f"Counter={counter} convergence={convergence}  timetaken={time.time()-t1}")
        print(f"Counter={counter}  convergence={convergence} sumBefore={sumBefore} sumAfter={sumAfter}")
    return convergenceArray
    

def checkerBoard(x):
    print("Starting checkerboard")
    mask = np.indices((n,n,n)).sum(axis=0)%2
    
    convergence=epsilon+1
    counter=0

    while(convergence>=epsilon):
        sumBefore = np.sum(x)
        dummyWhite = np.copy(x)
        dummyWhite[mask==1]=0
        dummyWhite = (np.roll(dummyWhite,1,axis=0)+np.roll(dummyWhite,1,axis=1)+np.roll(dummyWhite,1,axis=2)+np.roll(dummyWhite,-1,axis=0)+np.roll(dummyWhite,-1,axis=1)+np.roll(dummyWhite,-1,axis=2)+rho)/6
        setBoundaries3D(dummyWhite)
        dummyBlack = np.copy(dummyWhite)
        dummyBlack = (np.roll(dummyBlack,1,axis=0)+np.roll(dummyBlack,1,axis=1)+np.roll(dummyBlack,1,axis=2)+np.roll(dummyBlack,-1,axis=0)+np.roll(dummyBlack,-1,axis=1)+np.roll(dummyBlack,-1,axis=2)+rho)/6
        setBoundaries3D(dummyBlack)

        x = np.add(dummyWhite,dummyBlack)
        
        #setBoundaries3D(x)
        sumAfter = np.sum(x)
        convergence =abs(sumAfter-sumBefore)
        counter+=1
        print(f"Counter={counter}  convergence={convergence} sumBefore={sumBefore} sumAfter={sumAfter} ")

    # while(convergence>=epsilon):
    #     sumBefore= np.sum(x)

    #     for i in range(1,n-1,1):
    #         for j in range(1,n-1,1):
    #             for k in range(1,n-1,1):
    #                 if(mask[i,j,k]==1):

    #                     nn=calculateNearestNeighbours(x,i,j,k)
    #                     x[i,j,k] = (nn+rho[i,j,k])/6
          
    #     for i in range(1,n-1,1):
    #         for j in range(1,n-1,1):
    #             for k in range(1,n-1,1):
    #                 if(mask[i,j,k]==0):
    #                     nn=calculateNearestNeighbours(x,i,j,k)
    #                     x[i,j,k] = (nn+rho[i,j,k])/6
        
    #     sumAfter =np.sum(x)
    #     convergence=abs(sumAfter-sumBefore)
    #     print(f"Counter={counter}  convergence={convergence} sumBefore={sumBefore} sumAfter={sumAfter}")
    #     counter+=1


t1=time.time()
my = myUpdate(x)
#normal=normalGaussSeidel(x)
#checkerBoard(x)
t2=time.time()
print(f"Time takne for everything={t2-t1}s")


# x=[[1,3,5,6],
# [5,1,2,6],
# [9,5,15,7],
# [1,6,8,6]]
# x= np.array(x)
# n=len(x)
# mask = np.indices((n,n)).sum(axis=0)%2
# test = np.copy(x)
# print("Initial")
# print(test)
# test[mask==1]=0
# print("masked")
# print(test)
# test= (np.roll(test,1,axis=0)+np.roll(test,1,axis=1)+np.roll(test,-1,axis=1)+np.roll(test,-1,axis=0))/4
# print("rolling white")
# print(test)
# dummy = np.copy(test)
# dummy =  (np.roll(dummy,1,axis=0)+np.roll(dummy,1,axis=1)+np.roll(dummy,-1,axis=1)+np.roll(dummy,-1,axis=0))/4
# print("rolling black")
# print(dummy)

# dummy= np.add(dummy,test)
# print("final")
# print(dummy)
'''
I have array1 which is the original array              1 ARRAY
then i basically roll it? without changing the black nodes  
then I roll the black nodes with the new updated white nodes.


have one original array
then another array where all the white nodes are 0
then solve for white by rolling? 

then set all the blacks to 0? in the new updated white
then solve for blacks
'''