import numpy as np 
import time
import matplotlib.pyplot as plt

# A poisson class that will create an object which initialise variables shared by all poisson objects
class poisson(object):


    def __init__(self,n,epsilon=1,dx=1,minW=1,maxW=2):
        '''
        Constructor for poisson object. Initialises arrays and variables shared by all poisson object.
        Paramters:
        ----------
        n: Type integer
           Size of the model

        epsilon: Type float
                 Precision for a system to converge. 
        dx: Type float/Integer
            the time step used for lapalce and gradient operators.
        minW: Type float
              minimum value of omega for SOR method.
        maxW: Type float
              maximum value of omega for SOR method.
        '''
        self.n=n
        self.epsilon=epsilon
        self.dx=dx

        self.rho = np.zeros((n,n,n))        #rho array
        self.lattice=np.zeros((n,n,n))      #potential array
        self.nextLattice=np.zeros((n,n,n))  #a copy of potential array, used by jacobi method
        self.minW=minW
        self.maxW=maxW
   

    def setBoundaries3D(self,Array):
        '''
        Sets the borders of a 3D array to 0.
        Parameters:
        ----------
        Array: Type 3D array
        '''
        Array[0,:,:]= 0
        Array[:,0,:]= 0
        Array[:,:,0]= 0
        Array[-1,:,:]= 0
        Array[:,-1,:]= 0
        Array[:,:,-1]= 0
    
    def setBoundaries2D(self,Array):
        '''
        Sets the borders of a 2D array to 0.
        Parameters:
        ----------
        Array: Type 2D array
        '''
        Array[0,:]=0
        Array[-1,:]=0
        Array[:,0]=0
        Array[:,-1]=0 
    
    def setPointCharge2D(self):
        '''
        Initialises the rho array to a point charge in 2D. 
        '''
        mid=int(self.n/2)
        self.rho= np.zeros((self.n,self.n))
        self.rho[mid,mid]=1
    
    def setPointCharge3D(self):
        '''
        Initialises the rho array to a point charge in 3D. 
        '''
        mid=int(self.n/2)
        self.rho = np.zeros((self.n,self.n,self.n))
        self.rho[mid,mid,mid]=1
    
    def setChargedWire3D(self):
        '''
        Initialises the rho array to have a charged wire along the z axis.
        Boundaries are still 0.
        '''
        mid=int(self.n/2)
        for z in range(1,self.n-1):
            self.rho[mid,mid,z]=1
    
    def calculateNearestNeighbours3D(self,i,j,k):
        '''
        Calculates the num of the 6 nearest neighbours of the index i,j,k in the lattice.
        Parameters:
        -----------
        i,j,k: Type integer
               x,y,z axis index in the lattice
        '''
        n=self.n
        return self.lattice[(i+1)%n,j,k]+self.lattice[(i-1)%n,j,k]+self.lattice[i,(j-1)%n,k]+self.lattice[i,(j+1)%n,k]+self.lattice[i,j,(k-1)%n]+self.lattice[i,j,(k+1)%n]

    def dependentNeighbours(self,i,j,k):
        '''
        Calculates the depenedent neighbours of a point in the lattice.
        Dependet neighbours are neighbours that has already been iterated over/already vistited.
        Parameters:
        -----------
        i,j,k: Type integer
               x,y,z axis index in the lattice
        Returns:
        --------
        Sum of the depenedent neighbours
        '''
        n=self.n
        return self.lattice[(i-1)%n,j,k]+self.lattice[i,(j-1)%n,k]+self.lattice[i,j,(k-1)%n] 


    def dependentNeighbours2D(self,i,j,Array):
        '''
        Calculates the depenedent neighbours of a point in the Array.
        Dependet neighbours are neighbours that has already been iterated over/already vistited.
        Parameters:
        -----------
        Array: Type 2D Array
        i,j: Type integer
               x,y axis index in the lattice
        
        Returns:
        --------
        Sum of the depenedent neighbours
        '''
        n=self.n
        return Array[(i-1)%n,j]+Array[i,(j-1)%n]
    
    def getPotential(self):
        '''
        Returns a slice of the vector potential.
        '''
        return self.lattice[:,:,int(self.n/2)]


    def jacobiUpdate(self):
        '''
        Update method of the jacobi. 
        Using np roll to update my lattice in one go. Carries out np roll in 6 directions to carry out the jacobi update
        in one line. 
        SumBefore is calculated before the update method and sumAfter is calculated after the update method.
        Aftering taking the abs difference, we set sumBefore to have the same value as sumAfter
        '''
        print("Starting jacobi update method")

        #to allow the first while loop iteration            
        convergence = self.epsilon+1 
        #to track the convergence step 
        counter=0 
        sumBefore = np.sum(self.lattice)
        self.setBoundaries3D(self.lattice)
        while(convergence>=self.epsilon):
            t1=time.time()
            
            #update method of jacobi
            #new values are dependent on the values at previous index. 
            #The update method has to refer to the previous lattice, which should remain unchanged
            #then after its update, copy it. But using np.roll allows us to bypass that copying part.
            self.lattice= (self.rho+(
            np.roll(self.lattice,1,axis=0)+np.roll(self.lattice,-1,axis=0)
            +np.roll(self.lattice,1,axis=1)+np.roll(self.lattice,-1,axis=1)+
            np.roll(self.lattice,1,axis=2)+np.roll(self.lattice,-1,axis=2)))/6

            #setting boundaries to 0
            self.setBoundaries3D(self.lattice)
            #calculating the sum of the lattice then finding convergence
            sumAfter=np.sum(self.lattice)
            convergence = abs(sumAfter-sumBefore)

            #updating the variable sumBefore for the next iteration if convergence condition not met
            sumBefore=sumAfter
            counter+=1

            if(counter%100==0):
                print(f"Counter={counter} convergence={convergence} time={time.time()-t1}s")

        print("finished jacobi method")

    
    def gaussSeidelUpdate_roll(self):
        '''
        Update method for the gauss-seidel which precomutes the values in the update method which 
        are not dependent or has not changed. 
        Jacobi takes in (i-1,n+1) and (i+1,n). where n is the step. So we compute all the (i+1,n) using np.roll
        then calculate the depepndent values as we update the model using for loops.
        '''
        print("Starting gauss roll")

        #to pass the first while loop check.
        convergence = self.epsilon+1
        counter=0
        sumBefore = np.sum(self.lattice)
        #checking if convergence is below the precision if so, we stop else, we repeat the steps.
        while(convergence>=self.epsilon):
            t1=time.time()
            #precompute the neighbours that are from the same time step using np.roll
            independentNeighbours = np.roll(self.lattice,-1,axis=0)+np.roll(self.lattice,-1,axis=1)+np.roll(self.lattice,-1,axis=2)

            for i in range(1,self.n-1):
                for j in range(1,self.n-1):
                    for k in range(1,self.n-1):
                        #loop over each index and update the lattice.
                        self.lattice[i,j,k] = (independentNeighbours[i,j,k]+self.dependentNeighbours(i,j,k)+self.rho[i,j,k])/6

            #compute teh sum of the lattice
            sumAfter=np.sum(self.lattice)
            #compute the convergence then set teh value of sumBefore to sumAfter
            convergence = abs(sumAfter-sumBefore)
            counter+=1
            sumBefore=sumAfter

            if(counter%100==0):
                print(f"Counter={counter} convergence={convergence} time={time.time()-t1}s")

        print(f"Counter={counter} convergence={convergence}")

        print("Finished gauss roll")

    def gaussSeidel_CheckerBoard(self):
        '''
        Gauss-seidel udpate but using checkerBoard method. This is not equivalent to gauss-seidel but is still similar. 
        It first splits the model to black and white nodes, similar to a chess board. It then first updates all the black nodes,
        using the white nodes. Then it sets all the white nodes to 0, then use this newly updated black nodes to udpate the white nodes.
        We then concatenate these two array whilst setting the boundaries to 0 after each update. 
        '''
        print("STARTING CHECKER Board")

        #mask for a 3d array. Creates a checkerboard, splitting model to white and black nodes.
        mask = np.indices((self.n,self.n,self.n)).sum(axis=0)%2
        #to pass the first convergence check
        convergence=self.epsilon+1
        counter=0
        sumBefore = np.sum(self.lattice)

        #if convergence less than the precision, stop updating else keep udpating until condition met.
        while(convergence>=self.epsilon):
            t1=time.time()

            #first start by setting all the black nodes to 0
            self.lattice[mask==0]=0

            #use the white nodes to update the black nodes
            self.lattice = (np.roll(self.lattice,1,axis=0)+np.roll(self.lattice,1,axis=1)+np.roll(self.lattice,1,axis=2)+np.roll(self.lattice,-1,axis=0)+np.roll(self.lattice,-1,axis=1)+np.roll(self.lattice,-1,axis=2)+self.rho)/6
            self.setBoundaries3D(self.lattice)

            #then set all the nodes to 0
            self.lattice[mask==1]=0

            #we then use all the updated black nodes, to calculate the white nodes.
            self.nextLattice = (np.roll(self.lattice,1,axis=0)+np.roll(self.lattice,1,axis=1)+np.roll(self.lattice,1,axis=2)+np.roll(self.lattice,-1,axis=0)+np.roll(self.lattice,-1,axis=1)+np.roll(self.lattice,-1,axis=2)+self.rho)/6
            self.setBoundaries3D(self.nextLattice)
            #we then concatenate the white node array and teh black node array
            self.lattice = np.add(self.lattice,self.nextLattice)

            #calculate the sum of the new lattice and calculate the convergence
            sumAfter = np.sum(self.lattice)
            convergence= abs(sumAfter-sumBefore)
            counter+=1
            sumBefore = sumAfter
            if(counter%100==0):
                print(f"Counter={counter} convergence={convergence} time={time.time()-t1}s")

        print(f"Counter={counter}  convergence={convergence}")

        print("Finished CHECKER Board")


    def generate_SOR_CheckerBoard(self):
        '''
        Updating the model using the SOR method and the checkerboard method.
        First seperated the model by black and white nodes. Then set all black nodes to 0.
        Then calculate the black nodes by using the white nodes whilst following a similar appraoch to the SOR method. 
        Then use these new black nodes to calculate the white nodes.
        Then concatenate the array to return a newly updated model. 
        '''
        
        #different values of omega that we will run the udpate method on. 
        w=np.arange(self.minW,self.maxW,0.01)
        stableSweeps=[]

        #creating the mask for the model. Splitting it by black and white ndoes.
        mask = np.indices((self.n,self.n)).sum(axis=0)%2

        #carry out the update method for all values of omega.
        for x in w:
            t1=time.time()
            #for each value of omega, we reset the model to an empty array
            self.lattice=np.zeros((self.n,self.n))
            #set a point charge in 2d array. This is similar to taking a cut of a charged wire in a 3D model.
            self.setPointCharge2D()

            #to allow first while loop condition check.
            convergence=self.epsilon+1
            #sum of lattice at the start is always 0
            sumBefore=0
            counter=0

            print(f"STARTING point 3d sor for w={x}")

            while(convergence>self.epsilon):
                #SOR method requires you to add the original value of the lattice itself. So we save a copy of this lattice.
                original = np.copy(self.lattice)

                #set all black nodes to 0 then calculate the new values using the white nodes.
                self.lattice[mask==0]=0
                self.lattice =(1-x)*original+x*((np.roll(self.lattice,1,axis=0)+np.roll(self.lattice,1,axis=1)+np.roll(self.lattice,-1,axis=0)+np.roll(self.lattice,-1,axis=1)+self.rho)/4)
                self.setBoundaries2D(self.lattice)

                #then set all the white nodes to 0, then update it using the updated black nodes.
                self.lattice[mask==1]=0
                self.nextLattice = np.copy(self.lattice)
                self.nextLattice = (1-x)*(original)+x*((np.roll(self.nextLattice,1,axis=0)+np.roll(self.nextLattice,1,axis=1)+np.roll(self.nextLattice,-1,axis=0)+np.roll(self.nextLattice,-1,axis=1)+self.rho)/4)
                #we set the boundaries to 0 again
                self.setBoundaries2D(self.nextLattice)

                #we set all the black nodes to 0 for this copy of array because, we add the original array. 
                #adding this step makes it so we avoid double counting of the original array.
                self.nextLattice[mask==0]=0

                #concatenate the two arrays, the black and the white node arrays.
                self.lattice=np.add(self.lattice,self.nextLattice)
                #calculate the sum of the model and then the convergence
                sumAfter=np.sum(self.lattice)
                convergence=abs(sumAfter-sumBefore)

                counter+=1
                if(counter%100==0):
                    print(f"Counter={counter} convergence={convergence} sumBefore={sumBefore}  sumAfter={sumAfter}  ")
                sumBefore=sumAfter
            
            #store the number of steps required to converge. 
            stableSweeps.append(counter)
            print(f"Finsihed at counter={counter} at w={x} at time={time.time()-t1}s")
        
        #Plotting of the SOR data, the steps taken to converge for each value of omega.
        plt.plot(w,stableSweeps)
        plt.scatter(w,stableSweeps,s=5,color='k')
        plt.xlabel("w")
        plt.ylabel("sweeps")
        plt.title(f"Steps required to stabilise for different omega values. SOR_{self.n}")
        plt.savefig(f"figures/SOR/SOR_Checkerboard_{self.n}_{self.epsilon}.png")
        plt.show()
        combined=np.array((w,stableSweeps))
        np.savetxt(f"data/SOR_CHECKERBOARD_{self.n}_{self.epsilon}.dat",np.transpose(combined),fmt='%.8f')
        print(f"Minimum at {np.min(stableSweeps)}")


    def generate_SOR(self):
        '''
        This method uses a SOR update rule, which precomputes the independent parts of the gauss-seidel. 
        Uses a for loop to update each point in the model, turn wise turn.
        '''
        #different values of omega that we will run the udpate method on. 
        w=np.arange(self.minW,self.maxW,0.01)
        #array to store the number of steps required to converge for each omega.
        stableSweeps=[]
        n=self.n

        #calculate steps required to converge for each value of omega
        for x in w:
            t1=time.time()
            #intiliase the model and the rho array to 0 and point charge for each value of omega.
            self.lattice=np.zeros((n,n))
            self.setPointCharge2D()

            #to allow first pass of the while loop
            convergence=self.epsilon+1
            counter=0
            sumBefore=0
            print(f"STARTING W={x}")

            #carry out update until convergence is below epsilon/precision.
            while(convergence>self.epsilon):
                convergence=0
                #precompute all the independent variables in gauss-seidel
                independentNeighbours = np.roll(self.lattice,-1,axis=0)+np.roll(self.lattice,-1,axis=1)
                #then update the model by looping over each value
                for i in range(1,self.n-1):
                    for j in range(1,self.n-1):
                        old = self.lattice[i,j]
                        #SOR update method for one point. divinding by 4 because, 4 nearest neighbour for a 2D model.
                        self.lattice[i,j] = ((1-x)*self.lattice[i,j])+x*(self.dependentNeighbours2D(i,j,self.lattice)+independentNeighbours[i,j]+self.rho[i,j])/4
                        #compute the convergence step by step.
                        convergence+=abs(self.lattice[i,j]-old)
               
                if(counter%100==0):
                    print(f"counter={counter} convergence={convergence}")
                counter+=1

            
            print(f"Finsihed at counter={counter} at w={x} at time={time.time()-t1}s")
            #store the number of steps required to converge. 
            stableSweeps.append(counter)

        #Plotting of the SOR data, the steps taken to converge for each value of omega.
        plt.plot(w,stableSweeps)
        plt.xlabel("w")
        plt.ylabel("steps")
        plt.title(f"Steps required to stabilise for different omega values. SOR_{self.n}")
        plt.savefig(f"figures/SOR/SOR_{self.n}_{self.epsilon}.png")
        plt.show()
        np.savetxt(f"data/SOR/SOR_{self.n}_{self.epsilon}",np.transpose(np.array((w,stableSweeps))),fmt='%.8f')
        print(f"Minimum of {np.amin(stableSweeps)}")

  

    def distantCaculator(self,i,j,k):
        '''
        Return the distance from the centre of the model, for given indicies i,j and k.
        '''
        mid = int(self.n/2)
        return np.sqrt((mid-i)**2+(mid-j)**2+(mid-k)**2)

    def getPotential(self):
        '''
        Slices a 3D array at the half way point along the z-axis and return a 2d array.
        '''
        return self.lattice[:,:,int(self.n/2)]

    def getEx(self,i,j,k):
        '''
        Returns the electric field in the x direction of the index i,j,k
        '''
        gradient = self.lattice[(i+1)%self.n,j,k]- self.lattice[(i-1)%self.n,j,k]
        return -gradient/2

    def getEy(self,i,j,k):
        '''
        Returns the electric field in the y direction of the index i,j,k
        '''
        gradient = self.lattice[i,(j+1)%self.n,k]- self.lattice[i,(j-1)%self.n,k]
        return -gradient/2

    def getEz(self,i,j,k):
        '''
        Returns the electric field in the z direction of the index i,j,k
        '''
        gradient = self.lattice[i,j,(k+1)%self.n]- self.lattice[i,j,(k-1)%self.n]
        return -gradient/2

 
    def generateElectricData(self):
        '''
        Generates the electric field and potential data for a monopole then stores it in a data file.
        It also, calls the function to plot the vector electric field.
        '''
        x =[]
        y=[]
        z=[]
        potentialArray=[]
        ex=[]
        ey=[]
        ez=[]
        potential = self.getPotential()
        mid=int(self.n/2)
        distance=[]
        electricField = []
        for i in range(1,self.n-1):
            for j in range(1,self.n-1):
                if(i==mid and j==mid):
                    continue
                potentialArray.append(potential[i,j])
                x.append(i)
                y.append(j)
                z.append(mid)
                currentEx = self.getEx(i,j,mid)
                currentEy=self.getEy(i,j,mid)
                currentEz=self.getEz(i,j,mid)
                ex.append(currentEx)
                ey.append(currentEy)
                ez.append(currentEz)
                distance.append(self.distantCaculator(i,j,mid))
                electricField.append(np.sqrt(currentEx**2+currentEy**2+currentEz**2))
        

        withPosition = np.array((x,y,z,potentialArray,ex,ey,ez))
        withDistance = np.array((distance,potentialArray,electricField),dtype=float)
        print(withPosition.shape)
        print(withDistance.shape)
        np.savetxt(f"data/ElectricField/potentialData_{self.n}.dat",np.transpose(withPosition))
        np.savetxt(f"data/ElectricField/potentialDataVR_{self.n}.dat",np.transpose(withDistance))

        self.plotElectricField()


    
    def plotElectricField(self):
        '''
        Plots the vector electric field using plt.quiver
        The eleectric field is calculated once again here to avoid complications with meshgrid
        and this time we store the data in a 2D array of shape nxn.
        '''

        xGradient = np.zeros((self.n,self.n))
        yGradient = np.zeros((self.n,self.n))
        zGradient = np.zeros((self.n,self.n))
        mid =int(self.n/2)
        for i in range(1,self.n-1):
            for j in range(1,self.n-1):
                xGradient[i,j] = self.getEx(i,j,mid)
                yGradient[i,j] = self.getEy(i,j,mid)
                zGradient[i,j] = self.getEz(i,j,mid)

        #normalising the electric field
        for i in range(1,self.n-1):
            for j in range(1,self.n-1):
                norm = np.sqrt((xGradient[i,j]**2)+(yGradient[i,j]**2)+(zGradient[i,j]**2))
                xGradient[i,j]= xGradient[i,j]/norm
                yGradient[i,j]=yGradient[i,j]/norm
        
        #plotting the vector electric field using plt.quiver
        ranges = np.arange(0,self.n,1)
        x,y = np.meshgrid(ranges,ranges)
        plt.quiver(y,x,xGradient,yGradient,linewidth=0.5)
        plt.xlabel("X-axis")
        plt.ylabel("Y-axis")
        plt.title(f"Electric field for monopole n={self.n}")
        plt.savefig(f"figures/ElectricField/vectorElectricField_{self.n}.png")
        plt.show()


        #if the size of the mode is 100, we also plot a zoomed in version of the vector field
        if(self.n==100):
            ranges=np.arange(40,60)
            xGradient = np.zeros((20,20))
            yGradient = np.zeros((20,20))
            zGradient = np.zeros((20,20))
            mid =int(self.n/2)
            for i in range(20):
                for j in range(20):
                    xGradient[i,j] = self.getEx(i+40,j+40,mid)
                    yGradient[i,j] = self.getEy(i+40,j+40,mid)
                    zGradient[i,j] = self.getEz(i+40,j+40,mid)

            for i in range(20):
                for j in range(20):
                    norm = np.sqrt((xGradient[i,j]**2)+(yGradient[i,j]**2)+(zGradient[i,j]**2))
                    xGradient[i,j]= xGradient[i,j]/norm
                    yGradient[i,j]=yGradient[i,j]/norm
            x,y = np.meshgrid(ranges,ranges)
            plt.quiver(y,x,xGradient,yGradient,linewidth=0.5)
            plt.xlabel("X-axis")
            plt.ylabel("Y-axis")
            plt.title(f"Electric field for monopole n={self.n}")
            plt.savefig(f"figures/ElectricField/vectorElectricField_zoomed_{self.n}.png")
            plt.show()

    def generateMagneticData(self):
        '''
        Generates the electric field and potential data for a monopole then stores it in a data file.
        It also, calls the function to plot the vector electric field.
        '''
        x=[]
        y=[]
        potentialArray=[]
        bx=[]
        by=[]
        mid= int(self.n/2)
        distance=[]
        normalisedMagnetic=[]
        potential = self.getPotential()

        for i in range(1,self.n-1):
            for j in range(1,self.n-1):
                #avoid the middle point and boundaries as this results in error when taking the log.
                if(i==mid and j==mid):
                    continue 
                potentialArray.append(potential[i,j])
                x.append(i)
                y.append(j)
                #calculating the magnetic field in x and y direction, this follows the curl of vector potential.
                currentBx= (potential[i,(j+1)%self.n]-potential[i,(j-1)%self.n])/2
                currentBy = -(potential[(i+1)%self.n,j]-potential[(i-1)%self.n,j])/2
                bx.append(currentBx)
                by.append(currentBy)
                normalisedMagnetic.append(np.sqrt((currentBx**2)+(currentBy**2)))
                #calculate the distance to the
                distance.append(np.sqrt((mid-i)**2+(mid-j)**2))

        #store it in a data file
        withPosition = np.array((x,y,potentialArray,bx,by))
        withDistance = np.array((distance,potentialArray,normalisedMagnetic))
        np.savetxt(f"data/magneticField/potentialData_{self.n}.dat",np.transpose(withPosition))
        np.savetxt(f"data/magneticField/potentialDataVR_{self.n}.dat",np.transpose(withDistance))

        #call teh method to plot vector magnetic field.
        self.plot_MagneticField()

    
    def plot_MagneticField(self):
        '''
        Plots the vector magnetic field using plt.quiver
        The eleectric field is calculated once again here to avoid complications with meshgrid
        and this time we store the data in a 2D array of shape nxn.
        '''
        potential = self.getPotential()
        bx = np.zeros((self.n,self.n))
        by= np.zeros((self.n,self.n))


        for i in range(1,self.n-1):
            for j in range(1,self.n-1):
                bx[i,j] = (potential[i,(j+1)%self.n]-potential[i,(j-1)%self.n])/(2)
                by[i,j] = -(potential[(i+1)%self.n,j]-potential[(i-1)%self.n,j])/(2)

        #normalise the magnetic field
        for i in range(1,self.n-1):
            for j in range(1,self.n-1):
                norm = np.sqrt((bx[i,j]**2)+(by[i,j]**2))
                bx[i,j] /= norm
                by[i,j] /=norm
        
        #plot the vector field using plt.quiver
        ranges = np.arange(0,self.n,1)
        x,y = np.meshgrid(ranges,ranges)
        plt.quiver(y,x,bx,by,linewidth=0.5)
        plt.ylabel("y")
        plt.xlabel("x")
        plt.title("magnetic field for line of charge")
        plt.savefig(f"figures/magneticField/vectorMagneticField_{self.n}.png")
        plt.show()

        #if the size of the mode is 100, we also plot a zoomed in version of the vector field
        if(self.n==100):
            bx = np.zeros((20,20))
            by= np.zeros((20,20))


            for i in range(20):
                for j in range(20):
                    bx[i,j] = (potential[i+40,(j+40+1)%self.n]-potential[i+40,(j+40-1)%self.n])/(2)
                    by[i,j] = -(potential[(i+1+40)%self.n,j+40]-potential[(i+40-1)%self.n,j+40])/(2)

            for i in range(20):
                for j in range(20):
                    norm = np.sqrt((bx[i,j]**2)+(by[i,j]**2))
                    bx[i,j] /= norm
                    by[i,j] /=norm

            ranges = np.arange(40,60,1)

            x,y = np.meshgrid(ranges,ranges)
            plt.quiver(y,x,bx,by,linewidth=0.5)
            plt.ylabel("y")
            plt.xlabel("x")
            plt.title("Magnetic field for line of charge")
            plt.savefig(f"figures/magneticField/vectorMagneticField_zoomed_{self.n}.png")
            plt.show()
  