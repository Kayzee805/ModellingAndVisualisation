import numpy as np 
import time
import matplotlib.pyplot as plt

class poisson(object):
    def __init__(self,int n,double epsilon=1, int dx=1,double minW=1,double maxW=2):
        self.n=n
        self.epsilon=epsilon
        self.dx=dx
        self.rho = np.zeros((n,n,n))
        self.lattice=np.zeros((n,n,n))
        self.nextLattice=np.zeros((n,n,n))
        self.minW=minW
        self.maxW=maxW

    

    
    def getPotential(self):
        return self.lattice[:,:,int(self.n/2)]

    def setBoundaries3D(self,Array):
        #set boundaries to 0?
        Array[0,:,:]= 0
        Array[:,0,:]= 0
        Array[:,:,0]= 0
        Array[-1,:,:]= 0
        Array[:,-1,:]= 0
        Array[:,:,-1]= 0
    
    def setBoundaries2D(self,Array):
        Array[0,:]=0
        Array[-1,:]=0
        Array[:,0]=0
        Array[:,-1]=0 
    
    def setPointCharge2D(self):
        mid=int(self.n/2)
        self.rho= np.zeros((self.n,self.n))
        self.rho[mid,mid]=1
    
    def setPointCharge3D(self):
        mid=int(self.n/2)
        self.rho = np.zeros((self.n,self.n,self.n))
        self.rho[mid,mid,mid]=1
    
    def setChargedWire3D(self):
        mid=int(self.n/2)
        for z in range(1,self.n-1):
            self.rho[mid,mid,z]=1
    
    def calculateNearestNeighbours3D(self,i,j,k):
        #find the 6 neighbours, so +- for all i,j,k
        n=self.n
        return self.lattice[(i+1)%n,j,k]+self.lattice[(i-1)%n,j,k]+self.lattice[i,(j-1)%n,k]+self.lattice[i,(j+1)%n,k]+self.lattice[i,j,(k-1)%n]+self.lattice[i,j,(k+1)%n]

    def dependentNeighbours(self,i,j,k):
        n=self.n
        return self.lattice[(i-1)%n,j,k]+self.lattice[i,(j-1)%n,k]+self.lattice[i,j,(k-1)%n]


    def dependentNeighbours2D(self,i,j,ARRAY):
        n=self.n
        return ARRAY[(i-1)%n,j]+ARRAY[i,(j-1)%n]

    def jacobiUpdate(self):
        print("Starting jacobi")
        convergence = self.epsilon+1
        counter=0
        sumBefore = np.sum(self.lattice)
        while(convergence>=self.epsilon):
            t1=time.time()
            self.setBoundaries3D(self.lattice)
            
            self.nextLattice= (self.rho+(
            np.roll(self.lattice,1,axis=0)+np.roll(self.lattice,-1,axis=0)
            +np.roll(self.lattice,1,axis=1)+np.roll(self.lattice,-1,axis=1)+
            np.roll(self.lattice,1,axis=2)+np.roll(self.lattice,-1,axis=2)))/6

            self.setBoundaries3D(self.nextLattice)
            sumAfter=np.sum(self.nextLattice)
            convergence = abs(sumAfter-sumBefore)
            self.lattice = np.copy(self.nextLattice)
            sumBefore=sumAfter
            counter+=1

            if(counter%100==0):
                print(f"Counter={counter} convergence={convergence} time={time.time()-t1}s")

        print("finished jacobi")

    
    def gaussSeidelUpdate_roll(self):
        print("Starting gauss roll")
        convergence = self.epsilon+1
        counter=0
        sumBefore = np.sum(self.lattice)
        while(convergence>=self.epsilon):
            t1=time.time()
            independentNeighbours = np.roll(self.lattice,-1,axis=0)+np.roll(self.lattice,-1,axis=1)+np.roll(self.lattice,-1,axis=2)

            for i in range(1,self.n-1):
                for j in range(1,self.n-1):
                    for k in range(1,self.n-1):
                        self.lattice[i,j,k] = (independentNeighbours[i,j,k]+self.dependentNeighbours(i,j,k)+self.rho[i,j,k])/6
            sumAfter=np.sum(self.lattice)
            convergence = abs(sumAfter-sumBefore)
            counter+=1
            
            if(counter%100==0):
                print(f"Counter={counter} convergence={convergence} time={time.time()-t1}s sumBefore={sumBefore} sumAfter={sumAfter}")
            sumBefore=sumAfter

        print(f"Counter={counter} convergence={convergence} time={time.time()-t1}s")

        print("Finished gauss roll")

    def gaussSeidel_CheckerBoard(self):
        print("STARTING CHECKER Board")
        mask = np.indices((self.n,self.n,self.n)).sum(axis=0)%2
        convergence=self.epsilon+1
        counter=0
        sumBefore = np.sum(self.lattice)
        while(convergence>=self.epsilon):
            t1=time.time()
            self.lattice[mask==0]=0

            self.lattice = (np.roll(self.lattice,1,axis=0)+np.roll(self.lattice,1,axis=1)+np.roll(self.lattice,1,axis=2)+np.roll(self.lattice,-1,axis=0)+np.roll(self.lattice,-1,axis=1)+np.roll(self.lattice,-1,axis=2)+self.rho)/6
            self.setBoundaries3D(self.lattice)
            self.lattice[mask==1]=0

            self.nextLattice = (np.roll(self.lattice,1,axis=0)+np.roll(self.lattice,1,axis=1)+np.roll(self.lattice,1,axis=2)+np.roll(self.lattice,-1,axis=0)+np.roll(self.lattice,-1,axis=1)+np.roll(self.lattice,-1,axis=2)+self.rho)/6
            self.setBoundaries3D(self.nextLattice)
            self.lattice = np.add(self.lattice,self.nextLattice)

            sumAfter = np.sum(self.lattice)
            convergence= abs(sumAfter-sumBefore)
            counter+=1
            if(counter%100==0):
                print(f"Counter={counter} convergence={convergence} time={time.time()-t1}s sumBefore={sumBefore} sumAfter={sumAfter}")
            sumBefore = sumAfter

        print(f"Counter={counter}  convergence={convergence}")

        print("Finished CHECKER Board")

    def generate_SOR_Point(self):
        w=np.arange(1,2,0.01)
        stableSweeps=[]
        mask = np.indices((self.n,self.n)).sum(axis=0)%2
        print("starting sor point for checkerboard")
        for x in w:
            t1=time.time()
            self.lattice=np.zeros((self.n,self.n))
            self.setPointCharge2D()

            convergence=self.epsilon+1
            sumBefore=0
            counter=0
            print(f"STARTING point 3d sor for w={x}")

            while(convergence>self.epsilon):
    
                original = np.copy(self.lattice)
                self.lattice[mask==0]=0
                self.lattice =(1-x)*original+x*((np.roll(self.lattice,1,axis=0)+np.roll(self.lattice,1,axis=1)+np.roll(self.lattice,-1,axis=0)+np.roll(self.lattice,-1,axis=1)+self.rho)/4)
                self.setBoundaries2D(self.lattice)
                self.lattice[mask==1]=0
                self.nextLattice = np.copy(self.lattice)

                self.nextLattice = (1-x)*(original)+x*((np.roll(self.nextLattice,1,axis=0)+np.roll(self.nextLattice,1,axis=1)+np.roll(self.nextLattice,-1,axis=0)+np.roll(self.nextLattice,-1,axis=1)+self.rho)/4)
                self.setBoundaries2D(self.nextLattice)
                self.nextLattice[mask==0]=0

                self.lattice=np.add(self.lattice,self.nextLattice)
                sumAfter=np.sum(self.lattice)
                convergence=abs(sumAfter-sumBefore)

                counter+=1
                if(counter%1000==0):
                    print(f"Counter={counter} convergence={convergence} sumBefore={sumBefore}  sumAfter={sumAfter}  ")
                sumBefore=sumAfter
            stableSweeps.append(counter)
            print(f"Finsihed at counter={counter} at w={x} at time={time.time()-t1}s")
        
        plt.plot(w,stableSweeps)
        plt.scatter(w,stableSweeps,s=5,color='k')
        plt.xlabel("w")
        plt.ylabel("sweeps")
        plt.title(f"Steps required to stabilise for different omega values. SOR_{self.n}")
        plt.savefig(f"figures/SOR__{self.n}.png")
        plt.show()
        combined=np.array((w,stableSweeps))
        np.savetxt(f"data/SOR_{self.n}.dat",np.transpose(combined),fmt='%.4f')
        print(f"Minimum at {np.min(stableSweeps)}")


    def generate_SOR(self):
        w=np.arange(1,1.99,0.01)
        stableSweeps=[]
        cdef int n=self.n
        cdef double x,t1,convergence,sumAfter,sumBefore
        cdef int i,j,counter
        print(f"Generating normal 3 loops SOR")
        for x in w:
            t1=time.time()
            mid=int(self.n/2)-1
            self.lattice=np.zeros((n,n))
            self.rho=np.zeros((self.n,self.n))
            self.rho[mid,mid]=1
            convergence=self.epsilon+1
            counter=0
            print(f"STARTING W={x}")
            sumBefore=0
            while(convergence>self.epsilon):
                convergence=0
                independentNeighbours = np.roll(self.lattice,-1,axis=0)+np.roll(self.lattice,-1,axis=1)
                for i in range(1,self.n-1):
                    for j in range(1,self.n-1):
                        old = self.lattice[i,j]
                        self.lattice[i,j] = ((1-x)*self.lattice[i,j])+x*(self.dependentNeighbours2D(i,j,self.lattice)+independentNeighbours[i,j]+self.rho[i,j])/4
                        #self.lattice[i,j]= (1-x)*old +x*(self.lattice[(i+1)%n,j]+self.lattice[(i-1)%n,j]+self.lattice[i,(j-1)%n]+self.lattice[i,(j+1)%n]+self.rho[i,j])/4
                        convergence+=abs(self.lattice[i,j]-old)
               
                if(counter%100==0):
                    print(f"counter={counter} convergence={convergence}")
                counter+=1

            
            print(f"Finsihed at counter={counter} at w={x} at time={time.time()-t1}s")
            stableSweeps.append(counter)

        plt.plot(w,stableSweeps)
        plt.scatter(w,stableSweeps,s=5,color='k')
        plt.xlabel("w")
        plt.ylabel("sweeps")
        plt.title(f"Steps required to stabilise for different omega values. SOR_{self.n}")
        plt.savefig(f"figures/SOR__{self.n}.png")
        plt.show()
        combined=np.array((w,stableSweeps))
        np.savetxt(f"data/SOR_{self.n}.dat",np.transpose(combined),fmt='%.4f')
        print(f"Minimum at {np.min(stableSweeps)}")


    def distantCaculator(self,i,j,k):
        mid = int(self.n/2)
        return np.sqrt((mid-i)**2+(mid-j)**2+(mid-k)**2)


    def generateMagneticData(self):
        #x,y, potential,bx,by
        #distance, potential, magnetic field
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
                if(i==mid and j==mid):
                    continue 
                potentialArray.append(potential[i,j])
                x.append(i)
                y.append(j)
                currentBx= (potential[i,(j+1)%self.n]-potential[i,(j-1)%self.n])/2
                currentBy = -(potential[(i+1)%self.n,j]-potential[(i-1)%self.n,j])/2
                bx.append(currentBx)
                by.append(currentBy)
                normalisedMagnetic.append(np.sqrt((currentBx**2)+(currentBy**2)))
                distance.append(np.sqrt((mid-i)**2+(mid-j)**2))

        withPosition = np.array((x,y,potentialArray,bx,by))
        withDistance = np.array((distance,potentialArray,normalisedMagnetic))
        np.savetxt(f"data/potentialDataMagnetic_{self.n}.dat",np.transpose(withPosition))
        np.savetxt(f"data/potentialDataVRMagnetic_{self.n}.dat",np.transpose(withDistance))

        plt.scatter(np.log2(distance),np.log2(potentialArray),s=3,marker='x')
        plt.title("Potential for charged wire")
        plt.xlabel("log(distance)")
        plt.ylabel("log(potential)")
        plt.savefig(f"figures/magneticFieldPotential_{self.n}")
        plt.show()

        plt.scatter(np.log2(distance),np.log2(normalisedMagnetic),s=3,marker='x')
        plt.title("Magnetic field for charged wire")
        plt.xlabel("log(distance)")
        plt.ylabel("log(magnetic)")
        plt.savefig(f"figures/magneticFieldPotentialVR_{self.n}.png")
        plt.show()

        self.plot_MagneticField()


    def generateElectricData(self):
        #x,y,z, potential, ex, ey, ez
        #distance to centre, potential and electricfield
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
        np.savetxt(f"data/potentialData_{self.n}.dat",np.transpose(withPosition))
        np.savetxt(f"data/potentialDataVR_{self.n}.dat",np.transpose(withDistance))
       # np.savetxt(f"data/potemntialDataVR_LOG{self.n}.dat",np.log2(np.transpose(withdistance)))
        plt.scatter(np.log2(distance),np.log2(potentialArray),s=3,marker='x')
        plt.xlabel("distance")
        plt.ylabel("potential")
        plt.title("log(distance) vs log(potential)")
        plt.savefig(f"figures/potentialVR_{self.n}.png")
        plt.show()
        plt.scatter(np.log2(distance),np.log2(electricField),s=3,marker='x')
        plt.xlabel("distance")
        plt.ylabel("eField")
        plt.title("log(distance) vs log(ElectricField)")
        plt.savefig(f"figures/electricFieldVR_{self.n}.png")
        plt.show()
        self.plotElectricField()

        distance=np.log2(distance)
        potentialArray=np.log2(potentialArray)
        electricField = np.log2(electricField)
      
        
        newDistance = distance[(distance>0.5) & (distance<1.5)]
        newPotential=potentialArray[(distance>0.5) & (distance<1.5)]
        print(f"Using potential {newPotential[0]}")
        xfit,xin = np.polyfit(newDistance,newPotential,1)
        print(f"Potential Fit = xfit={xfit}  and xin={xin}")

        newDistance = distance[(distance>0.5) & (distance<2)]
        newElectric = electricField[(distance>0.5) & (distance<2)]
        xfit,xin = np.polyfit(newDistance,newElectric,1)
        print(f"ElectricField  Fit = xfit={xfit}  and xin={xin}")



        #take the fit of log(dsitance) and log(potential)
    
    def plotElectricField(self):
        # allArray = np.loadtxt(f"data/potentialData_{self.n}.dat")
        # ex = allArray[:,4]
        # ey = allArray[:,5]
        # ez = allArray[:,6]
        # xx,yy = np.meshgrid(ex,ey)
        # xGradient = np.zeros((self.n,self.n))
        # yGradient = np.zeros((self.n,self.n))
        # zGradient = np.zeros((self.n,self.n)) #because we now that ez is constant?
        # for i in range(self.n):
        #     for j in range(self.n):

        xGradient = np.zeros((self.n,self.n))
        yGradient = np.zeros((self.n,self.n))
        zGradient = np.zeros((self.n,self.n))
        mid =int(self.n/2)
        for i in range(1,self.n-1):
            for j in range(1,self.n-1):
                xGradient[i,j] = self.getEx(i,j,mid)
                yGradient[i,j] = self.getEy(i,j,mid)
                zGradient[i,j] = self.getEz(i,j,mid)

        for i in range(1,self.n-1):
            for j in range(1,self.n-1):
                norm = np.sqrt((xGradient[i,j]**2)+(yGradient[i,j]**2)+(zGradient[i,j]**2))
                xGradient[i,j]= xGradient[i,j]/norm
                yGradient[i,j]=yGradient[i,j]/norm
        
        ranges = np.arange(0,self.n,1)
        x,y = np.meshgrid(ranges,ranges)
        plt.quiver(y,x,xGradient,yGradient,linewidth=0.5)
        plt.xlabel("X-axis")
        plt.ylabel("Y-axis")
        plt.title(f"Electric field for monopole n={self.n}")
        plt.savefig(f"figures/ElectricFieldPointCharge_{self.n}.png")
        plt.show()

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
            plt.savefig(f"figures/ElectricFieldPointChargeZoomed_{self.n}.png")
            plt.show()


    def getPotential(self):
        return self.lattice[:,:,int(self.n/2)]

    def getEx(self,i,j,k):
        gradient = self.lattice[(i+1)%self.n,j,k]- self.lattice[(i-1)%self.n,j,k]
        return -gradient/2

    def getEy(self,i,j,k):
        gradient = self.lattice[i,(j+1)%self.n,k]- self.lattice[i,(j-1)%self.n,k]
        return -gradient/2

    def getEz(self,i,j,k):
        gradient = self.lattice[i,j,(k+1)%self.n]- self.lattice[i,j,(k-1)%self.n]
        return -gradient/2
    def getElectricField(self,i,j,k):
        gradient =-(self.lattice[i,j,(k+1)%self.n]- self.lattice[i,j,(k-1)%self.n]+self.lattice[i,(j+1)%self.n,k]- self.lattice[i,(j-1)%self.n,k]+self.lattice[(i+1)%self.n,j,k]- self.lattice[(i-1)%self.n,j,k])/2
        return gradient
        
    
    def plot_MagneticField(self):
        #x,y,potential,bx,by
        #distance, potential, magnetic field
        potential = self.getPotential()
        bx = np.zeros((self.n,self.n))
        by= np.zeros((self.n,self.n))


        for i in range(1,self.n-1):
            for j in range(1,self.n-1):
                bx[i,j] = (potential[i,(j+1)%self.n]-potential[i,(j-1)%self.n])/(2)
                by[i,j] = -(potential[(i+1)%self.n,j]-potential[(i-1)%self.n,j])/(2)
        for i in range(1,self.n-1):
            for j in range(1,self.n-1):
                norm = np.sqrt((bx[i,j]**2)+(by[i,j]**2))
                bx[i,j] /= norm
                by[i,j] /=norm
        
        
        ranges = np.arange(0,self.n,1)

        x,y = np.meshgrid(ranges,ranges)
        plt.quiver(y,x,bx,by,linewidth=0.5)
        plt.ylabel("y")
        plt.xlabel("x")
        plt.title("magnetic field for line of charge")
        plt.savefig(f"figures/MagneticField_{self.n}.png")
        plt.show()


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
            plt.savefig(f"figures/MagneticFieldZoomed_{self.n}.png")
            plt.show()
  

  
def sor(n,epsilon):
    model = poisson(n,epsilon)
    model.generate_SOR()
    #model.generate_SOR3D()
  #  model.generate_SOR_Point()
  #  model.generate_SOR_Point3D()