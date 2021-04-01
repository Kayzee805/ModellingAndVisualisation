# Modelling and Visualisation in Physics: Check Point 3: Poisson's equation

## Initial run
Run the file using the following command from the 'CP3/cahnHilliard/' directory:
```
python main.py N epsilon
```
Where N: size of the lattice.
epsilon: Precision a lattice must meet to converge.

## User input
A user input will be taken after the initial run.

---

### Task

__0: Generate Data__\

Two user input will be taken. 
### First to select the model to generate data for.
0: Generate data for the monopole.
1: Generate data for a charged wire.
2: Generate data for SOR.

### Second to select the type of update method to use.
0: Jacobi update method
1: Gauss-seidel update method

Note: option is given within the main method of "main.py" to choose if user wants to enable the checkerboard method.

---

__1: Plot existing data__\
Plot the graphs and calculate the fits for the potential and electric field against distance for the monopole and the potential and magnetic field against distance for the charged wire.

---

## Data file columns

### Electric field (the monopole)
The column order for all potentialData_...dat files under "data/electricField" directory has the form:
*X,Y,Z,Potential,Ex,Ey,Ez*. 
Where X,Y,Z are the coordinates in the lattice and Ex,Ey,Ez are the electric field in the respective direction.

The column order for all potentialDataVr_....dat files under "data/electricField" directory has the form:
*Distance,Potential,ElectricField*
Where distance is the distance to the monopole and ElectricField is the normalised electricField at each point in the lattice.

### Magnetic field (charged wire)
The column order for all potentialData_...dat files under "data/magneticField" directory has the form:
*X,Y,Potential,Bx,By*. 
Where X,Y are the coordinates in the lattice and Bx,By are the magnetic field in the respective direction.

The column order for all potentialDataVr_....dat files under "data/magneticField" directory has the form:
*Distance,Potential,MagneticField*
Where distance is the distance to the monopole and MagneticField is the normalised electricField at each point in the lattice.


### SOR (2D cut of charged wire)
The coloumn order for all SOR_....dat files under "data/SOR" directory has the form:
*omega, convergenceStep*
Where omega ranges from 1->2 and convergenceStep is the number of steps/updates required for the system to converge.