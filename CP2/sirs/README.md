# Modelling and Visualisation in Physics: Check Point 2: SIRS

## Install dependencies
```
pip install -r requirements.txt
```

## Initial run
Run the file using the following command from the 'CP2/sirs/' directory:
```
python main.py N p1 p2 p3
```
Where N: Lattice/system length. Creates a NxN lattice model.
p1: Probability of a susceptible site to be infected
p2: Probability of a infected site to recover
p3: Probability of a recovered site to be susceptible.

## User input
A user input will be taken after the initial run.

---

__1: Animation__\
Input 1 to start an animation of the system for 10,000 sweeps. The lattice is updated each sweep.

__2: Plot Data__\
Input 2 to plot all data. This should only be ran if data already exists if not, data should be generated for each task first. Figures are saved in the figures directory.

__3: Task 3__\
Input 3 to generate data for Task 3 of the SIRS check point. All data are saved in the data directory. 

__4: Task 4__\
Input 4 to generate data for Task 4 of the SIRS check point. 

__5: Task5__\
Input 5 to generate data for Task 5.

Further input will be required at this point. Note: the program will run for 5 simulation so can take up alot of time. Which is why a cython code has been provided.
Running it in parallel is recommended to generate the data.
#### Input 0 to generate data for task 5 part one. 50x50 lattice model.
#### Input 1 to generate data for task 5 part two. 100x100 lattice model.

---


## cython
A cython code has been provided to generate the data quicker. However, as it is not required proper documentation has not been provided yet. However, the methods have the same funcitonality
as the original pyhton code. The cython code can be ran from the cython directory.

---

__Setup__\

run the follwing in the cython directory:
```
python setup.py build_ext --inplace
```

then the each method can be ran using:
```
python runCython.py
```
Where a method can be called using cythonCode.task4(size,sweeps) to run the method in task4. This can be done for any methods in the cythonCode file.

---