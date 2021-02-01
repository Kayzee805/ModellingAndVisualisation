# Modelling and Visualisation in Physics: Check Point 1


## Initial run
Run the file using the following command:
```
python main.py N T

or just

main.py N T
```
Where N: Lattice size and T: Temperature of the system. 
For the animation Temperature is constant whereas, for generating data, Temperature will range
from 1 -> T and spread over 21 points. 

## User input
A user input will be taken after the initial run.

---

___Generating Data___\
input 0 to generate data for the Glauber dynamics.

input 1 to generate data for the Kawasaki dynamics.

---

___Animation___\
input 2 to start an animation of the system. Where further user input will be needed. \
0 for Glauber animation.\
1 for Kawasaki animation.

---

___Plotting___\
input 3 to plot and save the data generated.

## Data files

### Glauber Data

The data order for glauber for each column follows as:\
_Temperature, Average absolute Magnetisation, Average Energy, Specific Heat, Specific Heat error and Susceptibility_

---

### Kawasaki Data

The data order for Kawasaki for each column follows as:\
_Temperature, Average Energy, Specific Heat and Specific heat error_

---


