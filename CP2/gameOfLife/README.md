# Modelling and Visualisation in Physics: Check Point 2: Game Of Life

## Install dependencies
```
pip install -r requirements.txt
```

## Initial run
Run the file using the following command from the 'CP2/gameOfLife/' directory:
```
python main.py N sweeps
```
Where N: Lattice/system length. Creates a NxN lattice model.
sweeps: Number of updates for the system. One update per sweep.

## User input
A user input will be taken after the initial run.

---

### Initialisation type

0: Random initialisation
1: Blinker initialisaiton
2: Glider initialisation

---

### Task

__1: Animation__\
Animates the system for the given input and initailisation type.

__2: Generate Histogram data__\
Generates the data for the histogram

__3: Generate centre of mass data__\
Recommended to carry out for the glider initialisation as the task required.
Generates the data for centre of mass and outputs the x and y velocity of the glider.

__4: Plot data__\
Plots all the data, the histogram, centre of mass for x and y components.

---

