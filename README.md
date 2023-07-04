# Final Project of Marco Mendívil Carboni

This repository contains the code of my final project at university. The objective of the project is to simulate with Langevin (and Brownian) dynamics a coarse-grained model of a polymer and investigate some of its equilibrium and dynamical properties.

## Contents

| filename | description |
|-|-|
| `forces.md` | derivation of the forces from the model's potentials |
| `units.md` | units and parameters that I'm using to model DNA |
| `Alex_pol.h` | functions to calculate the Alexander polynomial |
| `analyse-confinement.c` | C program to analyse confinement simulations |
| `analyse-equilibrium.c` | C program to analyse equilibrium simulations |
| `analyse-translocation.c` | C program to analyse translocation simulations |
| `check-performance.cu` | CUDA program to check the performance |
| `simulate-confinement.cu` | CUDA program to make confinement simulations |
| `simulate-equilibrium.cu` | CUDA program to make equilibrium simulations |
| `simulate-translocation.cu` | CUDA program to make translocation simulations |
| `visualize-confinement.tcl` | Tcl script to visualize confinement simulations in vmd |
| `visualize-equilibrium.tcl` | Tcl script to visualize equilibrium simulations in vmd |
| `visualize-translocation.tcl` | Tcl script to visualize translocation simulations in vmd |

## Model

The model consists simply of $N$ beads (point particles) upon which act four forces:
- An ellastic force between adjacent beads
- A bending force between adjacent bonds
- External forces:
    - a pulling force in the equilibrium program
    - the force exerted by the pore and wall in the translocation program
    - the force exerted by the spherical wall in the confinement program
- A excluded volume force between all beads (with a minimum length scale of self-interaction of 2)

## Equilibrium Simulation

The equilibrium simulation program has only one argument: the relative path to the directory of the simulation `sim_dir`. First it will try to read a file with the name `sim_dir/parameters.dat`. In this file the user must specify:
- `T` the temperature (in Kelvin)
- `eta` the viscosity of the medium
- `N` the number of beads
- `f_pull` the pulling force
- `IC_type` the type of initial condition:
    - `f` file (`sim_dir/initial-condition.bin`)
    - `r` random
    - `l` linear 
- `meas_freq` the measure frequency at which the positions will be saved
- `sim_time` the desired simulation time

Afterwards it will create the initial condition specified unless there is a file with the name `sim_dir/checkpoint.bin` in which case it will resume the previous simulation which generated it. The program generates three output files:
- The `sim_dir/initial-condition.gro` file has the coordinates of each bead in the initial condition, it is only created the first time you run this simulation
- The `sim_dir/trajectory-file-X.trr` (with `X` a file index) has the trajectory (only positions)
- The `sim_dir/checkpoint.bin` file is a checkpoint to continue the simulation from that point

## Equilibrium Analysis

The equilibrium analysis program also takes as argument the name of the simulation. It calculates the gyration radius and the end-to-end distance (and vector) for every frame and saves them in the text file `sim_dir/time-series.dat`. Then it reads this file and computes the mean and variance of this quantities. It also makes an estimation of the standard deviation of their sample mean with the "blocking" method. If this method isn't successful it warns the user. It also estimates the thermalization time with the Standard Marginal Error Rule. All this is written into the file `sim_dir/statistics.dat`. With the positions at the end of each trajectory file it calculates as well the Alexander polynomial after closing in a reasonable way the chain. The alexander polynomial is converted into the knot type (when posible) and written into the file `sim_dir/knot-type.dat`. Finally the histograms of the distributions of the end-to-end distance and the distance to the center of mass are computed and saved in `sim_dir/rpd-histogram.dat`.

## Translocation Simulation

The translocation simulation program has only one argument: the relative path to the directory of the simulation `sim_dir`. First it will try to read a file with the name `sim_dir/adjustable-parameters.dat`. In this file the user must specify:
- `T` the temperature (in Kelvin)
- `N` the number of beads
- `V` the voltage difference (in mV)

Afterwards it will search in the simulation directory how many translocations have been started before and set the simulation index `sim_idx` appropriately. Since it only looks at how many previous translocations there are if some file is deleted the next time the program is executed it will overwrite previous simulations. The program generates three output files:
- The `sim_dir/initial-configuration-$sim_idx.gro` files have the coordinates of each bead in the initial configuration of each translocation
- The `sim_dir/pore-borders-times-$sim_idx.gro` files contain the times at which each basepair enters the pore for the first time and exits for the last time
- The `sim_dir/trajectory-positions-$sim_idx.gro` files have the trajectory positions which are saved with a fixed frequency

## Translocation Analysis

The translocation analysis program also takes as argument the name of the simulation. After reading the parameters it will search in the simulation directory how many translocations have finished. If some file is deleted the program will fail to read the next translocations and terminate. For each simulation it reads the initial configuration and the the pore borders times. It finds the knot type of each IC and saves it if it's different from the unknot. It also reads how many unsuccessful simulations there have been and the total translocation time $\tau$. It calculates the mean and standard deviation of the translocation time and writes it in the file `sim_dir/analysis-statistics.dat` along with the knot type list and the success rate. Finally it also calculates the mean waiting times and saves them in `sim_dir/analysis-waiting-times.dat`.

## Confinement Simulation

The confinement simulation program has only one argument: the relative path to the directory of the simulation `sim_dir`. First it will try to read a file with the name `sim_dir/adjustable-parameters.dat`. In this file the user must specify:
- `T` the temperature (in Kelvin)
- `N` the number of beads
- `R` the radius of the sphere (in program units)

Afterwards it will search in the simulation directory how many simulations have been started before and set the simulation index `sim_idx` appropriately. The program generates two output files:
- The `sim_dir/initial-configuration-$sim_idx.gro` files have the coordinates of each bead in the initial configuration of each simulation
- The `sim_dir/trajectory-positions-$sim_idx.gro` files have the trajectory positions which are saved with a fixed frequency

## Confinement Analysis

The confinement analysis program also takes as argument the name of the simulation. After reading the parameters it will search in the simulation directory how many simulations have finished. For each simulation it calculates the gyration radius for every frame and saves them in the text file `sim_dir/Rg2-time-series-$sim_idx.dat`. Then it reads this file, computes the mean and variance of this quantity and estimates the thermalization time. With the positions at the end of each trajectory it calculates as well the Alexander polynomial and converts it into the knot type (when posible). All this is written into the file `sim_dir/analysis-statistics.dat`. Then the histograms of the distributions of the distance to the center of the sphere and the distance to the center of mass are computed and saved in `sim_dir/analysis-rpd-histogram.dat`. Finally the average over all simulations of an order parameter is calculated and written to `sim_dir/analysis-order-parameter.dat`.

## Visualization

The equilibrium visualization script has two arguments: the simulation directory and an index. First it loads the initial configuration in `.gro` format since VMD needs to know the structure (elements) of the molecule to represent it. Then it creates the bonds between the beads, sets the beads radius and sets the representation style. Finally it loads the appropriate trajectory file. The translocation visualization script works in the same way but the index now represents distinct simulations and it also draws the pore at the end. The confinement visualization script works the same as the last one but at the end it draws the sphere.

Actually in this model each monomer does not represent a single atom but rather a lot of them. The ideal thing would be to have custom beads, but this doesn't seem possible so we use an arbitrary atom. If in the future different monomers make up our polymer we will asign different atoms to each one.
