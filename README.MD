# On a kinetic opinion formation model for pre-election polling

## Prerequisites

- Python: 3.7 or better
    - numpy: 1.18 or better
- [Fenics](https://fenicsproject.org/): 2019.1 or better

## Description

This code is used to run an opinion formation model for pre-election polling. It is split into two parts.

### Direct Simulation Monte Carlo (Bird's Scheme) to simulate the Boltzmann-type system of equations

We have a population with M individuals and then we choose individuals from the population non-exclusively, and perform the interaction rules, depending on the species those individuals fell into. We distribute these individuals using data from previous elections and population data. This can be changed depending on region or data. This code is also able to be edited to remove the dependence of the provided data to be used in other ways with little editting of the code.
    
Each time step in our simulation corresponds to M interactions and we run the simulation for T time steps averaging over the last 0.5% of time steps. It is encouraged to average the results of eight or more realisations to give a good approximation of the long-term density. We discretise the spatial domain $\I = [-1,1]$ in N sub-intervals. 

### Finite Element Method to simulate the Fokker-Planck system of equations

We simulate the system of Fokker-Planck equations derived in Düring Wright, the limit of which is a good approximation of the Boltzmann-type equation. We do this using Fenics.