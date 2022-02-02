# Molecular_Dynamics_Simulations
 Simulations in C language of Molecular Dynamics based on "Understanding Molecular Simulation From Algorithms to Applications" by Frenkel and Smit.

## Simulacao.c
It is a molecular dynamic simulation using Verlet Integration. The parameters - number of particles (N), density (rho), initial temperature (t0), time-step (dt) and end time (tmax) - are declared inside the code, there is no need to use the file "parameters.txt". This program creates a file named "saida.txt" with the potential energy (U), the kinetic energy (K) and total energy (E) computed for each step.

## Graficos.py
This programs is used to plot a graph of the energy of the sistem. According to Frenkel's book, they must have the following behaviour: 
- the total energy must remain constant; 
- the kinetic and potential energy may vary initially, but they must oscillate around their equilibrium value around the end;

Any other behaviour means there is an error in the simulation.
