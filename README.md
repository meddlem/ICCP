# ICCP
This program is a Molecular Dynamics simulation, of an Argon gas in the NVT 
ensemble, based on the velocity verlet algorithm. 

module structure:
- constants contains all constants used througout the simulation + definitions
  for double/integer precision.

- initialize: contains subroutines used for generating initial positions in an 
  FCC lattice, and initial velocities sampled from a maxwell-Boltzmann 
  distribution.

- interactions: contains the subroutine used to calculate interaction forces, 
  potential energies between particles, based on a neighbor list. 
  Also contains the subroutine for generating the neighbor list.

- main_functions: contains additional functions used in the main program 

- plotroutines contains all subroutines used to plot particles with plplot, 
  and a subroutines for creating simple lineplots through gnuplot


