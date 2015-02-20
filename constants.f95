module constants
  implicit none
  ! contains all constants used in the program

  ! units: eps=1, sigma=1, m=1

  ! dt = timestep, 
  ! rc = lennard-jones force cutoff length
  ! steps = number of timesteps
  ! N = number of particles in simulation
  ! up_nbrs_list = number of iterations between each update of neighbor list
  ! n_bins = number of bins used for determining pair correlation function
  ! meas_start = number of timesteps before measurements start

  ! NOTE: IF YOU MAKE CHANGES HERE RECOMPILE ALL MODULES: "make -B" 
  
  real(8), parameter :: dt = 0.001d0 
  real(8), parameter :: rc = 2.5d0
  real(8), parameter :: rm = 3.3d0
  real(8), parameter :: pi = 4d0*atan(1d0) 
  
  integer, parameter :: steps = 10000
  integer, parameter :: N = 6**3*4
  integer, parameter :: n_bins = 120
  integer, parameter :: up_nbrs_list = 25
  integer, parameter :: meas_start = 1000
  integer, parameter :: n_meas = steps + 1 - meas_start
  
  logical, parameter :: prtplt = .false.

end module
