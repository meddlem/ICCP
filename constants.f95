module constants
  implicit none
  ! module contains all constants used in the program
  ! units: eps=1, sigma=1, m=1
  
  ! dp: compiler specific kind for double precision float
  ! lng: compiler specific kind for long integer
  ! dt: timestep 
  ! rc: lennard-jones force cutoff length
  ! rm: cutoff distance for neighbor list
  ! steps: number of timesteps
  ! N: number of particles in simulation
  ! n_bins: number of bins used for determining pair correlation function
  ! up_nbrs_list: number of iterations between each update of neighbor list
  ! meas_start: number of timesteps before measurements start
  ! n_avg: number of timesteps to average over for computing stddev
  ! n_meas: total number of measurements
  ! prtplt: determines if particle positions are plotted during iteration

  ! NOTE: IF YOU MAKE ANY CHANGES HERE RECOMPILE ALL MODULES: "make -B" 
  integer, parameter :: dp = selected_real_kind(15,307)
  integer, parameter :: lng = selected_int_kind(8)

  real(dp), parameter :: dt = 0.0004_dp 
  real(dp), parameter :: rc = 2.5_dp
  real(dp), parameter :: rm = 3.3_dp
  real(dp), parameter :: pi = 4._dp*atan(1._dp) 
  
  integer, parameter :: steps = 5000
  integer, parameter :: N = 6**3*4
  integer, parameter :: n_bins = 120
  integer, parameter :: up_nbrs_list = 25
  integer, parameter :: meas_start = 1000
  integer, parameter :: n_avg = 100 
  integer, parameter :: n_meas = steps + 1 - meas_start
  
  logical, parameter :: prtplt = .false.
end module
