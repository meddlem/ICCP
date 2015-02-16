module constants 
  implicit none

  ! model parameters (constants):
  ! dt = timestep, rc = LJ potential cutoff length, T_init = inital temp, &
  ! m = mass, rho = number density, units: eps=1, sigma=1, m=1, &
  ! steps = number of timesteps, N = number of particles, in FCC lattice 
  real(8), parameter :: dt = 0.001d0, rc = 2.5d0, T_init = 1.585d0, & 
    rho = 0.65d0,  pi = 4d0*atan(1d0)
  integer, parameter :: steps = 3600, N = 6**3*4
  logical, parameter :: prtplt = .false.
  real(8), parameter :: L = (N/rho)**(1./3.)
  
  ! axis labels in final plot:
  character(10), parameter :: xlabel = "time", ylabel = "T", title = "plot", & 
    title1 = ""   

end module 
