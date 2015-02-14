program main
  ! modules to use 
  use plplot ! plotting library 
  use Plotroutines ! module for plotting particles
  use Inits ! module for initializing model
  use Interactions ! module for calculating the interaction forces
  implicit none

  ! model parameters (constants):
  ! dt = timestep, rc = LJ potential cutoff length, Tinit = inital temp, &
  ! m = mass, rho = number density, units: eps=1, sigma=1, m=1, &
  ! steps = number of timesteps, N = number of particles, in FCC lattice 
  real(8), parameter :: dt = 0.001d0, rc = 2.5d0, Tinit = 1d0, rho = 0.5d0, & 
    eps = 1d0, sigma = 1d0, m = 1d0  
  integer, parameter :: steps = 1000, N = 6**3*4
  
  ! axis labels in plot:
  character(10), parameter :: xlabel = "time", ylabel = "T", label = "plot"   
  
  ! declare variables:
  real(8) :: x(steps+1), xrange(2), yrange(2), T(steps+1), E(steps+1), EV, &
    r(N,3), p(N,3), F(N,3), L = (N/rho)**(1./3.) 
  integer :: i, starttime, endtime
  
  ! initialize r, p, F, and 3d plot:
  call InitR(r,L,N)
  call InitP(p,Tinit,N)
  call Force(F,EV,r,rc,L)
  call ParticlePlotinit(-0.1d0*L,1.1d0*L) 
  
  ! calculate temp (T) and total energy (E)
  T(1) = sum(p**2)/N
  E(1) = EV + 0.5d0*sum(p**2)

  ! print energy, momentum, temp
  print *, "Psum t=0 : ", sum(p,1) 
  print *, "E t=0 : ", E(1)/N 
  print *, "T t=0 : ", T(1)
  call system_clock(starttime)
  
  ! time integration using the "velocity Verlet" algorithm: 
  do i = 1,steps
    ! plot particle positions
    call ParticlePlot(r) 
    
    r = r + p*dt + 0.5d0*F*(dt**2) !update positions
    r = r - floor(r/L)*L ! enforce PBC
    p = p + 0.5d0*F*dt ! update momentum (1/2)
    call Force(F,EV,r,rc,L) ! update force to next timestep
    p = p + 0.5d0*F*dt ! update momentum (2/2)
    
    ! calculate temp, total energy at timestep i
    T(i+1) = sum(p**2)/N
    E(i+1) = EV + 0.5d0*N*T(i+1)
  enddo
  call system_clock(endtime)
  call plend()

  print *, "runtime =", endtime-starttime
  print *, "Psum final: ", sum(p,1)
  print *, "E final: ", E(steps+1)/N 
  print *, "T final: ", T(steps+1)

  ! variables for final plot
  x = (/(i,i=0, steps)/)
  xrange = [0, steps]
  yrange = [0d0, maxval(E/N)*1.1d0]
  
  call LinePlot(x,E/N,xrange,yrange,xlabel,ylabel,label)
end program main 
