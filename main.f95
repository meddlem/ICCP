program main
  ! modules to use 
  use plplot !plotting library 
  use Plotroutines !module for plotting particles
  use Inits !module for initializing model
  use Interactions !module for calculating the interaction forces
  implicit none

  ! model parameters (constants):
  ! dt = timestep, rc = potential cutoff length, Tinit ~ inital temp, m = mass, L=length, units: eps=1, sigma=1, m=1
  real(8), parameter :: dt = 0.001d0, rc = 2.5d0, Tinit = 1d0, L = 14d0, eps = 1d0, sigma = 1d0, m = 1d0  
  integer, parameter :: steps = 1000, N = 7**3*4 !number of particles, must fit the FCC lattice 
  character(10), parameter :: xlabel = "time", ylabel = "T", label = "plot" !axislabels
  ! declare variables
  real(8) :: x(steps+1), xrange(2), yrange(2), T(steps+1), E(steps+1), EV, r(N,3), v(N,3), F(N,3) 
  integer :: i

  ! initialize r, v, F
  call InitCell(r,L,N)
  call InitVel(v,Tinit,N)
  call Force(F,EV,r,rc,L)

  ! initialize 3Dplot
  call ParticlePlotinit(-0.1d0*L,1.1d0*L) 

  T(1) = sum(v**2)/N ! calculate T
  E(1) = EV + 0.5d0*sum(v**2) ! calculate the total energy

  print *, "Psum t=0: ", sum(v,1) !test if total momentum is constant 
  print *, "Initial E: ", E(1)/N 
  print *, "T initial: ", T(1)

  ! time integration using velocity Verlet algorithm: 
  do i = 1,steps
    call ParticlePlot(r) !calls plot points

    r = r + v*dt + 0.5d0*F*(dt**2) !update positions
    r = r - floor(r/L)*L !enforce PBC
    v = v + 0.5d0*F*dt !update velocity (1/2)
    call Force(F,EV,r,rc,L) !update force to next timestep
    v = v + 0.5d0*F*dt !update velocity (2/2)

    T(i+1) = sum(v**2)/N ! calculate the temperatue at timestep i
    E(i+1) = EV + 0.5d0*N*T(i+1) ! calculate the total energy
  enddo

  call plend()

  print *, "T final: ", T(steps+1)
  print *, "E final: ", E(steps+1)/N 
  print *, "Psum final: ", sum(v,1)

  x = (/(i,i=0, steps)/)
  xrange = [0, steps]
  yrange = [0d0, maxval(T)]
  
  call LinePlot(x,T,xrange,yrange,xlabel,ylabel,label)
end program main 
