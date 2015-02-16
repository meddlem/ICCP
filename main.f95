program main
  ! modules to use 
  use plplot ! plotting library 
  use fortplot
  use Plotroutines ! module for plotting particles
  use Inits ! module for initializing model
  use Interactions ! module for calculating the interaction forces
  implicit none

  ! model parameters (constants):
  ! dt = timestep, rc = LJ potential cutoff length, Tinit = inital temp, &
  ! m = mass, rho = number density, units: eps=1, sigma=1, m=1, &
  ! steps = number of timesteps, N = number of particles, in FCC lattice 
  real(8), parameter :: dt = 0.001d0, rc = 2.5d0, Tinit = 1d0, rho = 0.25d0, & 
    eps = 1d0, sigma = 1d0, m = 1d0  
  integer, parameter :: steps = 1000, N = 6**3*4
  logical, parameter :: prtplt = .true.

  ! axis labels in plot:
  character(10), parameter :: xlabel = "time", ylabel = "T", title = "plot", & 
    title1 = ""   
  
  ! declare variables:
  real(8) :: x(steps+1), T(steps+1), E(steps+1), EV, &
    Ek, r(N,3), p(N,3), F(N,3), L = (N/rho)**(1./3.) 
  integer :: i, starttime, endtime
  
  ! variables for final plot
  x = (/(i,i=0, steps)/)
  
  ! initialize r, p, F, and 3d plot:
  call InitR(r,L,N)
  call InitP(p,Tinit,N)
  call Force(F,EV,r,rc,L)
  if(prtplt .eqv. .true.) then
    call ParticlePlotinit(-0.1d0*L,1.1d0*L) 
  endif

  ! calculate temp (T) and total energy (E)
  Ek = 0.5d0*sum(p**2)
  T(1) = 2d0/3d0*Ek/N
  E(1) = EV + Ek

  ! print energy, momentum, temp
  print *, "Psum t=0 : ", sum(p,1) 
  print *, "E t=0 : ", E(1)/N 
  print *, "T t=0 : ", T(1)
  
  call system_clock(starttime)
  ! time integration using the "velocity Verlet" algorithm: 
  do i = 1,steps
    ! plot particle positions
    if((mod(i,3)==0) .and. (prtplt .eqv. .true.)) then 
      call ParticlePlot(r) 
    endif

    r = r + p*dt + 0.5d0*F*(dt**2) !update positions
    r = r - floor(r/L)*L ! enforce PBC
    p = p + 0.5d0*F*dt ! update momentum (1/2)
    call Force(F,EV,r,rc,L) ! update force to next timestep
    p = p + 0.5d0*F*dt ! update momentum (2/2)
    
    ! calculate temp, total energy at timestep i
    Ek = 0.5d0*sum(p**2)
    T(i+1) = 2d0/3d0*Ek/N
    E(i+1) = EV + Ek
  enddo
  
  call system_clock(endtime)
  if(prtplt .eqv. .true.) then
    call plend()
  endif 
  
  print *, "runtime =", (endtime-starttime)/1000, "s"
  print *, "Psum final: ", sum(p,1)
  print *, "E final: ", E(steps+1)/N 
  print *, "T final: ", T(steps+1)
  
  call gnulineplot(x,T,xlabel,ylabel,title1,title)
end program main 
