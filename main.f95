program main
  ! modules to use 
  use plplot ! plotting library 
  use main_functions
  use fortplot
  use Plotroutines ! module for plotting particles
  use Inits ! module for initializing model
  use Interactions ! module for calculating the interaction forces
  implicit none

  ! model parameters (constants):
  ! dt = timestep, rc = LJ potential cutoff length, T_init = inital temp, &
  ! m = mass, rho = number density, units: eps=1, sigma=1, m=1, &
  ! steps = number of timesteps, N = number of particles, in FCC lattice 
  real(8), parameter :: dt = 0.001d0, rc = 2.5d0, T_init = 1.55d0, rho = 0.45d0 
  integer, parameter :: steps = 3600, N = 6**3*4
  logical, parameter :: prtplt = .false.

  ! axis labels in plot:
  character(10), parameter :: xlabel = "time", ylabel = "T", title = "plot", & 
    title1 = ""   
  
  ! declare variables: (tidy up here)
  real(8) :: x(steps+1), T(steps+1), E(steps+1), EV(steps+1), virial(steps+1), &
    p_init(N,3), cvv(steps+1), r(N,3), p(N,3), F(N,3), L = (N/rho)**(1./3.), &
    pressure 
  integer :: i, start_time, end_time
  
  ! initialize r, p, F, and 3d plot:
  call InitR(r,L,N)
  call InitP(p,T_init,N)
  call Force(F,EV(1),virial(1),r,rc,L)
  if(prtplt .eqv. .true.) then
    call ParticlePlotinit(-0.1d0*L,1.1d0*L) 
  endif

  ! calculate initial temp, velocity correlation, energy
  call measure(E(1),Ev(1),T(1),T_init,p,p_init,cvv(1),r)
  ! p_init = p
  ! cvv(1) = sum(p_init**2)/N
  ! Ek = 0.5d0*sum(p**2)
  ! T(1) = 2d0/3d0*Ek/N
  ! E(1) = EV + Ek

  ! print energy, momentum, temp
  print *, "T t=0 : ", T(1)
  
  call system_clock(start_time)
  ! time integration using the "velocity Verlet" algorithm: 
  do i = 1,steps
    ! plot particle positions
    if((mod(i,3)==0) .and. (prtplt .eqv. .true.)) then 
      call ParticlePlot(r) 
    endif

    r = r + p*dt + 0.5d0*F*(dt**2) !update positions
    r = r - floor(r/L)*L ! enforce PBC
    p = p + 0.5d0*F*dt ! update momentum (1/2)
    call Force(F,EV(i+1),virial(i+1),r,rc,L) ! update force to next timestep
    p = p + 0.5d0*F*dt ! update momentum (2/2)
    
    call measure(E(i+1),Ev(i+1),T(i+1),T_init,p,p_init,cvv(i+1),r)
    call rescale(p,T(i+1),T_init)
    ! calculate temp, total energy at timestep i
!    Ek = 0.5d0*sum(p**2)
!    T(i+1) = 2d0/3d0*Ek/N
!    E(i+1) = EV + Ek
!    cvv(i+1) = sum(p*p_init)/N
    ! rescaling p: 
!    lambda = sqrt(1 + 0.025*(T_init/T(i+1)-1d0))
!    p = p*lambda
  enddo
  
  call system_clock(end_time)
  if(prtplt .eqv. .true.) then
    call plend()
  endif 
  
  ! some extra calculations
  cvv = cvv/cvv(1)
  pressure = 1d0 + 1d0/(3*N*T_init)*sum(virial(steps-1500:steps+1))/1501d0
  
  ! processing the results  
  if (maxval(abs(sum(p,1)))>1d-8) then
    print *, "warning: momentum not conserved"
  endif 

  print *, "runtime =", (end_time-start_time)/1000, "s"
  print *, "pressure =", pressure
  print *, "T final: ", T(steps+1)
  
  ! final plot
  x = dt*(/(i,i=0, steps)/)
  
  call gnulineplot(x,T,xlabel,ylabel,title1,title)
end program main 
