program main
  ! modules to use 
  use plplot ! plotting library 
  use main_functions ! additional functions and subroutines
  use plotroutines ! module for plotting particles
  use inits ! initialize the model
  use interactions ! calculating the interaction forces and energies
  implicit none

  ! model parameters (constants):
  ! dt = timestep, rc = LJ potential cutoff length, T_init = inital temp, &
  ! m = mass, rho = number density, units: eps=1, sigma=1, m=1, kB=1 &
  ! steps = number of timesteps, N = number of particles, in FCC lattice 
  real(8), parameter :: dt = 0.001d0, rc = 2.5d0, T_init = 0.9d0, & 
    rho = 0.55d0,  pi = 4d0*atan(1d0)
  integer, parameter :: steps = 3600, N = 6**3*4
  logical, parameter :: prtplt = .false.
  real(8), parameter :: L = (N/rho)**(1d0/3d0)
  
  ! declare variables: 
  real(8) :: r(N,3), r_init(N,3), p(N,3), p_init(N,3), F(N,3), T(steps+1), &
    E(steps+1), U(steps+1), virial(steps+1), cvv(steps+1), eq_pres, x(steps+1)
  integer :: i, start_time, end_time
  
  ! initialize the model:
  call InitR(r_init,L,N)
  call InitP(p_init,T_init,N)
  call Force(F,U(1),virial(1),r_init,rc,L)
  if(prtplt .eqv. .true.) then
    call ParticlePlotinit(-0.1d0*L,1.1d0*L) 
  endif
  p = p_init 
  r = r_init
  
  ! measure initial temp, velocity correlation, energy
  call measure(E(1),U(1),T(1),p,p_init,cvv(1),r,rc,rho)
  
  ! time integration using the "velocity Verlet" algorithm: 
  print *, "starting simulation: ", steps, "iterations"
  call system_clock(start_time)
  
  do i = 1,steps
    ! plot particle positions
    if((mod(i,3)==0) .and. (prtplt .eqv. .true.)) then 
      call ParticlePlot(r) 
    endif

    r = r + p*dt + 0.5d0*F*(dt**2) !update positions
    r = r - floor(r/L)*L ! enforce PBC on positions
    p = p + 0.5d0*F*dt ! update momentum (1/2)
    call Force(F,U(i+1),virial(i+1),r,rc,L) ! update force to next timestep
    p = p + 0.5d0*F*dt ! update momentum (2/2)

    ! calculate total, potential energy, temperature, correlation (v) 
    call measure(E(i+1),U(i+1),T(i+1),p,p_init,cvv(i+1),r,rc,rho)
    ! rescale the momentum to keep T fixed
    call rescale(p,T(i+1),T_init)
  enddo
  
  call system_clock(end_time)
  call plend()
  
  ! further calculations
  cvv = cvv/cvv(1) ! normalize cvv
  eq_pres = pressure(virial,rc,T_init,rho,N) ! calculate pressure 
  
  ! processing the results  
  if (maxval(abs(sum(p,1)))>1d-8) then
    print *, "warning: momentum not conserved"
  endif 

  print *, "runtime =", (end_time-start_time)/1000, "s"
  print *, "rho =", rho
  print *, "pressure =", eq_pres
  ! print *, "heat capacity =", heat_cap(E,T_init)
  print *, "T final: ", T(steps+1)
  
  ! final plot
  x = dt*(/(i,i=0, steps)/)
  call gnulineplot(x,T,"time","T","","",1)
  call gnulineplot(x,E/N,"time","E","","",2)
end program main 
