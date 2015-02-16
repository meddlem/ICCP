program main
  ! modules to use 
  use constants ! contains model constants
  use plplot ! plotting library 
  use main_functions ! additional functions and subroutines
  use fortplot ! interface to gnuplot
  use Plotroutines ! module for plotting particles
  use Inits ! initialize the model
  use Interactions ! calculating the interaction forces and energies
  implicit none

!  ! model parameters (constants):
!  ! dt = timestep, rc = LJ potential cutoff length, T_init = inital temp, &
!  ! m = mass, rho = number density, units: eps=1, sigma=1, m=1, &
!  ! steps = number of timesteps, N = number of particles, in FCC lattice 
!  real(8), parameter :: dt = 0.001d0, rc = 2.5d0, T_init = 1.585d0, & 
!    rho = 0.65d0,  pi = 4d0*atan(1d0)
!  integer, parameter :: steps = 3600, N = 6**3*4
!  logical, parameter :: prtplt = .false.
!
!  ! axis labels in final plot:
!  character(10), parameter :: xlabel = "time", ylabel = "T", title = "plot", & 
!    title1 = ""   
!  
  ! declare variables: (tidy up here)
  real(8) :: r(N,3), p(N,3), p_init(N,3), F(N,3), T(steps+1), E(steps+1), &
    EV(steps+1), virial(steps+1), cvv(steps+1), eq_pres, x(steps+1)
  integer :: i, start_time, end_time
  
  ! initialize r, p, F, and 3d plot:
  call InitR(r,L,N)
  call InitP(p,T_init,N)
  call Force(F,EV(1),virial(1),r,rc,L)
  if(prtplt .eqv. .true.) then
    call ParticlePlotinit(-0.1d0*L,1.1d0*L) 
  endif

  ! calculate initial temp, velocity correlation, energy
  call measure(E(1),Ev(1),T(1),p,p_init,cvv(1),r)
  
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

    ! calculate total, potential energy, temperature, correlation (v) 
    call measure(E(i+1),Ev(i+1),T(i+1),p,p_init,cvv(i+1),r)
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
  print *, "pressure =", eq_pres
  print *, "T final: ", T(steps+1)
  
  ! final plot
  x = dt*(/(i,i=0, steps)/)
  call gnulineplot(x,T,xlabel,ylabel,title1,title)
end program main 
