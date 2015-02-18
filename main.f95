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
  ! steps = number of timesteps, N = number of particles, in FCC lattice, &
  ! up_nbrs_list = number of iterations between updates of neighbor list
  real(8), parameter :: dt = 0.001d0, rc = 2.5d0, rm = 3.3d0, T_init = 1d0, & 
    rho = 0.88d0,  pi = 4d0*atan(1d0)
  integer, parameter :: steps = 3600, N = 6**3*4, up_nbrs_list = 25
  logical, parameter :: prtplt = .false.
  real(8), parameter :: L = (N/rho)**(1d0/3d0)
  
  ! declare variables: 
  ! r = array containing position vectors, r_init = r(t=0), &
  ! p = array containing momentum vectors, p_init = p(t=0), &
  ! F = array w. total force vectors, T = array w. temps at every timestep, &
  ! E = array w. energies at every timestep, U = (..) pot. energy at (..), &
  ! virial = (..) virial at (..), cvv = (..) momentum corr. coef. at (..), &
  ! eq_pres = pressure in equilibrium, nbrs_list = (..) neighbor pairs, &
  ! n_nbrs = total number of neighbor pairs 
  real(8) :: r(N,3), r_init(N,3), p(N,3), p_init(N,3), F(N,3), T(steps+1), &
    E(steps+1), U(steps+1), virial(steps+1), cvv(steps+1), eq_pres, &
    x(1000), g(1000)
  integer :: i, start_time, end_time, nbrs_list(N*(N-1)/2,2), n_nbrs, &
    bin(1000)
  ! initialize the model:
  call init_r(r_init,L,N)
  call init_p(p_init,T_init,N)
  call make_nbrs_list(nbrs_list,n_nbrs,r,rm,L,bin)
  call force(F,U(1),virial(1),r_init,rc,L,nbrs_list,n_nbrs)
  if(prtplt .eqv. .true.) then
    call ParticlePlotinit(-0.1d0*L,1.1d0*L) 
  endif
  p = p_init 
  r = r_init
  
  ! measure initial temp, velocity correlation, energy
  call measure(E(1),U(1),T(1),p,p_init,cvv(1),r,rc,rho)
  bin(:) = 0 !remember to remove this
  
  ! time integration using the "velocity Verlet" algorithm: 
  print '(A,I5,A)', "starting simulation: ", steps, " iterations"
  call system_clock(start_time)
  do i = 1,steps
    if((mod(i,3)==0) .and. (prtplt .eqv. .true.)) then 
      ! plot particle positions
      call ParticlePlot(r) 
    endif
    if(mod(i,up_nbrs_list)==0) then
      ! update list of neighboring particles
      call make_nbrs_list(nbrs_list,n_nbrs,r,rm,L,bin)
    endif
    
    r = r + p*dt + 0.5d0*F*(dt**2) !update positions
    r = r - floor(r/L)*L ! enforce PBC on positions
    p = p + 0.5d0*F*dt ! update momentum (1/2)
    call force(F,U(i+1),virial(i+1),r,rc,L,nbrs_list,n_nbrs) ! update force
    p = p + 0.5d0*F*dt ! update momentum (2/2)

    ! calculate energy, potential energy, temp, correlation 
    call measure(E(i+1),U(i+1),T(i+1),p,p_init,cvv(i+1),r,rc,rho)
    ! rescale the momentum to keep T fixed
    call rescale(p,T(i+1),T_init)
  enddo
  
  call system_clock(end_time)
  call plend()
  
  ! further calculations
  cvv = cvv/cvv(1) ! normalize velocity correlation
  U = U/N ! normalize potential energy
  eq_pres = pressure(virial,rc,T_init,rho,N) ! calculate pressure 
    ! processing the results  
  if (maxval(abs(sum(p,1)))>1d-8) then
    print *, "warning: momentum not conserved"
  endif 
  print '(A,I4,A)', "runtime = ", (end_time-start_time)/1000, " s"
  print *, "rho =", rho
  print *, "equilibrium pressure =", eq_pres
  print *, "heat capacity =", heat_cap(E,T_init,N)
  print *, "T final: ", T(steps+1)
  print *, "U equilibrium =", U(steps+1)
  call rdf(g,bin,rm)
  ! generate final plots
  x = dt*(/(i,i=0, 999)/)
  ! call gnulineplot(x,T,"time","T","","",1)
  call gnulineplot(x,g,"time","E/N","","",2)
end program main 
