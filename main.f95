program main
  ! modules to use 
  use constants 
  use inits ! initialize the model
  use main_functions ! additional functions and subroutines
  use plotroutines ! module for plotting particles
  use interactions ! calculating the interaction forces and energies
  use plplot, only : plend ! plotting library 
  implicit none

  ! declare variables: 
  ! r: (N,3) array containing position vectors, r_init = r(t=0)
  ! p: (N,3) array containing momentum vectors, p_init = p(t=0)
  ! F: (N,3) array containing total force vectors
  ! T: array w. temps at every timestep 
  ! E: array w. energies at every timestep 
  ! U: (..) pot. energy at (..)
  ! virial: (..) virial at (..)
  ! cvv: (..) momentum corr. coefficient at (..)
  ! eq_pres: pressure in equilibrium, 
  ! nbrs_list: array containing a list of neighbor pairs
  ! n_nbrs: total number pairs in nbrs_list 
  ! g: radial distribution function
  ! bin: bin containing pair seperations
  ! D: diffusion constant times 6t 
  ! start, end_time: record runtime of simulation

  real(8) :: r(N,3), r_init(N,3), p(N,3), p_init(N,3), F(N,3), T(steps+1), &
    E(steps+1), U(steps+1), virial(steps+1), cvv(steps+1), eq_pres, &
    x_axis(n_bins-1), t_axis(steps+1), g(n_bins), D(steps+1), L, rho, T_init
  integer :: i, start_time, end_time, nbrs_list(N*(N-1)/2,2), n_nbrs, &
    bin(n_bins), tmp_bin(n_bins)
  
  ! get userinput 
  write(*,'(A)',advance='no') "number density = " 
  read(*,*) rho
  write(*,'(A)',advance='no') "target temperature = " 
  read(*,*) T_init

  ! initialize needed vars 
  L = (N/rho)**(1d0/3d0)
  t_axis = dt*(/(i,i=0, steps)/)
  x_axis = rm*(/(i,i=0, n_bins-2)/)/n_bins
  bin = 0
  
  ! initialize the model
  call init_r(r_init,L)
  call init_p(p_init,T_init)
  call make_nbrs_list(nbrs_list,n_nbrs,tmp_bin,r,L)
  call force(F,U(1),virial(1),r_init,rho,L,nbrs_list,n_nbrs)
  
  if(prtplt .eqv. .true.) then
    call particle_plot_init(-0.1d0*L,1.1d0*L) 
  endif
  
  p = p_init 
  r = r_init
  
  ! measure initial temp, velocity correlation, energy
  call measure(E(1),U(1),D(1),T(1),cvv(1),p,p_init,r,r_init)

  ! time integration using the "velocity Verlet" algorithm: 
  print '(A,I5,A)', "starting simulation: ", steps, " iterations"
  call system_clock(start_time)
  do i = 1,steps
    ! plot particle positions
    if((mod(i,5)==0) .and. (prtplt .eqv. .true.)) then 
      call particle_plot(r) 
    endif
  
    ! update list of neighboring particles
    if(mod(i,up_nbrs_list)==0) then
      call make_nbrs_list(nbrs_list,n_nbrs,tmp_bin,r,L)
      if(i>meas_start) then
        bin = bin + tmp_bin
      endif
    endif
    
    r = r + p*dt + 0.5d0*F*(dt**2) !update positions
    r = r - floor(r/L)*L ! enforce PBC on positions
    p = p + 0.5d0*F*dt ! update momentum (1/2)
    call force(F,U(i+1),virial(i+1),r,rho,L,nbrs_list,n_nbrs) ! update force
    p = p + 0.5d0*F*dt ! update momentum (2/2)

    call measure(E(i+1),U(i+1),D(i+1),T(i+1),cvv(i+1),p,p_init,r,r_init)
    call rescale(p,T(i+1),T_init)
  enddo
  
  call system_clock(end_time)
  call plend()
  
  ! further calculations
  U = U/N ! normalize potential energy
  D = D/(6d0*t_axis) ! "normalize" diffusion constant 
  g = radial_df(bin,rho) 
  cvv = cvv/cvv(1) ! normalize velocity correlation
  eq_pres = pressure(virial,T_init,rho)  
  
  ! checks 
  if (maxval(abs(sum(p,1)))>1d-8) then
    print *, "warning: momentum was not conserved"
  endif

  ! processing the results  
  print '(A,I4,A)', " runtime = ", (end_time-start_time)/1000, " s"
  print *, "equilibrium pressure =", eq_pres
  print *, "heat capacity =", heat_cap(E,T_init)
  print *, "T final =", T(steps+1)
  print *, "U equilibrium =", U(steps+1)
  print *, "D =", sum(D(meas_start:meas_start+n_meas))/n_meas 
  
  ! generate final plots
  call gnu_line_plot(t_axis,D,"time","D","","",1)
  call gnu_line_plot(x_axis,g(1:n_bins-1),"r","g","","",2)
end program main 
