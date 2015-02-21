program main
  ! modules to use 
  use constants 
  use initialize 
  use main_functions 
  use plotroutines 
  use interactions 
  use plplot, only : plend  
  implicit none

  ! declare variables: 
  ! r: (N,3) array containing position vectors, r_init = r(t=0)
  ! p: (N,3) array containing momentum vectors, p_init = p(t=0)
  ! F: (N,3) array containing total force vectors
  ! T: array w. temps at every timestep 
  ! E: array w. energies at every timestep 
  ! U: (..) pot. energy at (..)
  ! virial: (..) virial at (..)
  ! cvv: (..) velocity correlation coefficient at (..)
  ! eq_pres: pressure in equilibrium, 
  ! nbrs_list: array containing a list of neighbor pairs
  ! n_nbrs: total number pairs in nbrs_list 
  ! g: radial distribution function
  ! bin: bin containing pair seperations
  ! D: diffusion constant times 6t 
  ! start, end_time: record runtime of simulation
  real(dp), allocatable :: r(:,:), r_init(:,:), p(:,:), p_init(:,:), F(:,:), &
    T(:), E(:), U(:), virial(:), cvv(:), x_axis(:), t_axis(:), r_squared(:)
  integer, allocatable :: nbrs_list(:,:), fold_count(:,:)
  real(dp) :: eq_pres, g(n_bins), L, rho, T_init, D
  integer :: i, start_time, end_time, n_nbrs, bin(n_bins), tmp_bin(n_bins) 
  
  ! allocate large arrays
  allocate(r(N,3),r_init(N,3),p(N,3),p_init(N,3),F(N,3),fold_count(N,3))
  allocate(T(steps+1),E(steps+1),U(steps+1),virial(steps+1),cvv(steps+1))
  allocate(x_axis(n_bins-1),t_axis(steps+1),r_squared(steps+1))
  allocate(nbrs_list(N*(N-1)/2,2))
  
  ! get required userinput 
  write(*,'(A)',advance='no') "number density = " 
  read(*,*) rho
  write(*,'(A)',advance='no') "target temperature = " 
  read(*,*) T_init

  ! initialize needed vars 
  L = (N/rho)**(1._dp/3._dp)
  t_axis = dt*(/(i,i=0, steps)/)
  x_axis = rm*(/(i,i=0, n_bins-1)/)/n_bins
  bin = 0
  fold_count = 0
  
  ! initialize the model
  call init_r(r_init,L)
  call init_p(p_init,T_init)
  call make_nbrs_list(nbrs_list,n_nbrs,tmp_bin,r,L)
  call force(F,U(1),virial(1),r_init,rho,L,nbrs_list,n_nbrs)
  if(prtplt .eqv. .true.) then
    call particle_plot_init(-0.1_dp*L,1.1_dp*L) 
  endif
  p = p_init 
  r = r_init
  
  ! measure initial temp, velocity correlation, energy
  call measure(E(1),U(1),r_squared(1),T(1),cvv(1),p,p_init,r,r_init,fold_count,L)

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
    
    r = r + p*dt + 0.5_dp*F*(dt**2) ! update positions
    call fold(r,fold_count,L)
    p = p + 0.5_dp*F*dt ! update momentum (1/2)
    call force(F,U(i+1),virial(i+1),r,rho,L,nbrs_list,n_nbrs) ! update force
    p = p + 0.5_dp*F*dt ! update momentum (2/2)
    
    call measure(E(i+1),U(i+1),r_squared(i+1),T(i+1),cvv(i+1),p,p_init,r,&
      r_init,fold_count,L)
    call rescale(p,T(i+1),T_init)
  enddo
  
  call system_clock(end_time)
  call plend()
  
  ! further calculations
  U = U/N ! normalize potential energy
  D = diff_const(r_squared) ! calculate diffusion constant 
  g = radial_df(bin,rho) 
  cvv = cvv/cvv(1) ! normalize velocity correlation
  eq_pres = pressure(virial,T_init,rho)  
  
  ! check 
  if (maxval(abs(sum(p,1)))>1d-8) then
    print *, "warning: momentum was not conserved"
  endif

  ! processing the results  
  print '(A,I4,A)', " runtime = ", (end_time-start_time)/1000, " s"
  print *, "equilibrium pressure =", eq_pres
  print *, "heat capacity =", heat_cap(E,T_init)
  print *, "T final =", T(steps+1)
  print *, "U equilibrium =", U(steps+1)
  print *, "D =", D 
  
  ! generate final plots
  call gnu_line_plot(t_axis,r_squared,"time","x^2","","",1)
  call gnu_line_plot(x_axis,g,"r","g","","",2)
end program main 
