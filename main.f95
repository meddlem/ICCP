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
  real(dp), allocatable :: r(:,:), r_0(:,:), p(:,:), p_0(:,:), F(:,:), &
    T(:), E(:), U(:), virial(:), cvv(:), x_axis(:), t_axis(:), r_squared(:)
  integer, allocatable :: nbrs_list(:,:), fold_count(:,:)
  real(dp) :: eq_pres, err_p, eq_U, err_U, err_heat, heat_c, g(n_bins), L, rho, T_init, &
    D, U_tmp, T_tmp, virial_tmp
  integer :: i, j, start_time, end_time, n_nbrs, bin(n_bins), tmp_bin(n_bins) 
  
  ! allocate large arrays
  allocate(r(N,3),r_0(N,3),p(N,3),p_0(N,3),F(N,3))
  allocate(T(n_meas),E(n_meas),U(n_meas),virial(n_meas),cvv(n_meas))
  allocate(x_axis(n_bins-1),t_axis(n_meas),r_squared(n_meas))
  allocate(nbrs_list(N*(N-1)/2,2),fold_count(N,3))
  
  ! get required userinput 
  write(*,'(A)',advance='no') "number density = " 
  read(*,*) rho
  write(*,'(A)',advance='no') "target temperature = " 
  read(*,*) T_init

  ! initialize needed vars 
  L = (N/rho)**(1._dp/3._dp)
  t_axis = dt*(/(i,i=0, n_meas-1)/)
  x_axis = rm*(/(i,i=0, n_bins-1)/)/n_bins
  bin = 0
  r_0 = r
  p_0 = p
  fold_count = 0

  ! initialize the model
  call init_r(r,L)
  call init_p(p,T_init)
  call make_nbrs_list(nbrs_list,n_nbrs,tmp_bin,r,L)
  call force(F,U_tmp,virial_tmp,r,rho,L,nbrs_list,n_nbrs)
  if(prtplt .eqv. .true.) then
    call particle_plot_init(-0.1_dp*L,1.1_dp*L) 
  endif
  
  ! time integration using the "velocity Verlet" algorithm: 
  print '(A,I5,A)', "starting simulation: ", steps, " iterations"
  call system_clock(start_time)

  do i = 1,steps
    ! plot particle positions
    if ((mod(i,5)==0) .and. (prtplt .eqv. .true.)) then 
      call particle_plot(r) 
    endif
    ! update neighbor list
    if (mod(i,up_nbrs_list) == 0) then
      call make_nbrs_list(nbrs_list,n_nbrs,tmp_bin,r,L)
      if(i>=meas_start) then
        bin = bin + tmp_bin
      endif
    endif
     
    r = r + p*dt + 0.5_dp*F*(dt**2) ! update positions
    call fold(r,fold_count,L) ! enforce PBC
    p = p + 0.5_dp*F*dt ! update momentum (1/2)
    call force(F,U_tmp,virial_tmp,r,rho,L,nbrs_list,n_nbrs) ! update force
    p = p + 0.5_dp*F*dt ! update momentum (2/2)

    if (i == meas_start) then
      r_0 = r
      p_0 = p
      fold_count = 0
    endif

    if (i < meas_start) then
      call measure(E(1),U_tmp,r_squared(1),T_tmp,cvv(1),p,p_0,r,r_0,fold_count,L)
    else
      j = i + 1 - meas_start
      call measure(E(j),U_tmp,r_squared(j),T_tmp,cvv(j),p,p_0,r,r_0,fold_count,L)
      U(j) = U_tmp
      virial(j) = virial_tmp
      T(j) = T_tmp
    endif
    
    if ((rescale_T .eqv. .true.) .or. (i<meas_start)) then
      call rescale(p,T_tmp,T_init)
    endif
  enddo
  
  call system_clock(end_time)
  call plend()
  
  ! further calculations
  D = diff_const(r_squared) ! calculate diffusion constant 
  g = radial_df(bin,rho) 
  cvv = cvv/cvv(1) ! normalize velocity correlation
  
  call heat_cap(heat_c,err_heat,E,T(n_meas))
  call pressure(eq_pres,err_p,virial,T(n_meas),rho)  
  call pot_energy(eq_U,err_U,U)
  
  ! check 
  if (maxval(abs(sum(p,1)))>1d-8) then
    print *, "warning: momentum was not conserved"
  endif

  ! output the results  
  print '(A,I4,A)', " runtime = ", (end_time-start_time)/1000, " s"
  print *, "equilibrium pressure =", eq_pres
  print *, "err p =", err_p
  print *, "heat capacity =", heat_c
  print *, "err hc =", err_heat
  print *, "T final =", T(n_meas)
  print *, "U equilibrium =", eq_U
  print *, "err U =", err_U
  print *, "D =", D 

  print *, "test E: ", blk_var(E/N)
  print *, "test U: ", blk_var(U/N)
  
  ! generate final plots
  call gnu_line_plot(t_axis,virial/N,"time","X","","",1)
  call gnu_line_plot(x_axis,g,"r","g","","",2)
end program main 
