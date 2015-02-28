program main
  ! modules to use 
  use constants 
  use io
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
  ! x_axis, t_axis: used for plot and correlation function
  ! r_2: mean squared displacement from start of measurement
  ! nbrs_list: array containing a list of neighbor pairs
  ! fold_count: array with number of coord shifts for each particle
  ! eq_pres: pressure in equilibrium 
  ! eq_U pot energy per particle during measurement (avg)
  ! Cv: specific heat per particle, at constant volume
  ! g: radial distribution function
  ! L: side of volume 
  ! rho: number density
  ! D: diffusion constant 
  ! T_init: initial temp set by user
  ! n_nbrs: total number pairs in nbrs_list 
  ! bin: bin containing pair seperations
  ! start, end_time, runtime: record length of simulation

  real(dp), allocatable :: r(:,:), r_0(:,:), p(:,:), p_0(:,:), F(:,:), T(:), &
    E(:), U(:), virial(:), cvv(:), x_axis(:), t_axis(:), r_2(:)
  integer, allocatable :: nbrs_list(:,:), fold_count(:,:)
  real(dp) :: eq_pres, err_p, eq_U, err_U, err_Cv, Cv, g(n_bins), L, &
    rho, D, err_D, T_init, T_tmp, mean_T, err_T, r_2_fit(n_meas), offset
  integer :: i, j, start_time, end_time, runtime, n_nbrs, bin(n_bins), &
    tmp_bin(n_bins) 
  
  ! allocate large arrays
  allocate(r(N,3), r_0(N,3), p(N,3), p_0(N,3) ,F(N,3) ,fold_count(N,3), &
    T(n_meas), E(n_meas), U(n_meas), virial(n_meas), cvv(n_meas), &
    t_axis(n_meas), r_2(n_meas), x_axis(n_bins-1), nbrs_list(N*(N-1)/2,2))
  
  ! get required userinput 
  call user_in(rho,T_init)  

  ! initialize a few needed variables 
  L = (N/rho)**(1._dp/3._dp)
  t_axis = dt*(/(i,i=0, n_meas-1)/)
  x_axis = rm*(/(i,i=0, n_bins-1)/)/n_bins
  fold_count = 0
  bin = 0
  j = 1

  ! initialize the model
  call init_r(r,L)
  call init_p(p,T_init)
  call make_nbrs_list(nbrs_list,n_nbrs,tmp_bin,r,L)
  call force(F,U(j),virial(j),r,rho,L,nbrs_list,n_nbrs)
  if(prtplt .eqv. .true.) call particle_plot_init(-0.1_dp*L,1.1_dp*L)
  
  call system_clock(start_time) 
  ! central part of simulation
  do i = 1,steps
    ! plot particle positions
    if ((mod(i,5) == 0) .and. (prtplt .eqv. .true.)) call particle_plot(r) 
    ! update neighbor list
    if (mod(i,up_nbrs_list) == 0) then
      call make_nbrs_list(nbrs_list,n_nbrs,tmp_bin,r,L)
      if (i>=m_start) bin = bin + tmp_bin
    endif
     
    ! time integration using the "velocity Verlet" algorithm: 
    r = r + p*dt + 0.5_dp*F*(dt**2) ! update positions
    call fold(r,fold_count,L) ! enforce PBC
    p = p + 0.5_dp*F*dt ! update momentum (1/2)
    call force(F,U(j),virial(j),r,rho,L,nbrs_list,n_nbrs) ! update force
    p = p + 0.5_dp*F*dt ! update momentum (2/2)

    T_tmp = meas_T(p)

    if (i >= m_start) then
      if (i == m_start) then
        p_0 = p; r_0 = r ! initialize measurement vars
        fold_count = 0 ! reset
      endif
      
      j = i + 1 - m_start
      call measure(E(j),U(j),r_2(j),cvv(j),p,p_0,r,r_0,T_tmp,fold_count,L)
      T(j) = T_tmp
    endif
    
    if ((rescale_T .eqv. .true.).or.(i < m_start)) call rescale(p,T_tmp,T_init)
  enddo
  
  call system_clock(end_time)
  call plend()
  
  ! further calculations
  cvv = cvv/cvv(1) ! normalize velocity correlation
  runtime = (end_time - start_time)/1000
  call radial_df(g,bin,rho) 
  call mean_temp(mean_T,err_T,T)
  call diff_const(D,err_D,offset,r_2,t_axis) ! calc D using linear regression fit
  call heat_cap(Cv,err_Cv,E,T,rescale_T)
  call pressure(eq_pres,err_p,virial,mean_T,rho)  
  call pot_energy(eq_U,err_U,U)
  
  ! check 
  call f_check(p)

  ! print measurement results  
  call results_out(runtime, rho, T_init, eq_pres, err_p, Cv, err_Cv, mean_T, &
    err_T, eq_U, err_U, D, err_D)
  
  ! generate final plots
  r_2_fit = 6._dp*D*dt*(/(i,i=0, n_meas-1)/) + offset
  call gnu_line_plot(t_axis,r_2,"t","<r^2>","measurement","",1,r_2_fit,"fit")
  call gnu_line_plot(x_axis,g,"r","g(r)","","",2)
  call gnu_line_plot(t_axis,cvv,"t","Cvv","","",3)
end program main 
