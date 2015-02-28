module main_functions
  use constants
  implicit none
  private
  integer, parameter :: m = n_meas
  public :: measure, fold, rescale, pressure, heat_cap, radial_df, &
    diff_const, pot_energy, meas_T, mean_temp 
  
contains
  pure subroutine measure(E, U, r_2, cvv, p, p_0, r, r_0, T, fold_count, L)
    ! calculates energy, Temperature, velocity correlation, diffusion coeff 
    real(dp), intent(out) :: E, cvv, r_2
    real(dp), intent(in) :: U, T, p(:,:), p_0(:,:), r(:,:), r_0(:,:), L
    integer, intent(in) :: fold_count(:,:)
    real(dp) :: Ek, r_tmp(N,3)
    
    r_tmp = unfold(r,fold_count,L)
    Ek = (3._dp/2._dp)*(N-1)*T

    r_2 = sum((r_tmp-r_0)**2)/N 
    E = U + Ek
    cvv = sum(p*p_0)/N
  end subroutine 

  pure function meas_T(p)
    real(dp), intent(in) :: p(:,:)
    real(dp) :: meas_T
    
    meas_T = sum(p**2)/(3._dp*(N-1))
  end function

  pure subroutine rescale(p, T, T_tgt)
    ! rescales momentum to keep temperature fixed, method by berendsen et al.
    real(dp), intent(inout) :: p(:,:)
    real(dp), intent(in) :: T_tgt, T
    real(dp) :: lambda

    lambda = sqrt(1._dp + 0.0025_dp*(T_tgt/T-1._dp))
    p = p*lambda
  end subroutine  

  pure subroutine fold(r, fold_count, L)
    ! enforce periodic BC on positions, track of number of shifts
    real(dp), intent(inout) :: r(:,:)
    real(dp), intent(in) :: L
    integer, intent(inout) :: fold_count(:,:)

    fold_count = fold_count + floor(r/L)
    r = r - floor(r/L)*L
  end subroutine

  pure function unfold(r, fold_count, L)
    real(dp), intent(in) :: r(:,:), L
    integer, intent(in) :: fold_count(:,:)
    real(dp) :: unfold(N,3)
    
    unfold = r + fold_count*L 
  end function

  pure subroutine radial_df(g, bin, rho)     
    ! calculates radial distribution function from binned particle seperations
    real(dp), intent(out) :: g(n_bins)
    integer, intent(in) :: bin(:) 
    real(dp), intent(in) :: rho
    real(dp) :: rs(n_bins), delta_r
    integer :: i, n_s

    n_s = n_meas/up_nbrs_list ! number of timesteps where we binned distances
    delta_r = rm/n_bins ! bin size
    rs = rm*(/(i,i=0,n_bins-1)/)/(n_bins)   

    g = 2._dp/(rho*n_s*(N-1))*real(bin,kind=dp)/(4._dp*pi*delta_r*rs**2)
  end subroutine 
  
  pure subroutine pressure(eq_pres, err_p, virial, mean_T, rho)
    ! calculates equilibrium pressure from virials
    real(dp), intent(out) :: eq_pres, err_p
    real(dp), intent(in) :: virial(:), mean_T, rho
    real(dp) :: c1, c2
    
    ! contribution due to virial (c1), long range correction(c2)
    c1 = 1._dp/(3._dp*N*mean_T)*sum(virial)/m
    c2 = ((16._dp*pi*rho)/mean_T)*(2._dp/(9._dp*(rc**9))-1._dp/(3._dp*(rc**3))) 
    
    eq_pres = 1._dp + c1 + c2
    err_p = 1._dp/(3._dp*N*mean_T)*std_err(virial) ! calc error
  end subroutine 

  pure subroutine heat_cap(Cv, err_Cv, E, T, rescale_T)
    ! calculates head capacity from measured energies
    real(dp), intent(out) :: Cv, err_Cv
    logical, intent(in) :: rescale_T
    real(dp), intent(in) :: E(:), T(:)
    real(dp) :: mu_E, mu_E_2, mu_E_4, mean_T, err_T, mu_T_2, mu_T_4
    
    call mean_temp(mean_T,err_T,T)
      
    if (rescale_T .eqv. .true.) then
      ! calculate heat capacity per particle, NVT ensemble 
      mu_E = sum(E)/m ! 1st moment
      mu_E_2 = sum((E-mu_E)**2)/m ! 2nd moment
      mu_E_4 = sum(block_avg(E-mu_E)**4)/n_blocks**2 

      Cv = 1._dp/(N*mean_T**2)*mu_E_2
      err_Cv = 1._dp/(N*mean_T**2)*sqrt(mu_E_4)
    else 
      ! in NVE ensemble use the lebowitz formula instead
      mu_T_2 = sum((T-mean_T)**2)/m
      mu_T_4 = sum(block_avg(T-mean_T)**4)/n_blocks**2

      Cv = 1._dp/(3._dp/2._dp - mu_T_2/(mean_T**2))
      err_Cv = sqrt(mu_T_4) ! this is incorrect.. 
    endif
  end subroutine

  pure subroutine pot_energy(eq_U, err_U, U)
    real(dp), intent(out) :: eq_U, err_U
    real(dp), intent(in) :: U(:)
    
    eq_U = sum(U/N)/n_meas
    err_U = std_err(U/N)
  end subroutine
  
  pure subroutine diff_const(D, err_D, offset, r_2, t)
    real(dp), intent(out) :: D, err_D, offset
    real(dp), intent(in) :: r_2(:), t(:)
    real(dp) :: slope, mu_t, mu_r_2, ss_rr, ss_tt, ss_rt, s
    ! calculate diffusion constant from slope of <r^2>, using linear regression
    ! see also: http://mathworld.wolfram.com/LeastSquaresFitting.html

    mu_t = sum(t)/m
    mu_r_2 = sum(r_2)/m

    ss_rr = sum((r_2 - mu_r_2)**2)
    ss_tt = sum((t - mu_t)**2)
    ss_rt = sum((t - mu_t)*(r_2 - mu_r_2))
    
    slope = ss_rt/ss_tt
    offset = mu_r_2 - slope*mu_t
    s = sqrt((ss_rr - slope*ss_rt)/(m-2))

    D = slope/6._dp 
    err_D = s/(sqrt(ss_tt)*6._dp) ! error in slope calc
  end subroutine 

  pure subroutine mean_temp(T_mean, err_T, T)
    ! calculate mean temperature and its error during measurement
    real(dp), intent(out) :: T_mean, err_T
    real(dp), intent(in) :: T(:)
    
    T_mean = sum(T)/m
    err_T = std_err(T)
  end subroutine
  
  pure function std_err(A)
    ! calculates std error of A from blocked data
    real(dp), intent(in) :: A(:)
    real(dp) :: std_err, sigma_blocks_2, Avg(n_blocks)

    Avg = block_avg(A)
  
    sigma_blocks_2 = sum((Avg - sum(Avg)/n_blocks)**2)/(n_blocks-1)
    std_err = sqrt(sigma_blocks_2/n_blocks)
  end function

  pure function block_avg(A)
    ! returns array containing block average of A
    real(dp), intent(in) :: A(:)
    real(dp) :: block_avg(n_blocks)
    integer :: j

    do j = 0,(n_blocks-1)
      block_avg(j+1) = sum(A(n_avg*j+1:n_avg*(j+1)))/n_avg
    enddo
  end function
end module
