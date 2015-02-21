module main_functions
  use constants
  implicit none
  private
  public :: measure, fold, rescale, pressure, heat_cap, radial_df, diff_const

contains
  pure subroutine measure(E,U,r_squared,T,cvv,p,p_init,r,r_init,fold_count,L)
    ! calculates energy, Temperature, velocity correlation, diffusion coeff 
    real(dp), intent(out) :: E, T, cvv, r_squared
    real(dp), intent(in) :: U, p(:,:), p_init(:,:), r(:,:), r_init(:,:), L
    integer, intent(in) :: fold_count(:,:)
    real(dp) :: Ek, r_tmp(N,3)
    
    r_tmp = unfold(r,fold_count,L)

    Ek = 0.5_dp*sum(p**2)
    r_squared = sum((r_tmp-r_init)**2)/N 
    T = 2._dp/3._dp*Ek/N
    E = U + Ek
    cvv = sum(p*p_init)/N
  end subroutine 
  
  pure subroutine fold(r,fold_count,L)
    ! enforce periodic BC on positions, track of number of shifts
    real(dp), intent(inout) :: r(:,:)
    real(dp), intent(in) :: L
    integer, intent(inout) :: fold_count(:,:)

    fold_count = fold_count + floor(r/L)
    r = r -floor(r/L)*L
  end subroutine

  pure function unfold(r,fold_count,L)
    real(dp), intent(in) :: r(:,:), L
    integer, intent(in) :: fold_count(:,:)
    real(dp) :: unfold(N,3)
    
    unfold = r + fold_count*L 
  end function

  pure subroutine rescale(p,T,T_tgt)
    ! rescales momentum to keep temperature fixed, method by berendsen et al.
    real(dp), intent(inout) :: p(:,:)
    real(dp), intent(in) :: T_tgt, T
    real(dp) :: lambda

    lambda = sqrt(1._dp + 0.025_dp*(T_tgt/T-1._dp))
    p = p*lambda
  end subroutine 

  pure function radial_df(bin,rho)     
    ! calculates radial distribution function from binned particle seperations
    integer, intent(in) :: bin(:) 
    real(dp), intent(in) :: rho
    real(dp) :: radial_df(n_bins), rs(n_bins), delta_r
    integer :: i, m

    m = n_meas/up_nbrs_list ! number of timesteps where we binned distances
    delta_r = rm/n_bins ! bin size
    rs = rm*(/(i,i=0,n_bins-1)/)/(n_bins)   

    radial_df = 2._dp/(rho*m*(N-1))*real(bin,kind=dp)/(4._dp*pi*delta_r*rs**2)
  end function 
  
  pure function pressure(virial,T_tgt,rho)
    ! calculates equilibrium pressure from virials
    real(dp), intent(in) :: virial(:), T_tgt, rho
    real(dp) :: pressure, c1, c2
    integer :: m, s
    
    s = meas_start
    m = n_meas

    ! correction due to virial theorem(c1), long range correction(c2)
    c1 = 1._dp/(3._dp*N*T_tgt)*sum(virial(s:s+m))/m
    c2 = ((16._dp*pi*rho)/T_tgt)*(2._dp/(9._dp*(rc**9)) - 1._dp/(3._dp*(rc**3))) 
    
    pressure = 1._dp + c1 + c2
  end function 

  pure function heat_cap(E,T_tgt)
    ! calculates head capacity from measured energies
    real(dp), intent(in) :: E(:), T_tgt
    real(dp) :: heat_cap, sigma_E_2
    integer :: m, s
    
    s = meas_start
    m = n_meas

    ! calculate heat capacity, NVT ensemble 
    ! sigma_u_2 = sum((U(s:s+m) - sum(U(s:s+m)/m))**2)/m
    sigma_E_2 = sum((E(s:s+m)-sum(E(s:s+m))/m)**2)/m
    ! heat_cap = (3._dp/2._dp)*N/(1 - (2._dp/3._dp)*sigma_u_2/(N*T_tgt**2))
    heat_cap = 1._dp/(T_tgt**2)*sigma_E_2
  end function
  
  pure function diff_const(r_squared)
    real(dp), intent(in) :: r_squared(:)
    real(dp) :: diff_const
    integer :: m
    
    m = n_meas
    
    diff_const = r_squared(3)/(6._dp*1) !strictly incorrect, you need te slope
  end function 
end module
