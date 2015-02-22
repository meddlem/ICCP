module main_functions
  use constants
  implicit none
  private
  integer, parameter :: m = n_meas
  public :: measure, fold, rescale, pressure, heat_cap, radial_df, diff_const, &
    blk_var 
  
contains
  pure subroutine measure(E,U,r_squared,T,cvv,p,p_0,r,r_0,fold_count,L)
    ! calculates energy, Temperature, velocity correlation, diffusion coeff 
    real(dp), intent(out) :: E, T, cvv, r_squared
    real(dp), intent(in) :: U, p(:,:), p_0(:,:), r(:,:), r_0(:,:), L
    integer, intent(in) :: fold_count(:,:)
    real(dp) :: Ek, r_tmp(N,3)
    
    r_tmp = unfold(r,fold_count,L)

    Ek = 0.5_dp*sum(p**2)
    r_squared = sum((r_tmp-r_0)**2)/N 
    T = 2._dp/3._dp*Ek/N
    E = U + Ek
    cvv = sum(p*p_0)/N
  end subroutine 
 
  pure subroutine rescale(p,T,T_tgt)
    ! rescales momentum to keep temperature fixed, method by berendsen et al.
    real(dp), intent(inout) :: p(:,:)
    real(dp), intent(in) :: T_tgt, T
    real(dp) :: lambda

    lambda = sqrt(1._dp + 0.025_dp*(T_tgt/T-1._dp))
    p = p*lambda
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

  pure function radial_df(bin,rho)     
    ! calculates radial distribution function from binned particle seperations
    integer, intent(in) :: bin(:) 
    real(dp), intent(in) :: rho
    real(dp) :: radial_df(n_bins), rs(n_bins), delta_r
    integer :: i, n_s

    n_s = n_meas/up_nbrs_list ! number of timesteps where we binned distances
    delta_r = rm/n_bins ! bin size
    rs = rm*(/(i,i=0,n_bins-1)/)/(n_bins)   

    radial_df = 2._dp/(rho*n_s*(N-1))*real(bin,kind=dp)/(4._dp*pi*delta_r*rs**2)
  end function 
  
  pure function pressure(virial,T_tgt,rho)
    ! calculates equilibrium pressure from virials
    real(dp), intent(in) :: virial(:), T_tgt, rho
    real(dp) :: pressure, c1, c2
    
    ! contribution due to virial (c1), long range correction(c2)
    c1 = 1._dp/(3._dp*N*T_tgt)*sum(virial)/m
    c2 = ((16._dp*pi*rho)/T_tgt)*(2._dp/(9._dp*(rc**9)) - 1._dp/(3._dp*(rc**3))) 
    
    pressure = 1._dp + c1 + c2
    ! sigma_p_2 = 

  end function 

  pure function heat_cap(E,T_tgt)
    ! calculates head capacity from measured energies
    real(dp), intent(in) :: E(:), T_tgt
    real(dp) :: heat_cap, sigma_E_2
    
    ! calculate heat capacity, NVT ensemble 
    ! sigma_u_2 = sum((U(s:s+m) - sum(U(s:s+m)/m))**2)/m
    sigma_E_2 = sum((E(1:m)-sum(E(1:m))/m)**2)/m
    ! heat_cap = (3._dp/2._dp)*N/(1 - (2._dp/3._dp)*sigma_u_2/(N*T_tgt**2))
    heat_cap = 1._dp/(T_tgt**2)*sigma_E_2
  end function
  
  pure function diff_const(r_squared)
    real(dp), intent(in) :: r_squared(:)
    real(dp) :: diff_const

    ! calculate diffusion constant from slope of <r^2>
    ! not exactly the proper way to do it though
    diff_const = (r_squared(m) - r_squared(1))/(6._dp*m*dt) 
  end function 
  
  pure function blk_var(X)
    ! calculates the block variance of the input measurement
    real(dp), intent(in) :: X(:)
    real(dp) :: blk_var, Avg(n_blocks)
    integer :: j

    do j = 0,(n_blocks-1)
      Avg(j+1) = sum(X(n_avg*j+1:n_avg*(j+1)))/n_avg
    enddo

    blk_var = sum((Avg - sum(Avg)/n_blocks)**2)/n_blocks
  end function
end module
