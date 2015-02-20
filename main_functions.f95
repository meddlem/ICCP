module main_functions
  use constants
  implicit none
  private
  public :: measure, rescale, pressure, heat_cap, radial_df

contains
  pure subroutine measure(E,U,D,T,cvv,p,p_init,r,r_init)
    ! calculates energy, Temperature, velocity correlation, diffusion coeff 
    real(8), intent(in) :: U, p(:,:), p_init(:,:), r(:,:), r_init(:,:)
    real(8), intent(out) :: E, T, cvv, D
    real(8) :: Ek
    
    Ek = 0.5d0*sum(p**2)
    D = sum((r-r_init)**2)/N 
    T = 2d0/3d0*Ek/N
    E = U + Ek
    cvv = sum(p*p_init)/N
  end subroutine 

  pure subroutine rescale(p,T,T_tgt)
    ! rescales momentum to keep temperature fixed, method by berendsen et al.
    real(8), intent(in) :: T_tgt, T
    real(8), intent(inout) :: p(:,:)
    real(8) :: lambda

    lambda = sqrt(1d0 + 0.025d0*(T_tgt/T-1d0))
    p = p*lambda
  end subroutine 

  pure function radial_df(bin,rho)     
    ! calculates radial distribution function from bin
    integer, intent(in) :: bin(:) 
    real(8), intent(in) :: rho
    real(8) :: radial_df(n_bins)
    integer :: i, m
    real(8) :: rs(n_bins), delta_r
    
    m = n_meas/up_nbrs_list ! number of time points where we binned distances
    delta_r = rm/n_bins ! bin size
    rs = rm*(/(i,i=0,n_bins-1)/)/(n_bins)   

    radial_df = 2d0/(rho*m*(N-1))*real(bin,kind=8)/(4d0*pi*delta_r*rs**2)
  end function 
  
  pure real(8) function pressure(virial,T_tgt,rho)
    ! calculates equilibrium pressure from virials
    real(8), intent(in) :: virial(:), T_tgt, rho
    integer :: m, s
    real(8) :: c1, c2

    s = meas_start
    m = n_meas

    ! correction due to virial theorem(c1), long range correction(c2)
    c1 = 1d0/(3d0*N*T_tgt)*sum(virial(s:s+m))/m
    c2 = ((16d0*pi*rho)/T_tgt)*(2d0/(9d0*(rc**9)) - 1d0/(3d0*(rc**3))) 
    
    pressure = 1d0 + c1 + c2
  end function 

  pure real(8) function heat_cap(E,T_tgt)
    real(8), intent(in) :: E(:), T_tgt
    integer :: m, s
    real(8) :: sigma_E_2
    
    s = meas_start
    m = n_meas

    ! calculate heat capacity, check if this is correct wrt number of steps etc
    ! sigma_u_2 = sum((U(s:s+m) - sum(U(s:s+m)/m))**2)/m
    sigma_E_2 = sum((E(s:s+m)-sum(E(s:s+m))/m)**2)/m
    ! dit is maar 1 manier, best voor constante T, niet voor constante E 
    ! heat_cap = (3d0/2d0)*N/(1 - (2d0/3d0)*sigma_u_2/(N*T_tgt**2))
    ! in NVT ensemble
    heat_cap = 1d0/(T_tgt**2)*sigma_E_2
  end function
end module
