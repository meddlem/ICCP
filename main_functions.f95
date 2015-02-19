module main_functions
  implicit none
  private
  public :: measure, rescale, pressure, heat_cap, radial_df
  real(8), parameter :: pi = 4d0*atan(1d0)

contains
  subroutine measure(E,U,D,T,p,p_init,Cvv,r,r_init,rc,rho)
    ! calculates energy, potential energy, Temperature, velocity correlation
    real(8), intent(in) :: p(:,:), p_init(:,:), r(:,:), r_init(:,:), rc, rho
    real(8), intent(out) :: E, T, Cvv, D
    real(8), intent(inout) :: U
    real(8) :: Ek
    integer :: N
    
    N = size(r,1)
    ! apply long range correction to potential energy
    U = U + 8d0*pi*N*rho*(1d0/(9d0*(rc**9))-1d0/(3d0*(rc**3)))
    
    Ek = 0.5d0*sum(p**2)
    D = sum((r-r_init)**2)/N 
    T = 2d0/3d0*Ek/N
    E = U + Ek
    Cvv = sum(p*p_init)/N
  end subroutine 

  subroutine rescale(p,T,T_tgt)
    ! rescales momentum to keep temperature fixed, method by berendsen et al.
    real(8), intent(in) :: T_tgt, T
    real(8), intent(inout) :: p(:,:)
    real(8) :: lambda

    lambda = sqrt(1d0 + 0.025d0*(T_tgt/T-1d0))
    p = p*lambda
  end subroutine 

  real(8) function pressure(virial,rc,T_tgt,rho,N)
    real(8), intent(in) :: virial(:), rc, T_tgt, rho
    integer, intent(in) :: N
    real(8) ::  c1, c2
    integer :: M, steps

    steps = size(virial,1) - 1
    M = steps/2 

    ! correction due to virial theorem(c1), long range correction(c2)
    c1 = 1d0/(3d0*N*T_tgt)*sum(virial(steps-M:steps+1))/(M + 1)
    c2 = ((16d0*pi*rho)/T_tgt)*(2d0/(9d0*(rc**9)) - 1d0/(3d0*(rc**3))) 
    
    pressure = 1d0 + c1 + c2
  end function 

  real(8) function heat_cap(U,T_tgt,N)
    real(8), intent(in) :: U(:), T_tgt
    integer, intent(in) :: N
    real(8) :: sigma_u_2
    integer :: M, steps
    steps = size(U)
    M = steps/2 
    ! calculate heat capacity, check if this is correct wrt number of steps etc
    ! heat_cap = sum((E(M:steps+1)-sum(E(M:steps+1))/steps)**2)/(M*(T_tgt**2))
    sigma_u_2 = sum((U(M:steps+1) - sum(U(M:steps+1)/M))**2)/M
    ! dit is maar 1 manier, best voor constante T, niet voor constante E 
    heat_cap = (3d0/2d0)*N/(1 - (2d0/3d0)*sigma_u_2/(N*T_tgt**2))
  end function

  subroutine radial_df(g,bin,n_bins,steps,rho,rm)     
    ! calculates radial distribution function from bin
    integer, intent(in) :: bin(:), n_bins, steps
    real(8), intent(in) :: rho, rm
    real(8), intent(out) :: g(n_bins)
    integer :: i
    real(8) :: rs(n_bins)
    
    rs = rm*(/(i,i=0,n_bins-1)/)/(n_bins)   
    g = 1d0/(steps/2)*real(bin,kind=8)/(4d0*pi*rho*rs**2)
  end subroutine 
end module
