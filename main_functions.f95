module main_functions
  implicit none
  private
  public :: measure, rescale, pressure, heat_cap
  real(8), parameter :: pi = 4d0*atan(1d0)

contains
  subroutine measure(E,U,T,p,p_init,Cvv,r,rc,rho)
    ! calculates energy, potential energy, Temperature, velocity correlation
    real(8), intent(in) :: p(:,:), p_init(:,:), r(:,:), rc, rho
    real(8), intent(out) :: E, T, Cvv
    real(8), intent(inout) :: U
    real(8) :: Ek
    integer :: N
    
    N = size(r,1)
    Ek = 0.5d0*sum(p**2)
    ! apply long range correction to potential energy
    U = U + 8d0*pi*rho*(1d0/(9d0*(rc**9))-1d0/(3d0*(rc**3)))

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
    c1 = 1d0/(3d0*N*T_tgt)*sum(virial(steps-M:steps+1))/(M+1)
    c2 = ((16d0*pi*rho)/T_tgt)*(2d0/(9d0*(rc**9)) - 1d0/(3d0*(rc**3))) 
    
    pressure = 1d0 + c1 + c2
  end function 

  real(8) function heat_cap(E,T_tgt)
    real(8), intent(in) :: E(:), T_tgt
    integer :: steps
    steps = size(E)
    
    heat_cap = sum(E-sum(E)/steps)/(steps*(T_tgt**2))
  end function
end module
