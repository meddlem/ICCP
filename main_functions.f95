module main_functions
  implicit none
  private
  public :: measure, rescale

contains
  subroutine measure(E,Ev,T,p,p_init,Cvv,r)
    ! calculates energy, potential energy, Temperature, velocity correlation
    real(8), intent(in) :: Ev, p(:,:), p_init(:,:), r(:,:)
    real(8), intent(out) :: E, T, Cvv
    real(8) :: Ek
    integer :: N
    
    N = size(r,1)
    Ek = 0.5d0*sum(p**2)
    
    T = 2d0/3d0*Ek/N
    E = EV + Ek
    Cvv = sum(p*p_init)/N
  end subroutine 

  subroutine rescale(p,T,T_init)
    ! rescales momentum to keep temperature fixed, method by berendsen et al.
    real(8), intent(in) :: T_init, T
    real(8), intent(inout) :: p(:,:)
    real(8) :: lambda
    lambda = sqrt(1d0 + 0.025d0*(T_init/T-1d0))
    p = p*lambda
  end subroutine 

end module
