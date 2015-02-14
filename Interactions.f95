module Interactions
  use omp_lib
  implicit none
  private 
  public Force

contains
  subroutine Force(F,EV,r,rc,L)
    real(8), intent(in) :: rc, L, r(:,:)
    real(8), intent(out) :: F(:,:), EV
    real(8) :: d, dr(3) 
    real(8), allocatable :: FMAT(:,:,:)
    integer :: i, j, N

    ! initialize, allocate large array
    N = size(r,1)
    EV = 0d0
    allocate(FMAT(N,N,3))

    !$omp parallel do private(dr, d, i, j)
    do i = 1,N
      j = 1
      do while (j<i)
        dr = r(i,:) - r(j,:) 
        dr = dr - nint(dr/L,kind=8)*L ! implement PBC
        d = sqrt(sum(dr**2d0)) 

        if (d>(sqrt(3d0)*L)) then !check
          print *, "warning: reduce dt", d
          stop
        endif

        if (d<rc) then ! only particle pairs with d<rc contribute to E, F
          FMAT(i,j,:) = 48d0*dr*(1d0/(d**14d0)-0.5d0/(d**8d0))
          EV = EV + 4d0*(1d0/(d**12d0)-1d0/(d**6d0))
        else
          FMAT(i,j,:) = 0d0
        endif

        FMAT(j,i,:) = - FMAT(i,j,:) !use N3
        j = j+1
      enddo
    enddo
    !$omp end parallel do
    F = sum(FMAT,2) !calculate total force vector on particle i
    deallocate(FMAT)
  end subroutine Force
end module 
