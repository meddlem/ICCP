module interactions
  use omp_lib
  implicit none
  private 
  public Force

contains
  subroutine Force(F,EV,virial,r,rc,L)
    ! calculates net force on each particle, total potential energy, &
    ! and the virial coefficient
    real(8), intent(in) :: rc, L, r(:,:)
    real(8), intent(out) :: F(:,:), EV, virial
    real(8) :: d, dr(3)
    real(8), allocatable :: FMAT(:,:,:), VMAT(:,:), f_dot_dr(:,:)
    integer :: i, j, N

    ! initialize, allocate large array
    N = size(r,1)
    allocate(FMAT(N,N,3),VMAT(N,N),f_dot_dr(N,N))
    FMAT(:,:,:) = 0d0
    VMAT(:,:) = 0d0
    f_dot_dr(:,:) = 0d0

    !$omp parallel do private(d,dr)
    do i = 1,N
      j = 1
      do while (j<i)
        dr = r(i,:) - r(j,:) 
        dr = dr - nint(dr/L,kind=8)*L ! implement PBC
        d = sqrt(sum(dr**2)) 

        if (d>(sqrt(3d0)*L)) then !check
          print *, "warning: reduce dt", d
          stop
        endif

        if (d<rc) then ! only particle pairs with d<rc contribute to E, F
          FMAT(i,j,:) = 48d0*dr*(1d0/(d**14)-0.5d0/(d**8))
          f_dot_dr(i,j) = 48d0*(1d0/(d**12)-0.5d0/(d**6))
          VMAT(i,j) = 4d0*(1d0/(d**12)-1d0/(d**6))
        endif

        FMAT(j,i,:) = - FMAT(i,j,:) ! use Newton3
        j = j+1
      enddo
    enddo
    !$omp end parallel do
    
    ! calculate total force vectors, interaction energy
    F = sum(FMAT,2) 
    virial = sum(f_dot_dr) 
    EV = sum(VMAT) 
    deallocate(FMAT,VMAT)
  end subroutine Force
end module 
