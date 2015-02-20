module interactions
  use constants
  use omp_lib
  implicit none
  private 
  public :: force, make_nbrs_list

contains
  subroutine force(F,U,virial,r,rho,L,nbrs_list,n_nbrs)
    ! calculates net force on each particle, total potential energy, &
    ! and the virial
    real(8), intent(out) :: F(:,:), U, virial
    real(8), intent(in) :: L, r(:,:), rho
    integer, intent(in) :: nbrs_list(:,:), n_nbrs
    real(8), allocatable :: FMAT(:,:,:), VMAT(:,:), f_dot_dr(:,:)
    real(8) :: d, dr(3)
    integer :: i, j, k

    ! initialize, allocate large array
    allocate(FMAT(N,N,3),VMAT(N,N),f_dot_dr(N,N))
    FMAT = 0d0
    VMAT = 0d0
    f_dot_dr = 0d0

    !$omp parallel do private(d,dr)
    do k = 1,n_nbrs
      ! get all neighboring particles from the list
      i = nbrs_list(k,1)
      j = nbrs_list(k,2)
      
      dr = r(i,:) - r(j,:) 
      dr = dr - nint(dr/L,kind=8)*L ! implement PBC
      d = sqrt(sum(dr**2)) 

      if (d<rc) then ! only particle pairs with d<rc contribute to E, F
        FMAT(i,j,:) = 48d0*dr*(1d0/(d**14)-0.5d0/(d**8))
        f_dot_dr(i,j) = 48d0*(1d0/(d**12)-0.5d0/(d**6))
        VMAT(i,j) = 4d0*(1d0/(d**12)-1d0/(d**6))
      endif

      FMAT(j,i,:) = - FMAT(i,j,:) !use newton 3
    enddo
    !$omp end parallel do
    
    ! calculate total force vectors, pot. energy, virial
    F = sum(FMAT,2) 
    virial = sum(f_dot_dr) 
    U = sum(VMAT) 
    ! apply long range correction to U
    U = U + 8d0*pi*N*rho*(1d0/(9d0*(rc**9)) - 1d0/(3d0*(rc**3))) 
    deallocate(FMAT,VMAT,f_dot_dr)
  end subroutine

  subroutine make_nbrs_list(nbrs_list,n_nbrs,bin,r,L)
    ! creates a list of all particles j within distance rm of particle i
    integer, intent(out) :: nbrs_list(:,:), n_nbrs, bin(:)
    real(8), intent(in) :: r(:,:), L
    real(8) :: dr(3), d
    integer :: i, j, k
    
    ! initialize variables
    nbrs_list = 0
    k = 0
    bin = 0 

    do i = 1,N
      do j = i+1,N
        dr = r(i,:) - r(j,:)
        dr = dr - nint(dr/L,kind=8)*L
        d = sqrt(sum(dr**2))
        
        if (d>(sqrt(3d0)*L)) then !check
          print *, "warning: reduce dt", d
          stop
        elseif (d<rm) then
          k = k+1
          nbrs_list(k,:) = [i, j]
          bin(int((n_bins-1)*d/rm)+1) = bin(int((n_bins-1)*d/rm)+1) + 1
        endif
      enddo
    enddo

    n_nbrs = k !get total # of nn pairs
  end subroutine
end module 
