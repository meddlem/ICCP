module interactions
  use omp_lib
  implicit none
  private 
  public :: force, make_nbrs_list

contains
  subroutine force(F,U,virial,r,rc,L,nbrs_list,n_nbrs)
    ! calculates net force on each particle, total potential energy, &
    ! and the virial
    real(8), intent(in) :: rc, L, r(:,:)
    real(8), intent(out) :: F(:,:), U, virial
    integer, intent(in) :: nbrs_list(:,:), n_nbrs
    real(8) :: d, dr(3)
    real(8), allocatable :: FMAT(:,:,:), VMAT(:,:), f_dot_dr(:,:)
    integer :: i, j, k, N

    ! initialize, allocate large array
    N = size(r,1)
    allocate(FMAT(N,N,3),VMAT(N,N),f_dot_dr(N,N))
    FMAT(:,:,:) = 0d0
    VMAT(:,:) = 0d0
    f_dot_dr(:,:) = 0d0

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
    
    ! calculate total force vectors, interaction energy
    F = sum(FMAT,2) 
    virial = sum(f_dot_dr) 
    U = sum(VMAT) 
    deallocate(FMAT,VMAT,f_dot_dr)
  end subroutine

  subroutine make_nbrs_list(nbrs_list,n_nbrs,r,rm,L,bin)
    ! creates a list of all particles j within range rm of particle i
    real(8), intent(in) :: r(:,:), L, rm
    integer, intent(out) :: nbrs_list(:,:), n_nbrs
    integer, intent(inout) :: bin(:)
    real(8) :: dr(3), d
    integer :: i, j, k, N, n_bins
    
    N = size(r,1)
    nbrs_list(N*(N-1)/2,2) = 0
    n_bins = size(bin)
    k = 0 

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
