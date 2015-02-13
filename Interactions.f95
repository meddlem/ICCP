module Interactions
      implicit none
      private 
      public Force

contains

subroutine Force(F,EV,r,rc,L)
      implicit none

      real(8), intent(in) :: rc, L, r(:,:)
      real(8), intent(out) :: F(:,:), EV
      integer :: N, i, j
      real(8) :: d, dr(3), FMAT(size(r,1),size(r,1),3)

      N = size(r,1)
      F(:,:) = 0d0
      FMAT(:,:,:) = 0d0
      EV = 0d0

      do i = 1,N
        j = 1
        do while (j<i)
                dr = r(i,:) - r(j,:) 
                dr = dr - nint(dr/L,kind=8)*L ! implement PBC
                d = dsqrt(sum(dr**2d0)) 
                
                if (d>(dsqrt(3d0)*L)) then !check
                        print *, "broken, reduce dt", d
                endif

                if (d<rc) then ! only particle pairs with d<rc contribute
                        FMAT(i,j,:) = 48d0*dr*(1d0/(d**14d0)-0.5d0/(d**8d0))
                        EV = EV + LJ(d)
                endif
                
                FMAT(j,i,:) = - FMAT(i,j,:) !use N3
                j = j+1
        enddo
      enddo

      F = sum(FMAT,2) !calculate total force vector on particle i

      contains 
        
      real(8) function LJ(d)
                real(8), intent(in) :: d
                
                LJ = 4d0*(1d0/(d**12d0)-1d0/(d**6d0))
      end function LJ

end subroutine Force

end module 
