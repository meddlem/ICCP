module Force
      implicit none
      private 
      public ForceMatrix

contains

subroutine ForceMatrix(r,F,m,rc,L)
      implicit none
      real(8), intent(in) :: m, rc, L, r(:,:)
      real(8), intent(out) :: F(:,:)
      integer :: N, i, j
      real(8) :: d, dr(3), FMAT(size(r,1),size(r,1),3)

      N = size(r,1)
      F(:,:) = 0.
      FMAT(:,:,:) = 0.

      do i = 1,N
        j = 1
        do while (j<i)
                dr = r(i,:) - r(j,:) ! probably goes wrong somewhere here 
                dr = dr - nint(dr/L)*L ! implement PBC: shortest possible dx,dy,dz between particles
                d = sqrt(sum(dr**2.)) 
                
                if (d>(3*L**2.)) then !check 
                        print *, "broken, reduce dt", d
                endif

                if (d<rc) then ! only particle pairs with d<rc contribute
                        FMAT(i,j,:) = m*dr*(1./(d**14.)-0.5/(d**8.))
                endif
                
                FMAT(j,i,:) = - FMAT(i,j,:) !use N3
                j = j+1
        enddo
      enddo

      F = sum(FMAT,2) !calculate total force vector on particle i
!
!      contains 
!        
!      real(8) function LJ(d)
!                real(8), intent(in) :: d
!                
!                LJ = 1./(d**12)-1./(d**6)
!      end function LJ

end subroutine ForceMatrix

end module 
