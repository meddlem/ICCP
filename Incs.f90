module Incs
      implicit none
      private
      public :: Vinc, Rinc
contains

subroutine Rinc(r,v,dt)
        ! this subroutine calculates particle positions
        implicit none
        real(8), intent(in) :: dt, v(:,:)
        real(8), intent(inout) :: r(:,:)

        r = r + v*dt !computes positions at next timestep  
end subroutine Rinc 

subroutine Vinc(r,v,dt,m,rc)
        ! calcs the particles velocity for next timestep 
        use LJforce ! dont forget about the module bro
        implicit none
        integer :: N, i, j, k
        real(8), intent(in) :: m, dt, r(:,:)
        real(8), intent(inout) :: v(:,:)
        real(8) :: rc, d, sigma, F(size(v,1)), deltaV(size(v,1),size(v,2)), dir(size(v,1),size(v,2))
        !using r at timestep n+1, v at n: compute v at timestep n+1
        sigma = 0.1 !move this 
        N = size(r,1)
        do i=1,N
                do j=1,N
                        if (i/=j) then ! no force by particle on itself
                                
                                d = sqrt(sum((r(i,:)-r(j,:))**2))      
                                if (d<rc) then
                                        !only particle pairs with d below cutoff contribute
                                        F(j) = LJ(d,sigma) !calculate magnitude of the force on i by j, LF(d,sigma)
                                        dir(j,:) =  (r(i,:)-r(j,:))/d !calculate direction of the (central) force
                                end if

                        else
                                F(j) = 0
                        end if

                enddo ! we got dat..  
                
                do k = 1,3
                        deltaV(i,k) = 1/m*sum(F(:)*dir(:,k)) !calculates change in velocity vector for particle j, based on all
                        ! interaction forces within cutoff distance
                enddo 
        enddo
                
        v = v + deltaV*dt 

end subroutine Vinc

end module


