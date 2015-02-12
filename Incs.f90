module Incs
      implicit none
      private
      public :: Vinc, Rinc
contains

subroutine Rinc(r,v,dt,L)
        ! this subroutine calculates particle positions
        implicit none
        real(8), intent(in) :: dt, L, v(:,:)
        real(8), intent(inout) :: r(:,:)
        
        r = r + v*dt !computes positions at next timestep  
        !if (abs(r(i,1))>L) then
                
        ! implement the BC here
end subroutine Rinc 

subroutine Vinc(r,v,dt,alpha,rc,sigma)
        ! calcs the particles velocity for next timestep 
        use functions ! dont forget about the module bro
        implicit none
        integer :: N, i, j, k
        real(8), intent(in) :: alpha, sigma, dt, r(:,:)
        real(8), intent(inout) :: v(:,:)
        real(8) :: rc, d, F(size(v,1)), deltaV(size(v,1),size(v,2)), dir(size(v,1),size(v,2))
        !using r at timestep n+1, v at n: compute v at timestep n+1
        
        N = size(r,1)
        ! this can be optimized by cell list method, but also take into account newtons 3rd law: you only need to compute the force
        ! between 2 particles once : F(i,j) = - F(j,i)
        do i=1,N
                F(:) = 0 ! initialize F for particle i
                dir(:,:) = 0 ! initialize dir for particle i 
                
                do j=1,N
                        ! implement BC here 

                        if (i/=j) then ! no force by particle on itself
                                
                                d = sqrt(sum((r(i,:)-r(j,:))**2))      
                                
                                if (d<rc) then ! only particle pairs with d below cutoff contribute
                                        F(j) = LJ(d,sigma) !calculate magnitude of the force on i by j, LF(d,sigma)
                                        dir(j,:) = (r(i,:)-r(j,:))/d !calculate direction of the (central) force
                                end if
                        end if
                enddo  
                
                do k = 1,3
                        deltaV(i,k) = alpha*sum(F(:)*dir(:,k)) !calculates change in velocity vector for particle j, based on all
                        ! interaction forces within cutoff distance
                enddo 
        enddo
                
        v = v + deltaV*dt ! update all velocity vectors  

end subroutine Vinc

end module
