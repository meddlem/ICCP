module Incs
      implicit none
      private
      public :: Vinc, Rinc

contains
        ! to do: create better time integration algorithm, see verlets paper
subroutine Rinc(r,v,dt,L)
        ! this subroutine calculates particle positions
        implicit none
        real(8), intent(in) :: dt, L, v(:,:)
        real(8), intent(inout) :: r(:,:)
        
        r = r+v*dt !computes positions at next timestep  
        r = r-floor(r/L)*L ! implements PBC for position vectors

end subroutine Rinc 

subroutine Vinc(r,v,dt,alpha,rc,sigma,L)
        ! calcs the particles velocity for next timestep 
        implicit none
        integer :: N, i, j
        real(8), intent(in) :: alpha, sigma, L, dt, r(:,:)
        real(8), intent(inout) :: v(:,:)
        real(8) :: rc, d, dr(3), F(size(v,1),size(v,1),3) 
        !using r at timestep n+1, v at n: compute v at timestep n+1
        
        N = size(r,1)
        F(:,:,:) = 0 ! initialize F for particle j on i

        ! the following can be optimized by cell list method
        do i=1,N
                j = 1
                do while (j<i)
                        
                        dr(:) = r(i,:) - r(j,:) 
                        dr(:) = dr(:) - nint(dr(:)/L)*L ! implement PBC: shortest possible dx,dy,dz between particles
                        d = sqrt(sum(dr(:)**2)) 

                        if ((d<rc) .AND. (d>0.)) then ! only particle pairs with d below cutoff contribute
                                ! use closest direction vector (considering PBC) for calculating the force of j on i
                                
                                F(i,j,:) = alpha*dr(:)*(2*(sigma**12)/(d**14)-(sigma**6)/(d**8))
                        endif
                        
                        F(j,i,:) = - F(i,j,:) !use N3
                        j = j+1

                enddo
       enddo

       v = v+sum(F,2)*dt ! update all velocity vectors, using force of particles j/=i on particle i  

end subroutine Vinc

end module
