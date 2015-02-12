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
        r(:,:) = r(:,:) - floor(r(:,:)/L)*L ! implements PBC 
end subroutine Rinc 

subroutine Vinc(r,v,dt,alpha,rc,sigma,L)
        ! calcs the particles velocity for next timestep 
        use functions ! dont forget about the module bro
        implicit none
        integer :: N, i, j
        real(8), intent(in) :: alpha, sigma, L, dt, r(:,:)
        real(8), intent(inout) :: v(:,:)
        real(8) :: rc, d, dr(3), deltaV(size(v,1),3), F(size(v,1),3)
        !using r at timestep n+1, v at n: compute v at timestep n+1
        
        N = size(r,1)
        ! this can be optimized by cell list method, but also take into account newtons 3rd law: you only need to compute the force
        ! between 2 particles once : F(i,j) = - F(j,i)
        do i=1,N
                
                F(:,:) = 0 ! initialize F for particle j on i 
                
                do j=1,N
                        if (i/=j) then ! no force by particle on itself
                                dr(:) = r(i,:) - r(j,:) 
                                dr(:) = dr(:) - nint(dr(:)/L)*L
                                
                                d = sqrt(sum(dr(:)**2))

                                if ((d<rc) .AND. (d>0.)) then ! only particle pairs with d below cutoff contribute
                                        ! use closest direction vector (considering PBC) for calculating the force
                                        
                                        F(j,:) = alpha*dr(:)*(2*(sigma**12)/(d**14)-(sigma**6)/(d**8))
                                        !calculate magnitude of the force on i by j, LF(d,sigma)
                                endif
                        endif
                enddo  
                
                deltaV(i,:) = sum(F(:,:),1) !calculates change in velocity vector for particle j, based on all
                ! interaction forces within cutoff distance
        enddo
                
       v = v+deltaV*dt ! update all velocity vectors  

end subroutine Vinc

end module
