program main
      use maxwell 
      use LJforce
      implicit none
      ! model parameters, constants
      real, parameter :: dt = 1, a = 1, rc = 5 ! defines: timestep, lattice constant (initial), cutoff length
      integer, parameter :: N=4000 !number of particles, must be multiple of 4 for now 
      real :: r(N,3), v(N,3) ! initialize r (position) and v (velocity) vectors
      real :: beta = 1.0, m = 1.0 ! some parameters controlling temp, particle mass, in LJ potential we only need sigma eps can be
      !absorbed into m
      
       ! best to start bins at t=0 and then update every iteration? 

      call InitCell(r,a,N)
      call InitVel(v,m,beta,N)
      print *, r(4,2) 
      call Rinc(r,v,dt,N)
      print *, r(4,2)
      call Vinc(r,v,dt,m,rc,N)
      call Rinc(r,v,dt,N)
      print *, r(6,2)

end program main 

subroutine Rinc(r,v,dt,N)
        implicit none
        integer :: N 
        real, intent(in) :: dt, v(N,3)
        real, intent(inout) :: r(N,3)

        r = r + v*dt !computes positions at next timestep  
end subroutine Rinc 

subroutine Vinc(r,v,dt,m,rc,N)
        use LJforce ! dont forget about the module bro
        implicit none
        integer :: N, i, j, k
        real, intent(in) :: m, dt, r(N,3)
        real, intent(inout) :: v(N,3)
        real :: rc,d, F(N), deltaV(N,3), dir(N,3)
        !using r at timestep n+1, v at n: compute v at timestep n+1
        
        do i=1,N
                do j=1,N
                        if (i/=j) then ! no force by particle on itself
                                
                                d = sqrt(sum((r(i,:)-r(j,:))**2))      
                                if (d<rc) then
                                        !only particle pairs with d below cutoff contribute
                                        F(j) = LJ(d,0.001) !calculate magnitude of the force on i by j, LF(d,sigma)
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

subroutine InitCell(r,a,N) 
           implicit none 

           real, intent(in) :: a
           integer, intent(in) :: N
           real, intent(out) :: r(N,3)
           integer :: i, j, k, atom, S, M
           real :: unitcell(4,3) = 0
           
           ! define unit cell of the FCC lattice (only nonzero coordinates) 
           
           unitcell(2,1:2) = a/sqrt(2.0)
           unitcell(3,2:3) = a/sqrt(2.0)
           unitcell(4,1:3:2) = a/sqrt(2.0)
           
           M = int(N**(1.0/3.0)/4.0) !nr of cell shifts in any direction
           S = 0
           
           ! shifts the unit cell in steps of a (in x,y,z) to form an FCC lattice, with N "atoms"

           do i=0,M-1
                do j=0,M-1
                        do k=0,M-1
                                do atom=1,4
                                        r(atom+S,:)=unitcell(atom,:) + a*real([i,j,k])
                                enddo

                                S = S+4
                        enddo
                enddo 
           enddo
end subroutine InitCell  

subroutine InitVel(v,m,beta,N)
        use maxwell !dont forget this here ! 
        implicit none

        real,intent(in) :: m, beta
        integer, intent(in) :: N
        real, intent(out) :: v(N,3)
        integer :: i, j
        ! pick velocity components randomly from maxwell boltzmann velocity distribution
        do i=1,N
                do j=1,3
                        v(i,j) = MB(m,beta) 
                enddo
        enddo

end subroutine InitVel 

