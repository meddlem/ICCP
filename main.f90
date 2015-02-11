program main
      use maxwell 
      use LJforce
      use plplot
      implicit none
      ! model parameters, constants
      real(8), parameter :: dt = 1., a = 1., rc = 5. ! defines: timestep, lattice constant (initial), potential cutoff length
      integer, parameter :: N = 4000 !number of particles, multiple of 4
      real(8), parameter :: beta = 1., m = 1., L = 10. ! beta ~ 1/T, m = mass, L=length 
      
      !some variables
      real(8) :: r(N,3), v(N,3) ! declare position and velocity vectors 
      integer :: i
      
      ! initialize r and v
      call InitCell(r,a,N)
      call InitVel(v,m,beta,N)
      
      ! plotting stuff 
      ! redefining colors
      call plscol0(0,255,255,255)
      call plscol0(1,255,0,0)
      call plscol0(2,0,255,0)
      call plscol0(3,0,0,255)
      call plscol0(4,255,0,255)
      call plscol0(5,0,255,255)
      call plscol0(6,255,255,0)
      call plscol0(7,0,0,0)
      call plscol0(8,255,70,0)
      call plscol0(9,128,128,128)

      call plsdev("xcairo")
      call plinit()
      call pladv(0)
      !define viewport, world coords for the edges
      call plvpor(0d0, 1d0, 0d0, 0.9_plflt)
      call plwind(-1d0,1d0,-1d0,1.5_plflt)
      call plw3d(1d0, 1d0, 1.2_plflt, 0d0, L, 0d0, L, 0d0, L, 20d0, 45d0)

      call plspause(.false.)

      !where the magic happens: 
      do i = 1,25
                call plot_points(r,N) !calls plot points
                call Rinc(r,v,dt,N)  !calc particle positions
                call Vinc(r,v,dt,m,rc,N) !calc velocities
      enddo

      call plend()

end program main 

subroutine Rinc(r,v,dt,N)
        ! this subroutine calculates particle positions
        implicit none
        integer :: N 
        real(8), intent(in) :: dt, v(N,3)
        real(8), intent(inout) :: r(N,3)

        r = r + v*dt !computes positions at next timestep  
end subroutine Rinc 

subroutine Vinc(r,v,dt,m,rc,N)
        ! calcs the particles velocity for next timestep 
        use LJforce ! dont forget about the module bro
        implicit none
        integer :: N, i, j, k
        real(8), intent(in) :: m, dt, r(N,3)
        real(8), intent(inout) :: v(N,3)
        real(8) :: rc, d, sigma, F(N), deltaV(N,3), dir(N,3)
        !using r at timestep n+1, v at n: compute v at timestep n+1
        sigma = 0.1
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

subroutine InitCell(r,a,N) 
        ! gives initial positions based on FCC lattice
           implicit none 

           real(8), intent(in) :: a
           integer, intent(in) :: N
           real(8), intent(out) :: r(N,3)
           integer :: i, j, k, atom, S, M
           real(8) :: unitcell(4,3) = 0
           
           ! define unit cell of the FCC lattice (only nonzero coordinates) 
           
           unitcell(2,1:2) = a/sqrt(2.0)
           unitcell(3,2:3) = a/sqrt(2.0)
           unitcell(4,1:3:2) = a/sqrt(2.0)
           
           M = int((N/4)**(1.0/3.0)) !nr of cell shifts in any direction
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
        ! gives initial velocites based on maxwell-boltzmann dist
        use maxwell !dont forget this here ! 
        implicit none

        real(8),intent(in) :: m, beta
        integer, intent(in) :: N
        real(8), intent(out) :: v(N,3)
        integer :: i, j
        ! pick velocity components randomly from maxwell boltzmann velocity distribution
        do i=1,N
                do j=1,3
                        v(i,j) = MB(m,beta) 
                enddo
        enddo

end subroutine InitVel 

subroutine plot_points(r,N)
        !plots all particle position
        use plplot
        real(8), intent(in) :: r(N,N)

        call plclear()
        call plcol0(1) !axis color
        call plbox3("bnstu", "x", 0d0, 0, "bnstu", "y", 0d0, 0, "bcdmnstuv", "z", 0d0, 0) !plots the axes etc
        call plcol0(3) !point color
        call plpoin3(r(:,1), r(:,2), r(:,3), 4) !this plots the points
        call plflush()
end subroutine plot_points 
