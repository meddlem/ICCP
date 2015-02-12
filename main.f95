program main
      ! modules to use 
      ! use functions
      use plplot !plotting library 
      use plotpoints !module for plotting particles
      use Inits !module for initializing model
      use Incs !module for calculating iterations

      implicit none
      ! model parameters (constants):
      ! dt = timestep, rc = potential cutoff length, beta ~ 1/T, m = mass, L=length, eps and sigma belong to LJ potential 
      real(8), parameter :: dt = 0.3, rc = 1., sigma = 0.01, beta = 1., m = 1., eps = 1., L = 10. 
      integer, parameter :: N = 864 !number of particles, multiple of 4
      real(8), parameter :: a = L/((N/4)**(1./3.)), alpha = 24. * eps / m 
      ! declare variables
      real(8) :: T, r(N,3), v(N,3) ! declare position and velocity vectors 
      integer :: i

      
      ! initialize r and v
      call InitCell(r,a,N)
      call InitVel(v,m,beta,N)
      print *, "vsum t=0:", sum(v) !test if total momentum is constant 
      
      ! initialize plot
      call plotinit(-0.1*L,1.1*L) 

      !where the magic happens, also there is a memory leak somewhere here: 
      do i = 1,50
                T = sum(v**2)/N ! calculate the temperatue at timestep i
                call plot_points(r) !calls plot points
                call Vinc(r,v,dt,alpha,rc,sigma,L) !calc velocities
                call Rinc(r,v,dt,L)  !calc particle positions
                print *, T
      enddo
      print *, "vsum final:", sum(v) !test if total momentum is constant
      call plend()

end program main 

subroutine plotinit(xmin,xmax)
        use plplot !library for plotting
        implicit none

        real(8), intent(in) :: xmin, xmax

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
        ! define viewport, world coords for the edges
        call plvpor(0d0, 1d0, 0d0, 0.9_plflt)
        call plwind(-1d0,1d0,-1d0,1.5_plflt)
        call plw3d(1d0, 1d0, 1.2_plflt, xmin, xmax, xmin, xmax, xmin, xmax, 20d0, 45d0)

        call plspause(.false.)

end subroutine plotinit 
