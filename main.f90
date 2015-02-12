program main
      ! modules to use 
      ! use functions
      use plplot !plotting library 
      use plotpoints !module for plotting particles
      use Inits !module for initializing model
      use Incs !module for calculating iterations

      implicit none
      ! model parameters, constants
      real(8), parameter :: dt = 0.3, a = 1., rc = 5., sigma = 0.01 ! defines: timestep, lattice constant (initial), potential cutoff length
      integer, parameter :: N = 4000 !number of particles, multiple of 4
      real(8), parameter :: beta = 1., m = 1., eps = 1., L = 10. ! beta ~ 1/T, m = mass, L=length 
      real(8), parameter :: alpha = 24 * eps / m ! combine m and epsilon (of the LJ potential) into one constant 

      ! some variables
      real(8) :: r(N,3), v(N,3) ! declare position and velocity vectors 
      integer :: i
      
      ! initialize r and v
      call InitCell(r,a,N)
      call InitVel(v,m,beta,N)
      
      ! initialize plot
      call plotinit(0d0,L) 

      !where the magic happens: 
      do i = 1,50
                call plot_points(r) !calls plot points
                call Rinc(r,v,dt,L)  !calc particle positions
                call Vinc(r,v,dt,alpha,rc,sigma) !calc velocities
      enddo

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
