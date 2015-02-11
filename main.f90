program main
      ! modules to use 
      use maxwell !gives random velocity occording to MB distribution
      use LJforce !function that gives force based on LJ potential
      use plplot !code for plotting particles
      use Inits !includes code for initializing model
      use Incs !code for calculating timesteps

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
                call Rinc(r,v,dt)  !calc particle positions
                call Vinc(r,v,dt,m,rc) !calc velocities
      enddo

      call plend()

end program main 

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
