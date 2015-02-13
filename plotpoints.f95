module plotpoints
      implicit none
      private
      public :: Pplot, Pplotinit

contains

subroutine Pplot(r)
        !plots all particle position
        use plplot
        real(8), intent(in) :: r(:,:) !r(N,3) 

        call plclear()
        call plcol0(1) !axis color
        call plbox3("bnstu", "x", 0d0, 0, "bnstu", "y", 0d0, 0, "bcdmnstuv", "z", 0d0, 0) !plots the axes etc
        
        call plcol0(3) !point color
        call plpoin3(r(:,1), r(:,2), r(:,3), 4) !this plots the points
        call plflush()
end subroutine Pplot 

subroutine Pplotinit(xmin,xmax)
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
        call plvpor(0d0, 1d0, 0d0, 0.9d0)
        call plwind(-1d0,1d0,-1d0,1.5d0)
        call plw3d(1d0, 1d0, 1.2d0, xmin, xmax, xmin, xmax, xmin, xmax, 20d0, 45d0)

        call plspause(.false.)

end subroutine Pplotinit 

end module 
