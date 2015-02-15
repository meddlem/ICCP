module Plotroutines
  use plplot
  implicit none
  private
  public :: ParticlePlot, ParticlePlotinit, LinePlot

contains
  subroutine ParticlePlotinit(xmin,xmax)
    real(8), intent(in) :: xmin, xmax

    ! plotting stuff 
    call Colors()
    call plsdev("xcairo")
    call plinit()
    call pladv(0)

    ! define viewport, world coords for the edges
    call plvpor(0d0, 1d0, 0d0, 0.9d0)
    call plwind(-1d0,1d0,-1d0,1.5d0)
    call plw3d(1d0, 1d0, 1.2d0, xmin, xmax, xmin, xmax, xmin, xmax, 20d0, 45d0)

    call plspause(.false.)
    
  end subroutine ParticlePlotinit  

  subroutine ParticlePlot(r)
    !plots all particle position
    real(8), intent(in) :: r(:,:) !r(N,3) 
    
    call plclear()
    call plcol0(1) !axis color
    call plbox3("bnstu", "x", 0d0, 0, "bnstu", "y", 0d0, 0, "bcdmnstuv", "z", &
      0d0, 0) !plots the axes etc
    call plcol0(3) !point color
    call plpoin3(r(:,1), r(:,2), r(:,3), 4) !this plots the points
    call plflush()
  end subroutine ParticlePlot 

  subroutine Lineplot(x,y,xrange,yrange,xlabel,ylabel,label)
    real(8), intent(in) :: xrange(2), yrange(2), x(:), y(:)
    character(*), intent(in) :: xlabel, ylabel, label 

    call plparseopts(PL_PARSE_FULL)
    call Colors() 
    call plsdev("xcairo")
    call plinit()

    !call plcol0(3)
    call plenv(xrange(1),xrange(2),yrange(1),yrange(2),0,0)
    call pllab(xlabel,ylabel,label) 

    !call plcol0(1)
    call plline(x,y)
    call plspause(.true.) 
    !call plcol(2)

    call plend()
  end subroutine LinePlot 

  subroutine Colors()
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
  end subroutine Colors

end module 
