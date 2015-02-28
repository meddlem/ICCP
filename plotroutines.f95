module plotroutines
  use constants
  use plplot
  implicit none
  private
  public :: particle_plot, particle_plot_init, gnu_line_plot 

contains
  subroutine gnu_line_plot(x,y1,xlabel,ylabel,label1,title,plot_no,y2,label2)
    real(dp), intent(in) :: x(:), y1(:)
    real(dp), intent(in), optional :: y2(:)
    character(*), intent(in) :: xlabel, ylabel, label1, title
    character(*), intent(in), optional :: label2
    integer, intent(in) :: plot_no
    character(1024) :: filename
    integer :: i, ret, m
    real(dp) :: xmin, xmax, xrange(2), ymin, ymax, yrange(2)
    
    if (size(y1)/=size(x)) print *, "error, arguments must be same size"
    
    if (present(y2)) then
      if (size(y1)/=size(y2)) then 
        print *, "error, arguments must be same size"
        return
      endif
    endif
    
    m = size(x)
    xmin = minval(x)
    xmax = maxval(x)
    ymin = minval(y1)
    ymax = maxval(y1)

    xrange = [0._dp, xmax+(xmax-xmin)*0.1_dp]
    yrange = [ymin-(ymax-ymin)*0.1_dp, ymax+(ymax-ymin)*0.1_dp]

    open(10,access = 'sequential',file = 'xydata.dat')
    
    do i=1,m
      if (present(y2)) then
        write(10,*) x(i),y1(i),y2(i) ! write datapoints to file
      else
        write(10,*) x(i),y1(i) ! write datapoints to file
      endif
    enddo
    
    close(10,status = 'keep')
    
    ! create gnuplot command file
    write(filename,'(A,I1,A)') 'set output "plot',plot_no,'.eps"'
    open(10,access = 'sequential',file = 'gplot.txt')
    
    ! set output terminal  
    ! write(10,*) 'set term pngcairo size 640,480 enhanced font "Verdana,10"'
    write(10,*) 'set term epscairo size 13cm,9cm font "Verdana,15"'
    write(10,*) filename
    ! set line color definitions
    write(10,*) &
      'set style line 1 lt 1 lc rgb "#ff0000" lw 2 #red'
    write(10,*) &
      'set style line 2 lt 1 lc rgb "#0000ff" lw 2 #blue'
    ! axes 
    write(10,*) 'set style line 11 lc rgb "#808080" lt 1'
    ! write(10,*) 'set border 3 back ls 11'
    write(10,*) 'set border 31 back ls 11'
    write(10,*) 'set tics nomirror scale 0.75'
    write(10,*) 'set key right center'
    write(10,*) 'set mxtics 2'
    write(10,*) 'set mytics 2'
    ! grid 
    write(10,*) 'set style line 12 lc rgb "#808080" lt 0 lw 1'
    write(10,*) 'set grid back ls 12'
    ! plotrange
    write(10,*) 'set xrange [',xrange(1),':',xrange(2),']'
    write(10,*) 'set yrange [',yrange(1),':',yrange(2),']'
    ! plot labels
    write(10,*) 'set title "'//TRIM(title)//'"'
    write(10,*) 'set xlabel '//'"'//TRIM(xlabel)//'"'
    write(10,*) 'set ylabel '//'"'//TRIM(ylabel)//'"'
    
    if (m>0) then
      write(10,*) 'plot "xydata.dat" using 1:2 with line ls 10 t "", \'
      write(10,*) &
      ' "xydata.dat" using 1:2 with line ls 1 t "'//TRIM(label1)//'", \'
      if (present(y2)) then
        write(10,*) &
        ' "xydata.dat" using 1:3 with line ls 2 t "'//TRIM(label2)//'"'
      endif
    endif
    
    close(10,status = 'keep')

    ! now call gnuplot and plot the curves
    call system('gnuplot gplot.txt',ret)
    call system('rm gplot.txt',ret)
    call system('rm xydata.dat',ret)
  end subroutine 

  subroutine particle_plot_init(xmin,xmax)
    real(dp), intent(in) :: xmin, xmax
    
    call Colors()
    call plsdev("xcairo")
    call plinit()
    call pladv(0)
    
    ! define viewport, world coords for the edges
    call plvpor(0._dp, 1._dp, 0._dp, 0.9_dp)
    call plwind(-1._dp,1._dp,-1._dp,1.5_dp)
    call plw3d(1._dp, 1._dp, 1.2_dp, xmin, xmax, xmin, xmax, xmin, xmax, &
      20._dp, 45._dp)
    call plspause(.false.)
  end subroutine
 
  subroutine particle_plot(r)
    ! plots all particle position
    real(dp), intent(in) :: r(:,:) 
   
    call plclear()
    call plcol0(1) !axis color
    call plbox3("bnstu", "x", 0._dp, 0, "bnstu", "y", 0._dp, 0, "bcdmnstuv", "z", &
      0._dp, 0) !plots the axes etc
    call plcol0(3) !point color
    call plpoin3(r(:,1), r(:,2), r(:,3), 4) !this plots the points
    call plflush()
  end subroutine
  
  subroutine Colors()
    ! redefining plplot colors
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
  end subroutine
end module 
