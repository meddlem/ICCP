module fortplot
  implicit none
  private
  public :: gnulineplot

  contains
    subroutine gnulineplot(x,y,xlabel,ylabel,title1,title)
      real(8), intent(in) :: x(:), y(:)
      character(*), intent(in) :: xlabel, ylabel, title1, title
      integer :: i, ret, n
      real(8) :: xmin, xmax, ymin, ymax
      
      n = size(x)
      xmin = minval(x)
      xmax = maxval(x)
      ymin = minval(y)
      ymax = maxval(y)
      
      ! write data to files
      open(10,access = 'sequential',file = 'xydata.dat')
      do i=1,n
        write(10,*) x(i),y(i)
      enddo
      close(10,status = 'keep')

      ! create gnuplot command file
      open(10,access = 'sequential',file = 'gplot.txt')
      ! write(10,*) 'set output "plot.eps"'
      write(10,*) 'set term png font "Fira Mono" 13'
      write(10,*) 'set output "plot.png"'
      ! write(10,*) 'set style line 10 lt 1 lc rgbcolor "#ffffff" lw 15 #thick white'
      write(10,*) &
        'set style line 11 lt 1 lc rgbcolor "#000000" lw 1 #black'
      write(10,*) &
        'set style line 1 lt 1 lc rgbcolor "#ff0000" lw 2 #red'
      ! write(10,*) 'set style line 2 lt 1 lc rgbcolor "#0000ff" lw 4  #blue'
      write(10,*) 'set border 3 #black'
      write(10,*) 'set xtics nomirror'
      write(10,*) 'set ytics nomirror'
      write(10,*) 'set xrange [',0d0,':',&
        xmax+(xmax-xmin)*0.1D0,']'
      write(10,*) 'set yrange [',0d0,':',&
        ymax+(ymax-ymin)*0.1D0,']'
      write(10,*) 'set key right center'
      write(10,*) 'set title "'//TRIM(title)//'"'
      write(10,*) 'set xlabel '//'"'//TRIM(xlabel)//'"'
      write(10,*) 'set ylabel '//'"'//TRIM(ylabel)//'"'
      
      if (n>0) then
        write(10,*) 'plot "xydata.dat" using 1:2 with line ls 10 t "", \'
        write(10,*) &
          '     "xydata.dat" using 1:2 with line ls 1 t "'//TRIM(title1)//'", \'
      endif 

      close(10,status = 'keep')

      ! now call gnuplot and plot the curve
      ret = system('gnuplot gplot.txt')
      ret = system('rm gplot.txt')
      ret = system('rm xydata.dat')
    end subroutine gnulineplot 
end module
