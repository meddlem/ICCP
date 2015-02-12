module plotpoints
      implicit none
      private
      public :: plot_points

contains

subroutine plot_points(r)
        !plots all particle position
        use plplot
        real(8), intent(in) :: r(:,:) !r(N,3) 

        call plclear()
        call plcol0(1) !axis color
        call plbox3("bnstu", "x", 0d0, 0, "bnstu", "y", 0d0, 0, "bcdmnstuv", "z", 0d0, 0) !plots the axes etc
        
        call plcol0(3) !point color
        call plpoin3(r(:,1), r(:,2), r(:,3), 4) !this plots the points
        call plflush()
end subroutine plot_points 

end module 
