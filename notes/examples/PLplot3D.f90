program PLplot3d
  implicit none

  call plsdev("xcairo")
  call plinit()

  call pladv(0)
  call plvpor(0d0, 1d0, 0d0, 1d0)
  call plwind(-1d0, 1d0, -2d0 / 3, 4d0 / 3)
  call plw3d(1d0, 1d0, 1d0, xmin, xmax, ymin, ymax, &
             zmin, zmax, 45d0, -45d0)
  call plspause(.false.)
  call plend()

contains

  subroutine plot_points(xyz)
    real(8), intent(in) :: xyz(:, :)

    call plclear()
    call plcol0(1)
    call plbox3("bnstu", "x", 0d0, 0, "bnstu", "y", &
                  0d0, 0, "bcnmstuv", "z", 0d0, 0)
    call plcol0(2)
    call plpoin3(xyz(1, :), xyz(2, :), xyz(3, :), 4)
    call plflush()
  end subroutine

end program PLplot3d
