function box_muller(xs) result(zs)
  implicit none
  real(8) :: pi = 4*atan(1d0), zs(2)
  real(8), intent(in) :: xs(2)

  zs(1) = sqrt(-2d0*log(xs(1)))*cos(2*pi*xs(2))
  zs(2) = sqrt(-2d0*log(xs(1)))*sin(2*pi*xs(2))
end function box_muller
