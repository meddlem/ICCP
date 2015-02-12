module functions
      implicit none
      private
      public :: MB, LJ
contains
 
real(8) function MB (m,beta)
        ! gives random velocity (component) based on maxwell boltzmann
        ! distribution
        real(8), intent(in) :: m, beta
        real(8) :: u, v, s, std 
        integer :: i = 0
        s = 1 !initialize s
        std = sqrt(1/(m*beta)) !define std of velocity distribution
        ! generate normal dist number with std as above, using polar box-muller
        do while (s>=1.0)
                i = i + 1
                u = 2.0 * rand() - 1.0
                v = 2.0 * rand() - 1.0
                s = u**2 + v**2
        end do

        MB = u * std * sqrt((- 2.0 * log(s))/s) !could also include v*.. 
end function MB

real(8) function LJ (r,sigma)
        ! gives interaction force between two particle from lennart jones 
        ! potential 
        real(8), intent(in) :: r, sigma
        LJ = (2*(sigma**12)/(r**13)-(sigma**6)/(r**7))
end function

end module
