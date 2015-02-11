module LJforce
        implicit none 
        private
        public :: LJ

contains
        real(8) function LJ (r,sigma)
                real(8), intent(in) :: r, sigma
                LJ = (2*(sigma**12)/(r**13)-(sigma**6)/(r**7))
        end function

end module
        
