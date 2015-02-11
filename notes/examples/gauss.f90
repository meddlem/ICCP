module Gauss

    implicit none
    private

    public smart_sum

contains

    integer function smart_sum(a, b) result(total)

        integer, intent(in) :: a, b

        total = (b * (b + 1)) / 2 - (a * (a - 1)) / 2

    end function

end module