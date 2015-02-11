program SumsBeforeGauss

    implicit none

    integer :: first, last

    print *, "First:"
    read *, first
    print *, "Last:"
    read *, last

    print *, "Total:", stupid_sum(first, last)

contains

    integer function stupid_sum(a, b) result(total)

        integer, intent(in) :: a, b

        integer :: i

        total = 0
        do i = a, b
            total = total + i
        end do

    end function

end program
