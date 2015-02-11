program SumsBeforeGauss

    implicit none

    integer :: first, last, total

    print *, "First:"
    read *, first
    print *, "Last:"
    read *, last

    call stupid_sum(first, last, total)

    print *, "Total:", total

contains

    subroutine stupid_sum(a, b, total)

        integer, intent(in) :: a, b
        integer, intent(out) :: total

        integer :: i

        total = 0
        do i = a, b
            total = total + i
        end do

    end subroutine

end program
