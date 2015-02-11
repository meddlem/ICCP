program SumsBeforeGauss

    implicit none

    integer :: first, last, total, i

    print *, "First:"
    read *, first
    print *, "Last:"
    read *, last

    total = 0
    do i = first, last
        total = total + i
    end do

    print *, "Total:", total

end program