program Sums

    use Gauss

    implicit none

    integer :: first, last

    print *, "First:"
    read *, first
    print *, "Last:"
    read *, last

    print *, "Total:", smart_sum(first, last)

end program