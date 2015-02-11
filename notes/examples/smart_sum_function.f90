! Gauss was a genius! 1 + 2 + 3 + ... + n = n*(n + 1) / 2
integer function smart_sum(a, b) result(total)

    integer, intent(in) :: a, b

    total = (b * (b + 1)) / 2 - (a * (a - 1)) / 2

end function
