program MyProg

    implicit none

    real(8), parameter :: PI = 4 * atan(1d0)
    integer, parameter :: N = 20

    real(8) :: radius, surface, circumference
    integer :: i, fibonacci(N)

    print *, "Give the radius"
    read *, radius
    surface = PI * radius**2
    circumference = 2 * PI * radius
    print *, "The surface of the circle is", surface
    print *, "and the circumference is", circumference

    fibonacci(1) = 0
    fibonacci(2) = 1
    do i = 3, N
        fibonacci(i) = fibonacci(i - 1) + fibonacci(i - 2)
    end do
    print *, "Fibonacci sequence:"
    print "(4I6)", fibonacci
    ! The ratio of neighboring Fibonacci numbers
    ! converges to the golden ratio.
    print *, "Golden ratio:", &
             1d0*fibonacci(N)/fibonacci(N - 1)

end program
