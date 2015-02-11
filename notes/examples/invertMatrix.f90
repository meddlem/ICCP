! Calculates the inverse of a real symmetric matrix.
!
! a - On input, a real symmetric matrix. Only the upper
!     triangular part is accessed. On output, A is 
!     overwritten by its inverse.
subroutine matinvrs(a)

    real(8), intent(inout) :: a(:, :)

    integer, external :: ilaenv

    integer              :: i, j, n, nb, lwork, info
    integer              :: ipiv(size(a, 1))
    real(8), allocatable :: work(:)

    n = size(a, 1)

    ! Calculate the optimal size of the workspace array.
    nb = ilaenv(1, "DSYTRI", "U", n, -1, -1, -1)
    lwork = n * nb

    allocate (work(lwork))

    ! Invert the matrix.
    call dsytrf("U", n, a, n, ipiv, work, lwork, info)
    if (info /= 0) stop "error in call to dsytrf"
    call dsytri("U", n, a, n, ipiv, work, info)
    if (info /= 0) stop "error in call to dsytri"

    deallocate (work)

    ! Copy the upper triangular part of A to the lower.
    do j = 1, n - 1
        do i = j + 1, n
            a(i, j) = a(j, i)
        end do
    end do

end subroutine
