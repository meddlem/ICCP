! Calculates the eigenvalues and, optionally, the 
! eigenvectors of a real symmetric matrix.
!
! jobz - If jobz = "N", compute eigenvalues only, 
!        if jobz = "V", compute full eigensystem.
! a    - On input, a real symmetric matrix. Only the upper 
!        diagonal part is accessed.
! w    - On output, contains the eigenvalues of A.
! v    - On output, if jobz = "V", the columns are the
!        eigenvectors of A.
subroutine mateigrs(jobz, a, w, v)

    implicit none

    character, intent(in) :: jobz
    real(8), intent(in)   :: a(:, :)
    real(8), intent(out)  :: w(:), v(:, :)

    integer              :: n, lwork, info
    real(8), allocatable :: work(:)

    n = size(a, 1)
    v = a

    ! Query the optimal size of the workspace.
    allocate (work(1))
    lwork = -1
    call dsyev(jobz, "U", n, v, n, w, work, lwork, info)
    lwork = int(work(1))
    deallocate (work)

    allocate (work(lwork))

    ! Diagonalize the matrix.
    call dsyev(jobz, "U", n, v, n, w, work, lwork, info)
    if (info /= 0) stop "error in call to dsyev"

    deallocate (work)

end subroutine
