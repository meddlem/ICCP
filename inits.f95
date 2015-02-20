module inits
  use constants
  implicit none
  private
  public :: init_p, init_r

contains
  subroutine init_r(r,L) 
    ! gives initial positions based on FCC lattice
    real(dp), intent(out) :: r(N,3)
    real(dp), intent(in) :: L
    integer :: i, j, k, atom, S, M
    real(dp) :: a, unitcell(4,3)

    unitcell(:,:) = 0._dp !initialize unitcell 
    a = L/((N/4)**(1._dp/3._dp)) ! calculate lattice constant

    ! define unit cell of the FCC lattice (only nonzero coordinates) 
    unitcell(2,1:2) = a/sqrt(2._dp)
    unitcell(3,2:3) = a/sqrt(2._dp)
    unitcell(4,1:3:2) = a/sqrt(2._dp)

    M = int((N/4)**(1./3.)) ! # unit cell shifts in any direction
    S = 0 

    ! shifts the unit cell in steps of a (in x,y,z) to form an FCC lattice,&
    ! with N "atoms"
    do i = 0,M-1
      do j = 0,M-1
        do k = 0,M-1
          do atom = 1,4
            r(atom+S,:) = unitcell(atom,:) + a*real([i,j,k],kind=8)
          enddo
          S = S+4
        enddo
      enddo 
    enddo
  end subroutine  

  subroutine init_p(p,T_tgt)
    ! gives initial momenta based on maxwell-boltzmann dist
    real(dp), intent(out) :: p(N,3)
    real(dp), intent(in) :: T_tgt
    integer :: i, j
    real(dp) :: Pavg(3)

    call init_random_seed()

    ! pick momentum components from MB distribution
    do i=1,N
      do j=1,3
        p(i,j) = MB(T_tgt) 
      enddo
    enddo
    
    Pavg = sum(p,1)/N
   
    do i =1,3
      p(:,i) = p(:,i) - Pavg(i) ! make sure total momentum = 0
    enddo 
  end subroutine 

  ! initialize random seed, taken from ICCP github
  subroutine init_random_seed()
    integer, allocatable :: seed(:)
    integer :: i, m, un, istat, dtime(8), pid, t(2), s
    integer(8) :: count, tms

    call random_seed(size = m)
    allocate(seed(m))
    open(newunit=un, file="/dev/urandom", access="stream",&
      form="unformatted", action="read", status="old", &
      iostat=istat)
    if (istat == 0) then
      read(un) seed
      close(un)
    else
      call system_clock(count)
      if (count /= 0) then
        t = transfer(count, t)
      else
        call date_and_time(values=dtime)
        tms = (dtime(1) - 1970)*365_8 * 24 * 60 * 60 * 1000 &
          + dtime(2) * 31_8 * 24 * 60 * 60 * 1000 &
          + dtime(3) * 24 * 60 * 60 * 60 * 1000 &
          + dtime(5) * 60 * 60 * 1000 &
          + dtime(6) * 60 * 1000 + dtime(7) * 1000 &
          + dtime(8)
        t = transfer(tms, t)
      end if
      s = ieor(t(1), t(2))
      pid = getpid() + 1099279 ! Add a prime
      s = ieor(s, pid)
      if (m >= 3) then
        seed(1) = t(1) + 36269
        seed(2) = t(2) + 72551
        seed(3) = pid
        if (m > 3) then
          seed(4:) = s + 37 * (/ (i, i = 0, m - 4) /)
        end if
      else
        seed = s + 37 * (/ (i, i = 0, m - 1 ) /)
      end if
    end if
    call random_seed(put=seed)
  end subroutine

  real(dp) function MB(T_tgt)
    ! gives random momentum (component) based on maxwell boltzmann distribution
    real(dp), intent(in) :: T_tgt
    real(dp) :: u(2), std

    std = sqrt(T_tgt) 
    ! generate normal dist number using box-muller
    call random_number(u)
    MB = std*sqrt(-2._dp*log(u(1)))*cos(2*pi*u(2))
  end function
end module 
