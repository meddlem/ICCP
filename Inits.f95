module Inits
        implicit none
        private
        public :: InitVel, InitCell

contains

subroutine InitCell(r,L,N) 
        ! gives initial positions based on FCC lattice
        implicit none 

        real(8), intent(in) :: L
        integer, intent(in) :: N
        real(8), intent(out) :: r(N,3)
        integer :: i, j, k, atom, S, M
        real(8) :: a, unitcell(4,3)
        
        unitcell(:,:) = 0d0 !initialize unitcell 
        a = L/((N/4)**(1d0/3d0)) ! calculate lattice constant
           
        ! define unit cell of the FCC lattice (only nonzero coordinates) 
        unitcell(2,1:2) = a/sqrt(2d0)
        unitcell(3,2:3) = a/sqrt(2d0)
        unitcell(4,1:3:2) = a/sqrt(2d0)
           
        M = int((N/4)**(1./3.)) !nr of cell shifts in any direction
        S = 0
   
        ! shifts the unit cell in steps of a (in x,y,z) to form an FCC lattice, with N "atoms"

        do i=0,M-1
                do j=0,M-1
                        do k=0,M-1
                                do atom=1,4
                                        r(atom+S,:) = unitcell(atom,:) + a*real([i,j,k],kind=8)
                                enddo

                                S = S+4
                        enddo
                enddo 
        enddo
end subroutine InitCell  

subroutine InitVel(v,beta,N)
        ! gives initial velocites based on maxwell-boltzmann dist
        implicit none

        real(8),intent(in) :: beta
        integer, intent(in) :: N
        real(8), intent(out) :: v(N,3)
        integer :: i, j
        real(8) :: Vavg(3)
        
        call init_random_seed()

        ! pick velocity components from MB velocity distribution
        do i=1,N
                do j=1,3
                        v(i,j) = MB(beta) 
                enddo
        enddo
        
        Vavg = sum(v,1)/N

        do i = 1,3
                v(:,i) = v(:,i) - Vavg(i) ! make sure total momentum = 0
        enddo 
        
contains

        real(8) function MB (beta)
                ! gives random velocity (component) based on maxwell boltzmann
                ! distribution
                real(8), intent(in) :: beta
                real(8) :: u(2), std, pi

                pi = 4*atan(1d0) 

                std = sqrt(1d0/beta) !define std of velocity distribution
                ! generate normal dist number with std as above, 
                ! using box-muller

                call random_number(u)

                MB = std*sqrt(-2d0*log(u(1)))*cos(2*pi*u(2))
        end function MB

! initialize random seed, taken from ICCP github
subroutine init_random_seed()
  implicit none
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid, t(2), s
  integer(8) :: count, tms

  call random_seed(size = n)
  allocate(seed(n))
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
        call date_and_time(values=dt)
        tms = (dt(1) - 1970)*365_8 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24 * 60 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
        t = transfer(tms, t)
     end if
     s = ieor(t(1), t(2))
     pid = getpid() + 1099279 ! Add a prime
     s = ieor(s, pid)
     if (n >= 3) then
        seed(1) = t(1) + 36269
        seed(2) = t(2) + 72551
        seed(3) = pid
        if (n > 3) then
           seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
        end if
     else
        seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
     end if
  end if
  call random_seed(put=seed)
end subroutine init_random_seed

end subroutine InitVel 

end module 
