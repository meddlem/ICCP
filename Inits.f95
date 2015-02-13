module Inits
        implicit none
        private
        public :: InitVel, InitCell

contains

subroutine InitCell(r,a,N) 
        ! gives initial positions based on FCC lattice
           implicit none 

           real(8), intent(in) :: a
           integer, intent(in) :: N
           real(8), intent(out) :: r(N,3)
           integer :: i, j, k, atom, S, M
           real(8) :: unitcell(4,3) = 0
           
           ! define unit cell of the FCC lattice (only nonzero coordinates) 
           
           unitcell(2,1:2) = a/dsqrt(2d0)
           unitcell(3,2:3) = a/dsqrt(2d0)
           unitcell(4,1:3:2) = a/dsqrt(2d0)
           
           M = int((N/4)**(1./3.)) !nr of cell shifts in any direction
           S = 0
           
           ! shifts the unit cell in steps of a (in x,y,z) to form an FCC lattice, with N "atoms"

           do i=0,M-1
                do j=0,M-1
                        do k=0,M-1
                                do atom=1,4
                                        r(atom+S,:)=unitcell(atom,:) + a*real([i,j,k])
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
        ! pick velocity components randomly from maxwell boltzmann velocity distribution
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
                real(8) :: u, v, s, std 

                s = 1 !initialize s
                std = dsqrt(1/beta) !define std of velocity distribution
                ! generate normal dist number with std as above, 
                ! using polar box-muller
                do while (s>=1d0)
                        u = 2d0 * rand() - 1d0
                        v = 2d0 * rand() - 1d0
                        s = u**2 + v**2
                end do

                MB = u * std * dsqrt((-2d0 * log(s))/s) 
        end function MB

end subroutine InitVel 

end module 
