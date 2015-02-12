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
           
           unitcell(2,1:2) = a/sqrt(2.0)
           unitcell(3,2:3) = a/sqrt(2.0)
           unitcell(4,1:3:2) = a/sqrt(2.0)
           
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

subroutine InitVel(v,m,beta,N)
        ! gives initial velocites based on maxwell-boltzmann dist
        use functions !dont forget this here ! 
        implicit none

        real(8),intent(in) :: m, beta
        integer, intent(in) :: N
        real(8), intent(out) :: v(N,3)
        integer :: i, j
        ! pick velocity components randomly from maxwell boltzmann velocity distribution
        do i=1,N
                do j=1,3
                        v(i,j) = MB(m,beta) 
                enddo
        enddo

end subroutine InitVel 

end module 
