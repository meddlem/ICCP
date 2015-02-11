program FCC
      implicit none 
     ! private
     ! public FCCcell
     FCCcell(10.0,1.0) 

contains
      real function FCCcell (L,a) result(lattice)
              real, intent(in):: L,a
              integer, parameter :: N=10 
              integer :: i, j, k, S, atom 
              ! defines datastructure for saving particle coordinates under a particle label 
              type :: particle
                      real :: coords(1:3)=0 
              end type particle

              type(particle), dimension(3*N*4) :: lattice
              type(particle), dimension(4) :: unitcell

              ! creates unit cell for FCC lattice  (only the nonzero coords)
              unitcell(2)%coords(1) = a/sqrt(2.0)
              unitcell(2)%coords(2) = a/sqrt(2.0)
              unitcell(3)%coords(2) = a/sqrt(2.0)
              unitcell(3)%coords(3) = a/sqrt(2.0)
              unitcell(4)%coords(1) = a/sqrt(2.0)
              unitcell(4)%coords(3) = a/sqrt(2.0)
              
              !generate an L*L*L FCC lattice by translation of the unit cell

              do i=0,N
                        do j=0,N
                                do k=0,N
                                        S=i+j+k
                                        do atom=1,4
                                                lattice(atom+4*S)%coords(1)=unitcell(atom)%coords(1)+a*real(i) !translated x pos 
                                                lattice(atom+4*S)%coords(2)=unitcell(atom)%coords(2)+a*real(j) ! y pos
                                                lattice(atom+4*S)%coords(3)=unitcell(atom)%coords(3)+a*real(k) ! z pos
                                        end do
                                end do
                        end do
              end do 
      end function

end program
