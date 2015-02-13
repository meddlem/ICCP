program main
! modules to use 
use plplot !plotting library 
use plotpoints !module for plotting particles
use Inits !module for initializing model
use Interactions !module for calculating the interaction forces

implicit none

! model parameters (constants):
! dt = timestep, rc = potential cutoff length, beta ~ 1/inital temp, m = mass, L=length, units: eps=1, sigma=1, m=1
real(8), parameter :: dt = 0.0005d0, rc = 2.5d0, beta = 1d0, L = 14d0, eps = 1d0, sigma = 1d0, m = 1d0  
integer, parameter :: N = 8**3*4 !number of particles, must fit the FCC lattice 
! declare variables
real(8) :: T, E, EV, r(N,3), v(N,3), F(N,3) 
integer :: i

! initialize r, v, F
call InitCell(r,L,N)
call InitVel(v,beta,N)
call Force(F,EV,r,rc,L)

! initialize 3Dplot
call Pplotinit(-0.1d0*L,1.1d0*L) 

T = sum(v**2)/N ! calculate T
E = EV + 0.5d0*sum(v**2) ! calculate the total energy

print *, "Psum t=0: ", sum(v) !test if total momentum is constant 
print *, "Initial E: ", E/N 
print *, "T initial: ", T

! time integration using velocity Verlet algorithm: 
do i = 1,1000
        !T = sum(v**2)/N ! calculate the temperatue at timestep i
        ! print *, T
        call Pplot(r) !calls plot points
       
        r = r + v*dt + 0.5d0*F*(dt**2) !update positions
        r = r - floor(r/L)*L !enforce PBC
        v = v + 0.5d0*F*dt !update velocity (1/2)
        call Force(F,EV,r,rc,L) !update force to next timestep
        v = v + 0.5d0*F*dt !update velocity (2/2)

        E = E + EV + 0.5d0*sum(v**2) ! calculate the total energy
enddo

T = sum(v**2)/N ! calculate total energy, should be constant

print *, "T final: ", T
print *, "E mean: ", E/(1000d0*N) 
print *, "Psum final: ", sum(v)

call plend()

end program main 
