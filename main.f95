program main
! modules to use 
! use functions
use plplot !plotting library 
use plotpoints !module for plotting particles
use Inits !module for initializing model
use Force !module for calculating iterations

implicit none
! model parameters (constants):
! dt = timestep, rc = potential cutoff length, beta ~ 1/T, m = mass, L=length, units: eps=1 and sigma=1
real(8), parameter :: dt = 0.00032, rc = 2.5, beta = 1., m = 1., L = 10. 
integer, parameter :: N = 6**3*4 !number of particles, must fit the FCC lattice 
real(8), parameter :: a = L/((N/4)**(1./3.)) ! lattice constant, strenght of potential (effective) 
! declare variables
real(8) :: T, r(N,3), v(N,3), F(N,3) ! declare position and velocity vectors 
integer :: i

! initialize r, v, F
call InitCell(r,a,N)
call InitVel(v,m,beta,N)
call ForceMatrix(r,F,m,rc,L)

!print *, "Psum t=0:", sum(v) !test if total momentum is constant 

! initialize plot
call plotinit(-0.1*L,1.1*L) 

T = sum(v**2)/N ! calculate T
print *, "T initial", T

! time integration using velocity Verlet algorithm: 
do i = 1,1200
       ! T = sum(v**2)/N ! calculate the temperatue at timestep i
      !  print *, T
        if (modulo(i,10)==0) then
                call plot_points(r) !calls plot points
        endif
        r = r + v*dt + 0.5*F*(dt**2) !update positions
        r = r - floor(r/L)*L !enforce PBC
        v = v + 0.5*F*dt !update velocity (1/2)
        call ForceMatrix(r,F,m,rc,L) !update force to next timestep
        v = v + 0.5*F*dt !update velocity (2/2)

enddo

T = sum(v**2)/N ! calculate total energy, should be constant
!print *, "T final:", T

print *, "Psum final:", sum(v) !test if total momentum is constant
call plend()

end program main 

subroutine plotinit(xmin,xmax)
        use plplot !library for plotting
        implicit none

        real(8), intent(in) :: xmin, xmax

        ! plotting stuff 
        ! redefining colors
        call plscol0(0,255,255,255)
        call plscol0(1,255,0,0)
        call plscol0(2,0,255,0)
        call plscol0(3,0,0,255)
        call plscol0(4,255,0,255)
        call plscol0(5,0,255,255)
        call plscol0(6,255,255,0)
        call plscol0(7,0,0,0)
        call plscol0(8,255,70,0)
        call plscol0(9,128,128,128)

        call plsdev("xcairo")
        call plinit()
        call pladv(0)
        ! define viewport, world coords for the edges
        call plvpor(0d0, 1d0, 0d0, 0.9_plflt)
        call plwind(-1d0,1d0,-1d0,1.5_plflt)
        call plw3d(1d0, 1d0, 1.2_plflt, xmin, xmax, xmin, xmax, xmin, xmax, 20d0, 45d0)

        call plspause(.false.)

end subroutine plotinit 
