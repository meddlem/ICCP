module io
  use constants
  implicit none
  private
  public :: user_in, results_out, f_check
contains

  subroutine user_in(rho,T_init)
    real(dp), intent(out) :: rho, T_init
  
    write(*,'(A)',advance='no') "number density = " 
    read(*,*) rho
    write(*,'(A)',advance='no') "target temperature = " 
    read(*,*) T_init
  end subroutine

  subroutine results_out(runtime,eq_pres,err_p,heat_cap,err_heat, &
      T, eq_U,err_U,D)
    real(dp), intent(in) :: eq_pres, err_p, heat_cap, err_heat, T, eq_U, &
      err_U, D
    integer, intent(in) :: runtime
      
    print '(A,I4,A)', " runtime = ", runtime, " s"
    print *, "equilibrium pressure =", eq_pres
    print *, "err p =", err_p
    print *, "heat capacity =", heat_cap
    print *, "err hc =", err_heat
    print *, "T final =", T
    print *, "U equilibrium =", eq_U
    print *, "err U =", err_U
    print *, "D =", D 
  end subroutine

  subroutine f_check(p)
    real(dp), intent(in) :: p(:,:)
    
    if (maxval(abs(sum(p,1)))>1d-8) then
      print *, "warning: momentum was not conserved"
    endif
  end subroutine
end module
