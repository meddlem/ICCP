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
    print '(A,I5,A)', "starting simulation: ", steps, " iterations"
  end subroutine

  subroutine results_out(runtime,eq_pres,err_p,heat_cap,err_heat, &
      mean_T, err_T, eq_U,err_U,D,err_D)
    real(dp), intent(in) :: eq_pres, err_p, heat_cap, err_heat, mean_T, &
      err_T, eq_U, err_U, D, err_D
    integer, intent(in) :: runtime
      
    print '(A,I4,A)', " runtime = ", runtime, " s"
    print '(A,F6.4,A,F6.4)', "equilibrium pressure =", eq_pres, "±", err_p
    print '(A,F7.2,A,F6.2)', "heat capacity =", heat_cap, "±", err_heat
    print '(A,F6.4,A,F6.4)', "T final =", mean_T, "±", err_T
    print '(A,F7.4,A,F7.4)', "U equilibrium =", eq_U, "±", err_U
    print '(A,F7.5,A,F7.5)', "D =", D, "±", err_D

  end subroutine

  subroutine f_check(p)
    real(dp), intent(in) :: p(:,:)
    
    if (maxval(abs(sum(p,1)))>1d-8) then
      print *, "warning: momentum was not conserved"
    endif
  end subroutine
end module
