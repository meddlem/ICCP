module io
  use constants
  implicit none
  private
  public :: user_in, results_out, f_check
contains

  subroutine user_in(rho,T_init)
    real(dp), intent(out) :: rho, T_init
  
    write(*,'(/,A,/)') '************ Input *************' 
    write(*,'(A)',advance='no') "Number density = " 
    read(*,*) rho
    write(*,'(A)',advance='no') "Target temperature = " 
    read(*,*) T_init
    write(*,'(A,I5,A)') "Starting simulation: ", steps, " iterations"
  end subroutine

  subroutine results_out(runtime, rho, T_init, eq_pres, err_p, heat_cap, &
      err_heat, mean_t, err_t, eq_u, err_u, d, err_d)
    real(dp), intent(in) :: rho, T_init, eq_pres, err_p, heat_cap, err_heat, &
      mean_T, err_T, eq_U, err_U, D, err_D
    integer, intent(in) :: runtime

    open(10,access = 'sequential',file = 'output.txt')
    
    write(10,'(/,A,/)') '*********** Summary ***********' 
    write(10,output_format) "Number density : ", rho
    write(10,output_format) "Initial temperature : ", T_init 
    
    write(10,'(/,A,/)') '*********** Output ************' 
    write(10,'(A,I6,A)') "Runtime : ", runtime, " s"
    write(10,output_format) "Equilibrium pressure : ", eq_pres, " ± ", err_p
    write(10,output_format) "Heat capacity : ", heat_cap, " ± ", err_heat
    write(10,output_format) "Simulation temperature : ", mean_T, " ± ", err_T
    write(10,output_format) "Potential energy : ", eq_U, " ± ", err_U
    write(10,output_format) "Diffusion constant : ", D, " ± ", err_D
    write(10,'(/,A,/)') '*******************************' 
    
    close(10,status = 'keep')
    
    call system('cat output.txt')
  end subroutine

  subroutine f_check(p)
    real(dp), intent(in) :: p(:,:)
    
    if (maxval(abs(sum(p,1)))>1d-8) then
      print *, "Warning: momentum was not conserved"
    endif
  end subroutine
end module
