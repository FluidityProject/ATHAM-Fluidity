module FW
  use FW_data
  implicit none
contains

  subroutine allocate_storage(number_of_tracers,n)
    integer :: n
!f2py integer, intent(hide), depend(number_of_tracers) :: n=shape(number_of_tracers,0)
    integer :: number_of_tracers(n)
    allocate(gas_tracers(number_of_tracers(1)))
    allocate(incompressible_tracers(number_of_tracers(2)))
    allocate(other_fields(number_of_tracers(3)))
  end subroutine allocate_storage

  subroutine finalize
    deallocate(gas_tracers)
    deallocate(incompressible_tracers)
    deallocate(other_fields)
  end subroutine finalize

  subroutine set_gas_tracers(new_val,old_val,source,n,i)
    integer :: n
!f2py integer, intent(hide), depend(new_val) :: n=shape(new_val,0)
    real, intent(in), dimension(n), target :: new_val, old_val, source
!f2py real, intent(inplace), dimension(n) :: new_val, old_val, source
    integer :: i
    gas_tracers(i)%new=>new_val
    gas_tracers(i)%old=>old_val
    gas_tracers(i)%source=>source
  end subroutine set_gas_tracers

  subroutine run_microphysics(current_time,dt)
    real, intent(in) :: current_time, dt
    interface
       subroutine microphysics_main(&
            t0,&
            t1,&
            t2,&
            time,timestep)
         use FW_data_type
         type(tracer), intent(inout), dimension(:) :: t0
         type(tracer), intent(inout), dimension(:) :: t1
         type(tracer), intent(inout), dimension(:) :: t2
         real, intent(in) :: time, timestep
       end subroutine microphysics_main
    end interface

    call microphysics_main(&
         gas_tracers,&
         incompressible_tracers,&
         other_fields,&
         current_time,dt)

  end subroutine run_microphysics

end module FW
