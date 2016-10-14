subroutine microphysics_main(current_time,dt,gt,it,of)

  use fw_data_type

  real, parameter :: L_v=2500000.0
  real, parameter, dimension(2) :: c_v=(/714.285714, 1368.5/),&
       c_p=(/1000.0, 1840.0/)
  real, parameter, dimension(2) :: c=(/4100.0,4100.0/),&
       rho_s=(/1024.0,1024.0/)

  type(tracer), dimension(:) :: gt, it, of
  real :: current_time, dt
  real, dimension(size(gt(1)%new),size(it)) :: q_s
  real, dimension(size(gt(1)%new)) :: CConN
  integer :: i

  CConN=ConN(of(1)%new,gt(1)%new,it,of(2)%new)
  gt(1)%source=CConN
  of(1)%source=L_v*CConN/cc_p(gt,it)
  
  contains

    function cc_p(gas,solid)
      type(tracer), intent(in), dimension(:) :: gas, solid
      real, dimension(size(gt(1)%new)) :: cc_p
      integer :: n

      cc_p=c_p(1)
      do n=1,size(gas)
         cc_p=cc_p+(c_p(n+1)-c_p(1))*gas(n)%new
      end do
      do n=1,size(solid)
         cc_p=cc_p+(c(n+1)-c_p(1))*solid(n)%new
      end do
    end function cc_p
      

    function p_sat(T)
      real, intent(in), dimension(:) :: T
      real, dimension(size(T)) :: p_sat
      integer :: n

      forall (n=1:size(T))
         p_sat(n)=611.2*exp(17.62*(max(T(n),273.0)-273.0)/(max(T(n),273.0)-30.0))
      end forall
    end function p_sat
      

    function q_sat(T,q_s,rho)
      real, intent(in), dimension(:) :: T,rho
      type(tracer), intent(in), dimension(:) :: q_s
      real, dimension(size(T)) :: q_sat,irhog
      integer :: n

      irhog=1.0/rho
      do n=1,size(q_s)
         irhog=irhog-q_s(n)%new/rho_s(n)
      end do
      q_sat=p_sat(T)/(R_v*T)*(irhog)

    end function q_sat

    function ConN(T,q_v,q_s,rho)
      real, intent(in), dimension(:) :: T,q_v,rho
      type(tracer), intent(in), dimension(:) :: q_s
      real, dimension(size(T)) :: ConN, qq_sat
      integer :: n
      
      qq_sat=q_sat(T,q_s,rho)

      forall (n=1:size(T))
         ConN(n)=max((q_v(n)-qq_sat(n))/dt,0.0)
      end forall

    end function ConN
      
end subroutine microphysics_main



subroutine microphysics_main_pointwise(current_time,dt,gt,it,of)

  use fw_data_type

  real, parameter :: L_v=2500000.0
  real, parameter, dimension(2) :: c_v=(/714.285714, 1368.5/),&
       c_p=(/1000.0, 1840.0/)
  real, parameter, dimension(2) :: c=(/4100.0,4100.0/),&
       rho_s=(/1024.0,1024.0/)

  real, dimension(:,:) :: gt, it, of
  real :: current_time, dt
  real :: CConN
  integer :: i

  CConN=ConN(of(1,1),gt(1,1),it,of(2,1))
  gt(1,3)=CConN
  of(1,3)=L_v*CConN/cc_p(gt,it)
  
  contains

    function cc_p(gas,solid)
      real, intent(in), dimension(:,:) :: gas, solid
      real :: cc_p
      integer :: n

      cc_p=c_p(1)
      do n=1,size(gas,1)
         cc_p=cc_p+(c_p(n+1)-c_p(1))*gas(n,1)
      end do
      do n=1,size(solid)
         cc_p=cc_p+(c(n+1)-c_p(1))*solid(n,1)
      end do
    end function cc_p
      

    function p_sat(T)
      real, intent(in):: T
      real  :: p_sat
      integer :: n
      p_sat=611.2*exp(17.62*(max(T,273.0)-273.0)/(max(T,273.0)-30.0))
    end function p_sat
      
    function q_sat(T,q_s,rho)
      real, intent(in) :: T,rho
      real, intent(in), dimension(:,:) :: q_s
      real  :: q_sat,irhog
      integer :: n

      irhog=1.0/rho
      do n=1,size(q_s)
         irhog=irhog-q_s(n,1)/rho_s(n)
      end do
      q_sat=p_sat(T)/(R_v*T)*(irhog)

    end function q_sat

    function ConN(T,q_v,q_s,rho)
      real, intent(in) :: T,q_v,rho
      real, intent(in), dimension(:,:) :: q_s
      real :: ConN, qq_sat
      integer :: n
      
      qq_sat=q_sat(T,q_s,rho)

      ConN=max((q_v-qq_sat)/dt,0.0)

    end function ConN
      
  end subroutine microphysics_main_pointwise
