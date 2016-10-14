import os
path=os.path.dirname(__file__) 

def MakeWrapperFiles(field_dict,call_str,pointwise):
    write_to_file(field_dict,call_str,pointwise)

def allocate_str(field_dict):
    s="""  subroutine allocate_storage(number_of_tracers,n)
    integer :: n
!f2py integer, intent(hide), depend(number_of_tracers) :: n=shape(number_of_tracers,0)
    integer :: number_of_tracers(n)
"""
    for n,k in enumerate(field_dict):
        s+="    if (allocated(%s)) deallocate(%s)\n"%(k,k)
        s+="    allocate(%s(number_of_tracers(%d)))\n"%(k,n+1)
    s+="  end subroutine allocate_storage\n\n"
    return s

def finalize_str(field_dict):
    s="  subroutine finalize\n"
    for k in field_dict:
        s+="deallocate(%s)\n"%k
    s+="  end subroutine finalize\n\n"
    return s

def set_field_str(fname):
    s="""  subroutine set_%s(i,new_val,n,old_val,source,m)
    integer :: m,n
!f2py integer, intent(hide), depend(new_val) :: n=shape(new_val,0)
!f2py integer, intent(hide), depend(new_val) :: m=shape(source,0)
    real, intent(in), dimension(n), target :: new_val, old_val
    real, intent(in), dimension(n), target, optional ::source
!f2py real, intent(inplace), dimension(n) :: new_val, old_val
!f2py real, intent(inplace), dimension(n), optional :: source
    integer :: i
    %s(i)%%new=>new_val
    %s(i)%%old=>old_val
    print*, present(source), m
    if (present(source) .and. m==n)&
       %s(i)%%source=>source
  end subroutine set_%s
"""%(fname,fname,fname,fname,fname)
    return s

def run_str(field_dict,call_string):
    s="""subroutine run_microphysics(current_time,dt)
    real, intent(in) :: current_time, dt
    interface
       subroutine %s(time,timestep"""%call_string
    for n,k in enumerate(field_dict):
        s+=',&\n            t%d'%n
    s+=')\n'
    s+='         use FW_data_type\n'
    s+='        real, intent(in) :: time, timestep\n'
    for n,k in enumerate(field_dict):
        s+='         type(basic_scalar), intent(inout), dimension(:) :: t%d\n'%n
    s+='       end subroutine %s\n'%call_string
    s+="""    end interface

    call %s(current_time,dt"""%call_string
    for k in field_dict:
        s+=',&\n         %s'%k
    s+=')\n\n'
    s+='  end subroutine run_microphysics\n\n'
    return s

def run_str_pointwise(field_dict,call_string):
    s="""subroutine run_microphysics(current_time,dt)
    real, intent(in) :: current_time, dt
    integer :: i,j\n

"""
    for n,k in enumerate(field_dict):
        s+='    real, dimension(size(%s),3) :: tracer%d\n'%(k,n)
    s+="""    interface\n
       subroutine %s(time,timestep"""%call_string
    for n,k in enumerate(field_dict):
        s+=',&\n            t%d'%n
    s+=')\n'
    s+='         use FW_data_type\n'
    s+='        real, intent(in) :: time, timestep\n'
    for n,k in enumerate(field_dict):
        s+='         real, intent(inout), dimension(:,:) :: t%d\n'%n
    s+='       end subroutine %s\n'%call_string
    s+="    end interface\n"
    s+="    do i=1, size(%s(0)%%new)\n"%(field_dict.keys()[0])

    for n,k in enumerate(field_dict):
        s+='    do j=1,size(%s)\n'%k
        s+='      tracer%d(j,1)=%s(j)%%new(i)\n'%(n,k)
        s+='      tracer%d(j,2)=%s(j)%%old(i)\n'%(n,k)
        s+='      if (associated(%s(j)%%source))&\n        tracer%d(j,3)=%s(j)%%source(i)\n'%(k,n,k)
        s+='    end do\n\n'
    s+="    call %s(current_time,dt"%call_string
    for k in range(len(field_dict)):
        s+=',&\n     tracer%d'%n
    s+=')\n\n'
    for n,k in enumerate(field_dict):
        s+='    do j=1,size(%s)\n'%k
        s+='      %s(j)%%new(i)=tracer%d(j,1)\n'%(k,n)
        s+='      %s(j)%%old(i)=tracer%d(j,2)\n'%(k,n)
        s+='      if (associated(%s(j)%%source))&\n        %s(j)%%source(i)=tracer%d(j,3)\n'%(k,k,n)
        s+='     end do\n\n'
    s+='    end do\n\n'
    s+='  end subroutine run_microphysics\n\n'
    return s

    
def write_to_file(field_dict={},
                  call_string='',
                  pointwise=False,
                  dirname=path+'/src',
                  src_name='FW_auto',
                  data_name='FW_data'):
    
    f=open(dirname+'/'+src_name+'.F90','w')

    s="""module FW_auto
    use FW_data
  implicit none
contains
"""
    f.write(s)
    f.write(allocate_str(field_dict))
    f.write(finalize_str(field_dict))
    
    for k in field_dict:
        f.write(set_field_str(k))

    if pointwise:
        f.write(run_str_pointwise(field_dict,call_string))
    else:
        f.write(run_str(field_dict,call_string))
    f.write("end module FW_Auto\n")
    f.close()      

    f=open(dirname+'/'+data_name+'.F90','w')
    f.write("""module %s
  use FW_data_type
"""%data_name)
    for k in field_dict:
        f.write('    type(basic_scalar), dimension(:), allocatable :: %s\n'%k)
    f.write('end module %s\n'%data_name)
    f.close()    
