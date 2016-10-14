#include "fdebug.h"

module hadapt_extrude_radiation
  !!< Extrude a given 2D mesh to a full 3D mesh.
  !!< The layer heights are specified by a sizing function
  !!< which can be arbitrary python.
  use elements
  use fields
  use spud
  use quadrature
  use global_parameters
  use fields_data_types
  use sparse_tools
  use vtk_interfaces
  use linked_lists
  use hadapt_combine_meshes
  use halos
  
  implicit none

  private
  public :: extrude_mesh_from_surface

  contains

  subroutine extrude_mesh_from_surface(h_mesh, option_path, out_mesh)
    !!< The horizontal 2D mesh.
    !!< Note: this must be linear.
    type(vector_field), intent(inout) :: h_mesh
    !!< options to be set for out_mesh,
    !!< at the moment: /name, and under from_mesh/extrude/:
    !!< height, sizing_function optionally top_surface_id and bottom_surface_id
    character(len=*), intent(in) :: option_path
    !!< The full extruded 3D mesh.
    type(vector_field), intent(out) :: out_mesh

    character(len=FIELD_NAME_LEN):: mesh_name, file_name  
    type(quadrature_type) :: quad
    type(element_type) :: full_shape
    type(vector_field) :: constant_z_mesh
    type(vector_field), dimension(:), allocatable :: z_meshes
    logical :: height_from_python, height_is_constant, surface_is_constant
    real:: height
    integer:: h_dim, column, quadrature_degree

    logical :: sigma_layers
    integer :: number_sigma_layers
    
    character(len=PYTHON_FUNC_LEN) :: height_function
    integer :: n_regions, r, dim, number_levels
    integer, dimension(:), allocatable :: region_ids
    logical :: apply_region_ids, constant_z_mesh_initialised
    integer, dimension(node_count(h_mesh)) :: visited
    logical, dimension(node_count(h_mesh)) :: column_visited
    real, dimension(h_mesh%dim) :: position
 
    ewrite(1,*) 'In extrude_mesh_from_surface'   

    allocate(z_meshes(node_count(h_mesh)))
    
    dim=h_mesh%dim

    call add_nelist(h_mesh%mesh)
    
    n_regions = option_count(trim(option_path)//'/from_mesh/extrude/regions')
    if(n_regions==0) then
      ewrite(-1,*) "Since r13369 it has been possible to extrude using different parameters"
      ewrite(-1,*) "in each region id of a horizontal mesh.  This means that the extrusion parameters"
      ewrite(-1,*) "must be included under geometry/mesh::MeshNameHere/from_mesh/extrude/regions."
      ewrite(-1,*) "This is necessary even if you want to use the same parameters"
      ewrite(-1,*) "for the whole mesh (using a single regions option and no region_ids)."
      ewrite(-1,*) "I've been told to extrude but have found no regions options."
      ewrite(-1,*) "Have you updated your flml since r13369?"
      FLExit("No regions options found under extrude.")
    elseif(n_regions<0) then
      FLAbort("Negative number of regions options found under extrude.")
    end if
    apply_region_ids = (n_regions>1)
    visited = 0 ! a little debugging check - can be removed later
   
    surface_is_constant=(maxval(h_mesh%val(dim,:))==minval(h_mesh%val(dim,:)))
    column_visited = .false.
    
    do r = 0, n_regions-1
      
      constant_z_mesh_initialised = .false.
      
      
      call get_extrusion_options(option_path, r, apply_region_ids, region_ids, &
                                 height_is_constant, height, height_from_python, height_function, &
                                 file_name, number_levels)
				       
      ! create a 1d vertical mesh under each surface node
      do column=1, node_count(h_mesh)
        position=h_mesh%val(:,column)
      
        ! decide if this column needs visiting...
        if(skip_column_extrude(h_mesh%mesh, column, &
                              apply_region_ids, column_visited(column), region_ids, &
                              visited_count = visited(column))) cycle
        
        if(height_is_constant.and.surface_is_constant) then
          if (.not. constant_z_mesh_initialised) then
            call compute_z_nodes_wrap_rad(dim, constant_z_mesh, position,	 &
                            height_is_constant, height, height_from_python, height_function, 	 &
                            number_levels)
            constant_z_mesh_initialised = .true.
          end if
          call get_previous_z_nodes(z_meshes(column), constant_z_mesh)
        else
          call compute_z_nodes_wrap_rad(dim, z_meshes(column), position,	 &
                            height_is_constant, height, height_from_python, height_function, 	 &
                            number_levels)
        end if
      end do
      
      if(apply_region_ids) deallocate(region_ids)
      
      if (constant_z_mesh_initialised) then
        call deallocate(constant_z_mesh)
      end if
    
    end do
    
#ifdef DDEBUG
    if(apply_region_ids) then
      ewrite(2,*) "Maximum number of times a node was visited: ", maxval(visited)
      ewrite(2,*) "Minimum number of times a node was visited: ", minval(visited)
      if(.not.isparallel()) then
        assert(minval(visited)>0)
      end if
    end if
#endif
      
    ! Now the tiresome business of making a shape function.
    h_dim = mesh_dim(h_mesh)
    call get_option("/geometry/quadrature/degree", quadrature_degree)
    quad = make_quadrature(vertices=h_dim + 2, dim=h_dim + 1, degree=quadrature_degree)
    full_shape = make_element_shape(vertices=h_dim + 2, dim=h_dim + 1, degree=1, quad=quad)
    call deallocate(quad)

    call get_option(trim(option_path)//'/name', mesh_name)

    ! combine the 1d vertical meshes into a full mesh
    call combine_z_meshes(h_mesh, z_meshes, out_mesh, full_shape, mesh_name, option_path, sigma_layers)

    do column=1, node_count(h_mesh)
      if (.not. node_owned(h_mesh, column)) cycle
      call deallocate(z_meshes(column))
    end do
    call deallocate(full_shape)
    deallocate(z_meshes)
    
  end subroutine extrude_mesh_from_surface

  subroutine get_extrusion_options(option_path, region_index, apply_region_ids, region_ids, &
                                   height_is_constant, height, height_from_python, height_function, & 
                                   file_name, number_levels)

    character(len=*), intent(in) :: option_path
    integer, intent(in) :: region_index
    logical, intent(in) :: apply_region_ids
    
    integer, dimension(:), allocatable :: region_ids
    
    logical, intent(out) :: height_is_constant, height_from_python
    real, intent(out) :: height
    character(len=PYTHON_FUNC_LEN), intent(out) :: height_function

    character(len=FIELD_NAME_LEN), intent(out) :: file_name
    
    integer, intent(out) :: number_levels
    
    integer, dimension(2) :: shape_option
    integer :: stat

    if(apply_region_ids) then
      shape_option=option_shape(trim(option_path)//"/from_mesh/extrude/regions["//int2str(region_index)//"]/region_ids")
      allocate(region_ids(1:shape_option(1)))
      call get_option(trim(option_path)//"/from_mesh/extrude/regions["//int2str(region_index)//"]/region_ids", region_ids)
    end if

    ! get the extrusion options
    height_from_python=.false.
    call get_option(trim(option_path)//&
                    '/from_mesh/extrude/regions['//int2str(region_index)//&
                    ']/top_height/constant', &
                      height, stat=stat)

    if (stat==0) then
      height_is_constant = .true.
    else
      height_is_constant = .false.
      call get_option(trim(option_path)//&
                      '/from_mesh/extrude/regions['//int2str(region_index)//&
                      ']/top_height/python', &
                       height_function, stat=stat)
      if (stat==0) height_from_python = .true.
      if (stat /= 0) then
        FLAbort("Unknown way of specifying bottom height function in mesh extrusion")
      end if
    end if
  
    call get_option(trim(option_path)//&
                    '/from_mesh/extrude/regions['//int2str(region_index)//&
                    ']/vertical_levels', &
                    number_levels, default=2)
  
  end subroutine get_extrusion_options

  subroutine compute_z_nodes_wrap_rad(dim, z_mesh, xy, height_is_constant, height, 	&
  				      height_from_python, height_function, number_levels)

    integer, intent(in) :: dim
    type(vector_field), intent(out) :: z_mesh
    real, dimension(dim), intent(in) :: xy
    logical, intent(in) :: height_is_constant, height_from_python
    real, intent(in) :: height
    character(len=*), intent(in) :: height_function
    integer, intent(in) :: number_levels

    real, dimension(1) :: tmp_height
    real, dimension(dim, 1) :: tmp_pos
    real :: lheight, depth
    
    if(height_is_constant) then
      lheight = height
    else 
      tmp_pos(:,1) = xy
      if (height_from_python) then
        call set_from_python_function(tmp_height, trim(height_function), tmp_pos, time=0.0)
        lheight = tmp_height(1)
      else
        FLAbort("Unknown way of specifying the bottom_height.")
      end if
    end if
    
    call compute_z_nodes_rad(dim, z_mesh, height_is_constant, lheight, xy, number_levels)
    
  end subroutine compute_z_nodes_wrap_rad

  subroutine get_previous_z_nodes(z_mesh, z_mesh_previous)
    type(vector_field), intent(inout) :: z_mesh, z_mesh_previous
    z_mesh = z_mesh_previous
    call incref(z_mesh)
  end subroutine get_previous_z_nodes

  subroutine compute_z_nodes_rad(dim, z_mesh, height_is_constant, height, xy, number_levels)
    !!< Figure out at what heights to put the layers.
    integer, intent(in) :: dim
    type(vector_field), intent(out) :: z_mesh
    real, intent(in):: height
    real, dimension(dim), intent(in):: xy
    logical, intent(in) :: height_is_constant
    ! to prevent infinitesimally thin bottom layer if sizing function
    ! is an integer mulitple of total height, the bottom layer needs
    ! to have at least this fraction of the layer height above it.
    ! The recommended value is 1e-3.
    integer, optional, intent(in) :: number_levels
    ! this is a safety gap:
    integer, parameter:: MAX_VERTICAL_NODES=1e6

    integer :: elements
    real :: constant_value

    type(rlist):: heights
    type(mesh_type) :: mesh
    type(element_type) :: oned_shape
    type(quadrature_type) :: oned_quad
    integer :: quadrature_degree
    integer :: ele
    integer, parameter :: loc=2
    integer :: node
    real, dimension(1:size(xy)+1):: xyz
    real :: delta_h, z, depth, tolerance
    character(len=PYTHON_FUNC_LEN) :: py_func
    
    tolerance = 1.	! Tolerance +/-1m

    call get_option("/geometry/quadrature/degree", quadrature_degree)
    oned_quad = make_quadrature(vertices=loc, dim=1, degree=quadrature_degree)
    oned_shape = make_element_shape(vertices=loc, dim=1, degree=1, quad=oned_quad)
    call deallocate(oned_quad)

    if (number_levels/=2) then
      depth=height-xy(size(xy))
      constant_value=depth/float(number_levels)
      py_func = " "
    else
      FLAbort("Need to supply a reasonable number of vertical levels!")
    end if

    ! first size(xy) coordinates remain fixed, 
    ! the last entry will be replaced with the appropriate height
    xyz(1:dim)=xy
    z=xyz(dim)
    node=2
    call insert(heights, z)
    
    do
      xyz(dim+1)=height
      if (height_is_constant) then
        delta_h = constant_value
      else
        FLAbort("Oops! The python interface to compute the vertical grid size does not exist yet!")
!        delta_h = get_delta_h(xyz, constant_value, py_func)
      endif 
      
      z=z+delta_h
      if (z > height+tolerance) then
        exit
      else if (z>=height-tolerance .and. z<=height+tolerance) then
        z=height
      endif
      
      call insert(heights, z)
      if (heights%length>MAX_VERTICAL_NODES) then
        ewrite(-1,*) "Check your extrude/sizing_function"
        FLExit("Maximum number of vertical layers reached")
      end if
    end do
    elements=heights%length-1

    ! Allocate and assign vertical mesh
    call allocate(mesh, elements+1, elements, oned_shape, "ZMesh")
    do ele=1,elements
      mesh%ndglno((ele-1) * loc + 1: ele*loc) = (/ele, ele+1/)
    end do
    call allocate(z_mesh, 1, mesh, "ZMeshCoordinates")
    call deallocate(mesh)
    call deallocate(oned_shape)

    call set(z_mesh, 1, (/0.0/))
    do node=1, elements+1
      call set(z_mesh, node,  (/ pop(heights) /))
    end do
    
    assert(oned_quad%refcount%count == 1)
    assert(oned_shape%refcount%count == 1)
    assert(z_mesh%refcount%count == 1)
    assert(mesh%refcount%count == 1)
      
  end subroutine compute_z_nodes_rad

  logical function skip_column_extrude(horizontal_mesh, column, &
                                       apply_region_ids, column_visited, region_ids, &
                                       visited_count)
    !!< this function decides if a column need extruding or not
    type(mesh_type), intent(in) :: horizontal_mesh
    integer, intent(in) :: column
    logical, intent(in) :: apply_region_ids
    logical, intent(inout) :: column_visited
    integer, dimension(:), intent(in) :: region_ids
    integer, intent(inout), optional :: visited_count
    
    integer, dimension(:), pointer :: eles
    logical :: node_in_region
    integer :: rs
    
    skip_column_extrude = .false.
    if(.not.node_owned(horizontal_mesh, column)) then
      skip_column_extrude = .true.
      return
    end if
    
    ! need to work out here if this column is in one of the current region ids!
    ! this is a bit arbitrary since nodes belong to multiple regions... therefore
    ! the extrusion height had better be continuous across region id boundaries!
    if(apply_region_ids) then
      if(column_visited) then
        skip_column_extrude = .true.
        return
      end if
      eles => node_neigh(horizontal_mesh, column)
      node_in_region = .false.
      region_id_loop: do rs = 1, size(region_ids)
        if(any(region_ids(rs)==horizontal_mesh%region_ids(eles))) then
          node_in_region = .true.
          exit region_id_loop
        end if
      end do region_id_loop
      if(.not.node_in_region) then
        skip_column_extrude = .true.
        return
      end if
      column_visited=.true.
      if(present(visited_count)) then
        visited_count = visited_count + 1
      end if
    end if
  
  end function skip_column_extrude

    
end module hadapt_extrude_radiation
