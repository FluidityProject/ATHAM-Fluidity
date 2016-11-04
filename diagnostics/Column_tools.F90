#include "fdebug.h"

module column_module
  !!< Extrude a given 2D mesh to a full 3D mesh.
  !!< The layer depths are specified by a sizing function
  !!< which can be arbitrary python.
  use elements
  use fields
  use spud
  use quadrature
  use fetools
  use vector_tools
  use state_module
  use parallel_tools
  use node_ownership
  use global_parameters
  use pickers_inquire
  use interpolation_module
  use sparse_tools
  use vtk_interfaces
  use linked_lists
  use populate_state_module
  use hadapt_combine_meshes
  use hadapt_advancing_front
  use parallel_fields, only : node_owned
  implicit none

  private
  
  public :: compute_z_nodes_simple, create_column_mesh, create_basis, &
  	    interpolate_to_column

  contains
  
  subroutine create_column_mesh (h_mesh, top, nlevel, out_mesh)
    real, intent(in) :: top
    integer, intent(in) :: nlevel
    type(vector_field), intent(inout) :: h_mesh, out_mesh
    
    integer :: i, column, h_dim, quadrature_degree
    type(vector_field), dimension(:), allocatable :: c_mesh
    type(quadrature_type) :: quad
    type(element_type) :: full_shape
    real :: depth, constant_sizing, min_bottom_layer_frac
    
    allocate(c_mesh(1:node_count(h_mesh)))
    
    do column=1,node_count(h_mesh)
    
      depth=top-node_val(h_mesh,h_mesh%dim,column)
      constant_sizing=depth/real(nlevel)
      min_bottom_layer_frac=1.e-3
    
      call compute_z_nodes_simple(c_mesh(column), depth, node_val(h_mesh, column), &
    			min_bottom_layer_frac, 'bottom_up', constant_sizing)

    enddo
      
    ! Now the tiresome business of making a shape function.
    h_dim = mesh_dim(h_mesh)
    call get_option("/geometry/quadrature/degree", quadrature_degree)
    quad = make_quadrature(vertices=h_dim+2, dim=h_dim+1, degree=quadrature_degree)
    full_shape = make_element_shape(vertices=h_dim+2, dim=h_dim+1, degree=1, quad=quad)
    call deallocate(quad)

    ! combine the 1d vertical meshes into a full mesh
    call combine_z_meshes_simple(h_mesh, c_mesh, out_mesh, full_shape, "ExtrudedMesh", "")
       
    do column=1, node_count(h_mesh)
      call deallocate(c_mesh(column))
    end do
    call deallocate(full_shape)
    deallocate(c_mesh)
    			 
  end subroutine create_column_mesh
  
  subroutine interpolate_to_column (field, node, X, s_mesh, c_mesh, column)
    integer, intent(in) :: node
    type(scalar_field), intent(in) :: field
    type(vector_field), intent(inout) :: X, c_mesh, s_mesh
    real, dimension(:), intent(inout) :: column
    
    integer :: i, nnodes, ele
    integer, dimension(:), allocatable :: map
    real, dimension(:), allocatable :: local_coord
    real, dimension(:,:), allocatable :: basis, coords
    
    nnodes=node_count(c_mesh)
    allocate(map(1:nnodes), coords(1:X%dim,nnodes))
    allocate(basis(1:nnodes,1:ele_loc(X,1)), local_coord(X%dim+1))
    
    do i = 1, X%dim-1
      coords(i,:)=s_mesh%val(i,node)
    enddo

    do i = 1, nnodes
      coords(X%dim,i)=c_mesh%val(1,i)
      print*, 'before picker inquire', coords(:,i)
      
      if (X%dim == 2) then
        call picker_inquire(X, coordx=coords(1,i), coordy=coords(2,i), ele=ele, &
	      local_coord=local_coord, global=.false.)
      else if (X%dim == 3) then
        call picker_inquire(X, coordx=coords(1,i), coordy=coords(2,i), coordz=coords(3,i), ele=ele, &
	      local_coord=local_coord, global=.false.)
      endif
      map(i) = ele
    enddo
      
    print*, 'before create basis', map
!    call create_basis(X, c_mesh, map, basis)
    
    ! Remap at reference height
!    call interpolate_field_column (field, column, map, basis)

!    print*, 'after interpolate', column
    deallocate(map, coords, basis, local_coord)
			 
  end subroutine interpolate_to_column

  subroutine compute_z_nodes_simple(z_mesh, depth, xy, min_bottom_layer_frac, &
                                         direction, sizing)
    !!< Figure out at what depths to put the layers.
    type(vector_field), intent(out) :: z_mesh
    real, intent(in):: depth
    real, dimension(:), intent(in):: xy
    ! to prevent infinitesimally thin bottom layer if sizing function
    ! is an integer mulitple of total depth, the bottom layer needs
    ! to have at least this fraction of the layer depth above it.
    ! The recommended value is 1e-3.
    real, intent(in) :: min_bottom_layer_frac
    real, optional, intent(in):: sizing
    character(len=*), intent(in) :: direction
    ! this is a safety gap:
    integer, parameter:: MAX_VERTICAL_NODES=1e6

    integer :: elements
    real :: constant_value, delta_h, z

    type(rlist):: depths
    type(mesh_type) :: mesh
    type(element_type) :: oned_shape
    type(quadrature_type) :: oned_quad
    integer :: quadrature_degree
    integer :: ele
    integer, parameter :: loc=2
    integer :: node
    integer :: list_size

    call get_option("/geometry/quadrature/degree", quadrature_degree)
    oned_quad = make_quadrature(vertices=loc, dim=1, degree=quadrature_degree)
    oned_shape = make_element_shape(vertices=loc, dim=1, degree=1, quad=oned_quad)
    call deallocate(oned_quad)

    constant_value=sizing
    delta_h = constant_value

    ! Start the mesh at z=0 and work down to z=-depth.
    z=0.0
    call insert(depths, z)
    do
      if (trim(direction)=='top_to_bottom') then
        z=z-delta_h
        if (depth > 0.0) then
          if (z<-depth+min_bottom_layer_frac*delta_h) exit
        else
          if (z>-depth+min_bottom_layer_frac*delta_h) exit
        end if
      else if (trim(direction)=='bottom_up') then
        z=z+delta_h
        if (z > depth+min_bottom_layer_frac*delta_h) then
          exit
        endif
      endif
      
      call insert(depths, z)
      if (depths%length>MAX_VERTICAL_NODES) then
        ewrite(-1,*) "Check your extrude/sizing_function"
        FLExit("Maximum number of vertical layers reached")
      end if
    end do
    elements=depths%length-1

    call allocate(mesh, elements+1, elements, oned_shape, "ZMesh")
    do ele=1,elements
      mesh%ndglno((ele-1)*loc+1:ele*loc) = (/ele,ele+1/)
    end do

    call allocate(z_mesh, 1, mesh, "ZMeshCoordinates")
    call deallocate(mesh)
    call deallocate(oned_shape)

    call set(z_mesh, 1, (/0.0/))
    do node=1, elements+1
      call set(z_mesh, node,  (/ pop(depths) /))
    end do
    
    assert(oned_quad%refcount%count == 1)
    assert(oned_shape%refcount%count == 1)
    assert(z_mesh%refcount%count == 1)
    assert(mesh%refcount%count == 1)
      
  end subroutine compute_z_nodes_simple
  
  subroutine create_basis (position, remapped_position, map, basis)
  
    type(vector_field), intent(inout) :: position, remapped_position
    real, dimension(node_count(remapped_position),ele_loc(position,1)), intent(out) :: basis
    integer, dimension(node_count(remapped_position)), intent(in) :: map
  
    integer, dimension(:), pointer :: node_list
    integer :: info, ni, nl, nj, ele, node  
    real, dimension(position%dim, ele_loc(position,1)) :: X_val
    real, dimension(position%dim) :: X_mean, X_new, phi
    real, dimension(position%dim,position%dim) :: phimat
    
    do node =1, size(map)
    
    ele=map(node)
    node_list=>ele_nodes(position,ele)

    X_new=node_val(remapped_position,node)

    X_val=ele_val(position,ele)
    X_mean=sum(X_val,2)/size(X_val,2)

    !Construct interpolation basis skipping nj
    do nj = 1, size(node_list)
      basis(node,:)=0.0
      
      nl=0
      do ni=1,size(X_val,2)
    	if (ni==nj) cycle
    	nl=nl+1
    	phimat(:,nl)=X_val(:,ni)-X_mean
      enddo
      phi=X_new-X_mean

      !Solve for basis coefficients
      call solve(phimat,phi,info)

      nl=0
      do ni=1,size(X_val,2)
    	if (ni==nj) cycle
    	nl=nl+1
    	basis(node,ni)=phi(nl)
      enddo

      !If it's good, no need to search further
      if (minval(basis(node,:)) == 0.) exit
    enddo
    
    enddo
  
  end subroutine create_basis

  subroutine interpolate_field_column (field, remapped_field, map, basis)
  
    type(scalar_field), intent(in) :: field
    real, dimension(:), intent(inout) :: remapped_field
    real, dimension(:,:), intent(in) :: basis
    integer, dimension(:), intent(in) :: map
  
    integer, dimension(field%mesh%shape%numbering%vertices) :: ele_vertices
    real, dimension(ele_loc(field,1)) :: T_val
    real :: T_mean, T_new
    integer :: ele, node
    

    !Find vertices on old mesh
    ele_vertices=local_vertices(field%mesh%shape)
    remapped_field=0.0
    do node=1,size(map)
      ele=map(node)

      T_val=ele_val(field,ele)
      T_mean=sum(T_val(ele_vertices))/size(ele_vertices)

      remapped_field(node)=T_mean+sum(basis(node,:)*(T_val(ele_vertices)-T_mean))
    enddo
    
  end subroutine interpolate_field_column

  subroutine combine_z_meshes_simple(h_mesh, z_meshes, out_mesh, full_shape, mesh_name, option_path)
  !! Given the h_mesh and a z_mesh under each node of it combines these
  !! into a full horiz+vertic. mesh
  type(vector_field), intent(inout):: h_mesh
  type(vector_field), dimension(:), intent(in):: z_meshes
  type(vector_field), intent(out):: out_mesh
  type(element_type), intent(in):: full_shape
  character(len=*), intent(in):: mesh_name, option_path
  
    type(csr_sparsity):: out_columns
    type(mesh_type):: mesh
    integer, dimension(:), allocatable:: no_hanging_nodes
    integer:: column, total_out_nodes, total_out_elements, z_elements, last_seen
    
    allocate(no_hanging_nodes(1:node_count(h_mesh)))
    no_hanging_nodes=0
    do column=1, size(z_meshes)
      if (node_owned(h_mesh, column)) then
        no_hanging_nodes(column)=ele_count(z_meshes(column))
      end if
    end do
    if (associated(h_mesh%mesh%halos)) then
      call halo_update(h_mesh%mesh%halos(2), no_hanging_nodes)
    end if
    
    total_out_nodes = 0
    total_out_elements = 0
    do column=1, size(z_meshes)
      z_elements = no_hanging_nodes(column)
      total_out_nodes = total_out_nodes + z_elements + 1
      assert(associated(h_mesh%mesh%adj_lists))
      assert(associated(h_mesh%mesh%adj_lists%nelist))
      total_out_elements = total_out_elements + z_elements * row_length(h_mesh%mesh%adj_lists%nelist, column)
    end do

    call allocate(mesh, total_out_nodes, total_out_elements, full_shape, mesh_name)
    ! allocate mapping between extruded nodes to surface node (column number)
    ! it lies under
    allocate(mesh%columns(total_out_nodes))
    mesh%columns = 0
    ! allocate mapping between extruded elements to surface elements in horizontal mesh
    ! it lies under
    allocate(mesh%element_columns(total_out_elements))
    mesh%element_columns = 0
    ! if the horizontal mesh has region ids these can be propagated down into the full
    ! mesh.  allocate space for this...
    if(associated(h_mesh%mesh%region_ids)) then
      allocate(mesh%region_ids(total_out_elements))
    end if
    
    call allocate(out_mesh, mesh_dim(h_mesh)+1, mesh, trim(mesh_name)//"Coordinate")
    call deallocate(mesh)
    
    out_mesh%mesh%option_path=option_path
    out_mesh%option_path=""
    out_mesh%mesh%periodic = mesh_periodic(h_mesh)

    last_seen = 0
    do column=1,node_count(h_mesh)
      if (node_owned(h_mesh, column)) then
        call append_to_structures(column, z_meshes(column), h_mesh, out_mesh, last_seen)
      else
        out_mesh%mesh%columns(last_seen+1:last_seen+no_hanging_nodes(column)+1) = column
        last_seen = last_seen + no_hanging_nodes(column)+1
      end if
    end do
    assert(all(out_mesh%mesh%columns>0))
      
    call create_columns_sparsity(out_columns, out_mesh%mesh)
    
    if (associated(h_mesh%mesh%halos)) then
      ! derive l2 node halo for the out_mesh
      call derive_extruded_l2_node_halo(h_mesh%mesh, out_mesh%mesh, out_columns)
      ! positions in the non-owned columns can now simply be halo-updated
      call halo_update(out_mesh)
    end if
      
    call generate_layered_mesh(out_mesh, h_mesh, direction='bottom_up')

    call deallocate(out_columns)
    
  end subroutine combine_z_meshes_simple
    
end module column_module
