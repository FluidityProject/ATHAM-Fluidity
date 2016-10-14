#include "fdebug.h"
!#define EXTRA_SPECIAL_GRADATION_DEBUGGING

module gradation_metric
!!< This module implements a gradation algorithm
!!< to ensure smooth mesh transitions from small
!!< to large edge lengths.
!!< The algorithm implemented is a modification of
!!< "Anisotropic mesh gradation control", Li et. al,
!!< 13th International Meshing Roundtable, 2004

  use fldebug
  use spud
  use vector_tools
  use sparse_tools, only: csr_sparsity, csr_matrix, CSR_INTEGER
  use unittest_tools
  use adjacency_lists
  use linked_lists
  use metric_tools
  use fields
  use vtk_interfaces
  use node_boundary
  use edge_length_module
  use field_derivatives
  use gradation_tools
  
  implicit none

  private
  public :: initialise_gradation_metric
  public :: form_gradation_metric
  public :: use_gradation_metric, gradation_initialised

  ! These are the DEFAULTS ONLY IF YOU DON'T CALL
  ! initialise_gradation_metric.
  ! initialise_gradation_metric changes them for real code.
  ! Bottom line: if you want to change whether gradation is used
  ! or the gradation constant, CHANGE THE VALUES in INITIALISE_GRADATION_METRIC
  logical :: use_gradation_metric = .false.
  logical :: gradation_initialised = .false.

  contains

  subroutine initialise_gradation_metric
  
    use_gradation_metric=have_option("/mesh_adaptivity/hr_adaptivity/enable_gradation")

    if (have_option("/mesh_adaptivity/hr_adaptivity/enable_gradation/gradation_parameter")) then
      call get_option("/mesh_adaptivity/hr_adaptivity/enable_gradation/gradation_parameter", gamma0)
    else
      gamma0 = 1.5
    end if

    ewrite(2,*) 'gradation: ', use_gradation_metric, gamma0
    
  end subroutine initialise_gradation_metric

  subroutine form_gradation_metric(positions, error_metric, noits)
    type(tensor_field), intent(inout) :: error_metric !!< The metric formed so far
    type(vector_field), intent(in) :: positions
    integer, intent(out), optional :: noits

    type(csr_matrix) :: nnlist !!< Node-node adjacency list
    type(csr_sparsity), pointer :: nn_sparsity !!< Node-node adjacency list sparsity
    type(elist) :: edgelist    !!< Linked list of edges
    type(mesh_type) :: mesh
    integer :: p, q !! the nodes
    real, dimension(error_metric%dim(1), error_metric%dim(2)) :: vec_P, vec_Q ! eigenvectors
    real, dimension(error_metric%dim(1)) :: val_P, val_Q ! the eigenvalues
    logical :: vals_changed_P, vals_changed_Q !!< have P or Q changed? If so need to reform
    logical :: vecs_changed_P, vecs_changed_Q !!< the metric and update any surrounding nodes.

    integer :: boundcount_P, boundcount_Q, expected_boundcount
    logical :: do_warp_directions

    integer :: dim, count

    integer :: global_its, end_marker ! count how many sweeps this involves

    type(scalar_field) :: edgelen
    integer, save :: adaptcnt = 0

#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
    integer :: stepcount
    type(scalar_field) :: nodefield
#endif
    logical :: debug_metric

    dim = error_metric%dim(1)
    debug_metric = have_option("/mesh_adaptivity/hr_adaptivity/debug/write_metric_stages")

    ewrite(2,*) "++: Applying gradation"

    mesh = error_metric%mesh
    domain_scale = domain_length_scale(positions)
    call initialise_boundcount(mesh, positions)
#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
    write(0,*) "domain_scale == ", domain_scale
#endif

    if (domain_is_2d()) then
      expected_boundcount = 1
    else
      expected_boundcount = 0
    end if

    !! Here I describe a convention I use in the nnlist.
    !! Normally the val array of the CSR matrix doesn't exist
    !! (as the value doesn't really matter). Here I set
    !! (i,j) to < 0 if (i,j) is NOT in the linked list of edges,
    !! and > 0 if it IS in the linked list of edges.
    !! |(i,j)| is the number of times it's been checked + 1.

    nn_sparsity => extract_nnlist(mesh)
    call allocate(nnlist, nn_sparsity, type=CSR_INTEGER)
    nnlist%ival = -1

    call construct_edge_list(mesh, nnlist, edgelist) 

    if (debug_metric) then
      call allocate(edgelen, error_metric%mesh, "Desired edge lengths")
    end if

    ! OK. So now we have the edge list.
#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
    stepcount = 0
    call allocate(nodefield, error_metric%mesh, "Node number")
    call get_node_field(error_metric%mesh, nodefield)
#endif

    end_marker = edgelist%length
    global_its = 0

    do while (edgelist%length /= 0)

      ! Count the number of sweeps through the mesh.
      end_marker = end_marker - 1
      if (end_marker == 0) then
        global_its = global_its + 1
        end_marker = edgelist%length
      end if

      call wrap_pop(nnlist, edgelist, p, q, count)                               ! fetch the nodes
      
#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
      write(0,*) "----------------------------------------------------"
      write(0,*) "stepcount == ", stepcount
      write(0,*) "(p, q) == (", p, ", ", q, ")" 
#endif

      call eigendecomposition_symmetric(node_val(error_metric, p), vec_P, val_P) ! decompose
      call eigendecomposition_symmetric(node_val(error_metric, q), vec_Q, val_Q)

#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
      call check_basis(vec_P)
      call check_basis(vec_Q)
      write(0,*) "Input P:"
      call write_matrix(node_val(error_metric, p), "P")
      call write_matrix(vec_P, "vec_P")
      call write_vector(val_P, "val_P")
      write(0,*) "Input Q:"
      call write_matrix(node_val(error_metric, q), "Q")
      call write_matrix(vec_Q, "vec_Q")
      call write_vector(val_Q, "val_Q")
#endif

      vecs_changed_P = .false. ; vecs_changed_Q = .false.

      boundcount_P = node_boundary_count(p)
      boundcount_Q = node_boundary_count(q)

       do_warp_directions = .false.

#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
      write(0,*) "boundcount_P == ", boundcount_P
      write(0,*) "boundcount_Q == ", boundcount_Q
      write(0,*) "do_warp_directions == ", do_warp_directions
#endif

      if (do_warp_directions) then
        call warp_directions(vec_P, val_P, vecs_changed_P, vec_Q, val_Q, vecs_changed_Q, distance(positions, p, q))
      end if

#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
      write(0,*) "P after warping:"
      call write_matrix(vec_P, "vec_P")
      write(0,*) "Q after warping:"
      call write_matrix(vec_Q, "vec_Q")
      write(0,*) "vecs_changed_P == ", vecs_changed_P, "vecs_changed_Q == ", vecs_changed_Q
      call check_basis(vec_P)
      call check_basis(vec_Q)
#endif


      vals_changed_P = .false. ; vals_changed_Q = .false.
      call reduce_edgelen(node_val(error_metric, p), vec_P, val_P, vals_changed_P, &
                          node_val(error_metric, q), vec_Q, val_Q, vals_changed_Q, &
                          positions, p, q)

#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
      write(0,*) "P after reducing:"
      call write_vector(val_P, "val_P")
      write(0,*) "Q after reducing:"
      call write_vector(val_Q, "val_Q")
      write(0,*) "vals_changed_P == ", vals_changed_P, "vals_changed_Q == ", vals_changed_Q
#endif

      if (vals_changed_P) then
        call eigenrecomposition(error_metric%val(:, :, p), vec_P, val_P)
        call tag_edges(nnlist, edgelist, p, q, count)
      end if
      if (vals_changed_Q) then
        call eigenrecomposition(error_metric%val(:, :, q), vec_Q, val_Q)
        call tag_edges(nnlist, edgelist, q, p, count)
      end if
      if (count <= max_rot_its) then ! honour directional changes for the first 4 sweeps
        if (vecs_changed_P) then
          call eigenrecomposition(error_metric%val(:, :, p), vec_P, val_P)
          call tag_edges(nnlist, edgelist, p, q, count)
        end if
        if (vecs_changed_Q) then
          call eigenrecomposition(error_metric%val(:, :, q), vec_Q, val_Q)
          call tag_edges(nnlist, edgelist, q, p, count)
        end if
      end if

#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
      call check_metric(error_metric)
      call get_edge_lengths(error_metric, edgelen)
      call vtk_write_fields("data/gradation_debug", stepcount, positions, positions%mesh, &
                            sfields=(/nodefield, edgelen/), tfields=(/error_metric/))
      stepcount = stepcount + 1
#endif
    end do

    call deallocate(nnlist)

    if (debug_metric) then
      call get_edge_lengths(error_metric, edgelen)
      call vtk_write_fields(trim("gradation_metric"), adaptcnt, positions, positions%mesh, &
                             sfields=(/edgelen/), tfields=(/error_metric/))
      call deallocate(edgelen)
      adaptcnt = adaptcnt + 1
    endif
#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
    call deallocate(nodefield)
#endif

  ewrite(2,*) "Finished gradation algorithm: global iterations == ", global_its
  if (present(noits)) noits = global_its
  end subroutine form_gradation_metric

end module gradation_metric
