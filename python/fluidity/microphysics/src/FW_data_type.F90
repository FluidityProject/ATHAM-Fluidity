module FW_data_type
  type basic_scalar
     sequence
     real, dimension(:), pointer :: new=>null()
     real, dimension(:), pointer :: old=>null()
     real, dimension(:), pointer :: source=>null()
  end type basic_scalar

  type basic_vector
     sequence
     real, dimension(:,:), pointer :: new=>null()
     real, dimension(:,:), pointer :: old=>null()
     real, dimension(:,:), pointer :: source=>null()
  end type basic_vector

  type basic_tensor
     sequence
     real, dimension(:,:,:), pointer :: new=>null()
     real, dimension(:,:,:), pointer :: old=>null()
     real, dimension(:,:,:), pointer :: source=>null()
  end type basic_tensor

end module FW_data_type

