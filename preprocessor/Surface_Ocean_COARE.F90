!    Copyright (C) 2007 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

module Surface_Ocean_COARE

use spud
use fldebug
use elements
use fields
use field_options
use state_module
use vector_tools
use boundary_conditions
use equation_of_state
use interpolation_module
use node_boundary
use node_ownership

  IMPLICIT NONE

  real, parameter :: tdk=273.16, dter=0.3, visw=1.e-6, tcw=0.6, be=0.026, rhow=1000., rich=0.65, zi=900., a=0.018, b=0.729, Beta=1.2, von=0.4, fdg=1.00, pi=3.141593

  private
  
  public :: get_COARE3_bc_values, set_COARE3_bc_value
  
  interface set_COARE3_bc_value
    module procedure set_COARE3_bc_value_scalar, set_COARE3_bc_value_vector
  end interface set_COARE3_bc_value

  contains

  subroutine get_COARE3_bc_values(state, bc_ind, bc_name, bc_path, surface_mesh, surface_element_list)
    implicit none

    type(state_type), intent(inout) :: state
    type(mesh_type), intent(inout) :: surface_mesh
    character(len=*), intent(in):: bc_path, bc_name
    integer, dimension(:), pointer, intent(in) :: surface_element_list
    integer, intent(in) :: bc_ind
        
    type(mesh_type), pointer :: mesh
    type(scalar_field) :: temperature, qc_p, qc_v, eosdensity
    type(scalar_field), target  :: dummyscalar
    type(vector_field), target  :: dummyvector
    type(scalar_field), pointer :: pt, pressure
    type(scalar_field), pointer :: sea_surface_t, sensible_flux, latent_flux, rain_flux
    type(vector_field), pointer :: momentum_stress, wind
    
    integer :: j, f, nfields, stat
    real :: ds, ps, cp, cv, pcom, hs, hl, rf, tem
     
    ewrite(1,*) 'Start get_COARE3_bc_values'
    
    ! Extract surface fields
    pressure=>extract_scalar_field(state,'Pressure')
    wind=>extract_vector_field(state,'Velocity')
    
    call get_thermo_variable(state, pt)
    
    mesh=>extract_mesh(state, trim(pt%mesh%name))

    ! get thermodynamic fields
    call allocate(qc_p,mesh,name="CP")
    call allocate(qc_v,mesh,name="CV")
    call allocate(temperature,mesh,name="Temperature")
    call allocate(eosdensity,mesh,name="EOSDensity")
    call allocate(dummyscalar,mesh,name="DummyScalar")
    call allocate(dummyvector,wind%dim,mesh,name="DummyVector")
    
    call zero(dummyscalar)

    call compressible_eos(state, density=eosdensity, temperature=temperature, qc_p=qc_p, qc_v=qc_v)

    sea_surface_t => extract_surface_field(pt, bc_name, "sea_surface_temperature")
    sensible_flux => extract_surface_field(pt, bc_name, "sensible_heat_flux",stat=stat)
    if (stat /= 0) sensible_flux => dummyscalar
    latent_flux => extract_surface_field(pt, bc_name, "latent_heat_flux",stat=stat)
    if (stat /= 0) latent_flux => dummyscalar
    rain_flux => extract_surface_field(pt, bc_name, "rain_heat_flux",stat=stat)
    if (stat /= 0) rain_flux => dummyscalar
    momentum_stress => extract_surface_field(wind, bc_name, "momentum_stress",stat=stat)
    if (stat /= 0) momentum_stress => dummyvector
    
    ! Run the COARE algorithm.
    call COARE_interface (state, mesh, surface_mesh, surface_element_list, bc_name, bc_path, 	   &
    			  wind, pressure, eosdensity, temperature, qc_p, qc_v,	  		   &
			  sea_surface_t, sensible_flux, latent_flux, rain_flux, 		   &
			  momentum_stress )

    ewrite_minmax(sensible_flux)
    ewrite_minmax(latent_flux)
    ewrite_minmax(rain_flux)
    ewrite_minmax(momentum_stress)
    
    call deallocate(dummyscalar)
    call deallocate(dummyvector)
    call deallocate(temperature)
    call deallocate(eosdensity)
    call deallocate(qc_p)
    call deallocate(qc_v)
    
    ewrite(1,*) 'End get_COARE3_bc_values'
  
  end subroutine get_COARE3_bc_values
  
  subroutine set_COARE3_bc_value_scalar(state, surface_field, bc_name, field_name, surface_mesh, surface_node_list)
    implicit none

    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: surface_field
    type(mesh_type), pointer, intent(in) :: surface_mesh
    integer, dimension(:), pointer, intent(in) :: surface_node_list
    character(len=*), intent(in):: field_name, bc_name
  
    integer :: j, node
    real :: ps, ds, cp, cv, hs, hl, rf, pcom, tem
    type(scalar_field), pointer :: thermal
    type(scalar_field), pointer :: sensible_flux, latent_flux, rain_flux  
    type(scalar_field) :: temperature, eosdensity, eospressure, qc_p, qc_v
    
    ewrite(1,*) 'Setting COARE boundary conditions for scalar: ', trim(field_name), ' ', trim(bc_name)

    call get_thermo_variable(state,thermal)
    
    call allocate(temperature,thermal%mesh,'Temperature')
    call allocate(eosdensity,thermal%mesh,'EOSDensity')
    call allocate(eospressure,thermal%mesh,'EOSPressure')
    call allocate(qc_p,thermal%mesh,'EOSCP')
    call allocate(qc_v,thermal%mesh,'EOSCV')
    
    call compressible_eos(state, pressure=eospressure, density=eosdensity, temperature=temperature, qc_p=qc_p, qc_v=qc_v)

    sensible_flux => extract_surface_field(thermal, bc_name, "sensible_heat_flux")
    latent_flux => extract_surface_field(thermal, bc_name, "latent_heat_flux")
    rain_flux => extract_surface_field(thermal, bc_name, "rain_heat_flux")

    assert(node_count(sensible_flux) == size(surface_node_list))

    do j=1,size(surface_node_list)
      node=surface_node_list(j)
      
      ps=node_val(eospressure,node)
      ds=node_val(eosdensity,node)
      tem=node_val(temperature,node)
      cp=node_val(qc_p,node)
      cv=node_val(qc_v,node)
      
      hs=node_val(sensible_flux,j)
      hl=node_val(latent_flux,j)
      rf=node_val(rain_flux,j)
      
      hs=hs/(ds*cp)
      hl=hl/(ds*cal_flv(tem))
      rf=rf/(ds*cal_flv(tem))
    
      if (trim(field_name) == "PotentialTemperature") then        
        pcom=(ps/1.e+05)**((cp-cv)/cp)
    	call set(surface_field, j, hs/pcom)
      else if (trim(field_name)=="ConservedPotentialTemperature") then        
        pcom=(ps/1.e+05)**((cp-cv)/cp)
    	call set(surface_field, j, ds*hs/pcom)
      else if (trim(field_name)=="Temperature") then
    	call set(surface_field, j, hs+hl+rf)
      else if (trim(field_name)=="TotalWaterQ" .or. trim(field_name)=="VapourWaterQ") then
    	call set(surface_field, j, hl-rf)
      end if
    enddo
   
    call deallocate(temperature)
    call deallocate(eosdensity)
    call deallocate(eospressure)
    call deallocate(qc_p)
    call deallocate(qc_v)
  
  end subroutine set_COARE3_bc_value_scalar
  
  subroutine set_COARE3_bc_value_vector(state, surface_field, bc_name, field_name, surface_mesh, surface_node_list)
    implicit none

    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: surface_field
    type(mesh_type), pointer, intent(in) :: surface_mesh
    integer, dimension(:), pointer, intent(in) :: surface_node_list
    character(len=*), intent(in):: field_name, bc_name
  
    integer :: i, j, node
    real :: ds    
    real, dimension(surface_field%dim) :: tau
    type(vector_field), pointer :: reference_wind, momentum_stress  
    type(scalar_field) :: eosdensity
    
    ewrite(1,*) 'Setting COARE boundary conditions for vector: ', trim(bc_name)

    reference_wind=>extract_vector_field(state,'Velocity')
    
    call allocate(eosdensity,reference_wind%mesh,'EOSDensity')
    
    call compressible_eos(state, density=eosdensity)

    momentum_stress => extract_surface_field(reference_wind, bc_name, "momentum_stress")
    
    assert(node_count(momentum_stress) == size(surface_node_list))

    if (trim(field_name)=="Velocity") then      
    
    do i=1,size(surface_node_list)
      node=surface_node_list(i)
      
      ds=node_val(eosdensity,node)
      tau=node_val(momentum_stress, i)/ds
	 
      call addto(surface_field, i, tau)
    enddo
        
    endif
    
    call deallocate(eosdensity)

  end subroutine set_COARE3_bc_value_vector
  
  subroutine COARE_interface (state, mesh, surface_mesh, surface_element_list, bc_name, bc_path, &
      wind, pressure, density, temperature, qc_p, qc_v, sea_surface_t, &
      sensible_flux, latent_flux, rain_flux, momentum_stress )
      
    implicit none
    
    type(state_type), intent(inout) :: state
    type(mesh_type), intent(inout) :: mesh, surface_mesh
    type(scalar_field), intent(in) :: pressure, density, temperature, qc_p, qc_v
    type(vector_field), intent(in) :: wind
    type(scalar_field), intent(inout) :: sea_surface_t, sensible_flux, latent_flux, rain_flux
    type(vector_field), intent(inout) :: momentum_stress
    character(len=*), intent(in):: bc_name, bc_path
    integer, dimension(:), pointer, intent(in) :: surface_element_list
   
    logical :: cool_skin, warm_layer, wave
    integer, dimension(:), pointer:: surface_nodes
    integer :: i, j, jcool, jwave, ele
    real :: lon, lat, dt, current_time, mps_to_mmph, rhow=1000.
    real :: sst, hs, hl, rf
    real :: zref, ps, ts, qs, ds, rs, rl, cp, cv, pcom, unorm
    real, dimension(momentum_stress%dim) :: us, tau
    real, dimension(node_count(surface_mesh),ele_loc(mesh,1)) :: basis
    integer, dimension(node_count(surface_mesh)) :: map

    type(mesh_type) :: X_surface_mesh
    type(scalar_field), pointer :: q_r, q_r_vel, q_w, sw, lw, sfield
    type(vector_field), pointer :: position, vfield
    type(scalar_field), target :: dummyscalar
    type(scalar_field) :: precipitation, surface_pressure, surface_short_wave, surface_long_wave, surface_precip, surface_qr, surface_qr_vel
    type(scalar_field) :: reference_height_temperature, reference_height_moisture, reference_height_density, reference_height_cp, reference_height_cv
    type(vector_field) :: reference_height_wind
    type(vector_field) :: virtual_bc_position
    integer :: f, nfields
    
    ewrite(1,*) 'Start COARE3_interface'
    
    mps_to_mmph=3.6e6
    zref=10.
    jcool=0
    jwave=0
      
    ! Microphysics properties
    call allocate(dummyscalar,mesh,name="DummyScalar")
    call zero(dummyscalar)
    
    if (has_scalar_field(state,"TotalWaterQ")) then
      q_w=>extract_scalar_field(state,"TotalWaterQ")
    else if (has_scalar_field(state,"VapourWaterQ")) then
      q_w=>extract_scalar_field(state,"VapourWaterQ") 
    else
      q_w=>dummyscalar
    end if
    if (has_scalar_field(state,"Qrain")) then
      q_r=>extract_scalar_field(state,"Qrain")
    else
      q_r=>dummyscalar
    end if
    if (has_scalar_field(state,"QrainSinkingVelocity")) then
      q_r_vel=>extract_scalar_field(state,"QrainSinkingVelocity")
    else
      q_r_vel=>dummyscalar
    end if
    
    ! Radiation fluxes
    if (has_scalar_field(state,"SWRadiationFlux")) then
      sw=>extract_scalar_field(state,"SWRadiationFlux")
    else
      sw=>dummyscalar
    end if
    if (has_scalar_field(state,"LWRadiationFlux")) then
      lw=>extract_scalar_field(state,"LWRadiationFlux")
    else
      lw=>dummyscalar
    end if
    
    ! Get precipitation
    call allocate (precipitation, mesh, name='Precipitation')
    call set(precipitation,q_r)
    call scale(precipitation,q_r_vel)
    call scale(precipitation,density)
    call scale(precipitation,1./rhow)
    call scale(precipitation,mps_to_mmph)

    ! Allocate surface and reference 
    call allocate(surface_precip,surface_mesh,name="SurfacePrecipitation")
    call allocate(surface_pressure,surface_mesh,name="SurfacePressure")
    call allocate(surface_short_wave,surface_mesh,name="SurfaceSW")
    call allocate(surface_long_wave,surface_mesh,name="SurfaceLW")
    call allocate(surface_qr,surface_mesh,name="SurfaceQRSinkingVel")
    call allocate(surface_qr_vel,surface_mesh,name="SurfaceQR")
    
    call allocate(reference_height_temperature,surface_mesh,name="ReferenceTemperature")
    call allocate(reference_height_density,surface_mesh,name="ReferenceDensity")
    call allocate(reference_height_moisture,surface_mesh,name="ReferenceMoisture")    
    call allocate(reference_height_cp,surface_mesh,name="ReferenceCP")
    call allocate(reference_height_cv,surface_mesh,name="ReferenceCV")        
    call allocate(reference_height_wind,wind%dim,surface_mesh,name="ReferenceWind")
    
    ! Remap surface fields
    call remap_field_to_surface(precipitation, surface_precip, surface_element_list)
    call remap_field_to_surface(pressure, surface_pressure, surface_element_list)
    call remap_field_to_surface(sw, surface_short_wave, surface_element_list)
    call remap_field_to_surface(lw, surface_long_wave, surface_element_list)
    call remap_field_to_surface(q_r, surface_qr, surface_element_list)
    call remap_field_to_surface(q_r_vel, surface_qr_vel, surface_element_list)

    ! We need to create a virtual boundary, 10m above the real one, for the flux algorithm
    position => extract_vector_field(state, "Coordinate")
    virtual_bc_position=get_coordinates_remapped_to_surface(position, surface_mesh, surface_element_list)
    do ele=1,node_count(surface_mesh)
      call addto (virtual_bc_position,virtual_bc_position%dim,ele,zref)      
    enddo
      
    !Find node ownership: in which elements do we find the new nodes
    map = get_element_mapping(position, virtual_bc_position)
    
    ! Create interpolation basis
    call create_surface_mesh(X_surface_mesh, surface_nodes, position%mesh, name='CoordinateSurfaceMesh')

    call linear_interpolation (temperature, position, reference_height_temperature, virtual_bc_position, map)
    call linear_interpolation (density, position, reference_height_density, virtual_bc_position, map)
    call linear_interpolation (q_w, position, reference_height_moisture, virtual_bc_position, map)
    call linear_interpolation (qc_p, position, reference_height_cp, virtual_bc_position, map)
    call linear_interpolation (qc_v, position, reference_height_cv, virtual_bc_position, map)
    call linear_interpolation (wind, position, reference_height_wind, virtual_bc_position, map)

    ! Get model options 	
    if (have_option(trim(bc_path)//"/cool_skin_calculation")) jcool=1
    if (have_option(trim(bc_path)//"/wave_calculation")) jwave=2
    
    call get_option('/physical_parameters/location_date/longitude', lon, default=0.0)
    call get_option('/physical_parameters/location_date/latitude', lat, default=0.0)
    call get_option('/timestepping/timestep',dt)
    call get_option('/timestepping/current_time',current_time)
    
    ! Call COARE and get surface fluxes: Careful with the units
    ! Loop over surface nodes
    do j=1,node_count(surface_mesh)
      if (have_option('/radiation_model')) then
    	rs=node_val(surface_short_wave,j)
    	rl=node_val(surface_long_wave,j)
      else
        rs=0.0
        rl=0.0
      endif

      us=node_val(reference_height_wind,j)	 	      !U units: m/s
      ts=node_val(reference_height_temperature,j)-273.16      !T units: C
      ds=node_val(reference_height_density,j)		      !D units: kg/m3
      qs=node_val(reference_height_moisture,j)*1000.	      !Q units: g/kg
      cp=node_val(reference_height_cp,j)
      cv=node_val(reference_height_cv,j)

      ps=node_val(surface_pressure,j)/100.	      !P units: mb
      sst=node_val(sea_surface_t,j)-273.16	      !SST units: C
      rf=node_val(surface_precip,j)
      hs=node_val(sensible_flux,j)
      hl=node_val(latent_flux,j)

      do i =1,momentum_stress%dim
        tau(i)=node_val(momentum_stress,i,j)
      enddo

      call compute_COARE_fluxes (dt,current_time,lon,lat,	      &
    				 zref,ps,us,ts,ds,qs,rs,rl,cp,cv,     &
        			 sst,tau,hs,hl,rf,		      &
        			 jcool,jwave)

      call set(sensible_flux,j,hs)
      call set(latent_flux,j,hl)
      call set(rain_flux,j,rf)
      call set(sea_surface_t,j,sst+273.16)

      do i=1,momentum_stress%dim     
    	call set(momentum_stress,i,j,tau(i))  
      enddo    
    end do
    
    call deallocate(dummyscalar)    
    call deallocate(precipitation)
    call deallocate(surface_pressure)
    call deallocate(surface_precip)
    call deallocate(surface_short_wave)
    call deallocate(surface_long_wave)
    call deallocate(surface_qr)
    call deallocate(surface_qr_vel)
    
    call deallocate(reference_height_temperature)
    call deallocate(reference_height_density)
    call deallocate(reference_height_moisture)
    call deallocate(reference_height_cp)
    call deallocate(reference_height_cv)
    call deallocate(reference_height_wind)
    call deallocate(virtual_bc_position)
    call deallocate(X_surface_mesh)
    
    ewrite(1,*) 'End COARE3_interface'

contains

    subroutine create_basis (position, remapped_position, surface_nodes, map, basis)
    
      type(vector_field), intent(inout) :: position, remapped_position
      integer, dimension(:), pointer, intent(in) :: surface_nodes
      real, dimension(node_count(remapped_position),ele_loc(position,1)), intent(out) :: basis
      integer, dimension(node_count(remapped_position)), intent(in) :: map
   
      integer, dimension(:), pointer :: node_list
      integer :: info, node, ni, nl, nj      
      real, dimension(position%dim, ele_loc(position,1)) :: X_val
      real, dimension(position%dim) :: X_mean, X_new, phi
      real, dimension(position%dim,position%dim) :: phimat
      
      do node=1,node_count(remapped_position)
      
        ele=map(node)
	node_list=>ele_nodes(position,ele)
	
	X_new=node_val(remapped_position,node)
	
	X_val=ele_val(position,ele)
	X_mean=sum(X_val,2)/size(X_val,2)
	
	!Construct interpolation basis skipping nj
	do nj = 1, size(node_list)
	
	  if (any(surface_nodes == node_list(nj))) cycle
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

    subroutine interpolate_reference_height_scalar (field, remapped_field, map, basis)
    
      type(scalar_field), intent(in) :: field
      type(scalar_field), intent(inout) :: remapped_field
      real, dimension(:,:), intent(in) :: basis
      integer, dimension(:), intent(in) :: map
   
      integer, dimension(field%mesh%shape%numbering%vertices) :: ele_vertices
      real, dimension(ele_loc(field,1)) :: T_val
      real :: T_mean, T_new
      integer :: ele, node
      
      call zero(remapped_field)

      !Find vertices on old mesh
      ele_vertices=local_vertices(field%mesh%shape)
      
      do node=1,node_count(remapped_field)

        ele=map(node)
	
	T_val=ele_val(field,ele)
	T_mean=sum(T_val(ele_vertices))/size(ele_vertices)
	
	!Interpolate scalar value
	T_new=T_mean+sum(basis(node,:)*(T_val(ele_vertices)-T_mean))
	call addto(remapped_field, node, T_new)
	      
      enddo
      
    end subroutine interpolate_reference_height_scalar  

    subroutine interpolate_reference_height_vector (field, remapped_field, map, basis)
    
      type(vector_field), intent(in) :: field
      type(vector_field), intent(inout) :: remapped_field
      real, dimension(:,:), intent(in) :: basis
      integer, dimension(:), intent(in) :: map
   
      integer, dimension(field%mesh%shape%numbering%vertices) :: ele_vertices
      real, dimension(field%dim,ele_loc(field,1)) :: T_val
      real, dimension(field%dim) :: T_mean, T_new
      integer :: node, ele, nd
      
      call zero(remapped_field)

      !Find vertices on old mesh
      ele_vertices=local_vertices(field%mesh%shape)
      
      do node=1,node_count(remapped_field)

        ele=map(node)
	
	T_val=ele_val(field,ele)
	do nd = 1, field%dim
	  T_mean(nd)=sum(T_val(nd,ele_vertices))/size(ele_vertices)
	enddo
	
	!Interpolate scalar value
	do nd = 1, field%dim
  	  T_new(nd)=T_mean(nd)+sum(basis(node,:)*(T_val(nd,ele_vertices)-T_mean(nd)))
	enddo
	call addto(remapped_field, node, T_new)
	      
      enddo
      
    end subroutine interpolate_reference_height_vector
      
  end subroutine COARE_interface

  subroutine compute_COARE_fluxes (dt,current_time,lon,lat,	&
  				   Z,P,U,T,D,Q,Rs,Rl,cp,cv,	&
				   Ts,Tau,Hs,Hl,RF,		&
				   jcool,jwave)

!*********** basic specifications  *****
!	zu=			height of wind measurement
!	zt=			height of air temperature measurement
!	zq=			height of air humidity measurement
!	ts_depth	depth of water temperature measurement
!	jwarm=		0=no warm layer calc, 1 =do warm layer calc
!	jcool=		0=no cool skin calc, 1=do cool skin calc
!       jwave=          0= Charnock, 1=Oost et al, 2=Taylor and Yelland
!
!***********   input data **************
!	YYYYMMHHMMSS=		date in toga coare format, Y2K version
!	u=			wind speed (m/s), height zu
!	us=			surface current (m/s)
!	ts=			bulk surface sea temp (cent)
!	t=			air temp (cent), height zt
!	qs=			sea surface sat specific humidity (g/kg)
!	q=			air specific humidity (g/kg), height zq
!	Rs=			downward solar flux (w/m^2)
!	Rl=			downward IR flux (w/m^2)
!	zi=			inversion height (m)
!	P=			air pressure (mb)
!	rain=		rain rate (mm/hr)
!	lon=		longitude (deg E=+)
!	lat=		latitude (deg N=+)
!
!********** output data  ***************
!	hsb=			sensible heat flux (w/m^2)
!	hlb=			latent heat flux (w/m^2)
!	RF=			rain heat flux(w/m^2)
!	wbar=	   		webb mean w (m/s)
!	tau=			stress (nt/m^2)
!	zo=			velocity roughness length (m)
!	zot			temperature roughness length (m)
!	zoq=			moisture roughness length (m)
!	L=			Monin_Obukhov stability length
!	usr=			turbulent friction velocity (m/s), including gustiness
!	tsr			temperature scaling parameter (K)
!	qsr			humidity scaling parameter (g/g)
!	dter=			cool skin temperature depression (K)
!	dqer=			cool skin humidity depression (g/g)
!	tkt=			cool skin thickness (m)
!	Cd=			velocity drag coefficient at zu, referenced to u
!	Ch=			heat transfer coefficient at zt
!	Ce=			moisture transfer coefficient at zq
!	Cdn_10=			10-m velocity drag coeeficient, including gustiness
!	Chn_10=			10-m heat transfer coeeficient, including gustiness
!	Cen_10=			10-m humidity transfer coeeficient, including gustiness
!
!  Compared to the original COARE 3.0 algorithm, the warm-layer correction has been
!  removed. This latter requires time integrals of the energy and momentum fluxes 
!  across the atmosphere/ocean interface and therefore assumes that the complete
!  flux history (at least from the night before) is known. This is obviously not
!  the case in general.
!
  implicit none 
  
  integer, intent(in) :: jcool, jwave
  real, intent(in) :: dt, current_time, lon, lat
  real, dimension(:), intent(in) :: U
  real, intent(in) :: P, T, Q, D, Z, Rl, Rs, cp, cv
  
  real, intent(inout) :: Ts, Hs, Hl, RF
  real, dimension(:), intent(inout) :: Tau
  
  ! internal variables
  integer :: jwarm
  real :: x(20), y(23)
  real :: umag, us, qs, hwave, twave, pcom
  real :: hsb, hlb, taub, zo, zot, zoq, L, usr, tsr, qsr, dter, dqer, tkt, wbar, Cd, Ch, Ce, Cdn_10, Chn_10, Cen_10, Wg, S

    ! Various initialisations
    jwarm=0
    us=0.
    dter=0.
    dqer=0.
    qs=qsee(Ts, P) 

    ! Winds 
    if(size(U)==2 .or. size(U)==1) then
      umag=U(1)
      taub=Tau(1)
    else
      umag=sqrt(U(1)**2.+U(2)**2.)
      taub=sqrt(Tau(1)**2.+Tau(2)**2.)
    endif
    
    ! Wave properties
    twave=b*umag 
    hwave=a*umag**2.*(1+.015*umag)
    
    ! input array for cor30a (x has 20 elements)
    x=(/umag, Us, Ts, T, qs, Q, D, Rs, Rl, RF, zi, P, Z, Z, Z, lat, lon, twave, hwave, 0./)
    
    if (jwarm/=0) call warm_layer(dt,current_time,jcool,x,Hs,Hl,RF,taub)
    
    ! Call modified LKB routine
    call cor30a(jcool,jwave,x,y) 
    
    ! Get outputs from cor30 (*: variables not included in original output)
    hsb=y(1)			!sensible heat flux W/m/m
    hlb=y(2)			!latent
    taub=y(3)			!stress
    zo=y(4)			!vel roughness*
    zot=y(5)			!temp "*
    zoq=y(6)			!hum  "*
    L=y(7)			!Ob Length*
    usr=y(8)			!ustar*
    tsr=y(9)			!tstar*
    qsr=y(10)			!qstar  [g/g]*
    dter=y(11)  		!cool skin delta T
    dqer=y(12)  		!  "   "     "   q*
    tkt=y(13)			!thickness of cool skin
    RF=y(14)			!rain heat flux
    wbar=y(15)  		!webb mean w	 
    Cd=y(16)			!drag @ zu*
    Ch=y(17)			!Heat transfer coeff*
    Ce=y(18)			!Dalton (humidity transfer coeff)*
    Cdn_10=y(19)		!neutral drag @ 10m [includes gustiness]*
    Chn_10=y(20)		!Heat transfer coeff @ 10m*
    Cen_10=y(21)		!Humidity transfer coeff @ 10m*
    Wg=y(22) 
    S=y(23)			!Average wind speed including gusts @ 10m*
    
    Ts=Ts-dter*jcool 

    ! Fluxes in raw flux form
    Hs=hsb	
    Hl=hlb	
    Tau(1)=taub*U(1)/umag
    Tau(2)=taub*U(2)/umag
    
    ! Fluxes in bulk form
!    pcom=(P/1000.)**((cp-cv)/cp)
!    Hs=D*cp*Ch*S*(Ts-T/pcom)
!    Hl=D*flv*Ce*S*(qs-Q)
!    Tau(1)=D*Cd*S*U(1)
!    Tau(2)=D*Cd*S*U(2)
     
end subroutine compute_COARE_fluxes

  subroutine warm_layer (dtime,time,jcool,x,Hs,Hl,RF,Tau)
  
    ! This routine has been kept from the original COARE algorithm
    ! but is not supposed to be called at the moment (jwarm=0).
  
    integer :: i,jcool
    real :: dtime, time, x(20), Hs, Hl, RF, Tau
    real :: umag, Us, Ts, T, qs, Q, D, Rs, Rl, rain, zinv, P, zu, zt, zq, lat, lon, twave, hwave
    real :: intime, chktime, newtime, loc, locx
    real :: rhoa, Rnl, Rns, ctd1, ctd2, grav, Le, visa, cpw, qcol_ac, tau_ac
    real :: Al, tsea, dsea, qjoule, qr_out, q_pwp, tk_pwp, fxp, dt_wrm, ts_depth

    umag=x(1) !wind speed (m/s)  at height zu (m)
    us=x(2) !surface current speed in the wind direction (m/s)
    ts=x(3) !bulk water temperature (C) if jcool=1, interface water T if jcool=0  
    t=x(4) !bulk air temperature (C), height zt
    Qs=x(5)/1000 !bulk water spec hum (g/kg) if jcool=1, ...
    Q=x(6)/1000 !bulk air spec hum (g/kg), height zq
    rhoa=x(7) !Air density
    Rs=x(8) !downward solar flux (W/m**2)
    Rl=x(9) !downard IR flux (W/m**2)
    rain=x(10) !rain rate (mm/hr)
    zinv=x(11) !PBL depth (m)
    P=x(12) !Atmos surface pressure (mb)
    zu=x(13) !wind speed measurement height (m)
    zt=x(14) !air T measurement height (m)
    zq=x(15) !air q measurement height (m)
    lat=x(16) !latitude (deg, N=+)
    lon=x(17) !longitude
    twave=x(18) !wave period (s)
    hwave=x(19) !wave height (m)    

!    if(version_af)then
      ts_depth=.05 !bulk water temperature sensor depth, ETL sea snake
!    else if(version_ah)then
!      ts_depth=6. !bulk water temperature sensor depth, ETL sea snake
!     endif

    qcol_ac=0.
    tau_ac=0.
   
    ! Get right time and location
    intime=time/24./3600.
    loc=(lon+7.5)/15 
    locx=loc     
    tsea=ts
    
    Al=2.1e-5*(tsea+3.2)**0.79 
    cpw=c_p_l
    grav=grv(lat) 
    Le=(2.501-.00237*tsea)*1e6 
    visa=1.326e-5*(1+6.542e-3*T+8.301e-6*T*T-4.84e-9*T*T*T) 

    ! Apply warm layer
    chktime=loc+intime 
    newtime=(chktime-24*floor(chktime/24))*3600 

    ctd1=sqrt(2*rich*cpw/(Al*grav*rhow)) 
    ctd2=sqrt(2*Al*grav/(rich*rhow))/(cpw**1.5) 

    Rnl=0.97*(5.67e-8*(ts-dter*real(jcool)+273.16)**4-Rl)
    Rns=0.945*Rs 	!oceanic albedo=0.055 daily average
    qr_out=Rnl+Hs+Hl+RF     !total cooling at surface
    q_pwp=fxp*Rns-qr_out    !tot heat abs in warm layer
    
    ! Check total energy threshold
    if (q_pwp > 50 .or. qcol_ac > 0) then !check for threshold
      ! increment momentum integral
      tau_ac=tau_ac+max(.002,Tau)*dtime       
    	
      ! increment energy integral (needs iterations)
      if ((qcol_ac+q_pwp*dtime) > 0) then  !check threshold for warm layer existence
    	do i=1,5 
    	  fxp=1 - (0.28*0.014*(1-exp(-tk_pwp/0.014)) + 0.27*0.357*(1-exp(-tk_pwp/0.357)) &
	  	+ 0.45*12.82*(1-exp(-tk_pwp/12.82))) / tk_pwp 
    	  qjoule=(fxp*Rns-qr_out)*dtime 
    	  if (qcol_ac+qjoule > 0) then
    	    tk_pwp=min(19.,ctd1*tau_ac/sqrt(qcol_ac+qjoule)) 
    	  endif 
    	enddo !  end i loop
      else  !warm layer wiped out
    	fxp=0.75 
    	tk_pwp=19.0
	qjoule=(fxp*Rns-qr_out)*dtime 
      endif 
      qcol_ac=qcol_ac+qjoule  

      if (qcol_ac > 0) dt_wrm=ctd2*(qcol_ac)**1.5/tau_ac 
    endif ! end threshold check
    
    if (tk_pwp < ts_depth) then
      dsea=dt_wrm 
    else
      dsea=dt_wrm*ts_depth/tk_pwp 
    endif 
    
    Ts=tsea+dsea 
    qs=qsee(Ts, P)
    
    !Update x array 
    x=(/umag, Us, Ts, T, qs, Q, rhoa, Rs, Rl, rain, zinv, P, zu, zt, zq, lat, lon, twave, hwave, 0./)

  end subroutine warm_layer
    
  subroutine cor30a(jcool,jwave,x,y)
  
  implicit none
  
  !version with shortened iteration  modified Rt and Rq
  !uses wave information wave period in s and wave ht in m
  !no wave, standard coare 2.6 charnock:  jwave=0 
  !Oost et al.  zo=50/2/pi L (u*/c)**4.5 if jwave=1
  !taylor and yelland  zo=1200 h*(L/h)**4.5 jwave=2

  integer i, nits
  integer jcool, jwave

  real x(20), y(23)
  real u,us,ts,t,Qs,Q,Rs,Rl,rain,P,zu,zt,zq,lat,lon,twave,hwave
  real grav,cpw,Le,cpa,rhoa,visa,Al,bigc,wetc,lwave,cwave,Rns,Rnl,du,dt,dq,qout,dels,qcol,alq,xlamx,alfac,bf,cc,cd10,ch10,charn,ct,ct10,dtmp,dwat,hl_webb,zinv
  real l10,ribcu,ribu,rr,ta,u10,ut,zet,zetu,zo10,zot10
  real hsb, hlb, tau, zo, zot, zoq, L, usr, tsr, qsr, dter, dqer, tkt, RF, wbar, Cd, Ch, Ce, Cdn_10, Chn_10, Cen_10, ug 

    u=x(1) !wind speed (m/s)  at height zu (m)
    us=x(2) !surface current speed in the wind direction (m/s)
    ts=x(3) !bulk water temperature (C) if jcool=1, interface water T if jcool=0  
    t=x(4) !bulk air temperature (C), height zt
    Qs=x(5)/1000 !bulk water spec hum (g/kg) if jcool=1, ...
    Q=x(6)/1000 !bulk air spec hum (g/kg), height zq
    rhoa=x(7) !Air density
    Rs=x(8) !downward solar flux (W/m**2)
    Rl=x(9) !downard IR flux (W/m**2)
    rain=x(10) !rain rate (mm/hr)
    zinv=x(11) !PBL depth (m)
    P=x(12) !Atmos surface pressure (mb)
    zu=x(13) !wind speed measurement height (m)
    zt=x(14) !air T measurement height (m)
    zq=x(15) !air q measurement height (m)
    lat=x(16) !latitude (deg, N=+)
    lon=x(17) !longitude
    twave=x(18) !wave period (s)
    hwave=x(19) !wave height (m)
  
    ! Atmosphere physical constants (gravity, Cp, latent heat, viscosity)
    grav=grv(lat) !9.82 
    cpa=c_p
    Le=(2.501-.00237*ts)*1e6 
    visa=1.326e-5*(1+6.542e-3*t+8.301e-6*t*t-4.84e-9*t*t*t) 

    ! Water surface physical constants (Cp, albedo, )
    cpw=c_p_l
    Al=2.1e-5*(ts+3.2)**0.79 
    bigc=16.*grav*cpw*(rhow*visw)**3/(tcw*tcw*rhoa*rhoa) 
    wetc=0.622*Le*Qs/((c_p-c_v)*(ts+tdk)**2.) 
    
    lwave=grav/2./pi*twave**2. 
    cwave=grav/2./pi*twave 
    
    Rns=Rs*.945 
    Rnl=0.97*(5.67e-8*(ts-0.3*jcool+tdk)**4-Rl) 
        
    ! first guess
    du=u-us 
    dt=ts-t-.0098*zt 
    dq=Qs-Q 
    ta=t+tdk 
    ug=.5 
    dter=0.3  
    dqer=wetc*dter 
    ut=sqrt(du*du+ug*ug) 
    u10=ut*log(10/1e-4)/log(zu/1e-4) 
    usr=.035*u10 
    zo10=0.011*usr*usr/grav+0.11*visa/usr 
    Cd10=(von/log(10/zo10))**2 
    Ch10=0.00115 
    Ct10=Ch10/sqrt(Cd10) 
    zot10=10/exp(von/Ct10) 
    Cd=(von/log(zu/zo10))**2 
    Ct=von/log(zt/zot10) 
    CC=von*Ct/Cd 
    Ribcu=-zu/zinv/.004/Beta**3 
    Ribu=-grav*zu/ta*((dt-dter*jcool)+.61*ta*dq)/ut**2 
    nits=3 
    
    if (Ribu < 0) then 
      zetu=CC*Ribu/(1+Ribu/Ribcu) 
    else 
      zetu=CC*Ribu*(1+27/9*Ribu/CC)
    endif 
    
    L10=zu/zetu 
    if (zetu > 50) then 
      nits=1 
    endif 
    
    usr=ut*von/(log(zu/zo10)-psiuo(zu/L10))
    tsr=-(dt-dter*jcool)*von*fdg/(log(zt/zot10)-psit_30(zt/L10)) 
    qsr=-(dq-wetc*dter*jcool)*von*fdg/(log(zq/zot10)-psit_30(zq/L10)) 
    tkt=.001
    
    charn=0.011 
    if (ut > 10) then
      charn=0.011+(ut-10)/(18-10)*(0.018-0.011) 
    endif 
    if (ut > 18) then
      charn=0.018 
    endif 
    
    ! bulk loop
    do i=1, nits 
      zet=von*grav*zu/ta*(tsr*(1+0.61*Q)+.61*ta*qsr)/(usr*usr)/(1+0.61*Q) 
      if (jwave .EQ. 0) zo=charn*usr*usr/grav+0.11*visa/usr  
!      if (jwave .EQ. 1) zo=50/2/pi*lwave*(usr/cwave)**4.5+0.11*visa/usr !Oost et al
      if (jwave .EQ. 2) zo=1200*hwave*(hwave/lwave)**4.5+0.11*visa/usr !Taylor and Yelland
      
      rr=zo*usr/visa 
      L=zu/zet 
      zoq=min(1.15e-4,5.5e-5/rr**.6) 
      zot=zoq 
      usr=ut*von/(log(zu/zo)-psiuo(zu/L)) 
      tsr=-(dt-dter*jcool)*von*fdg/(log(zt/zot)-psit_30(zt/L)) 
      qsr=-(dq-wetc*dter*jcool)*von*fdg/(log(zq/zoq)-psit_30(zq/L)) 
      Bf=-grav/ta*usr*(tsr+.61*ta*qsr) 
      if (Bf > 0) then
        ug=Beta*(Bf*zinv)**.333 
      else
        ug=.2 
      endif
    
      ut=sqrt(du*du+ug*ug) 
      Rnl=0.97*(5.67e-8*(ts-dter*jcool+tdk)**4-Rl) 
      hsb=-rhoa*cpa*usr*tsr 
      hlb=-rhoa*Le*usr*qsr 
      qout=Rnl+hsb+hlb 
      dels=Rns*(.065+11*tkt-6.6e-5/tkt*(1-exp(-tkt/8.0e-4))) ! Eq.16 Shortwave
      qcol=qout-dels 
      alq=Al*qcol+be*hlb*cpw/Le  ! Eq. 7 Buoy flux water

      if (alq > 0) then 
        xlamx=6/(1+(bigc*alq/usr**4)**.75)**.333    ! Eq 13 Saunders
        tkt=xlamx*visw/(sqrt(rhoa/rhow)*usr)    !Eq.11 Sub. thk
      else
        xlamx=6.0 
        tkt=min(.01,xlamx*visw/(sqrt(rhoa/rhow)*usr))   !Eq.11 Sub. thk
      endif 
     
      dter=qcol*tkt/tcw !  Eq.12 Cool skin
      dqer=wetc*dter
    enddo 
    
    tau=rhoa*usr*usr*du/ut
    hsb=-rhoa*cpa*usr*tsr 
    hlb=-rhoa*Le*usr*qsr 
     
    ! rain heat flux
    dwat=2.11e-5*((t+tdk)/tdk)**1.94 !! water vapour diffusivity
    dtmp=(1.+3.309e-3*t-1.44e-6*t*t)*0.02411/(rhoa*cpa)   !!heat diffusivity
    alfac=1/(1+(wetc*Le*dwat)/(cpa*dtmp))    !! wet bulb factor
    RF=rain*alfac*cpw*((ts-t-dter*jcool)+(Qs-Q-dqer*jcool)*Le/cpa)/3600 
    
    ! Webb et al. correection
    wbar=1.61*hlb/Le/(1+1.61*Q)/rhoa+hsb/rhoa/cpa/ta !formulation in hlb already includes webb
    !wbar=1.61*hlb/Le/rhoa+(1+1.61*Q)*hsb/rhoa/cpa/ta 
    hl_webb=rhoa*wbar*Q*Le 
    
    ! compute transfer coeffs relative to ut @meas. ht 
    Cd=tau/rhoa/ut/max(.1,du) 
    Ch=-usr*tsr/ut/(dt-dter*jcool) 
    Ce=-usr*qsr/(dq-dqer*jcool)/ut 
    
    ! 10-m neutral coeff realtive to ut
    Cdn_10=von*von/log(10/zo)/log(10/zo) 
    Chn_10=von*von*fdg/log(10/zo)/log(10/zot) 
    Cen_10=von*von*fdg/log(10/zo)/log(10/zoq) 
    
   ! Fill in output array
   y=(/hsb, hlb, tau, zo, zot, zoq, L, usr, tsr, qsr, dter, dqer, tkt, RF, wbar, Cd, Ch, Ce, Cdn_10, Chn_10, Cen_10, ug, ut /) 

  end subroutine cor30a
 
  real function grv(lat)
    real :: lat
    real :: x, gamma, c1, c2, c3, c4, phi
    gamma=9.7803267715
    c1=0.0052790414
    c2=0.0000232718
    c3=0.0000001262
    c4=0.0000000007
    
    phi=lat*pi/180
    x=sin(phi)
    grv=gamma*(1+(c1*x**2)+(c2*x**4)+(c3*x**6)+(c4*x**8))
    return
  end function grv

  real function psit_30(zet)
    real :: zet
    real :: x, psik, psic, f, c
    x=(1.-(15*zet))**.5 
    psik=2*log((1+x)/2) 
    x=(1.-(34.15*zet))**.3333 
    psic=1.5*log((1.+x+x*x)/3.)-sqrt(3.)*atan((1.+2.*x)/sqrt(3.))+4.*atan(1.)/sqrt(3.) 
    f=zet*zet/(1+zet*zet) 
    psit_30=(1-f)*psik+f*psic   
   
    if(zet>0)then 
      c=min(50.,.35*zet) 
      psit_30=-((1.+2./3.*zet)**1.5+.6667*(zet-14.28)/exp(c)+8.525)
   endif
   return
  end function psit_30

  real function psiuo(zet)
    real :: zet
    real :: x, psik, psic, f, c
    x=(1.-15.*zet)**.25 
    psik=2.*log((1.+x)/2.)+log((1.+x*x)/2.)-2.*atan(x)+2.*atan(1.) 
    x=(1.-10.15*zet)**.3333 
    psic=1.5*log((1.+x+x*x)/3.)-sqrt(3.)*atan((1.+2.*x)/sqrt(3.))+4.*atan(1.)/sqrt(3.) 
    f=zet*zet/(1+zet*zet) 
    psiuo=(1-f)*psik+f*psic                                                
    if(zet>0)then 
      c=min(50.,.35*zet) 
      psiuo=-((1+1.0*zet)**1.0+.667*(zet-14.28)/exp(c)+8.525)
    endif 
    return
  end function psiuo 

  real function qsat(y)
    real :: y(2)
    real :: x, p, es
    x=y(1) !temp
    p=y(2) !pressure
    es=6.112*exp(17.502*x/(x+241.0))*(1.0007+3.46e-6*p)
    qsat=es*622./(p-.378*es)
    return
  end function qsat

  real function qsee(ts,Pa)
    real :: ts,Pa
    real :: p, es, x
    x=ts
    p=Pa
    es=6.112*exp(17.502*x/(x+240.97))*.98*(1.0007+3.46e-6*p)
    qsee=es*621.97/(p-.378*es)
    return
  end function

end module Surface_Ocean_COARE
