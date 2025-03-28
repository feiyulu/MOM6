!> Interfaces for MOM6 ensembles and data assimilation.
module MOM_oda_ml_mod

! This file is part of MOM6. see LICENSE.md for the license.

! MOM infrastructure
use MOM_cpu_clock, only : cpu_clock_begin, cpu_clock_end, cpu_clock_id
use MOM_verticalGrid, only : verticalGrid_type
use gsw_mod_toolbox, only : gsw_ct_from_pt, gsw_sigma0
use netcdf, only : nf90_open, nf90_inq_varid, nf90_get_var, nf90_close, nf90_close
use netcdf, only : nf90_nowrite, nf90_noerr

implicit none ; private

public :: oda_ml_init, oda_ml_end, oda_ml_inference

! Data structure to save the ML configuration, input, and output data
type, public :: ocean_oda_ml_config ; private
    character(len=255)  :: filename
    real, dimension(16,66)  :: l1_weight
    real, dimension(16,16)  :: l2_weight, l3_weight
    real, dimension(16) :: l1_bias, l2_bias, l3_bias
    real, dimension(:), allocatable :: z_l
    real, dimension(:), allocatable :: z_i
    integer :: nk
end type ocean_oda_ml_config

type, public :: ocean_oda_ml_data
    integer :: nk
    real :: dyCu_left, dyCu_right, dxCv_south, dxCv_north, areacello
    real :: bathyT, bathyU_left, bathyU_right, bathyV_south, bathyV_north
    real :: mask2dT, OBCmaskCu_left, OBCmaskCu_right, OBCmaskCv_south, OBCmaskCv_north 
    !! Input features
    real :: SSH !<sea surface height (m) across ensembles
    real :: taux_left !<zonal wind stress
    real :: taux_right !<zonal wind stress
    real :: tauy_north !<meridional wind stress
    real :: tauy_south !<zonal wind stress
    real :: latent !<latent heat flux
    real :: sensible !<sensile heat flux
    real :: lw !<longwave radiation flux
    real :: sw !<shortwave radiation flux
    real, pointer, dimension(:) :: T=>NULL() !<layer potential temperature (degC) across ensembles
    real, pointer, dimension(:) :: S=>NULL() !<layer salinity (psu or g kg-1) across ensembles
    real, pointer, dimension(:) :: U_left=>NULL() !<layer zonal velocity (m s-1) across ensembles
    real, pointer, dimension(:) :: U_right=>NULL() !<layer zonal velocity (m s-1) across ensembles
    real, pointer, dimension(:) :: V_north=>NULL() !<layer meridional velocity (m s-1) across ensembles
    real, pointer, dimension(:) :: V_south=>NULL() !<layer meridional velocity (m s-1) across ensembles

    !! Output predictions
    real, pointer, dimension(:) :: T_inc=>NULL()
    real, pointer, dimension(:) :: S_inc=>NULL()

end type ocean_oda_ml_data

real :: PRHO_change = 0.03
real :: reference_depth = 10
!real :: value_for_control_depth = 1E6
!real :: value_for_control_PRHO = 999
!real :: value_for_control_oceanzvars = 1E5
real :: ReLU_zero = 0
real, dimension(15) :: target_sigmas = (/0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9/)
real, dimension(16) :: output_flux_sigmas = (/0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0/)
character(len=255)  :: danni_ANN_name = '/gpfs/f5/gfdl_sd/world-shared/Danni.Du/ECDA_data/ML/danni_ANN_M4_TzBzdivmldtaufluxUz_2003_2016_10epoch.nc'
real :: seconds_in_30_days = 3600*24*30

integer :: id_clock_ml_remapping
integer :: id_clock_ml_normalization
integer :: id_clock_ml_inference

#include <MOM_memory.h>

character(len=40)  :: mdl = "MOM_oda_ml" !< This module's name.

contains

    subroutine oda_ml_inference(ml_config,ml_data)
        type(ocean_oda_ml_config), pointer, intent(in) :: ml_config
        type(ocean_oda_ml_data), pointer, intent(in) :: ml_data
        
        real :: SA, PT, CT, PRHO ,tauamp
        real, dimension(:), allocatable :: PRHO_profile
        real :: PRHO_mld, PRHO_10m 
        real :: mld_depth
        integer :: zl_index_mld, zl_index10m, zl_index_3mld, right_index
        real, dimension(:), allocatable :: zl_to_sigma, zi_to_sigma
        real :: thetao, so, uo_left, uo_right, vo_south, vo_north, div, thetao_top, thetao_bottom, so_top,so_bottom
        real :: PRHO_top, PRHO_bottom, uo_right_top, uo_right_bottom, uo_left_top, uo_left_bottom
        real, dimension(15) :: thetao_zgrad_sigma, so_zgrad_sigma, PRHO_zgrad_sigma, div_sigma, output_DT_sigmas, uo_zgrad_sigma
        real, dimension(:), allocatable :: thetao_zgrad_profile, so_zgrad_profile, div_profile, PRHO_zgrad_profile, uo_zgrad_profile
        real, dimension(66) :: ANN_input
        real, dimension(:), allocatable :: output_DT_at_zl, output_flux_at_zi
        real, dimension(:), allocatable :: z_l
        real, dimension(16) :: l1_output, l2_output, l3_output
        integer :: zz, i
        real :: mask_Tuv
        
        ml_data%T_inc=0.0
        
        mask_Tuv = ml_data%mask2dT + ml_data%OBCmaskCu_left + ml_data%OBCmaskCu_right + ml_data%OBCmaskCv_south + ml_data%OBCmaskCv_north
        if (mask_Tuv < 5.0) then
            ml_data%T_inc=0.0
        else
            allocate(z_l(ml_config%nk),source=0.0)
            z_l = ml_config%z_l

            allocate(PRHO_profile(ml_data%nk),source=0.0)
            do zz  = 1, ml_data%nk
                SA = ml_data%S(zz)
                PT = ml_data%T(zz)
                CT = gsw_ct_from_pt(SA, PT)
                PRHO = gsw_sigma0(SA, CT)
                PRHO_profile(zz) = PRHO
            end do

            ! find 3 MLD
            ! first find the first index below 10m
            call find_right_index(z_l, reference_depth, ml_data%bathyT, z_l, zl_index10m)
            ! if zl_index10m is not found or it is 1, set mld_depth to be bathyT, so that the ml inference will not be done
            if (zl_index10m <= 1) then
                mld_depth = ml_data%bathyT
            else
                ! the 10m potential density
                call interpolate(z_l(zl_index10m-1),z_l(zl_index10m),PRHO_profile(zl_index10m-1), PRHO_profile(zl_index10m),reference_depth,PRHO_10m)
                ! the MLD potential density
                PRHO_mld = PRHO_10m + PRHO_change
                ! the first z_l index below MLD
                call find_right_index(PRHO_profile, PRHO_mld, ml_data%bathyT, z_l, zl_index_mld)
                ! if zl_index10m is not found or it is 1, set mld_depth to be bathyT, so that the ml inference will not be done
                if (zl_index_mld <= 1) then
                    mld_depth = ml_data%bathyT
                else
                    ! the MLD depth
                    call interpolate(PRHO_profile(zl_index_mld-1),PRHO_profile(zl_index_mld),z_l(zl_index_mld-1),z_l(zl_index_mld),PRHO_mld,mld_depth)
                    ! the MLD must be below 10m
                    if (mld_depth < 10) then
                        mld_depth = 10
                    end if
                end if
            end if

            ! the first z_l index below 3MLD
            call find_right_index(z_l(1:ml_config%nk-1), 3*mld_depth, ml_data%bathyT, z_l, zl_index_3mld)

            if (zl_index_3mld == 0) then ! if 3 mld not found
                ml_data%T_inc=0.0
            else
                if (z_l(zl_index_3mld+1) > ml_data%bathyT .OR. &
                    z_l(zl_index_3mld+1) > ml_data%bathyU_left .OR. &
                    z_l(zl_index_3mld+1) > ml_data%bathyU_right .OR. &
                    z_l(zl_index_3mld) > ml_data%bathyV_south .OR. &
                    z_l(zl_index_3mld) > ml_data%bathyV_north) then
                    ml_data%T_inc=0.0
                else ! if above bathy, then get the vertical profiles

                    zi_to_sigma = ml_config%z_i(2:zl_index_3mld + 1)/mld_depth
                    allocate(thetao_zgrad_profile(zl_index_3mld),source=0.0)
                    !allocate(so_zgrad_profile(zl_index_3mld),source=0.0)
                    allocate(PRHO_zgrad_profile(zl_index_3mld),source=0.0)
                    allocate(uo_zgrad_profile(zl_index_3mld),source=0.0)

                    do zz = 1, zl_index_3mld
                        thetao_top = ml_data%T(zz)
                        so_top = ml_data%S(zz)
                        CT = gsw_ct_from_pt(so_top,thetao_top)
                        PRHO_top = gsw_sigma0(so_top,CT)
                        uo_right_top = ml_data%U_right(zz)
                        uo_left_top = ml_data%U_left(zz)
                    
                        thetao_bottom = ml_data%T(zz+1)
                        so_bottom = ml_data%S(zz+1)
                        CT = gsw_ct_from_pt(so_bottom,thetao_bottom)
                        PRHO_bottom = gsw_sigma0(so_bottom,CT)
                        uo_right_bottom = ml_data%U_right(zz+1)
                        uo_left_bottom = ml_data%U_left(zz+1)

                        thetao_zgrad_profile(zz) = (thetao_top - thetao_bottom)/(z_l(zz+1) - z_l(zz))
                        !so_zgrad_profile(zz) = (so_top - so_bottom)/(z_l(zz+1) - z_l(zz))
                        PRHO_zgrad_profile(zz) = (PRHO_top - PRHO_bottom)/(z_l(zz+1) - z_l(zz))
                        uo_zgrad_profile(zz) = (uo_right_top-uo_right_bottom+uo_left_top-uo_left_bottom)/(2*(z_l(zz+1) - z_l(zz)))
                    end do

                    zl_to_sigma = z_l(1:zl_index_3mld)/mld_depth

                    allocate(div_profile(zl_index_3mld),source=0.0)

                    do zz = 1, zl_index_3mld
                        call compute_current_divergence(ml_data%U_left(zz)*ml_data%dyCu_left, ml_data%U_right(zz)*ml_data%dyCu_right, &
                                ml_data%V_south(zz)*ml_data%dxCv_south, ml_data%V_north(zz)*ml_data%dxCv_north, &
                                ml_data%areacello, div)
                        div_profile(zz) = div
                    end do 
            
            
                    ! interpolate values to target_sigmas
                    do i = 1, 15
                        call find_right_index_clean(zi_to_sigma, target_sigmas(i), right_index)
                        ! quality control done in the previous steps, so that right_index >=1 and right_index <= zl_index_3mld
                        if (right_index == 1) then
                            thetao_zgrad_sigma(i) = thetao_zgrad_profile(1)
                            !so_zgrad_sigma(i) = so_zgrad_profile(1)
                            PRHO_zgrad_sigma(i) = PRHO_zgrad_profile(1)
                            uo_zgrad_sigma(i) = uo_zgrad_profile(1)
                        else
                            call interpolate(zi_to_sigma(right_index-1),zi_to_sigma(right_index),thetao_zgrad_profile(right_index-1), &
                                    thetao_zgrad_profile(right_index),target_sigmas(i),thetao_zgrad_sigma(i))
                            !call interpolate(zi_to_sigma(right_index-1),zi_to_sigma(right_index),so_zgrad_profile(right_index-1), &
                                    !so_zgrad_profile(right_index),target_sigmas(i),so_zgrad_sigma(i))
                            call interpolate(zi_to_sigma(right_index-1),zi_to_sigma(right_index),PRHO_zgrad_profile(right_index-1), &
                                    PRHO_zgrad_profile(right_index),target_sigmas(i),PRHO_zgrad_sigma(i))
                            call interpolate(zi_to_sigma(right_index-1),zi_to_sigma(right_index),uo_zgrad_profile(right_index-1), &
                                    uo_zgrad_profile(right_index),target_sigmas(i),uo_zgrad_sigma(i))
                        end if

                        call find_right_index_clean(zl_to_sigma, target_sigmas(i), right_index)
                        if (right_index == 1) then
                            div_sigma(i) = div_profile(1)
                        else
                            call interpolate(zl_to_sigma(right_index-1),zl_to_sigma(right_index),div_profile(right_index-1),div_profile(right_index),target_sigmas(i),div_sigma(i))
                        end if
                    end do
                    
                    tauamp = sqrt(((ml_data%taux_left+ml_data%taux_right)/2)**2+((ml_data%tauy_south+ml_data%tauy_north)/2)**2)

                    ! subroutine(input,DA tendency)        
                    ANN_input(1:15) = thetao_zgrad_sigma*100
                    ANN_input(16:30) = PRHO_zgrad_sigma*100
                    ANN_input(31:45) = div_sigma*1E7
                    ANN_input(46) = mld_depth*0.1
                    ANN_input(47) = tauamp*100
                    ANN_input(48) = ml_data%latent*0.1
                    ANN_input(49) = ml_data%sensible*0.1
                    ANN_input(50) = ml_data%lw*0.1
                    ANN_input(51) = ml_data%sw*0.1
                    ANN_input(52:66) = uo_zgrad_sigma*1000

                    l1_output = max(ReLU_zero, matmul(ml_config%l1_weight, ANN_input) + ml_config%l1_bias)
                    l2_output = max(ReLU_zero, matmul(ml_config%l2_weight, l1_output) + ml_config%l2_bias)
                    l3_output = matmul(ml_config%l3_weight, l2_output) + ml_config%l3_bias
                    
                    ! l3_output is the predicted flux
                    output_DT_sigmas =  (l3_output(1:15)-l3_output(2:16))/(0.2*mld_depth)
                    
                    !allocate(output_flux_at_zi(zl_index_3mld+1))
                    !output_flux_at_zi(1) = l3_output(1)
                    !do zz = 1, zl_index_3mld
                        !call find_right_index_clean(output_flux_sigmas, zi_to_sigma(zz), right_index)
                        !if (right_index == 0) then
                            !output_flux_at_zi(zz+1) = 0.0
                        ! it is known that right_index > 1
                        !else
                            !call interpolate(output_flux_sigmas(right_index-1),output_flux_sigmas(right_index),l3_output(right_index-1),&
                                    !l3_output(right_index),zi_to_sigma(zz),output_flux_at_zi(zz+1))
                    
                        !end if       
                    !end do
                
                    allocate(output_DT_at_zl(zl_index_3mld))
                    do zz = 1, zl_index_3mld
                        call find_right_index_clean(target_sigmas, zl_to_sigma(zz), right_index)
                        if (right_index == 0) then
                            output_DT_at_zl(zz) = 0.0
                        else if (right_index == 1) then
                            output_DT_at_zl(zz) = output_DT_sigmas(1)
                        else
                            call interpolate(target_sigmas(right_index-1),target_sigmas(right_index),output_DT_sigmas(right_index-1),&
                                    output_DT_sigmas(right_index),zl_to_sigma(zz),output_DT_at_zl(zz))
                        end if
                    end do

                    ml_data%T_inc(1:zl_index_3mld)=output_DT_at_zl / seconds_in_30_days
                endif
            endif ! end if 3 mld exceeds total number of levels
        
        endif
        ml_data%S_inc=0.0

    end subroutine oda_ml_inference

    subroutine oda_ml_init(ml_config,ml_data,GV)
        type(ocean_oda_ml_config), pointer, intent(in) :: ml_config
        type(ocean_oda_ml_data), pointer, intent(in) :: ml_data
        type(verticalGrid_type), pointer, intent(in) :: GV   !< The ocean's vertical grid structure

        ! load the NN weights and biases
        call read_ANN_file(ml_config)

        allocate(ml_config%z_l(GV%ke), source=0.0)
        ml_config%z_l = GV%sLayer

        allocate(ml_config%z_i(GV%ke+1), source=0.0)
        ml_config%z_i = GV%sInterface

        ml_config%nk = GV%ke
        call init_oda_ml_features(ml_data,GV%ke)    

    end subroutine oda_ml_init

    subroutine oda_ml_end(ml_config,ml_data)
        type(ocean_oda_ml_config), pointer, intent(in) :: ml_config
        type(ocean_oda_ml_data), pointer, intent(in) :: ml_data

    end subroutine oda_ml_end

    subroutine init_oda_ml_features(ml_data,nk)
        type(ocean_oda_ml_data), pointer, intent(in) :: ml_data
        integer, intent(in) :: nk

        ml_data%nk = nk
        allocate(ml_data%T(nk),source=0.0)
        allocate(ml_data%S(nk),source=0.0)
        allocate(ml_data%U_left(nk),source=0.0)
        allocate(ml_data%U_right(nk),source=0.0)
        allocate(ml_data%V_north(nk),source=0.0)
        allocate(ml_data%V_south(nk),source=0.0)

        allocate(ml_data%T_inc(nk),source=0.0)
        allocate(ml_data%S_inc(nk),source=0.0)

    end subroutine init_oda_ml_features

    Subroutine read_ANN_file(ml_config)
        implicit none
        type(ocean_oda_ml_config), pointer, intent(in) :: ml_config

        ! character(len=*), intent(in) :: filename
        ! real, dimension(16,66), intent(out) :: l1_weight
        ! real, dimension(16,16), intent(out) :: l2_weight, l3_weight
        ! real, dimension(16), intent(out) :: l1_bias, l2_bias, l3_bias

        real, dimension(66,16)  :: l1_weight_temp
        real, dimension(16,16) :: l2_weight_temp, l3_weight_temp
        integer :: ncid, varid, retval
        character(len = 255) :: varname

        ml_config%filename = danni_ANN_name

        ! Open the NetCDF file
        retval = nf90_open(ml_config%filename, nf90_nowrite, ncid)
        if (retval /= nf90_noerr) then
        print *, 'Error: Unable to open file'
        stop
        endif

        ! Get the variable ID
        varname = 'l1_weight'
        retval = nf90_inq_varid(ncid, varname, varid)
        if (retval == nf90_noerr) then
        ! Read the dimension values
        retval = nf90_get_var(ncid, varid, l1_weight_temp)
        ml_config%l1_weight = transpose(l1_weight_temp)
        if (retval /= nf90_noerr) then
            print *, 'Error: Unable to get l1 weight values'
            stop
        endif
        else
        print *, 'Error: l1 weight variable not found'
        endif

        varname = 'l2_weight'
        retval = nf90_inq_varid(ncid, varname, varid)
        if (retval == nf90_noerr) then
        ! Read the dimension values
        retval = nf90_get_var(ncid, varid, l2_weight_temp)
        ml_config%l2_weight = transpose(l2_weight_temp)
        if (retval /= nf90_noerr) then
            print *, 'Error: Unable to get l2 weight values'
            stop
        endif
        else
        print *, 'Error: l2 weight variable not found'
        endif

        varname = 'l3_weight'
        retval = nf90_inq_varid(ncid, varname, varid)
        if (retval == nf90_noerr) then
        ! Read the dimension values
        retval = nf90_get_var(ncid, varid, l3_weight_temp)
        ml_config%l3_weight = transpose(l3_weight_temp)
        if (retval /= nf90_noerr) then
            print *, 'Error: Unable to get l3 weight values'
            stop
        endif
        else
        print *, 'Error: l3 weight variable not found'
        endif

        varname = 'l1_bias'
        retval = nf90_inq_varid(ncid, varname, varid)
        if (retval == nf90_noerr) then
        ! Read the dimension values
        retval = nf90_get_var(ncid, varid, ml_config%l1_bias)
        if (retval /= nf90_noerr) then
            print *, 'Error: Unable to get l1 bias values'
            stop
        endif
        else
        print *, 'Error: l1 bias variable not found'
        endif

        varname = 'l2_bias'
        retval = nf90_inq_varid(ncid, varname, varid)
        if (retval == nf90_noerr) then
        ! Read the dimension values
        retval = nf90_get_var(ncid, varid, ml_config%l2_bias)
        if (retval /= nf90_noerr) then
            print *, 'Error: Unable to get l2 bias values'
            stop
        endif
        else
        print *, 'Error: l2 bias variable not found'
        endif

        varname = 'l3_bias'
        retval = nf90_inq_varid(ncid, varname, varid)
        if (retval == nf90_noerr) then
        ! Read the dimension values
        retval = nf90_get_var(ncid, varid, ml_config%l3_bias)
        if (retval /= nf90_noerr) then
            print *, 'Error: Unable to get l3 bias values'
            stop
        endif
        else
        print *, 'Error: l3 bias variable not found'
        endif

        ! Close the NetCDF file
        retval = nf90_close(ncid)
        if (retval /= nf90_noerr) then
        print *, 'Error: Unable to close file'
        stop
        endif
    end subroutine read_ANN_file

    ! 1D linear interpolation; it is guaranteed that x1 <= thisx < x2
    subroutine interpolate(x1,x2,y1,y2,thisx,thisy)
        implicit none
        real, intent(in) :: x1, x2, y1, y2, thisx
        real, intent(out) :: thisy
        thisy = (thisx-x1)/(x2-x1)*(y2-y1) + y1
    end subroutine interpolate


    ! Subroutine to find the index. 1D array (index) is greater than the given value
    Subroutine find_right_index(array1d_for_indexing, value_for_indexing, bathy_for_control, array1d_for_depth, right_index)
        implicit none
        real, intent(in) :: array1d_for_indexing(:), array1d_for_depth(:)
        integer :: array1d_i
        integer, intent(out) :: right_index
        real, intent(in) :: value_for_indexing, bathy_for_control
        right_index = 0 ! if there's so such right index, it will be 0 
        do array1d_i = 1, size(array1d_for_indexing)
            if (array1d_for_depth(array1d_i) > bathy_for_control) then
                return
            else if (array1d_for_indexing(array1d_i) > value_for_indexing) then
                right_index = array1d_i
                return
            end if
        end do
    end subroutine find_right_index

    ! Subroutine to find the index. 1D array (index) is greater than the given value, no control value needed.
    Subroutine find_right_index_clean(array1d_for_indexing, value_for_indexing, right_index)
        implicit none
        real(8), intent(in) :: array1d_for_indexing(:)
        integer :: array1d_i
        integer, intent(out) :: right_index
        real(8), intent(in) :: value_for_indexing
        right_index = 0 ! if there's so such right index, it will be 0 
        do array1d_i = 1, size(array1d_for_indexing)
            if (array1d_for_indexing(array1d_i) > value_for_indexing) then
                right_index = array1d_i
                return
            end if
        end do
    end subroutine find_right_index_clean

    ! Subroutine to compute divergence
    Subroutine compute_current_divergence(uy_left, uy_right, vx_south, vy_north, area, div)
        implicit none
        real, intent(in) :: uy_left, uy_right, vx_south, vy_north, area
        real, intent(out) :: div

        div = (uy_right - uy_left + vy_north - vx_south) / area

    end subroutine compute_current_divergence

end module MOM_oda_ml_mod
