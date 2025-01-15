!> Interfaces for MOM6 ensembles and data assimilation.
module MOM_oda_ml_mod

! This file is part of MOM6. see LICENSE.md for the license.

! MOM infrastructure
use MOM_cpu_clock, only : cpu_clock_begin, cpu_clock_end, cpu_clock_id
! MOM Modules

implicit none ; private

public :: oda_ml_init, oda_ml_end, oda_ml_inference

! Data structure to save the ML configuration, input, and output data
type, public :: ocean_oda_ml_struct
    !! Normalization parameters

    !! Weights

    !! Input features
    real :: SSH !<sea surface height (m) across ensembles
    real :: taux !<zonal wind stress
    real :: tauy !<meridional wind stress
    real :: latent !<latent heat flux
    real :: sensible !<sensile heat flux
    real :: lw !<longwave radiation flux
    real :: sw !<shortwave radiation flux
    real :: MLD !<shortwave radiation flux
    real, pointer, dimension(:) :: T=>NULL() !<layer potential temperature (degC) across ensembles
    real, pointer, dimension(:) :: S=>NULL() !<layer salinity (psu or g kg-1) across ensembles
    ! real, pointer, dimension(:) :: U=>NULL() !<layer zonal velocity (m s-1) across ensembles
    ! real, pointer, dimension(:) :: V=>NULL() !<layer meridional velocity (m s-1) across ensembles
    real, pointer, dimension(:) :: Rho=>NULL() !<layer salinity (psu or g kg-1) across ensembles

    !! Output predictions
    real, pointer, dimension(:) :: T_inc=>NULL()
    real, pointer, dimension(:) :: S_inc=>NULL()
end type ocean_oda_ml_struct

integer :: id_clock_ml_remapping
integer :: id_clock_ml_normalization
integer :: id_clock_ml_inference

#include <MOM_memory.h>

character(len=40)  :: mdl = "MOM_oda_ml" !< This module's name.

contains

subroutine oda_ml_inference(ml_CS)
    type(ocean_oda_ml_struct), pointer, intent(in) :: ml_CS
    
end subroutine oda_ml_inference

subroutine oda_ml_init(ml_CS,nk)
    type(ocean_oda_ml_struct), pointer, intent(in) :: ml_CS
    integer, intent(in) :: nk

    call init_oda_ml_features(ml_CS,nk)
    
end subroutine oda_ml_init

subroutine oda_ml_end(ml_CS)
    type(ocean_oda_ml_struct), pointer, intent(in) :: ml_CS
    
end subroutine oda_ml_end

subroutine init_oda_ml_features(ml_CS,nk)
    type(ocean_oda_ml_struct), pointer, intent(in) :: ml_CS
    integer, intent(in) :: nk

    allocate(ml_CS%T(nk),source=0.0)
    allocate(ml_CS%S(nk),source=0.0)
    ! allocate(ml_CS%U(nk),source=0.0)
    ! allocate(ml_CS%V(nk),source=0.0)
    allocate(ml_CS%Rho(nk),source=0.0)

    allocate(ml_CS%T_inc(nk),source=0.0)
    allocate(ml_CS%S_inc(nk),source=0.0)

end subroutine init_oda_ml_features

end module MOM_oda_ml_mod
