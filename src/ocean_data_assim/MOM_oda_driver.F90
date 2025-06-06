!> Interfaces for MOM6 ensembles and data assimilation.
module MOM_oda_driver_mod

! This file is part of MOM6. see LICENSE.md for the license.

! MOM infrastructure
use MOM_coms, only : PE_here, num_PEs
use MOM_coms, only : set_PElist, set_rootPE, Get_PElist, broadcast
use MOM_domains, only : domain2d, global_field, get_domain_extent
use MOM_domains, only : pass_var, redistribute_array, broadcast_domain
use MOM_diag_mediator, only : register_diag_field, diag_axis_init, post_data
use MOM_diag_mediator, only : enable_averaging, disable_averaging
use MOM_diag_mediator, only : register_scalar_field
use MOM_ensemble_manager, only : get_ensemble_id, get_ensemble_size
use MOM_ensemble_manager, only : get_ensemble_pelist, get_ensemble_filter_pelist
use MOM_error_handler, only : stdout, stdlog, MOM_error
use MOM_forcing_type, only : forcing, mech_forcing
use MOM_io, only : SINGLE_FILE
use MOM_interp_infra, only : init_extern_field, get_external_field_info
use MOM_interp_infra, only : time_interp_extern
use MOM_interpolate, only : external_field
use MOM_remapping,    only : remappingSchemesDoc
use MOM_time_manager, only : time_type, real_to_time, get_date
use MOM_time_manager, only : operator(+), operator(>=), operator(/=)
use MOM_time_manager, only : operator(==), operator(<)
use MOM_cpu_clock, only : cpu_clock_begin, cpu_clock_end, cpu_clock_id
use MOM_horizontal_regridding, only : horiz_interp_and_extrap_tracer
use MOM_restart,              only : MOM_restart_CS, register_restart_field

! ODA Modules
use mpp_domains_mod, only : mpp_get_data_domain,mpp_get_compute_domain,mpp_get_global_domain
use ocean_da_types_mod, only : grid_type, ocean_profile_type
use ocean_da_types_mod, only : ensemble_control_struct, ocean_control_struct
use ocean_da_core_mod, only : ocean_da_core_init, get_profiles
use MOM_oda_ml_mod, only: oda_ml_init, oda_ml_end, oda_ml_inference
use MOM_oda_ml_mod, only: ocean_oda_ml_data, ocean_oda_ml_config
!This preprocessing directive enables the SPEAR online ensemble data assimilation
!configuration. Existing community based APIs for data assimilation are currently
!called offline for forecast applications using information read from a MOM6 state file.
!The SPEAR configuration (https://doi.org/10.1029/2020MS002149) calculated increments
!efficiently online. A community-based set of APIs should be implemented in place
!of the CPP directive when this is available.
#ifdef ENABLE_ECDA
use eakf_oda_mod, only : ensemble_filter
#endif
use kdtree, only : kd_root !# A kd-tree object using JEDI APIs
! MOM Modules
use MOM_io, only : slasher, MOM_read_data
use MOM_diag_mediator, only : diag_ctrl
use MOM_error_handler, only : FATAL, WARNING, MOM_error, MOM_mesg, is_root_pe
use MOM_get_input, only : get_MOM_input, directories
use MOM_grid, only : ocean_grid_type, MOM_grid_init
use MOM_grid_initialize, only : set_grid_metrics
use MOM_hor_index, only : hor_index_type, hor_index_init
use MOM_dyn_horgrid, only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
use MOM_transcribe_grid, only : copy_dyngrid_to_MOM_grid, copy_MOM_grid_to_dyngrid
use MOM_fixed_initialization, only : MOM_initialize_fixed, MOM_initialize_topography
use MOM_coord_initialization, only : MOM_initialize_coord
use MOM_file_parser, only : read_param, get_param, param_file_type
use MOM_string_functions, only : lowercase
use MOM_ALE, only : ALE_CS, ALE_initThicknessToCoord, ALE_init, ALE_updateVerticalGridType
use MOM_domains, only : MOM_domains_init, MOM_domain_type, clone_MOM_domain
use MOM_remapping, only : remapping_CS, initialize_remapping, remapping_core_h
use MOM_regridding, only : regridding_CS, initialize_regridding
use MOM_regridding, only : regridding_main, set_regrid_params
use MOM_unit_scaling, only : unit_scale_type, unit_scaling_init
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type, verticalGridInit

implicit none ; private

public :: init_oda, oda_end, set_prior_tracer, get_posterior_tracer
public :: set_analysis_time, oda, apply_oda_tracer_increments
public :: set_oda_restart_fields, init_oda_diags

!>@{ CPU time clock ID
integer :: id_clock_oda_init
integer :: id_clock_get_prior
integer :: id_clock_bias_adjustment
integer :: id_clock_ml_bias_correction
integer :: id_clock_ensemble_filter
integer :: id_clock_apply_increments
!>@}

#include <MOM_memory.h>

!> A structure with a pointer to a domain2d, to allow for the creation of arrays of pointers.
type :: ptr_mpp_domain
  type(domain2d), pointer :: mpp_domain => NULL() !< pointer to a domain2d
end type ptr_mpp_domain

!> A structure containing integer handles for bias adjustment of tracers
type :: INC_CS
  integer :: fldno = 0 !< The number of tracers
  type(external_field) :: T  !< The handle for the temperature file
  type(external_field) :: S  !< The handle for the salinity file
end type INC_CS

!> Control structure that contains a transpose of the ocean state across ensemble members.
type, public :: ODA_CS ; private
  type(ensemble_control_struct), pointer :: Ocean_prior=> NULL() !< ensemble ocean prior states in DA space
  real :: prior_ave_counter
  type(ensemble_control_struct), pointer :: Ocean_posterior=> NULL() !< ensemble ocean posterior states
                                                                  !! or increments to prior in DA space
  ! type(ensemble_control_struct), pointer :: Ocean_increment=> NULL() !< A separate structure for
  !                                                                 !! increment diagnostics
  type(ocean_control_struct), pointer :: Ocean_background_ave=> NULL() !< ocean averaged prior states in model space
  type(ocean_oda_ml_data), pointer :: ml_data => NULL()
  type(ocean_oda_ml_config), pointer :: ml_config => NULL()
  integer :: nk !< number of vertical layers used for DA
  type(ocean_grid_type), pointer :: Grid => NULL() !< MOM6 grid type and decomposition for the DA
  type(ocean_grid_type), pointer :: model_G => NULL() !< MOM6 grid type and decomposition for the model
  type(MOM_domain_type), pointer, dimension(:) :: domains => NULL() !< Pointer to mpp_domain objects
                                                                       !! for ensemble members
  type(verticalGrid_type), pointer :: GV => NULL() !< vertical grid for DA
  type(verticalGrid_type), pointer :: model_GV => NULL() !< vertical grid for DA
  type(unit_scale_type), pointer :: &
    US => NULL()    !< structure containing various unit conversion factors for DA

  type(domain2d), pointer :: mpp_domain => NULL() !< Pointer to a mpp domain object for DA
  type(grid_type), pointer :: oda_grid !< local tracer grid
  real, pointer, dimension(:,:,:) :: h => NULL() !<layer thicknesses [H ~> m or kg m-2] for DA
  real, pointer, dimension(:,:,:) :: T_tend => NULL() !<layer temperature tendency from DA [C T-1 ~> degC s-1]
  real, pointer, dimension(:,:,:) :: S_tend => NULL() !<layer salinity tendency from DA [S T-1 ~> ppt s-1]
  real, pointer, dimension(:,:,:) :: T_bc_tend => NULL() !< The layer temperature tendency due
                                                         !! to bias adjustment [C T-1 ~> degC s-1]
  real, pointer, dimension(:,:,:) :: S_bc_tend => NULL() !< The layer salinity tendency due
                                                         !! to bias adjustment [S T-1 ~> ppt s-1]
  real, pointer, dimension(:,:,:) :: T_ml_tend => NULL() !< The layer temperature tendency due
                                                         !! to ML bias adjustment [C T-1 ~> degC s-1]
  real, pointer, dimension(:,:,:) :: S_ml_tend => NULL() !< The layer salinity tendency due
                                                         !! to bias adjustment [S T-1 ~> ppt s-1]
  integer :: ni          !< global i-direction grid size
  integer :: nj          !< global j-direction grid size
  logical :: reentrant_x !< grid is reentrant in the x direction
  logical :: reentrant_y !< grid is reentrant in the y direction
  logical :: tripolar_N !< grid is folded at its north edge
  logical :: symmetric !< Values at C-grid locations are symmetric
  logical :: use_basin_mask !< If true, use a basin file to delineate weakly coupled ocean basins
  logical :: do_T_bias_adjustment !< If true, use spatio-temporally varying climatological tendency
                                !! adjustment for Temperature and Salinity
  logical :: do_S_bias_adjustment !< If true, use spatio-temporally varying climatological tendency
  real :: T_bias_adjustment_multiplier !< A scaling for the bias adjustment
  real :: S_bias_adjustment_multiplier !< A scaling for the bias adjustment
  logical :: do_T_ml_bias_adjustment !< If true, use machine learning-trained tendency
                                !! adjustment for Temperature and Salinity
  logical :: do_S_ml_bias_adjustment !< If true, use machine learning-trained tendency
  real :: T_ml_bias_adjustment_multiplier !< A scaling for the bias adjustment
  real :: S_ml_bias_adjustment_multiplier !< A scaling for the bias adjustment
  integer :: assim_method !< Method: NO_ASSIM,EAKF_ASSIM or OI_ASSIM
  integer :: ensemble_size !< Size of the ensemble
  integer :: ensemble_id = 0 !< id of the current ensemble member
  integer, pointer, dimension(:,:) :: ensemble_pelist !< PE list for ensemble members
  integer, pointer, dimension(:) :: filter_pelist !< PE list for ensemble members
  real :: assim_interval !< analysis interval [ T ~> s]
  real :: prior_interval !< prior accumulation interval [ T ~> s]
  real :: apply_interval !< apply increments interval [ T ~> s]
  ! Profiles local to the analysis domain
  type(ocean_profile_type), pointer :: Profiles => NULL() !< pointer to linked list of all available profiles
  type(ocean_profile_type), pointer :: CProfiles => NULL()!< pointer to linked list of current profiles
  type(kd_root), pointer :: kdroot => NULL() !< A structure for storing nearest neighbors
  type(ALE_CS), pointer :: ALE_CS=>NULL() !< ALE control structure for DA
  logical :: use_ALE_algorithm !< true is using ALE remapping
  type(regridding_CS) :: regridCS !< ALE control structure for regridding
  type(remapping_CS) :: remapCS !< ALE control structure for remapping
  type(time_type) :: Time !< Current Analysis time
  type(time_type) :: Prior_Time !< Current Prior time for time averaging
  type(time_type) :: Apply_Time !< Current Prior time for time averaging
  type(diag_ctrl), pointer :: diag_cs=> NULL() !<Pointer to diagnostics control structure
  type(INC_CS) :: INC_CS !< A Structure containing integer file handles for bias adjustment
  integer :: id_inc_t = -1 !< A diagnostic handle for the temperature climatological adjustment
  integer :: id_inc_s = -1 !< A diagnostic handle for the salinity climatological adjustment
  integer :: id_inc_ml_t = -1 !< A diagnostic handle for the temperature climatological adjustment
  integer :: id_inc_ml_s = -1 !< A diagnostic handle for the salinity climatological adjustment
  integer :: id_inc_t_z = -1 !< A diagnostic handle for the temperature climatological adjustment
  integer :: id_inc_s_z = -1 !< A diagnostic handle for the salinity climatological adjustment
  integer :: id_inc_ml_t_z = -1 !< A diagnostic handle for the temperature climatological adjustment
  integer :: id_inc_ml_s_z = -1 !< A diagnostic handle for the salinity climatological adjustment
  integer :: id_prior_count = -1
  integer :: id_prior_t = -1, id_prior_s = -1, id_prior_u = -1, id_prior_v = -1
  integer :: id_prior_t_z = -1, id_prior_s_z = -1, id_prior_u_z = -1, id_prior_v_z = -1
  integer :: id_prior_ssh = -1, id_prior_taux = -1, id_prior_tauy = -1
  integer :: id_prior_sw = -1, id_prior_lw = -1, id_prior_latent = -1, id_prior_sensible = -1
  integer :: answer_date    !< The vintage of the order of arithmetic and expressions in the
                            !! remapping invoked by the ODA driver.  Values below 20190101 recover
                            !! the answers from the end of 2018, while higher values use updated
                            !! and more robust forms of the same expressions.
end type ODA_CS


!>@{  DA parameters
integer, parameter :: NO_ASSIM = 0, OI_ASSIM=1, EAKF_ASSIM=2
!>@}
character(len=40)  :: mdl = "MOM_oda_driver" !< This module's name.

contains

!> initialize First_guess (prior) and Analysis grid
!! information for all ensemble members
subroutine init_oda(Time, G, GV, US, CS)

  type(time_type), intent(in) :: Time !< The current model time.
  type(ocean_grid_type), pointer, intent(in) :: G !< domain and grid information for ocean model
  type(verticalGrid_type), pointer, intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  ! type(diag_ctrl), target, intent(inout) :: diag_CS !< A pointer to a diagnostic control structure
  type(ODA_CS), pointer, intent(inout) :: CS  !< The DA control structure

! Local variables
  type(thermo_var_ptrs) :: tv_dummy
  type(dyn_horgrid_type), pointer :: dG=> NULL()
  type(hor_index_type), pointer :: HI=> NULL()
  type(directories) :: dirs

  type(grid_type), pointer :: T_grid => NULL() !< global tracer grid
  type(param_file_type) :: PF
  integer :: n
  ! integer :: isd, ied, jsd, jed
  ! integer :: is_oda, ie_oda, js_oda, je_oda
  ! integer :: isd_oda, ied_oda, jsd_oda, jed_oda
  integer, dimension(4) :: fld_sz
  character(len=32) :: assim_method
  integer :: npes_pm, ens_info(6)
  character(len=30) :: coord_mode
  character(len=200) :: inputdir, basin_file
  character(len=80) :: basin_var
  character(len=80) :: remap_scheme
  character(len=80) :: bias_correction_file, inc_file
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags.
  character(len=160) :: mesg
  integer :: yr, mon, day, hr, min, sec
  
  if (associated(CS)) call MOM_error(FATAL, 'Calling oda_init with associated control structure')
  allocate(CS)

  id_clock_get_prior=cpu_clock_id('(ODA getting prior)')
  id_clock_bias_adjustment=cpu_clock_id('(ODA getting bias correction)')
  id_clock_ml_bias_correction=cpu_clock_id('(ML inference of bias correction)')
  id_clock_ensemble_filter=cpu_clock_id('(ODA ensemble filter)')
  id_clock_apply_increments=cpu_clock_id('(ODA applying increments)')
  id_clock_oda_init=cpu_clock_id('(ODA initialization)')
  call cpu_clock_begin(id_clock_oda_init)

! Use ens1 parameters , this could be changed at a later time
! if it were desirable to have alternate parameters, e.g. for the grid
! for the analysis
  call get_MOM_input(PF,dirs,ensemble_num=0)
  call unit_scaling_init(PF, CS%US)

  call get_param(PF, mdl, "ASSIM_METHOD", assim_method,  &
       "String which determines the data assimilation method "//&
       "Valid methods are: \'EAKF\',\'OI\', and \'NO_ASSIM\'", default='NO_ASSIM')
  call get_param(PF, mdl, "ASSIM_INTERVAL", CS%assim_interval,  &
       "data assimilation update interval in hours",default=-1.0,units="hours",scale=3600.*US%s_to_T)
  call get_param(PF, mdl, "PRIOR_INTERVAL", CS%prior_interval,  &
       "prior averaging interval in hours",default=2.0,units="hours",scale=3600.*US%s_to_T)
  call get_param(PF, mdl, "APPLY_INTERVAL", CS%apply_interval,  &
       "increment application interval in hours",default=2.0,units="hours",scale=3600.*US%s_to_T)  
  if (CS%assim_interval < 0.) then
     call get_param(PF, mdl, "ASSIM_FREQUENCY", CS%assim_interval,  &
          "data assimilation update  in hours. This parameter name will \n"//&
          "be deprecated in the future. ASSIM_INTERVAL should be used instead.",default=-1.0, &
          units="hours",scale=3600.*US%s_to_T)
  endif

  call get_param(PF, mdl, "USE_REGRIDDING", CS%use_ALE_algorithm , &
                "If True, use the ALE algorithm (regridding/remapping).\n"//&
                "If False, use the layered isopycnal algorithm.", default=.false. )
  call get_param(PF, mdl, "REENTRANT_X", CS%reentrant_x, &
       "If true, the domain is zonally reentrant.", default=.true.)
  call get_param(PF, mdl, "REENTRANT_Y", CS%reentrant_y, &
       "If true, the domain is meridionally reentrant.", &
       default=.false.)
  call get_param(PF, mdl, "TRIPOLAR_N", CS%tripolar_N, &
       "Use tripolar connectivity at the northern edge of the "//&
       "domain.  With TRIPOLAR_N, NIGLOBAL must be even.", &
       default=.false.)
  call get_param(PF, mdl, "APPLY_TEMP_TENDENCY_ADJUSTMENT", CS%do_T_bias_adjustment, &
       "If true, add a spatio-temporally varying climatological adjustment "//&
       "to temperature.", &
       default=.false.)
  call get_param(PF, mdl, "APPLY_SALT_TENDENCY_ADJUSTMENT", CS%do_S_bias_adjustment, &
       "If true, add a spatio-temporally varying climatological adjustment "//&
       "to salinity.", &
       default=.false.)
  if (CS%do_T_bias_adjustment) then
    call get_param(PF, mdl, "TEMP_ADJUSTMENT_FACTOR", CS%T_bias_adjustment_multiplier, &
       "A multiplicative scaling factor for the climatological tracer tendency adjustment ", &
       units="nondim", default=0.0)
  endif
  if (CS%do_S_bias_adjustment) then
    call get_param(PF, mdl, "SALT_ADJUSTMENT_FACTOR", CS%S_bias_adjustment_multiplier, &
       "A multiplicative scaling factor for the climatological tracer tendency adjustment ", &
       units="nondim", default=0.0)
  endif
  write(mesg,*) 'OTA adjustment multiplier', CS%T_bias_adjustment_multiplier, CS%S_bias_adjustment_multiplier
  call MOM_mesg("ODA init: "//trim(mesg))
  call get_param(PF, mdl, "APPLY_ML_TEMP_TENDENCY_ADJUSTMENT", CS%do_T_ml_bias_adjustment, &
       "If true, add a machine learning-trained adjustment "//&
       "to temperature.", &
       default=.false.)
  call get_param(PF, mdl, "APPLY_ML_SALT_TENDENCY_ADJUSTMENT", CS%do_S_ml_bias_adjustment, &
       "If true, add a machine learning-trained adjustment "//&
       "to salinity.", &
       default=.false.)
  call get_param(PF, mdl, "ML_TEMP_ADJUSTMENT_FACTOR", CS%T_ml_bias_adjustment_multiplier, &
      "A multiplicative scaling factor for the machine learning tracer tendency adjustment ", &
      units="nondim", default=0.0)
  call get_param(PF, mdl, "ML_SALT_ADJUSTMENT_FACTOR", CS%S_ml_bias_adjustment_multiplier, &
      "A multiplicative scaling factor for the machine learning tracer tendency adjustment ", &
      units="nondim", default=0.0)
  write(mesg,*) 'ML adjustment multiplier', CS%T_ml_bias_adjustment_multiplier, CS%S_ml_bias_adjustment_multiplier
  call MOM_mesg("ODA init: "//trim(mesg))
  call get_param(PF, mdl, "USE_BASIN_MASK", CS%use_basin_mask, &
       "If true, add a basin mask to delineate weakly connected "//&
       "ocean basins for the purpose of data assimilation.", &
       default=.false.)

  call get_param(PF, mdl, "NIGLOBAL", CS%ni, &
       "The total number of thickness grid points in the "//&
       "x-direction in the physical domain.")
  call get_param(PF, mdl, "NJGLOBAL", CS%nj, &
       "The total number of thickness grid points in the "//&
       "y-direction in the physical domain.")
  call get_param(PF, mdl, "INPUTDIR", inputdir)
  call get_param(PF, mdl, "ODA_REMAPPING_SCHEME", remap_scheme, &
                 "This sets the reconstruction scheme used "//&
                 "for vertical remapping for all variables. "//&
                 "It can be one of the following schemes: "//&
                 trim(remappingSchemesDoc), default="PPM_H4")
  call get_param(PF, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)
  call get_param(PF, mdl, "ODA_ANSWER_DATE", CS%answer_date, &
               "The vintage of the order of arithmetic and expressions used by the ODA driver "//&
               "Values below 20190101 recover the answers from the end of 2018, while higher "//&
               "values use updated and more robust forms of the same expressions.", &
               default=default_answer_date, do_not_log=.not.GV%Boussinesq)
  if (.not.GV%Boussinesq) CS%answer_date = max(CS%answer_date, 20230701)
  inputdir = slasher(inputdir)

  select case(lowercase(trim(assim_method)))
    case('eakf')
      CS%assim_method = EAKF_ASSIM
    case('oi')
      CS%assim_method = OI_ASSIM
    case('no_assim')
      CS%assim_method = NO_ASSIM
    case default
      call MOM_error(FATAL, "Invalid assimilation method provided")
  end select

  ens_info = get_ensemble_size()
  CS%ensemble_size = ens_info(1)
  npes_pm=ens_info(3)
  CS%ensemble_id = get_ensemble_id()
  !! Switch to global pelist
  allocate(CS%ensemble_pelist(CS%ensemble_size,npes_pm))
  allocate(CS%filter_pelist(CS%ensemble_size*npes_pm))
  call get_ensemble_pelist(CS%ensemble_pelist, 'ocean')
  call get_ensemble_filter_pelist(CS%filter_pelist, 'ocean')

  call set_PElist(CS%filter_pelist)

  allocate(CS%domains(CS%ensemble_size))
  CS%domains(CS%ensemble_id)%mpp_domain => G%Domain%mpp_domain ! this should go away
  do n=1,CS%ensemble_size
    if (.not. associated(CS%domains(n)%mpp_domain)) allocate(CS%domains(n)%mpp_domain)
    call set_rootPE(CS%ensemble_pelist(n,1)) ! this line is not in Feiyu's version (needed?)
    call broadcast_domain(CS%domains(n)%mpp_domain)
  enddo
  call set_rootPE(CS%filter_pelist(1)) ! this line is not in Feiyu's version (needed?)
  CS%model_G => G
  CS%model_GV => GV
  allocate(CS%Grid)
  ! params NIHALO_ODA, NJHALO_ODA set the DA halo size
  call MOM_domains_init(CS%Grid%Domain, PF, param_suffix='_ODA', US=CS%US)
  allocate(HI)
  call hor_index_init(CS%Grid%Domain, HI, PF)
  call verticalGridInit( PF, CS%GV, CS%US )
  allocate(dG)
  call create_dyn_horgrid(dG, HI)
  call clone_MOM_domain(CS%Grid%Domain, dG%Domain,symmetric=.false.)
  call set_grid_metrics(dG, PF, CS%US)
  call MOM_initialize_topography(dG%bathyT, dG%max_depth, dG, PF, CS%US)
  call MOM_initialize_coord(CS%GV, CS%US, PF, tv_dummy, dG%max_depth)
  call ALE_init(PF, CS%GV, CS%US, dG%max_depth, CS%ALE_CS)
  call MOM_grid_init(CS%Grid, PF, global_indexing=.false.)
  call ALE_updateVerticalGridType(CS%ALE_CS, CS%GV)
  call copy_dyngrid_to_MOM_grid(dG, CS%Grid, CS%US)
  CS%mpp_domain => CS%Grid%Domain%mpp_domain
  CS%Grid%ke = CS%GV%ke
  CS%nk = CS%GV%ke
  ! initialize storage for prior and posterior
  if (.NOT. CS%assim_method == NO_ASSIM) then
    allocate(CS%Ocean_prior)
    call init_ocean_ensemble(CS%Ocean_prior,CS%Grid,CS%GV,CS%ensemble_size)
    allocate(CS%Ocean_posterior)
    call init_ocean_ensemble(CS%Ocean_posterior,CS%Grid,CS%GV,CS%ensemble_size)
  !   allocate(CS%Ocean_increment)
  !   call init_ocean_ensemble(CS%Ocean_increment,CS%Grid,CS%GV,CS%ensemble_size)
  endif

  call get_param(PF, 'oda_driver', "REGRIDDING_COORDINATE_MODE", coord_mode, &
       "Coordinate mode for vertical regridding.", &
       default="ZSTAR", fail_if_missing=.false.)
  call initialize_regridding(CS%regridCS, CS%GV, CS%US, dG%max_depth,PF,'oda_driver',coord_mode,'','')
  call initialize_remapping(CS%remapCS,remap_scheme,answer_date = CS%answer_date)
  call set_regrid_params(CS%regridCS, min_thickness=0.)
  
  allocate(CS%oda_grid)
  CS%oda_grid%x => CS%Grid%geolonT
  CS%oda_grid%y => CS%Grid%geolatT
  CS%oda_grid%bathyT => CS%Grid%bathyT

  if (CS%use_basin_mask) then
    call get_param(PF, 'oda_driver', "BASIN_FILE", basin_file, &
          "A file in which to find the basin masks.", default="basin.nc")
    basin_file = trim(inputdir) // trim(basin_file)
    call get_param(PF, 'oda_driver', "BASIN_VAR", basin_var, &
          "The basin mask variable in BASIN_FILE.", default="basin")
    ! Need different data domain indices for the ODA ensemble basin mask.
    ! call get_domain_extent(CS%Grid%Domain,is_oda,ie_oda,js_oda,je_oda,isd_oda,ied_oda,jsd_oda,jed_oda)
    ! allocate(CS%oda_grid%basin_mask(isd_oda:ied_oda,jsd_oda:jed_oda), source=0.0)
    allocate(CS%oda_grid%basin_mask(CS%Grid%isd:CS%Grid%ied,CS%Grid%jsd:CS%Grid%jed), source=0.0)
    call MOM_read_data(basin_file, basin_var, CS%oda_grid%basin_mask, CS%Grid%domain, timelevel=1)
  endif

  if (.not. associated(CS%h)) then
    allocate(CS%h(G%isd:G%ied,G%jsd:G%jed,CS%GV%ke), source=CS%GV%Angstrom_H)
    ! assign thicknesses
    call ALE_initThicknessToCoord(CS%ALE_CS, G, CS%GV, CS%h)
  endif

  !!  get global grid information from ocean model needed for ODA initialization
  call set_up_global_tgrid(T_grid, CS, G)
  call ocean_da_core_init(CS%mpp_domain, T_grid, CS%Profiles, Time)
  deallocate(T_grid)

  CS%Time = Time
  CS%Prior_Time = Time
  CS%Apply_Time = Time
  call get_date(Time, yr, mon, day, hr, min, sec)
  write(mesg,*)  'Model Time: ', yr, mon, day, hr, min, sec
  call MOM_mesg("ODA_INIT: "//trim(mesg))

  !! switch back to ensemble member pelist
  call set_PElist(CS%ensemble_pelist(CS%ensemble_id,:))

  allocate(CS%Ocean_background_ave)
  call init_ocean_background(CS,CS%Ocean_background_ave,G,CS%GV)
  CS%prior_ave_counter = CS%assim_interval / CS%prior_interval

  allocate(CS%T_tend(G%isd:G%ied,G%jsd:G%jed,CS%GV%ke), source=0.0)
  allocate(CS%S_tend(G%isd:G%ied,G%jsd:G%jed,CS%GV%ke), source=0.0)

  if (CS%do_T_bias_adjustment .or. CS%do_S_bias_adjustment) then
    call get_param(PF, mdl, "TEMP_SALT_ADJUSTMENT_FILE", bias_correction_file,  &
                "The name of the file containing temperature and salinity "//&
                "tendency adjustments", default='temp_salt_adjustment.nc')

    inc_file = trim(inputdir) // trim(bias_correction_file)
    CS%INC_CS%T = init_extern_field(inc_file, "temp_increment", &
          correct_leap_year_inconsistency=.true.,verbose=.true.,domain=G%Domain%mpp_domain)
    CS%INC_CS%S = init_extern_field(inc_file, "salt_increment", &
          correct_leap_year_inconsistency=.true.,verbose=.true.,domain=G%Domain%mpp_domain)
    call get_external_field_info(CS%INC_CS%T, size=fld_sz)
    CS%INC_CS%fldno = 2
    if (CS%nk /= fld_sz(3)) call MOM_error(FATAL,'Increment levels /= ODA levels')

    allocate(CS%T_bc_tend(G%isd:G%ied,G%jsd:G%jed,CS%GV%ke), source=0.0)
    allocate(CS%S_bc_tend(G%isd:G%ied,G%jsd:G%jed,CS%GV%ke), source=0.0)

  endif

  if (CS%do_T_ml_bias_adjustment .or. CS%do_S_ml_bias_adjustment) then

    allocate(CS%ml_data)
    allocate(CS%ml_config)
    call oda_ml_init(CS%ml_config, CS%ml_data, CS%GV)

    allocate(CS%T_ml_tend(G%isd:G%ied,G%jsd:G%jed,CS%GV%ke), source=0.0)
    allocate(CS%S_ml_tend(G%isd:G%ied,G%jsd:G%jed,CS%GV%ke), source=0.0)

  endif

  call cpu_clock_end(id_clock_oda_init)

!  if (CS%write_obs) then
!    temp_fid = open_profile_file("temp_"//trim(obs_file))
!    salt_fid = open_profile_file("salt_"//trim(obs_file))
!  end if

end subroutine init_oda

subroutine init_oda_diags(Time, US, diag_CS, CS)

  type(time_type), intent(in) :: Time !< The current model time.
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(diag_ctrl), target, intent(inout) :: diag_CS !< A pointer to a diagnostic control structure
  type(ODA_CS), pointer, intent(inout) :: CS  !< The DA control structure

  call cpu_clock_begin(id_clock_oda_init)

  ! set up diag variables for analysis increments
  CS%diag_CS => diag_CS

  CS%id_inc_t = register_diag_field('ocean_model', 'temp_increment', diag_CS%axesTL, &
      Time, 'Ocean potential temperature increments', 'degC', conversion=US%C_to_degC)
  CS%id_inc_s = register_diag_field('ocean_model', 'salt_increment', diag_CS%axesTL, &
      Time, 'Ocean salinity increments', 'psu', conversion=US%S_to_ppt)
  CS%id_inc_t_z = register_diag_field('ocean_model', 'temp_increment_z', diag_CS%axesTZ, &
      Time, 'Ocean potential temperature increments', 'degC', conversion=US%C_to_degC)
  CS%id_inc_s_z = register_diag_field('ocean_model', 'salt_increment_z', diag_CS%axesTZ, &
      Time, 'Ocean salinity increments', 'psu', conversion=US%S_to_ppt)

  CS%id_prior_count = register_scalar_field('ocean_model', 'prior_count', Time, diag_CS, 'Prior Count' )

  CS%id_prior_t = register_diag_field('ocean_model', 'thetao_prior', diag_CS%axesTL, &
      Time, 'Accumulated ocean potential temperature for DA/ML', 'degC', conversion=US%C_to_degC)
  CS%id_prior_s = register_diag_field('ocean_model', 'so_prior', diag_CS%axesTL, &
      Time, 'Accumulated ocean salinity for DA/ML', 'psu', conversion=US%S_to_ppt)
  CS%id_prior_ssh = register_diag_field('ocean_model', 'SSH_prior', diag_CS%axesT1, &
      Time, 'Accumulated Sea Surface Height for DA/ML', 'm', conversion=US%Z_to_m)

  CS%id_prior_u = register_diag_field('ocean_model', 'uo_prior', diag_CS%axesCuL, &
    Time, 'Accumulated ocean zonal velocity for DA/ML', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_prior_v = register_diag_field('ocean_model', 'vo_prior', diag_CS%axesCvL, &
    Time, 'Accumulated ocean meridional velocity for DA/ML', 'm s-1', conversion=US%L_T_to_m_s)

  CS%id_prior_t_z = register_diag_field('ocean_model', 'thetao_prior_z', diag_CS%axesTZ, &
    Time, 'Accumulated ocean potential temperature for DA/ML', 'degC', conversion=US%C_to_degC)
  CS%id_prior_s_z = register_diag_field('ocean_model', 'so_prior_z', diag_CS%axesTZ, &
    Time, 'Accumulated ocean salinity for DA/ML', 'psu', conversion=US%S_to_ppt)

  CS%id_prior_u_z = register_diag_field('ocean_model', 'uo_prior_z', diag_CS%axesCuZ, &
    Time, 'Accumulated ocean zonal velocity for DA/ML', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_prior_v_z = register_diag_field('ocean_model', 'vo_prior_z', diag_CS%axesCvZ, &
    Time, 'Accumulated ocean meridional velocity for DA/ML', 'm s-1', conversion=US%L_T_to_m_s)

  CS%id_prior_taux = register_diag_field('ocean_model', 'taux_prior', diag_CS%axesCu1, &
    Time, 'Accumulated zonal surface stress for DA/ML', 'Pa', conversion=US%RLZ_T2_to_Pa)
  CS%id_prior_tauy = register_diag_field('ocean_model', 'tauy_prior', diag_CS%axesCv1, &
    Time, 'Accumulated meridional surface stress for DA/ML', 'Pa', conversion=US%RLZ_T2_to_Pa)

  CS%id_prior_sw = register_diag_field('ocean_model', 'SW_prior', diag_CS%axesT1, &
    Time, 'Accumulated shortwave radiation flux into ocean for DA/ML', 'W m-2', conversion=US%QRZ_T_to_W_m2)
  CS%id_prior_lw = register_diag_field('ocean_model', 'LW_prior', diag_CS%axesT1, &
    Time, 'Accumulated longwave radiation flux into ocean for DA/ML', 'W m-2', conversion=US%QRZ_T_to_W_m2)

  CS%id_prior_latent = register_diag_field('ocean_model', 'latent_prior', diag_CS%axesT1, &
    Time, 'Accumulated latent heat flux into ocean for DA/ML', 'W m-2', conversion=US%QRZ_T_to_W_m2)
  CS%id_prior_sensible = register_diag_field('ocean_model', 'sensible_prior', diag_CS%axesT1, &
    Time, 'Accumulated sensible heat flux into ocean for DA/ML', 'W m-2', conversion=US%QRZ_T_to_W_m2)

  if (CS%do_T_ml_bias_adjustment .or. CS%do_S_ml_bias_adjustment) then
    if (CS%do_T_ml_bias_adjustment) then
      CS%id_inc_ml_t = register_diag_field('ocean_model', 'temp_ml_increment', diag_CS%axesTL, &
        Time, 'Ocean potential temperature increments predicted by ML', 'degC', conversion=US%C_to_degC)
      CS%id_inc_ml_t_z = register_diag_field('ocean_model', 'temp_ml_increment_z', diag_CS%axesTZ, &
      Time, 'Ocean potential temperature increments predicted by ML', 'degC', conversion=US%C_to_degC)
    endif

    if (CS%do_S_ml_bias_adjustment) then
      CS%id_inc_ml_s = register_diag_field('ocean_model', 'salt_ml_increment', diag_CS%axesTL, &
        Time, 'Ocean salinity increments predicted by ML', 'psu', conversion=US%S_to_ppt)
      CS%id_inc_ml_s_z = register_diag_field('ocean_model', 'salt_ml_increment_z', diag_CS%axesTZ, &
        Time, 'Ocean salinity increments predicted by ML', 'psu', conversion=US%S_to_ppt)
    endif
  endif

  call cpu_clock_end(id_clock_oda_init)

end subroutine init_oda_diags

!> Copy ensemble member tracers to ensemble vector.
subroutine set_prior_tracer(Time, G, GV, h, tv, model_u, model_v, model_ssh, fluxes, forces, CS)
  type(time_type), intent(in)    :: Time !< The current model time
  type(ocean_grid_type), pointer :: G !< domain and grid information for ocean model
  type(verticalGrid_type),               intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h   !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                 intent(in) :: tv   !< A structure pointing to various thermodynamic variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)   :: model_u
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)   :: model_v
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: model_ssh
  type(mech_forcing), intent(in) :: forces !< A structure with the driving mechanical forces
  type(forcing), intent(in) :: fluxes  !< A structure with pointers to themodynamic,
                                                     !! tracer and mass exchange forcing fields

  type(ODA_CS), pointer :: CS !< ocean DA control structure
  real, dimension(SZI_(G),SZJ_(G),CS%nk) :: T  ! Temperature on the analysis grid [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),CS%nk) :: S  ! Salinity on the analysis grid [S ~> ppt]
  real, dimension(SZIB_(G),SZJ_(G),CS%nk) :: U      !< zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),CS%nk) :: V      !< meridional velocity [L T-1 ~> m s-1]

  integer :: i, j, m
  integer :: isc, iec, jsc, jec, iscB, iecB, jscB, jecB
  real :: h_neglect, h_neglect_edge                 ! small thicknesses [H ~> m or kg m-2]
  integer :: isd, ied, jsd, jed
  character(len=160) :: mesg  ! The text of an error message
  integer :: yr, mon, day, hr, min, sec

  ! return if not time for averaging prio
  if (Time < CS%Prior_Time) return
  call cpu_clock_begin(id_clock_get_prior)

  if (.not. associated(CS%Grid)) call MOM_ERROR(FATAL,'ODA_CS ensemble horizontal grid not associated')
  if (.not. associated(CS%GV)) call MOM_ERROR(FATAL,'ODA_CS ensemble vertical grid not associated')

  if (CS%answer_date >= 20190101) then
    h_neglect = GV%H_subroundoff ; h_neglect_edge = GV%H_subroundoff
  elseif (GV%Boussinesq) then
    h_neglect = GV%m_to_H * 1.0e-30 ; h_neglect_edge = GV%m_to_H * 1.0e-10
  else
    h_neglect = GV%kg_m2_to_H * 1.0e-30 ; h_neglect_edge = GV%kg_m2_to_H * 1.0e-10
  endif

  isc  = G%isc ; iec  = G%iec  ; jsc  = G%jsc  ; jec  = G%jec
  iscB  = G%iscB ; iecB  = G%iecB  ; jscB  = G%jscB  ; jecB  = G%jecB

  ! remap temperature and salinity from the ensemble member to the analysis grid
  do j=jsc,jec ; do i=isc,iec
    call remapping_core_h(CS%remapCS, GV%ke, h(i,j,:), tv%T(i,j,:), &
         CS%nk, CS%h(i,j,:), T(i,j,:), h_neglect, h_neglect_edge)
    call remapping_core_h(CS%remapCS, GV%ke, h(i,j,:), tv%S(i,j,:), &
         CS%nk, CS%h(i,j,:), S(i,j,:), h_neglect, h_neglect_edge)
  enddo ; enddo
  ! remap U and V from the ensemble member to the analysis grid
  ! do j=jsc,jec ; do i=iscB,iecB
  !   call remapping_core_h(CS%remapCS, GV%ke, h(i,j,:), model_u(i,j,:), &
  !        CS%nk, CS%h(i,j,:), U(i,j,:), h_neglect, h_neglect_edge)
  ! enddo ; enddo
  ! do j=jscB,jecB ; do i=isc,iec
  !   call remapping_core_h(CS%remapCS, GV%ke, h(i,j,:), model_v(i,j,:), &
  !        CS%nk, CS%h(i,j,:), V(i,j,:), h_neglect, h_neglect_edge)
  ! enddo ; enddo
  
  ! cast ensemble members to the analysis domain
  if (CS%prior_ave_counter < 0.5) then
    CS%Ocean_background_ave%T = 0.0
    CS%Ocean_background_ave%S = 0.0
    ! CS%Ocean_background_ave%SSH = 0.0

    if (CS%do_T_ml_bias_adjustment .or. CS%do_S_ml_bias_adjustment) then
      CS%Ocean_background_ave%U = 0.0
      CS%Ocean_background_ave%V = 0.0
      CS%Ocean_background_ave%taux = 0.0
      CS%Ocean_background_ave%tauy = 0.0
      CS%Ocean_background_ave%latent = 0.0
      CS%Ocean_background_ave%sensible = 0.0
      CS%Ocean_background_ave%lw = 0.0
      CS%Ocean_background_ave%sw = 0.0
    endif

    call MOM_mesg("ODA Reset background accumulation")
  endif
  
  CS%Ocean_background_ave%T(isc:iec,jsc:jec,:) = CS%Ocean_background_ave%T(isc:iec,jsc:jec,:) + T(isc:iec,jsc:jec,:)
  CS%Ocean_background_ave%S(isc:iec,jsc:jec,:) = CS%Ocean_background_ave%S(isc:iec,jsc:jec,:) + S(isc:iec,jsc:jec,:)
  ! CS%Ocean_background_ave%SSH(isc:iec,jsc:jec) = CS%Ocean_background_ave%SSH(isc:iec,jsc:jec) + model_ssh(isc:iec,jsc:jec)

  if (CS%do_T_ml_bias_adjustment .or. CS%do_S_ml_bias_adjustment) then
    CS%Ocean_background_ave%U(iscB:iecB,jsc:jec,:) = CS%Ocean_background_ave%U(iscB:iecB,jsc:jec,:) + U(iscB:iecB,jsc:jec,:)
    CS%Ocean_background_ave%V(isc:iec,jscB:jecB,:) = CS%Ocean_background_ave%V(isc:iec,jscB:jecB,:) + V(isc:iec,jscB:jecB,:)
    CS%Ocean_background_ave%taux(iscB:iecB,jsc:jec) = CS%Ocean_background_ave%taux(iscB:iecB,jsc:jec) + forces%taux(iscB:iecB,jsc:jec)
    CS%Ocean_background_ave%tauy(isc:iec,jscB:jecB) = CS%Ocean_background_ave%tauy(isc:iec,jscB:jecB) + forces%tauy(isc:iec,jscB:jecB)
    CS%Ocean_background_ave%latent(isc:iec,jsc:jec) = CS%Ocean_background_ave%latent(isc:iec,jsc:jec) + fluxes%latent(isc:iec,jsc:jec)
    CS%Ocean_background_ave%sensible(isc:iec,jsc:jec) = CS%Ocean_background_ave%sensible(isc:iec,jsc:jec) + fluxes%sens(isc:iec,jsc:jec)
    CS%Ocean_background_ave%lw(isc:iec,jsc:jec) = CS%Ocean_background_ave%lw(isc:iec,jsc:jec) + fluxes%lw(isc:iec,jsc:jec)
    CS%Ocean_background_ave%sw(isc:iec,jsc:jec) = CS%Ocean_background_ave%sw(isc:iec,jsc:jec) + fluxes%sw(isc:iec,jsc:jec)
  endif

  call pass_var(CS%Ocean_background_ave%T,G%Domain)
  call pass_var(CS%Ocean_background_ave%S,G%Domain)
  ! call pass_var(CS%Ocean_background_ave%SSH,G%Domain)

  if (CS%do_T_ml_bias_adjustment .or. CS%do_S_ml_bias_adjustment) then
    call pass_var(CS%Ocean_background_ave%U,G%Domain)
    call pass_var(CS%Ocean_background_ave%V,G%Domain)
    call pass_var(CS%Ocean_background_ave%taux,G%Domain)
    call pass_var(CS%Ocean_background_ave%tauy,G%Domain)
    call pass_var(CS%Ocean_background_ave%latent,G%Domain)
    call pass_var(CS%Ocean_background_ave%sensible,G%Domain)
    call pass_var(CS%Ocean_background_ave%lw,G%Domain)
    call pass_var(CS%Ocean_background_ave%sw,G%Domain)
  endif

  CS%prior_ave_counter = CS%prior_ave_counter + 1.0

  call get_date(CS%Prior_Time, yr, mon, day, hr, min, sec)
  write(mesg,*) 'Count: ', INT(CS%prior_ave_counter),' Prior Time: ', yr, mon, day, hr, min, sec
  call MOM_mesg("ODA get_prior: "//trim(mesg))
  
  if (Time >= CS%Prior_Time) then
    ! increment the analysis time to the next step
    CS%Prior_Time = CS%Prior_Time + real_to_time(CS%US%T_to_s*(CS%prior_interval))
  endif
  if (CS%Prior_Time < Time) then
    call MOM_error(FATAL, " set_prior_time: " // &
         "prior averaging interval appears to be shorter than " // &
         "the model timestep")
  endif

  call enable_averaging(CS%prior_interval, CS%Prior_Time, CS%diag_CS)

  if (CS%id_prior_count > 0) call post_data(CS%id_prior_count, CS%prior_ave_counter, CS%diag_CS)
  if (CS%id_prior_t > 0) call post_data(CS%id_prior_t, tv%T, CS%diag_CS)
  if (CS%id_prior_s > 0) call post_data(CS%id_prior_s, tv%S, CS%diag_CS)
  if (CS%id_prior_u > 0) call post_data(CS%id_prior_u, model_u, CS%diag_CS)
  if (CS%id_prior_v > 0) call post_data(CS%id_prior_v, model_v, CS%diag_CS)
  if (CS%id_prior_t_z > 0) call post_data(CS%id_prior_t_z, T, CS%diag_CS)
  if (CS%id_prior_s_z > 0) call post_data(CS%id_prior_s_z, S, CS%diag_CS)
  if (CS%id_prior_u_z > 0) call post_data(CS%id_prior_u_z, U, CS%diag_CS)
  if (CS%id_prior_v_z > 0) call post_data(CS%id_prior_v_z, V, CS%diag_CS)
  if (CS%id_prior_ssh > 0) call post_data(CS%id_prior_ssh, model_ssh, CS%diag_CS)
  if (CS%id_prior_taux > 0) call post_data(CS%id_prior_taux, forces%taux, CS%diag_CS)
  if (CS%id_prior_tauy > 0) call post_data(CS%id_prior_tauy, forces%tauy, CS%diag_CS)
  if (CS%id_prior_latent > 0) call post_data(CS%id_prior_latent, fluxes%latent, CS%diag_CS)
  if (CS%id_prior_sensible > 0) call post_data(CS%id_prior_sensible, fluxes%sens, CS%diag_CS)
  if (CS%id_prior_lw > 0) call post_data(CS%id_prior_lw, fluxes%lw, CS%diag_CS)
  if (CS%id_prior_sw > 0) call post_data(CS%id_prior_sw, fluxes%sw, CS%diag_CS)
  
  call disable_averaging(CS%diag_CS)

  call cpu_clock_end(id_clock_get_prior)
  
  return

end subroutine set_prior_tracer

!> Returns posterior adjustments or full state
!!Note that only those PEs associated with an ensemble member receive data
subroutine get_posterior_tracer(Time, CS, increment)
  type(time_type), intent(in) :: Time !< the current model time
  type(ODA_CS), pointer :: CS !< ocean DA control structure
  logical, optional, intent(in) :: increment !< True if returning increment only

  integer :: m
  logical :: get_inc
  type(time_type) :: Time_Next

  ! return if not analysis time (retain pointers for h and tv)
  if (Time < CS%Time .or. CS%assim_method == NO_ASSIM) return

  !! switch to global pelist
  call set_PElist(CS%filter_pelist)
  call MOM_mesg('Getting posterior')

  !! Calculate and redistribute increments to CS%tv right after assimilation
  !! Retain CS%tv to calculate increments for IAU updates CS%tv_inc otherwise
  get_inc = .true.
  if (present(increment)) get_inc = increment

  if (get_inc) then
    CS%Ocean_posterior%T = CS%Ocean_posterior%T - CS%Ocean_prior%T
    CS%Ocean_posterior%S = CS%Ocean_posterior%S - CS%Ocean_prior%S
  endif
  ! It may be necessary to check whether the increment and ocean state have the
  ! same dimensionally rescaled units.
  do m=1,CS%ensemble_size
    if (get_inc) then
      call redistribute_array(CS%mpp_domain, CS%Ocean_posterior%T(:,:,:,m),&
           CS%domains(m)%mpp_domain, CS%T_tend, complete=.true.)
      call redistribute_array(CS%mpp_domain, CS%Ocean_posterior%S(:,:,:,m),&
           CS%domains(m)%mpp_domain, CS%S_tend, complete=.true.)
    else
      call redistribute_array(CS%mpp_domain, CS%Ocean_posterior%T(:,:,:,m),&
           CS%domains(m)%mpp_domain, CS%T_tend, complete=.true.)
      call redistribute_array(CS%mpp_domain, CS%Ocean_posterior%S(:,:,:,m),&
           CS%domains(m)%mpp_domain, CS%S_tend, complete=.true.)
    endif
  enddo

  !! switch back to ensemble member pelist
  call set_PElist(CS%ensemble_pelist(CS%ensemble_id,:))

  call pass_var(CS%T_tend,CS%domains(CS%ensemble_id))
  call pass_var(CS%S_tend,CS%domains(CS%ensemble_id))

  !convert to a tendency (degC or PSU per second)
  CS%T_tend = CS%T_tend / (CS%assim_interval)
  CS%S_tend = CS%S_tend / (CS%assim_interval)

end subroutine get_posterior_tracer

!> Gather observations and call ODA routines
subroutine oda(Time, CS)
  type(time_type), intent(in) :: Time !< the current model time
  type(oda_CS), pointer :: CS !< A pointer the ocean DA control structure
  integer :: m
  character(len=160) :: mesg  ! The text of an error message
  integer :: yr, mon, day, hr, min, sec
  integer :: isc, iec, jsc, jec

  if ( Time >= CS%Time ) then
    call cpu_clock_begin(id_clock_ensemble_filter)

    isc=CS%model_G%isc; iec=CS%model_G%iec; jsc=CS%model_G%jsc; jec=CS%model_G%jec

    call get_date(Time, yr, mon, day, hr, min, sec)
    write(mesg,*) 'Count: ', INT(CS%prior_ave_counter),' Model Time: ', yr, mon, day, hr, min, sec
    call MOM_mesg("ODA averaging prior: "//trim(mesg))

    CS%Ocean_background_ave%T = CS%Ocean_background_ave%T / (CS%prior_ave_counter)
    CS%Ocean_background_ave%S = CS%Ocean_background_ave%S / (CS%prior_ave_counter)
    ! CS%Ocean_background_ave%SSH = CS%Ocean_background_ave%SSH / (CS%prior_ave_counter)
    if (CS%do_T_ml_bias_adjustment .or. CS%do_S_ml_bias_adjustment) then
      CS%Ocean_background_ave%U = CS%Ocean_background_ave%U / (CS%prior_ave_counter)
      CS%Ocean_background_ave%V = CS%Ocean_background_ave%V / (CS%prior_ave_counter)
      CS%Ocean_background_ave%taux = CS%Ocean_background_ave%taux / (CS%prior_ave_counter)
      CS%Ocean_background_ave%tauy = CS%Ocean_background_ave%tauy / (CS%prior_ave_counter)
      CS%Ocean_background_ave%latent = CS%Ocean_background_ave%latent / (CS%prior_ave_counter)
      CS%Ocean_background_ave%sensible = CS%Ocean_background_ave%sensible / (CS%prior_ave_counter)
      CS%Ocean_background_ave%lw = CS%Ocean_background_ave%lw / (CS%prior_ave_counter)
      CS%Ocean_background_ave%sw = CS%Ocean_background_ave%sw / (CS%prior_ave_counter)
    endif

    call pass_var(CS%Ocean_background_ave%T,CS%model_G%Domain)
    call pass_var(CS%Ocean_background_ave%S,CS%model_G%Domain)
    ! call pass_var(CS%Ocean_background_ave%SSH,CS%model_G%Domain)
    if (CS%do_T_ml_bias_adjustment .or. CS%do_S_ml_bias_adjustment) then
      call pass_var(CS%Ocean_background_ave%U,CS%model_G%Domain)
      call pass_var(CS%Ocean_background_ave%V,CS%model_G%Domain)
      call pass_var(CS%Ocean_background_ave%taux,CS%model_G%Domain)
      call pass_var(CS%Ocean_background_ave%tauy,CS%model_G%Domain)
      call pass_var(CS%Ocean_background_ave%latent,CS%model_G%Domain)
      call pass_var(CS%Ocean_background_ave%sensible,CS%model_G%Domain)
      call pass_var(CS%Ocean_background_ave%lw,CS%model_G%Domain)
      call pass_var(CS%Ocean_background_ave%sw,CS%model_G%Domain)
    endif
    
    CS%prior_ave_counter = 0.0

    if (CS%do_T_bias_adjustment .or. CS%do_S_bias_adjustment) call get_bias_correction_tracer(Time, CS%US, CS)

    if (CS%do_T_ml_bias_adjustment .or. CS%do_S_ml_bias_adjustment) then

      call get_ML_bias_correction(Time, CS%US, CS)

      if (CS%do_T_ml_bias_adjustment) then
        CS%Ocean_background_ave%T(isc:iec,jsc:jec,:) = CS%Ocean_background_ave%T(isc:iec,jsc:jec,:) + &
          CS%T_ml_tend(isc:iec,jsc:jec,:) * CS%assim_interval
        call pass_var(CS%Ocean_background_ave%T, CS%model_G%Domain)
      endif

      if (CS%do_S_ml_bias_adjustment) then
        CS%Ocean_background_ave%S(isc:iec,jsc:jec,:) = CS%Ocean_background_ave%S(isc:iec,jsc:jec,:) + &
          CS%S_ml_tend(isc:iec,jsc:jec,:) * CS%assim_interval
        call pass_var(CS%Ocean_background_ave%S, CS%model_G%Domain)
      endif

    endif

    !! switch to global pelist
    if (.NOT. CS%assim_method == NO_ASSIM) then
      call set_PElist(CS%filter_pelist)
      call get_profiles(Time, CS%Profiles, CS%CProfiles)

      do m=1,CS%ensemble_size
        call redistribute_array(CS%domains(m)%mpp_domain, CS%Ocean_background_ave%T,&
            CS%mpp_domain, CS%Ocean_prior%T(:,:,:,m), complete=.true.)
        call redistribute_array(CS%domains(m)%mpp_domain, CS%Ocean_background_ave%S,&
            CS%mpp_domain, CS%Ocean_prior%S(:,:,:,m), complete=.true.)
        ! call redistribute_array(CS%domains(m)%mpp_domain, CS%Ocean_background_ave%SSH,&
        !     CS%mpp_domain, CS%Ocean_prior%SSH(:,:,m), complete=.true.)
      enddo

      do m=1,CS%ensemble_size
        call pass_var(CS%Ocean_prior%T(:,:,:,m),CS%Grid%domain)
        call pass_var(CS%Ocean_prior%S(:,:,:,m),CS%Grid%domain)
        ! call pass_var(CS%Ocean_prior%SSH(:,:,m),CS%Grid%domain)
      enddo

      #ifdef ENABLE_ECDA
        call ensemble_filter(CS%Ocean_prior, CS%Ocean_posterior, CS%CProfiles, CS%kdroot, CS%mpp_domain, CS%oda_grid)
      #endif
      !! switch back to ensemble member pelist
      call set_PElist(CS%ensemble_pelist(CS%ensemble_id,:))
    endif

    call get_posterior_tracer(Time, CS, increment=.true.)

    call cpu_clock_end(id_clock_ensemble_filter)

  endif

  return
end subroutine oda

subroutine get_bias_correction_tracer(Time, US, CS)
  type(time_type), intent(in) :: Time !< the current model time
  type(unit_scale_type), intent(in) :: US !< A dimensional unit scaling type
  type(ODA_CS), pointer :: CS !< ocean DA control structure

  ! Local variables
  real, allocatable, dimension(:,:,:) :: T_bias ! Estimated temperature tendency bias [C T-1 ~> degC s-1]
  real, allocatable, dimension(:,:,:) :: S_bias ! Estimated salinity tendency bias [S T-1 ~> ppt s-1]
  real, allocatable, dimension(:,:,:) :: valid_flag ! Valid value flag on the horizontal model grid
                                                    ! and input-file vertical levels [nondim]
  real, allocatable, dimension(:), target :: z_in       ! Cell center depths for input data [Z ~> m]
  real, allocatable, dimension(:), target :: z_edges_in ! Cell edge depths for input data [Z ~> m]
  real :: missing_value ! A value indicating that there is no valid input data at this point [CU ~> conc]
  integer, dimension(3) :: fld_sz
  integer :: i,j,k


  call cpu_clock_begin(id_clock_bias_adjustment)
  call horiz_interp_and_extrap_tracer(CS%INC_CS%T, Time, CS%model_G, T_bias, &
            valid_flag, z_in, z_edges_in, missing_value, scale=US%degC_to_C*US%s_to_T, spongeOngrid=.true.)
  call horiz_interp_and_extrap_tracer(CS%INC_CS%S, Time, CS%model_G, S_bias, &
            valid_flag, z_in, z_edges_in, missing_value, scale=US%ppt_to_S*US%s_to_T, spongeOngrid=.true.)

  ! This should be replaced to use mask_z instead of the following lines
  ! which are intended to zero land values using an arbitrary limit.
  fld_sz=shape(T_bias)
  do i=1,fld_sz(1)
    do j=1,fld_sz(2)
      do k=1,fld_sz(3)
! The following two lines are needed for backward compatibility (ANSWER_DATE< 20181231) 
! Need to discuss if we need to use the valid_flag instead.
        if (ABS(T_bias(i,j,k)) > 1.0E-3*US%degC_to_C) T_bias(i,j,k) = 0.0
        if (ABS(S_bias(i,j,k)) > 1.0E-3*US%ppt_to_S) S_bias(i,j,k) = 0.0
!        if (valid_flag(i,j,k)==0.) then
!          T_bias(i,j,k)=0.0
!          S_bias(i,j,k)=0.0
!        endif
      enddo
    enddo
  enddo

  CS%T_bc_tend = T_bias * CS%T_bias_adjustment_multiplier
  CS%S_bc_tend = S_bias * CS%S_bias_adjustment_multiplier

  call pass_var(CS%T_bc_tend, CS%domains(CS%ensemble_id))
  call pass_var(CS%S_bc_tend, CS%domains(CS%ensemble_id))

  call cpu_clock_end(id_clock_bias_adjustment)

end subroutine get_bias_correction_tracer

subroutine get_ML_bias_correction(Time, US, CS)
  type(time_type), intent(in) :: Time !< the current model time
  type(unit_scale_type), intent(in) :: US !< A dimensional unit scaling type
  type(ODA_CS), pointer :: CS !< ocean DA control structure

  ! Local variables
  integer :: isc, iec, jsc, jec
  integer :: i,j,k
  character(len=160) :: mesg  ! The text of an error message

  call cpu_clock_begin(id_clock_ml_bias_correction)

  call MOM_mesg('Doing ML inference')

  isc=CS%model_G%isc; iec=CS%model_G%iec; jsc=CS%model_G%jsc; jec=CS%model_G%jec
  !! Loop through all local gridpoints
  do j=jsc,jec ; do i=isc,iec

    if (CS%model_G%geolatT(i,j) > 60.0 .or. CS%model_G%geolatT(i,j) < -60.0) then
      CS%T_ml_tend(i,j,:) = 0.0
      CS%S_ml_tend(i,j,:) = 0.0
    else
      !! put local variables into ml_data
      CS%ml_data%T = CS%Ocean_background_ave%T(i,j,:)
      CS%ml_data%S = CS%Ocean_background_ave%S(i,j,:)
      CS%ml_data%U_left = CS%Ocean_background_ave%U(i-1,j,:)
      CS%ml_data%U_right = CS%Ocean_background_ave%U(i,j,:)
      CS%ml_data%V_north = CS%Ocean_background_ave%V(i,j,:)
      CS%ml_data%V_south = CS%Ocean_background_ave%V(i,j-1,:)
      CS%ml_data%latent = CS%Ocean_background_ave%latent(i,j)
      CS%ml_data%sensible = CS%Ocean_background_ave%sensible(i,j)
      CS%ml_data%lw = CS%Ocean_background_ave%lw(i,j)
      CS%ml_data%sw = CS%Ocean_background_ave%sw(i,j)
      CS%ml_data%taux_left = CS%Ocean_background_ave%taux(i-1,j)
      CS%ml_data%taux_right = CS%Ocean_background_ave%taux(i,j)
      CS%ml_data%tauy_north = CS%Ocean_background_ave%tauy(i,j)
      CS%ml_data%tauy_south = CS%Ocean_background_ave%tauy(i,j-1)
      
      CS%ml_data%dyCu_left = CS%model_G%dyCu(i-1,j)
      CS%ml_data%dyCu_right = CS%model_G%dyCu(i,j)
      CS%ml_data%dxCv_north = CS%model_G%dxCv(i,j)
      CS%ml_data%dxCv_south = CS%model_G%dxCv(i,j-1)
      CS%ml_data%areacello = CS%model_G%areaT(i,j)
      
      CS%ml_data%bathyT = CS%model_G%bathyT(i,j)
      CS%ml_data%bathyU_left = (CS%model_G%bathyT(i-1,j)+CS%model_G%bathyT(i,j))/2
      CS%ml_data%bathyU_right = (CS%model_G%bathyT(i,j)+CS%model_G%bathyT(i+1,j))/2
      CS%ml_data%bathyV_south = (CS%model_G%bathyT(i,j-1)+CS%model_G%bathyT(i,j))/2
      CS%ml_data%bathyV_north = (CS%model_G%bathyT(i,j)+CS%model_G%bathyT(i,j+1))/2

      CS%ml_data%mask2dT = CS%model_G%mask2dT(i,j)
      CS%ml_data%OBCmaskCu_left = CS%model_G%OBCmaskCu(i-1,j)
      CS%ml_data%OBCmaskCu_right = CS%model_G%OBCmaskCu(i,j)
      CS%ml_data%OBCmaskCv_south = CS%model_G%OBCmaskCv(i,j-1)
      CS%ml_data%OBCmaskCv_north = CS%model_G%OBCmaskCu(i,j)

      !! Call inference subroutine with the concatenated vector
      call oda_ml_inference(CS%ml_config, CS%ml_data)

      CS%T_ml_tend(i,j,:) = CS%ml_data%T_inc * CS%T_ml_bias_adjustment_multiplier
      CS%S_ml_tend(i,j,:) = CS%ml_data%S_inc * CS%S_ml_bias_adjustment_multiplier

      do k=1,CS%nk
        if (CS%T_ml_tend(i,j,k) > 1.0E-5*US%degC_to_C) CS%T_ml_tend(i,j,k) = 1.0E-5
        if (CS%T_ml_tend(i,j,k) < -1.0E-5*US%degC_to_C) CS%T_ml_tend(i,j,k) = -1.0E-5
      enddo

    endif
  enddo; enddo

  call pass_var(CS%T_ml_tend, CS%domains(CS%ensemble_id))
  call pass_var(CS%S_ml_tend, CS%domains(CS%ensemble_id))

  call cpu_clock_end(id_clock_ml_bias_correction)

end subroutine get_ML_bias_correction

!> Finalize DA module
subroutine oda_end(CS)
  type(ODA_CS), intent(inout) :: CS !< the ocean DA control structure

end subroutine oda_end

!> Initialize DA module
subroutine init_ocean_ensemble(CS,Grid,GV,ens_size)
  type(ensemble_control_struct), pointer :: CS !< Pointer to ODA control structure
  type(ocean_grid_type), pointer :: Grid !< Pointer to ocean analysis grid
  type(verticalGrid_type), pointer :: GV !< Pointer to DA vertical grid
  integer, intent(in) :: ens_size !< ensemble size

  integer :: isd, ied, jsd, jed, nk, isdB, iedB, jsdB, jedB
  character(len=160) :: mesg  ! The text of an error message
  
  nk=GV%ke
  isd  = Grid%isd ; ied  = Grid%ied  ; jsd  = Grid%jsd  ; jed  = Grid%jed
  isdB = Grid%isdB; iedB = Grid%iedB ; jsdB = Grid%jsdB ; jedB = Grid%jedB
  CS%ensemble_size=ens_size

  allocate(CS%T(isd:ied,jsd:jed,nk,ens_size),source=0.0)
  allocate(CS%S(isd:ied,jsd:jed,nk,ens_size),source=0.0)
  ! allocate(CS%SSH(isd:ied,jsd:jed,ens_size),source=0.0)

  return
end subroutine init_ocean_ensemble

!> Initialize background variables
subroutine init_ocean_background(CS,PriorCS,Grid,GV)
  type(ODA_CS), pointer, intent(inout) :: CS !< the DA control structure
  type(ocean_control_struct), pointer :: PriorCS !< Pointer to ODA control structure
  type(ocean_grid_type), pointer :: Grid !< Pointer to ocean analysis grid
  type(verticalGrid_type), pointer :: GV !< Pointer to DA vertical grid

  integer :: isd, ied, jsd, jed, nk, isdB, iedB, jsdB, jedB
  character(len=160) :: mesg  ! The text of an error message

  nk=GV%ke
  isd  = Grid%isd ; ied  = Grid%ied  ; jsd  = Grid%jsd  ; jed  = Grid%jed
  isdB = Grid%isdB; iedB = Grid%iedB ; jsdB = Grid%jsdB ; jedB = Grid%jedB

  allocate(PriorCS%T(isd:ied,jsd:jed,nk),source=0.0)
  allocate(PriorCS%S(isd:ied,jsd:jed,nk),source=0.0)
  ! allocate(PriorCS%SSH(isd:ied,jsd:jed),source=0.0)
  
  if (CS%do_T_ml_bias_adjustment .or. CS%do_S_ml_bias_adjustment) then
    allocate(PriorCS%U(isdB:iedB,jsd:jed,nk),source=0.0)
    allocate(PriorCS%V(isd:ied,jsdB:jedB,nk),source=0.0)
    allocate(PriorCS%taux(isdB:iedB,jsd:jed),source=0.0)
    allocate(PriorCS%tauy(isd:ied,jsdB:jedB),source=0.0)
    allocate(PriorCS%latent(isd:ied,jsd:jed),source=0.0)
    allocate(PriorCS%sensible(isd:ied,jsd:jed),source=0.0)
    allocate(PriorCS%lw(isd:ied,jsd:jed),source=0.0)
    allocate(PriorCS%sw(isd:ied,jsd:jed),source=0.0)
  endif

  return
end subroutine init_ocean_background

!> Set the next analysis time
subroutine set_analysis_time(Time,CS)
  type(time_type), intent(in) :: Time !< the current model time
  type(ODA_CS), pointer, intent(inout) :: CS !< the DA control structure

  character(len=160) :: mesg  ! The text of an error message
  integer :: yr, mon, day, hr, min, sec

  if (Time >= CS%Time) then
    ! increment the analysis time to the next step
    CS%Time = CS%Time + real_to_time(CS%US%T_to_s*(CS%assim_interval))

    call get_date(Time, yr, mon, day, hr, min, sec)
    write(mesg,*) 'Model Time: ', yr, mon, day, hr, min, sec
    call MOM_mesg("set_analysis_time: "//trim(mesg))
    call get_date(CS%time, yr, mon, day, hr, min, sec)
    write(mesg,*) 'Assimilation Time: ', yr, mon, day, hr, min, sec
    call MOM_mesg("set_analysis_time: "//trim(mesg))
  endif
  if (CS%Time < Time) then
    call MOM_error(FATAL, " set_analysis_time: " // &
         "assimilation interval appears to be shorter than " // &
         "the model timestep")
  endif
  return

end subroutine set_analysis_time

!> Apply increments to tracers
subroutine apply_oda_tracer_increments(Time, G, GV, tv, h, CS)
  ! real,                     intent(in)    :: dt !< The tracer timestep [T ~> s]
  ! type(time_type), intent(in)             :: Time_end !< Time at the end of the interval
  type(time_type), intent(in)             :: Time !< Time at the end of the interval
  type(ocean_grid_type),    intent(in)    :: G  !< ocean grid structure
  type(verticalGrid_type),  intent(in)    :: GV !< The ocean's vertical grid structure
  type(thermo_var_ptrs),    intent(inout) :: tv !< A structure pointing to various thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: h  !< layer thickness [H ~> m or kg m-2]
  type(ODA_CS), pointer                   :: CS !< the data assimilation structure

  !! local variables
  integer :: i, j
  integer :: isc, iec, jsc, jec
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: T_tend_inc !< an adjustment to the temperature
                                                    !! tendency [C T-1 -> degC s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: S_tend_inc !< an adjustment to the salinity
                                                    !! tendency [S T-1 -> ppt s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: T_ml_tend_inc !< an adjustment to the temperature
                                                    !! tendency [C T-1 -> degC s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: S_ml_tend_inc !< an adjustment to the salinity
                                                    !! tendency [S T-1 -> ppt s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(CS%Grid)) :: T_tend !< The temperature tendency adjustment from
                                                           !! DA [C T-1 ~> degC s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(CS%Grid)) :: S_tend !< The salinity tendency adjustment from DA
                                                          !! [S T-1 ~> ppt s-1]
  real :: h_neglect, h_neglect_edge                 ! small thicknesses [H ~> m or kg m-2]
  character(len=160) :: mesg  ! The text of an error message
  integer :: yr, mon, day, hr, min, sec

  if (.not. associated(CS)) return
  if (Time < CS%Apply_Time) return
  if (CS%assim_method == NO_ASSIM .and. (.not. CS%do_T_bias_adjustment) .and. (.not. CS%do_S_bias_adjustment) &
    .and. (.not. CS%do_T_ml_bias_adjustment) .and. (.not. CS%do_S_ml_bias_adjustment)) return

  call cpu_clock_begin(id_clock_apply_increments)

  call get_date(CS%Apply_Time, yr, mon, day, hr, min, sec)
  write(mesg,*) 'apply_int: ', INT(CS%apply_interval),' Apply Time: ', yr, mon, day, hr, min, sec
  call MOM_mesg("ODA applying increments: "//trim(mesg))

  T_tend_inc(:,:,:) = 0.0; S_tend_inc(:,:,:) = 0.0; T_tend(:,:,:) = 0.0; S_tend(:,:,:) = 0.0
  if (.NOT. CS%assim_method == NO_ASSIM) then
    T_tend = T_tend + CS%T_tend
    S_tend = S_tend + CS%S_tend
  endif
  if (CS%do_T_bias_adjustment ) then
    T_tend = T_tend + CS%T_bc_tend
  endif
  if (CS%do_S_bias_adjustment ) then
    S_tend = S_tend + CS%S_bc_tend
  endif
  if (CS%do_T_ml_bias_adjustment ) then
    T_tend = T_tend + CS%T_ml_tend
  endif
  if (CS%do_S_ml_bias_adjustment ) then
    S_tend = S_tend + CS%S_ml_tend
  endif

  if (CS%answer_date >= 20190101) then
    h_neglect = GV%H_subroundoff ; h_neglect_edge = GV%H_subroundoff
  elseif (GV%Boussinesq) then
    h_neglect = GV%m_to_H * 1.0e-30 ; h_neglect_edge = GV%m_to_H * 1.0e-10
  else
    h_neglect = GV%kg_m2_to_H * 1.0e-30 ; h_neglect_edge = GV%kg_m2_to_H * 1.0e-10
  endif

  isc=G%isc; iec=G%iec; jsc=G%jsc; jec=G%jec
  do j=jsc,jec; do i=isc,iec
    call remapping_core_h(CS%remapCS, CS%nk, CS%h(i,j,:), T_tend(i,j,:), &
         G%ke, h(i,j,:), T_tend_inc(i,j,:), h_neglect, h_neglect_edge)
    call remapping_core_h(CS%remapCS, CS%nk, CS%h(i,j,:), S_tend(i,j,:), &
         G%ke, h(i,j,:), S_tend_inc(i,j,:), h_neglect, h_neglect_edge)
  enddo; enddo

  call pass_var(T_tend_inc, G%Domain)
  call pass_var(S_tend_inc, G%Domain)

  if (CS%do_T_ml_bias_adjustment .or. CS%do_S_ml_bias_adjustment) then
    do j=jsc,jec; do i=isc,iec
      call remapping_core_h(CS%remapCS, CS%nk, CS%h(i,j,:), CS%T_ml_tend(i,j,:), &
          G%ke, h(i,j,:), T_ml_tend_inc(i,j,:), h_neglect, h_neglect_edge)
      call remapping_core_h(CS%remapCS, CS%nk, CS%h(i,j,:), CS%S_ml_tend(i,j,:), &
          G%ke, h(i,j,:), S_ml_tend_inc(i,j,:), h_neglect, h_neglect_edge)
    enddo; enddo
    call pass_var(T_ml_tend_inc, G%Domain)
    call pass_var(S_ml_tend_inc, G%Domain)
  endif

  tv%T(isc:iec,jsc:jec,:) = tv%T(isc:iec,jsc:jec,:) + T_tend_inc(isc:iec,jsc:jec,:)*CS%apply_interval
  tv%S(isc:iec,jsc:jec,:) = tv%S(isc:iec,jsc:jec,:) + S_tend_inc(isc:iec,jsc:jec,:)*CS%apply_interval

  call pass_var(tv%T, G%Domain)
  call pass_var(tv%S, G%Domain)

  if (Time >= CS%Apply_Time) then
    ! increment the analysis time to the next step
    CS%Apply_Time = CS%Apply_Time + real_to_time(CS%US%T_to_s*(CS%apply_interval))
  endif
  if (CS%Apply_Time < Time) then
    call MOM_error(FATAL, " set_apply_time: " // &
         "increment application interval appears to be shorter than " // &
         "the model timestep")
  endif

  call enable_averaging(CS%apply_interval, CS%Apply_Time, CS%diag_CS)
  if (CS%id_inc_t > 0) call post_data(CS%id_inc_t, T_tend_inc, CS%diag_CS)
  if (CS%id_inc_s > 0) call post_data(CS%id_inc_s, S_tend_inc, CS%diag_CS)
  if (CS%id_inc_t_z > 0) call post_data(CS%id_inc_t_z, T_tend, CS%diag_CS)
  if (CS%id_inc_s_z > 0) call post_data(CS%id_inc_s_z, S_tend, CS%diag_CS)
  if (CS%do_T_ml_bias_adjustment) then
    if (CS%id_inc_ml_t > 0) call post_data(CS%id_inc_ml_t, T_ml_tend_inc, CS%diag_CS)
    if (CS%id_inc_ml_t_z > 0) call post_data(CS%id_inc_ml_t_z, CS%T_ml_tend, CS%diag_CS)
  endif
  if (CS%do_S_ml_bias_adjustment) then
    if (CS%id_inc_ml_s > 0) call post_data(CS%id_inc_ml_s, S_ml_tend_inc, CS%diag_CS)
    if (CS%id_inc_ml_s_z > 0) call post_data(CS%id_inc_ml_s_z, CS%S_ml_tend, CS%diag_CS)
  endif
  call disable_averaging(CS%diag_CS)

  call cpu_clock_end(id_clock_apply_increments)

end subroutine apply_oda_tracer_increments

  subroutine set_up_global_tgrid(T_grid, CS, G)
    type(grid_type), pointer :: T_grid !< global tracer grid
    type(ODA_CS), pointer, intent(in) :: CS !< A pointer to DA control structure.
    type(ocean_grid_type), pointer :: G !< domain and grid information for ocean model

    ! local variables
    real, dimension(:,:), allocatable :: global2D, global2D_old
    integer :: i, j, k, ii, jj
    integer :: isd_oda,ied_oda,jsd_oda,jed_oda
    integer :: isc,iec,jsc,jec,isg,ieg,jsg,jeg
    character(len=160) :: mesg  ! The text of an error message
    !    get global grid information from ocean_model

    !if (associated(T_grid)) call MOM_error(FATAL,'MOM_oda_driver:set_up_global_tgrid called with associated T_grid')
    allocate(T_grid)

    T_grid%ni = CS%ni
    T_grid%nj = CS%nj
    T_grid%nk = CS%nk
    allocate(T_grid%x(CS%ni,CS%nj))
    allocate(T_grid%y(CS%ni,CS%nj))
    call global_field(CS%mpp_domain, CS%Grid%geolonT, T_grid%x)
    call global_field(CS%mpp_domain, CS%Grid%geolatT, T_grid%y)
    if (CS%use_basin_mask) then
      allocate(T_grid%basin_mask(CS%ni,CS%nj))
      call global_field(CS%mpp_domain, CS%oda_grid%basin_mask, T_grid%basin_mask)
    endif
    
    T_grid%bathyT => CS%Grid%bathyT
    
    if (.not. associated(T_grid%h)) then
      allocate(T_grid%h(CS%Grid%isd:CS%Grid%ied,CS%Grid%jsd:CS%Grid%jed,CS%nk), source=CS%GV%Angstrom_H)
      ! assign thicknesses
      call ALE_initThicknessToCoord(CS%ALE_CS, CS%Grid, CS%GV, T_grid%h)
      call pass_var(T_grid%h,CS%Grid%domain)
    endif

    allocate(T_grid%mask(CS%Grid%isd:CS%Grid%ied,CS%Grid%jsd:CS%Grid%jed,CS%nk), source=0.0)
    allocate(T_grid%z(CS%Grid%isd:CS%Grid%ied,CS%Grid%jsd:CS%Grid%jed,CS%nk), source=0.0)

    do k = 1, CS%nk
      do i=CS%Grid%isd,CS%Grid%ied ; do j=CS%Grid%jsd,CS%Grid%jed
        if ( T_grid%h(i,j,k) > 1 ) then
           T_grid%mask(i,j,k) = 1.0
        endif
      enddo; enddo
      if (k == 1) then
         T_grid%z(:,:,k) = &
            T_grid%h(:,:,k)/2
      else
         T_grid%z(:,:,k) = T_grid%z(:,:,k-1) + &
            (T_grid%h(:,:,k) + &
            T_grid%h(:,:,k-1))/2
      endif
    enddo

  end subroutine set_up_global_tgrid

  subroutine set_oda_restart_fields(US, CS, restart_CSp)
    ! type(verticalGrid_type),  intent(inout) :: GV         !< ocean vertical grid structure
    type(unit_scale_type),    intent(inout) :: US         !< A dimensional unit scaling type
    ! type(param_file_type),    intent(in) :: param_file    !< opened file for parsing to get parameters
    type(ODA_CS),             intent(in) :: CS            !< control structure set up by initialize_MOM
    type(MOM_restart_CS),     pointer    :: restart_CSp   !< pointer to the restart control
                                                          !! structure that will be used for MOM.
    ! Local variables

    if (associated(CS%Ocean_background_ave%T)) &
      call register_restart_field(CS%Ocean_background_ave%T, "Temp_prior", .false., restart_CSp, &
                                  "Accumulated Potential Temperature", "degC", conversion=US%C_to_degC)

    if (associated(CS%Ocean_background_ave%S)) &
      call register_restart_field(CS%Ocean_background_ave%S, "Salt_prior", .false., restart_CSp, &
                                  "Accumulated Salinity", "PPT", conversion=US%S_to_ppt)

    if (associated(CS%Ocean_background_ave%SSH)) &
      call register_restart_field(CS%Ocean_background_ave%SSH, "SSH_prior", .false., restart_CSp, &
                                  "Accumulated SSH", 'm', conversion=US%Z_to_m)
                                  
    if (CS%do_T_ml_bias_adjustment .or. CS%do_S_ml_bias_adjustment) then

      if (associated(CS%Ocean_background_ave%U)) &
        call register_restart_field(CS%Ocean_background_ave%U, "u_prior", .false., restart_CSp, &
                                    "Accumulated Zonal Velocity", "m s-1", conversion=US%L_T_to_m_s)

      if (associated(CS%Ocean_background_ave%V)) &
        call register_restart_field(CS%Ocean_background_ave%V, "v_prior", .false., restart_CSp, &
                                    "Accumulated Meridional Velocity", "m s-1", conversion=US%L_T_to_m_s)

      if (associated(CS%Ocean_background_ave%taux)) &
        call register_restart_field(CS%Ocean_background_ave%taux, "taux_prior", .false., restart_CSp, &
                                    "Accumulated zonal surface stress", 'Pa', conversion=US%RLZ_T2_to_Pa)

      if (associated(CS%Ocean_background_ave%tauy)) &
        call register_restart_field(CS%Ocean_background_ave%tauy, "tauy_prior", .false., restart_CSp, &
                                    "Accumulated meridional surface stress", 'Pa', conversion=US%RLZ_T2_to_Pa)

      if (associated(CS%Ocean_background_ave%latent)) &
        call register_restart_field(CS%Ocean_background_ave%latent, "latent_prior", .false., restart_CSp, &
                                    "Accumulated latent heat flux into ocean", 'W m-2', conversion=US%QRZ_T_to_W_m2)

      if (associated(CS%Ocean_background_ave%sensible)) &
        call register_restart_field(CS%Ocean_background_ave%sensible, "sensible_prior", .false., restart_CSp, &
                                    "Accumulated sensible heat flux into ocean", 'W m-2', conversion=US%QRZ_T_to_W_m2)
                                    
      if (associated(CS%Ocean_background_ave%lw)) &
        call register_restart_field(CS%Ocean_background_ave%lw, "LW_prior", .false., restart_CSp, &
                                    "Accumulated longwave radiation flux into ocean", 'W m-2', conversion=US%QRZ_T_to_W_m2)

      if (associated(CS%Ocean_background_ave%sw)) &
        call register_restart_field(CS%Ocean_background_ave%sw, "SW_prior", .false., restart_CSp, &
                                    "Accumulated shortwave radiation flux into ocean", 'W m-2', conversion=US%QRZ_T_to_W_m2)
    endif

  end subroutine set_oda_restart_fields

!> \namespace MOM_oda_driver_mod
!!
!! \section section_ODA The Ocean data assimilation (DA) and Ensemble Framework
!!
!! The DA framework implements ensemble capability in MOM6.   Currently, this framework
!! is enabled using the cpp directive ENSEMBLE_OCEAN.  The ensembles need to be generated
!! at the level of the calling routine for oda_init or above. The ensemble instances may
!! exist on overlapping or non-overlapping processors. The ensemble information is accessed
!! via the FMS ensemble manager. An independent PE layout is used to gather (prior) ensemble
!! member information where this information is stored in the ODA control structure.  This
!! module was developed in collaboration with Feiyu Lu and Tony Rosati in the GFDL prediction
!! group for use in their coupled ensemble framework. These interfaces should be suitable for
!! interfacing MOM6 to other data assimilation packages as well.

end module MOM_oda_driver_mod
