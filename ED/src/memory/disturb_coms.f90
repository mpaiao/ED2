!==========================================================================================!
!==========================================================================================!
!   module disturb_coms                                                                    !
!   This module contains variables used to control the disturbance rates.                  !
!   N.B.: Variables that are not parameters should be initialized in the subroutine        !
!         initialize_disturb_params in ed_params.f90, since some compilers don't actually  !
!         initialize variables in modules.                                                 !
!------------------------------------------------------------------------------------------!
module disturb_coms
   use ed_max_dims, only : str_len      & ! intent(in)
                         , maxgrds      & ! intent(in)
                         , n_pft        & ! intent(in)
                         , n_dist_types ! ! intent(in)
   implicit none


   !=======================================================================================!
   !=======================================================================================!
   !    General parameters.                                                                !
   !---------------------------------------------------------------------------------------!

   !------ Number of land use transitions in the input land use disturbance dataset. ------!
   integer, parameter :: num_lu_trans = 19

   !---------------------------------------------------------------------------------------!
   !     Maximum number of years in which land use can be applied.  In case the simulation !
   ! runs longer than this, the missing years will be filled with zeroes.  The first and   !
   ! last year of each is checked in landuse_init.  This variable is also used by the fire !
   ! ignition files.                                                                       !
   !---------------------------------------------------------------------------------------!
   integer, parameter :: max_lu_years = 2500 

   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     Namelist variables.                                                               !
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Fire model.  Possible values are:                                                 !
   !                   0. (ED-2.2 default) No fires.                                       !
   !                   1. (deprecated) Fire will be triggered with enough biomass and      !
   !                      integrated ground water depth less than a threshold.  Based on   !
   !                      ED-1, the threshold assumes that the soil depth is 1 m, so       !
   !                      deeper soils must be much drier to allow fires to happen. Kept   !
   !                      as an option for legacy, but strongly discouraged.               !
   !                   2. (ED-2.2 default) Fire will be triggered with enough biomass and  !
   !                      the total soil water at the top 50cm falls below a (relative)    !
   !                      threshold.                                                       !
   !                   3. (Beta) The EMBERFire model (Longo et al. in prep).  Empirical    !
   !                      approach loosely based on option 2, with the following           !
   !                      differences:                                                     !
   !                      - Fuel.  Fast/structural carbon + understory (h<2m) biomass.     !
   !                      - Fire ignition probably based on Nesterov index, similar to     !
   !                        to SPITFIRE (Thonicke et al. 2010)                             !
   !                      - Mortality.  Simple sigmoidal function of bark thickness.       !
   !                        (no account for fire intensity).                               !
   !                      - Combustion. A fraction of fuels is volatilised during the      !
   !                        disturbance, and the remainder is accumulated as litter.       !
   !                   4. (In development) FIRESTARTER model.  This is a more process-     !
   !                      -based fire model that builds on fire density, ignition and      !
   !                      termination from HESFIRE (LePage et al. 2015), and fire damage   !
   !                      and spread rates based on SPITFIRE (Thonicke et al. 2010).       !
   !                      By default this option is not going to run, as it is still in    !
   !                      development.                                                     !
   !---------------------------------------------------------------------------------------!
   integer :: include_fire
   !---------------------------------------------------------------------------------------!



   !----- Dimensionless parameter controlling speed of fire spread. -----------------------!
   real :: fire_parameter
   
   !----- Fractions of fast and structural carbon and nitrogen lost through combustion. ---!
   real :: fe_combusted_fast_c
   real :: fe_combusted_struct_c
   real :: fe_combusted_fast_n
   real :: fe_combusted_struct_n
   !---- Maximum height for non-grass cohort to be considered part of fuel. ---------------!
   real :: fuel_height_max
   !---- Flag: only allow burning in polygons with some anthropogenic activity? (T|F) -----!
   logical :: fe_anth_ignt_only

   !---------------------------------------------------------------------------------------!
   !     Anthropogenic disturbance.  1 means that anthropogenic disturbances will be       !
   ! included, whereas 0 means that it won't.                                              !
   !---------------------------------------------------------------------------------------!
   integer :: ianth_disturb
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Treefall disturbance                                                              !
   ! > 0. Usual disturbance rate, in 1/years, at which treefall gaps form.                 !
   ! = 0. No treefall disturbance;                                                         !
   ! < 0. Treefall will be added as a mortality rate (it will kill plants, but it won't    !
   !      create a new patch).                                                             !
   !---------------------------------------------------------------------------------------!
   real :: treefall_disturbance_rate
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     This is the time until we start knocking down trees.  Used only when              !
   ! treefall_disturbance_rate > 0.                                                        !
   !---------------------------------------------------------------------------------------!
   real :: time2canopy
   !---------------------------------------------------------------------------------------!

   !----- Minimum relative area required for a patch to be created or maintained. ---------!
   real :: min_patch_area 
   !---------------------------------------------------------------------------------------!


   !----- The prefix for land use disturbance rates. The path and prefix must be included. !
   character(len=str_len), dimension(maxgrds) :: lu_database 
   !----- File with plantation fraction.  If no file is available, leave it blank. --------!
   character(len=str_len), dimension(maxgrds) :: plantation_file
   !----- File with initial land use area scale.  If no file is available, leave it blank. !
   character(len=str_len), dimension(maxgrds) :: lu_rescale_file
   !----- The prefix for socio-economic indices. The path and prefix must be included. ----!
   character(len=str_len), dimension(maxgrds) :: sei_database 
   !----- The prefix for flash rate densities. The path and prefix must be included. ------!
   character(len=str_len), dimension(maxgrds) :: flash_database 
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Selective logging characteristics (when IANTH_DISTURB is set to 2).  Included in !
   ! ED2IN because selective logging is site-dependent.  Variable description.             !
   !                                                                                       !
   ! SL_SCALE           -- This flag assumes whether the simulation scale is local or      !
   !                       landscape.  This controls the recurrence of logging.            !
   !                       0.  Local. The simulation represents one logging unit.  Apply   !
   !                           logging only once every SL_NYRS                             !
   !                       1.  Landscape.  The simulation represents a landscape. Logging  !
   !                           occurs every year but it is restricted to patches with age  !
   !                           greater than or equal to SL_NYRS                            !
   ! SL_YR_FIRST        -- The first year to apply logging.  In case IANTH_DISTURB is 2 it !
   !                       must be a simulation year (i.e. between IYEARA and IYEARZ).     !
   ! SL_NYRS            -- This variable defines the logging cycle, in years (see variable !
   !                       SL_SCALE above)                                                 !
   ! SL_PFT             -- PFTs that can be harvested.                                     !
   ! SL_PROB_HARVEST    -- Logging intensity (one value for each PFT provided in SL_PFT).  !
   !                       Values should be between 0.0 and 1.0, with 0 meaning no         !
   !                       removal, and 1 removal of all trees needed to meet demands.     !
   ! SL_MINDBH_HARVEST  -- Minimum DBH for logging (one value for each PFT provided in     !
   !                       SL_PFT).                                                        !
   ! SL_BIOMASS_HARVEST -- Target biomass to be harvested in each cycle, in kgC/m2.  If    !
   !                       zero, then all trees that meet the minimum DBH and minimum      !
   !                       patch age will be logged.  In case you don't want logging to    !
   !                       occur, don't set this value to zero! Instead, set IANTH_DISTURB !
   !                       to zero.                                                        !
   !                                                                                       !
   ! The following variables are used when IANTH_DISTURB is 1 or 2.                        !
   !                                                                                       !
   ! SL_SKID_REL_AREA    -- area damaged by skid trails (relative to felled area).         !
   ! SL_SKID_S_GTHARV    -- survivorship of trees with DBH > MINDBH in skid trails.        !
   ! SL_SKID_S_LTHARV    -- survivorship of trees with DBH < MINDBH in skid trails.        !
   ! SL_FELLING_S_LTHARV -- survivorship of trees with DBH < MINDBH in felling gaps.       !
   !                                                                                       !
   ! Cropland variables, used when IANTH_DISTURB is 1 or 2.                                !
   !                                                                                       !
   ! CL_FSEEDS_HARVEST   -- fraction of seeds that is harvested.                           !
   ! CL_FSTORAGE_HARVEST -- fraction of non-structural carbon that is harvested.           !
   ! CL_FLEAF_HARVEST    -- fraction of leaves that is harvested in croplands.             !
   !---------------------------------------------------------------------------------------!
   integer                        :: sl_scale
   integer                        :: sl_yr_first
   integer                        :: sl_nyrs
   integer     , dimension(n_pft) :: sl_pft
   real(kind=4), dimension(n_pft) :: sl_prob_harvest
   real(kind=4), dimension(n_pft) :: sl_mindbh_harvest
   real(kind=4)                   :: sl_biomass_harvest
   real(kind=4)                   :: sl_skid_rel_area
   real(kind=4)                   :: sl_skid_s_gtharv
   real(kind=4)                   :: sl_skid_s_ltharv
   real(kind=4)                   :: sl_felling_s_ltharv
   real(kind=4)                   :: cl_fseeds_harvest
   real(kind=4)                   :: cl_fstorage_harvest
   real(kind=4)                   :: cl_fleaf_harvest
   !---------------------------------------------------------------------------------------!



   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    Patch dynamics variables, to be set in ed_params.f90.                              !
   !---------------------------------------------------------------------------------------!
   !----- Only trees above this height create a gap when they fall. -----------------------!
   real                          :: treefall_hite_threshold
   !----- Flag to decide whether or not to limit disturbance to patches with tall trees. --!
   logical                       :: does_hite_limit_tfpatch
   !---------------------------------------------------------------------------------------!
   !      Minimum age above which we disregard the disturbance type (land use) and assume  !
   ! old growth, thus allowing patch fusion to occur.                                      !
   !---------------------------------------------------------------------------------------!
   real, dimension(n_dist_types) :: min_oldgrowth
   !---------------------------------------------------------------------------------------!
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    Forestry variables, to be set in ed_params.f90.                                    !
   !---------------------------------------------------------------------------------------!
   !----- Set to 1 if to do forest harvesting. --------------------------------------------!
   integer :: forestry_on
   !----- Set to 1 if to do agriculture. --------------------------------------------------!
   integer :: agriculture_on
   !----- Earliest year at which plantations occur. ---------------------------------------!
   integer :: plantation_year
   !----- Number of years that a plantation requires to reach maturity. -------------------!
   real :: plantation_rotation
   !----- Years that a non-plantation patch requires to reach maturity. -------------------!
   real :: mature_harvest_age
   !----- Minimum plantation fraction to consider the site a plantation. ------------------!
   real :: min_plantation_frac
   !----- Minimum site biomass for even trying harvest. -----------------------------------!
   real :: min_harvest_biomass
   !---------------------------------------------------------------------------------------!
   !    Maximum distance to the current polygon that we still consider the file grid point !
   ! to be representative of the polygon for plantation fraction.                          !
   !---------------------------------------------------------------------------------------!
   real :: max_plantation_dist
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     Fire parameters.                                                                  !
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Fire may occur if total equivalent water depth (ground + underground) falls below !
   ! this threshold and include_fire is 1.  Units: meters.                                 !
   !---------------------------------------------------------------------------------------!
   real :: fire_dryness_threshold 

   !---------------------------------------------------------------------------------------!
   ! SM_FIRE        -- This is used only when INCLUDE_FIRE = 2 or 3, and it has different  !
   !                   meanings.  The sign here matters.                                   !
   !                   When INCLUDE_FIRE = 2:                                              !
   !                      >= 0. - Minimum relative soil moisture above dry air of the top  !
   !                              1m that will prevent fires to happen.                    !
   !                      <  0. - Minimum mean soil moisture potential in MPa of the top   !
   !                              1m that will prevent fires to happen.  The dry air soil  !
   !                              potential is defined as -3.1 MPa, so make sure SM_FIRE   !
   !                              is greater than this value.                              !
   !                   When INCLUDE_FIRE = 3, only positive values are allowed.  This is   !
   !                   the minimum water deficit, in kg/m2/30 days, to trigger fires.      !
   !---------------------------------------------------------------------------------------!
   real :: sm_fire

   !---------------------------------------------------------------------------------------!
   !     Depth to be compared with the soil average when include_fire is 2 or 3.           !
   ! Units: meters, and this must be negative, consistent with slz.                        !
   !---------------------------------------------------------------------------------------!
   real :: fire_smoist_depth         

   !----- k level of the deepest layer to be considered. ----------------------------------!
   integer :: k_fire_first
   !=======================================================================================!
   !=======================================================================================!



   !=======================================================================================!
   !=======================================================================================!
   !      Parameters for the FIRESTARTER model.                                            !
   !                                                                                       !
   ! Lower bound: minimum value for the variable in the model, even if the actual value    !
   !              falls below this limit.                                                  !
   ! Upper bound: maximum value for the variable in the model, even if the actual value    !
   !              goes above this limit.                                                   !
   !                                                                                       !
   ! Note that the parameters for HDI replace the parameters for GDP (as it is already     !
   ! normalised and accounts for other socio-economic factors).  Likewise, we use soil     !
   ! matric potential instead of soil moisture as the former is more descriptive of        !
   ! critical dryness than the latter.                                                     !
   !                                                                                       !
   ! Note: Parameters denoted with (*) should NOT be initialised through xml.              !
   !---------------------------------------------------------------------------------------!
   !----- Global parameters. --------------------------------------------------------------!
   real(kind=4) :: fh_grid       ! HESFIRE equivalent grid size                [       deg]
   real(kind=4) :: fh_f0001      ! Fraction of struct. C that is 1-hr fuel     [        --]
   real(kind=4) :: fh_f0010      ! Fraction of struct. C that is 10-hr fuel    [        --]
   real(kind=4) :: fh_f0100      ! Fraction of struct. C that is 100-hr fuel   [        --]
   real(kind=4) :: fh_f1000      ! Fraction of struct. C that is 1000-hr fuel  [        --]
   real(kind=4) :: fh_pcpg_ni0   ! Precipitation rate above which NI is reset  [   kg/m2/s]
   !----- Ignition parameters. ------------------------------------------------------------!
   real(kind=4) :: fi_cg_ignp    ! Cloud-to-ground ignition probability        [       ---]
   real(kind=4) :: fi_lu_ignd    ! Land use ignition density                   [    1/m2/s]
   real(kind=4) :: fi_sf_maxage  ! Max. age for secondary forests to be LU     [        yr]
   real(kind=4) :: fi_lu_exp     ! Land use exponent                           [       ---]
   real(kind=4) :: fi_lu_upr     ! Upper bound for land use (saturation point) [       ---]
   real(kind=4) :: fi_lu_off     ! (*) Offset for the land use ignition        [       ---]
   real(kind=4) :: fi_hdi_upr    ! Upper bound for HDI effect on termination   [       ---]
   real(kind=4) :: fi_hdi_exp    ! HDI shape parameter to modulate ignitions   [       ---]
   !----- Spread parameters. --------------------------------------------------------------!
   real(kind=4) :: fs_ba_frag    ! Burnt area fragmentation (min age to burn)  [        yr]
   real(kind=4) :: fs_rhv_lwr    ! Lower bound for relative humidity           [       ---]
   real(kind=4) :: fs_rhv_upr    ! Upper bound for relative humidity           [       ---]
   real(kind=4) :: fs_rhv_dti    ! (*) 1. / ( fs_rhv_upr - fs_rhv_lwr )        [       ---]
   real(kind=4) :: fs_rhv_exp    ! Exponent for relative humidity              [       ---]
   real(kind=4) :: fs_smpot_lwr  ! Lower bound for soil matric potential       [         m]
   real(kind=4) :: fs_smpot_upr  ! Upper bound for soil matric potential       [         m]
   real(kind=4) :: fs_smpot_dti  ! (*) 1. / ( fs_smpot_upr - fs_smpot_lwr )    [       1/m]
   real(kind=4) :: fs_smpot_exp  ! Exponent for soil matric potential          [       ---]
   real(kind=4) :: fs_temp_lwr   ! Lower bound for temperature                 [         K]
   real(kind=4) :: fs_temp_upr   ! Upper bound for temperature                 [         K]
   real(kind=4) :: fs_temp_dti   ! (*) 1. / ( fs_temp_upr - fs_temp_lwr )      [       1/K]
   real(kind=4) :: fs_temp_exp   ! Exponent for temperature                    [       ---]
   real(kind=4) :: fs_lbr_slp    ! Slope of the length-breadth ratio           [       ---]
   real(kind=4) :: fs_lbr_exp    ! Exponential factor for wind                 [       s/m]
   real(kind=4) :: fs_gw_infty   ! Value of g(W) at maximum wind speed         [       ---]
   !----- Fire intensity parameters. ------------------------------------------------------!
   real(kind=4) :: fx_a0001      ! Moist. sens. parameter (1-hr fuel)          [   1/degC2]
   real(kind=4) :: fx_a0010      ! Moist. sens. parameter (10-hr fuel)         [   1/degC2]
   real(kind=4) :: fx_a0100      ! Moist. sens. parameter (100-hr fuel)        [   1/degC2]
   real(kind=4) :: fx_rmfac      ! Relative factor for living fuel moisture    [       ---]
   real(kind=4) :: fx_tlh_slope  ! Slope for duration of lethal heating        [   m2 s/kg]
   real(kind=4) :: fx_tlc_slope  ! Slope for critical lethal heating time      [      1/cm]
   real(kind=4) :: fx_pmtau_di   ! Intercept for cambial damage mortality      [          ]
   real(kind=4) :: fx_pmtau_ds   ! Slope for cambial damage mortality          [          ]
   real(kind=4) :: fx_c0001_di   ! 1-hr fuel consump. factor: intercept/dry    [       ---]
   real(kind=4) :: fx_c0001_ds   ! 1-hr fuel consump. factor: slope/dry        [       ---]
   real(kind=4) :: fx_c0001_mi   ! 1-hr fuel consump. factor: intercept/moist  [       ---]
   real(kind=4) :: fx_c0001_ms   ! 1-hr fuel consump. factor: slope/moist      [       ---]
   real(kind=4) :: fx_c0010_di   ! 10-hr fuel consump. factor: intercept/dry   [       ---]
   real(kind=4) :: fx_c0010_ds   ! 10-hr fuel consump. factor: slope/dry       [       ---]
   real(kind=4) :: fx_c0010_mi   ! 10-hr fuel consump. factor: intercept/moist [       ---]
   real(kind=4) :: fx_c0010_ms   ! 10-hr fuel consump. factor: slope/moist     [       ---]
   real(kind=4) :: fx_c0100_di   ! 100-hr fuel consump. factor: intercept/dry  [       ---]
   real(kind=4) :: fx_c0100_ds   ! 100-hr fuel consump. factor: slope/dry      [       ---]
   real(kind=4) :: fx_c0100_mi   ! 100-hr fuel consump. factor: inter./moist   [       ---]
   real(kind=4) :: fx_c0100_ms   ! 100-hr fuel consump. factor: slope/moist    [       ---]
   real(kind=4) :: fx_c1000_di   ! 1000-hr fuel consump. factor: intercept/dry [       ---]
   real(kind=4) :: fx_c1000_ds   ! 1000-hr fuel consump. factor: slope/dry     [       ---]
   real(kind=4) :: fx_c1000_mi   ! 1000-hr fuel consump. factor: inter./moist  [       ---]
   real(kind=4) :: fx_c1000_ms   ! 1000-hr fuel consump. factor: slope/moist   [       ---]
   !----- Termination parameters. ---------------------------------------------------------!
   real(kind=4) :: ft_fint_lwr   ! Lower bound for fire intensity              [       W/m]
   real(kind=4) :: ft_fint_upr   ! Upper bound for fire intensity              [       W/m]
   real(kind=4) :: ft_fint_exp   ! Exponent for fire intensity                 [       ---]
   real(kind=4) :: ft_fint_dti   ! (*) 1. / ( ft_fint_upr - ft_fint_lwr )      [       m/W]
   real(kind=4) :: ft_frag_exp   ! Exponent for fragmentation effect           [       ---]
   real(kind=4) :: ft_lu_upr     ! Upper bound for land use effect             [       ---]
   real(kind=4) :: ft_lu_exp     ! Exponent for land use effect                [       ---]
   real(kind=4) :: ft_hdi_upr    ! Upper bound for HDI effect on termination   [       ---]
   real(kind=4) :: ft_hdi_exp    ! Exponent for HDI effect on termination      [       ---]
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     Constants used to define the maximum spread rate.  These are defined based on the !
   ! A18's revision of the R72 fire spread model, using "very dry" moisture conditions as  !
   ! defined by SB05.                                                                      !
   !                                                                                       !
   ! References:                                                                           !
   !                                                                                       !
   ! Andrews PL. 2018. The Rothermel surface fire spread model and associated develop-     !
   !    ments: A compre- hensive explanation. Gen. Tech. Rep. RMRS-GTR-371, U.S.           !
   !    Department of Agriculture, Forest Service, Rocky Mountain Research Station, Fort   !
   !    Collins, CO, U.S.A. https://www.fs.usda.gov/treesearch/pubs/55928 (A18).           !
   !                                                                                       !
   ! Rothermel RC. 1972. A mathematical model for predicting fire spread in wildland       !
   !    fuels. Res. Pap. INT- 115, U.S. Department of Agriculture, Intermountain Forest    !
   !    and Range Experiment Station, Ogden, UT, U. S. A.,                                 !
   !    https://www.fs.usda.gov/treesearch/pubs/32533 (R72).                               !
   !                                                                                       !
   ! Scott JH , Burgan RE. 2005. Standard fire behavior fuel models: a comprehensive set   !
   !    for use with Rothermel's surface fire spread model. Gen. Tech. Rep. RMRS-GTR-153,  !
   !    U.S. Department of Agriculture, Forest Service, Rocky Mountain Research Station,   !
   !    Fort Collins, CO, U.S.A. doi:10.2737/RMRS-GTR-153 (SB05).                          !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
   !      Dimensions (currently these are hard parameter.                                  !
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
   !----- Maximum number of SAV ratio breaks.  This has to be a fixed parameter. ----------!
   integer, parameter :: n_sbmax  = 20
   !----- Number of fuel statuses (dead/alive). -------------------------------------------!
   integer, parameter :: n_fst    = 2
   !----- Number of fuel classes. (1-hr, 10-hr, 100-hr, 1000-hr, herb, woody). ------------!
   integer, parameter :: n_fcl    = 6
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
   !----- Total mineral content.                   [     --]. -----------------------------!
   real(kind=4)                            :: fr_ST
   !----- Effective mineral content.               [     --]. -----------------------------!
   real(kind=4)                            :: fr_Se
   !----- Heat content                             [   J/kg]. -----------------------------!
   real(kind=4)                            :: fr_h
   !----- Particle density                         [  kg/m3]. -----------------------------!
   real(kind=4)                            :: fr_rhop
   !----- Number of SAV ratio bins. -------------------------------------------------------!
   integer                                 :: n_sbins
   !----- SAV ratio breaks for net load            [    1/m]. -----------------------------!
   real(kind=4), dimension(n_sbmax)        :: fr_sig_brks
   !----- Exponential factors for moisture of extinction (1=dead,2=live). -----------------!
   real(kind=4), dimension(2)              :: fr_epsil_ee
   !----- Linear Coefficients to obtain live fuel moisture of extinction. -----------------!
   real(kind=4), dimension(2)              :: fr_mxl_aa
   !----- Coefficients for moisture dampening coefficient. --------------------------------!
   real(kind=4), dimension(4)              :: fr_eta_m_aa
   !----- Power coefficients for mineral dampening coefficient. ---------------------------!
   real(kind=4), dimension(2)              :: fr_eta_s_uu
   !----- Coefficients for optimal packing ratio. -----------------------------------------!
   real(kind=4), dimension(2)              :: fr_beta_op_uu
   !----- Coefficients for maximum reaction velocity. -------------------------------------!
   real(kind=4), dimension(4)              :: fr_gamma_xx
   !----- Coefficients for ancillary variable used by optimim reaction velocity. ----------!
   real(kind=4), dimension(2)              :: fr_AA_uu
   !----- Coefficients for propagating flux ratio. ----------------------------------------!
   real(kind=4), dimension(6)              :: fr_xi_xx
   !----- Coefficients for corrected saturated wind speed. --------------------------------!
   real(kind=4), dimension(2)              :: fr_Umax_uu
   !----- Coefficients for intermediate parameters used by Rothermel's wind function. -----!
   real(kind=4), dimension(2)              :: fr_BB_uu
   real(kind=4), dimension(4)              :: fr_CC_xx
   real(kind=4), dimension(2)              :: fr_EE_ee
   !----- Coefficients for heat of pre-ignition. ------------------------------------------!
   real(kind=4), dimension(2)              :: fr_Qig_aa
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
   !      The following variables are fuel properties by class and status.  These are      !
   ! organised into matrices, following the original R72/A18 model.                        !
   !                                                                                       !
   !  Rows    -- Fuel classes: 1-hr, 10-hr, 100-hr, 1000-hr, herb, woody.                  !
   !  Columns -- Fuel status: dead and alive                                               !
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
   !------ Default SAV ratios              [     1/m]. ------------------------------------!
   real(kind=4), dimension(n_fst,n_fcl) :: fr_sigma_ij
   !------ Default "dry" moisture (used for estimating the maximum ROS) [kg_H2O/kg_dry]. --!
   real(kind=4), dimension(n_fst,n_fcl) :: fr_moist_ij
   !------ Default low heat content            [    J/kg]. --------------------------------!
   real(kind=4), dimension(n_fst,n_fcl) :: fr_hh_ij
   !------ Flag for dead/alive pools. (0 = alive; 1 = dead, set to real for convenience). -!
   real(kind=4), dimension(      n_fcl) :: fr_dead_j
   !------ Indices for sigma by size classes. ---------------------------------------------!
   integer     , dimension(n_fst,n_fcl) :: fr_sgclss_ij
   !------ Dummy value for effective sigma when fuel area is zero [1/m]. ------------------!
   real                                 :: fr_sigma_00
   !------ Dummy value for live-to-dead ratio. --------------------------------------------!
   real                                 :: fr_g_W_00
   !------ Fuel moisture extinction (dead fuels). -----------------------------------------!
   real                                 :: fr_Mxdead
   !------ Effective mineral content (currently the same for dead and live fuels). --------!
   real(kind=4), dimension(n_fst      ) :: fr_Se_i
   !------ Fuel depth [m]. ----------------------------------------------------------------!
   real                                 :: fr_depth
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     Variable type that contains the land use change disturbances.                     !
   !---------------------------------------------------------------------------------------!
   type lutime
      !------ The disturbance year. -------------------------------------------------------!
      integer :: landuse_year ! the year

      !------------------------------------------------------------------------------------!
      !    The land use information:  the current columns are:                             !
      !  ====== Disturbance rates ======                                                   !
      !  1 - Cropland to pasture                                           [         1/yr] !
      !  2 - Pasture to cropland                                           [         1/yr] !
      !  3 - Pasture to primary forest                                     [         1/yr] !
      !  4 - Primary forest to pasture                                     [         1/yr] !
      !  5 - Primary forest to cropland                                    [         1/yr] !
      !  6 - Cropland to primary forest                                    [         1/yr] !
      !  7 - Secondary forest to cropland                                  [         1/yr] !
      !  8 - Cropland to secondary forest                                  [         1/yr] !
      !  9 - Secondary forest to pasture                                   [         1/yr] !
      ! 10 - Pasture to secondary forest                                   [         1/yr] !
      ! 11 - Primary forest to secondary forest                            [         1/yr] !
      !  ====== Biomass to be harvested. ======                                            !
      ! 12 - Wood harvest on mature secondary forest land.                 [          kgC] !
      ! 13 - Wood harvest on mature secondary forest land.                 [       kgC/m²] !
      ! 14 - Wood harvest on primary forested land.                        [          kgC] !
      ! 15 - Wood harvest on primary forested land.                        [       kgC/m²] !
      ! 16 - Wood harvest on young secondary forest land.                  [          kgC] !
      ! 17 - Wood harvest on young secondary forest land.                  [       kgC/m²] !
      ! 18 - Wood harvest on primary non-forested land.                    [          kgC] !
      ! 19 - Wood harvest on primary non-forested land.                    [       kgC/m²] !
      !  ====== Special flags. ======                                                      !
      ! 12 - Secondary forest is harvested using the probability of harvesting when the    !
      !      DBH is above the minimum DBH.                                                 !
      ! 14 - Primary forest is harvested using the probability of harvesting when the DBH  !
      !      is above the minimum DBH.                                                     !
      !------------------------------------------------------------------------------------!
      real, dimension(num_lu_trans) :: landuse
   end type lutime
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     Variable type that contains the socio-economic index data, used to define fire    !
   ! ignition and suppression.  We currently assume one data entry per year.  This could   !
   ! be modified if monthly data become available.                                         !
   !---------------------------------------------------------------------------------------!
   type seitime
      integer :: sei_year ! the year                                            [       --]
      real    :: hdi      ! Human Development Index                             [       --]
      real    :: gdpc     ! GDP                                                 [ USD yr-1]
   end type seitime
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     Variable type that contains the flash rate density data, used to define fire      !
   ! ignition.  We currently assume that only lightning climatology data are available, so !
   ! we provide one data entry per month.  The structure is implemented so this could      !
   ! be easily modified if monthly data with interannual variability become available.     !
   !---------------------------------------------------------------------------------------!
   type flashtime
      integer :: flash_month ! the month                                        [       --]
      real    :: frd         ! Flash rate density (any flash)                   [   1/m2/s]
      real    :: c2g         ! Cloud-to-ground flash rate density               [   1/m2/s]
   end type flashtime
   !=======================================================================================!
   !=======================================================================================!
end module disturb_coms
!==========================================================================================!
!==========================================================================================!
