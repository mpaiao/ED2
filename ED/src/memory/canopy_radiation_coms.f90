!==========================================================================================!
!==========================================================================================!
!    This module contains several parameters used in the canopy radiation solver.          !
!                                                                                          !
!  IMPORTANT: Do not initialize non-parameters in their modules - not all compilers will   !
!             actually initialize them.  Instead, assign them at init_can_rad_params sub-  !
!             routine (ed_params.f90).                                                     !
!------------------------------------------------------------------------------------------!
module canopy_radiation_coms

   use ed_max_dims , only : n_pft
   implicit none 


   !---------------------------------------------------------------------------------------!
   ! ICANRAD -- Specifies how vertical canopy radiation is solved.  This variable sets     !
   !            both shortwave and longwave.                                               !
   !            0.  Two-stream model (Medvigy 2006), with the possibility to apply         !
   !                finite crown area to direct shortwave radiation.                       !
   !            1.  Multiple-scattering model (Zhao and Qualls 2005,2006), with the        !
   !                possibility to apply finite crown area to all radiation fluxes.        !
   !---------------------------------------------------------------------------------------!
   integer :: icanrad
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! IHRZRAD      -- Specifies how horizontal canopy radiation is solved.                  !
   !                 0.  Default ED-2.0: no horizontal patch shading.  All patches receive !
   !                     the same amount of light at the top.                              !
   !                 1.  A realized map of the plant community is built by randomly        !
   !                     assigning gaps associated with gaps (number of gaps proportional  !
   !                     to the patch area), and populating them with individuals,         !
   !                     respecting the cohort distribution in each patch.  The crown      !
   !                     closure index is calculated for the entire landscape and used     !
   !                     to change the amount of direct light reaching the top of the      !
   !                     canopy.  Patches are then split into 1-3 patches based on the     !
   !                     light condition (so expect simulations to be slower).  This       !
   !                     method is under development, suggestions on how to improve are    !
   !                     welcome.                                                          !
   !                 2.  Similar to option 1, except that height for trees with DBH >      !
   !                     DBH_crit are rescaled to calculate CCI.                           !
   !                 3.  Dummy horizontal canopy radiation.  This applies the same method  !
   !                     as 1 and 2 to split patches, but it does not change radiation     !
   !                     reaching the top of the canopy.  This is only useful to isolate   !
   !                     the effect of heterogeneous illumination from the patch count.    !
   !                 4.  Same as 0., but patch fusion takes into account the correction    !
   !                     for emergent trees.                                               !
   !---------------------------------------------------------------------------------------!
   integer :: ihrzrad
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     The following variables are temporary namelist variables used to control the      !
   ! radiation properties of leaves.                                                       !
   ! LTRANS_VIS   -- Leaf transmittance on visible.                                        !
   ! LTRANS_NIR   -- Leaf transmittance on near infrared.                                  !
   ! LREFLECT_VIS -- Leaf reflectance on visible.                                          !
   ! LREFLECT_NIR -- Leaf reflectance on near infrared.                                    !
   ! ORIENT_TREE  -- Leaf orientation parameter for tropical trees                         !
   ! ORIENT_GRASS -- Leaf orientation parameter for tropical grasses                       !
   ! CLUMP_TREE   -- Leaf clumping factor for tropical trees                               !
   ! CLUMP_GRASS  -- Leaf clumping factor for tropical grasses                             !
   !---------------------------------------------------------------------------------------!
   real :: ltrans_vis
   real :: ltrans_nir
   real :: lreflect_vis
   real :: lreflect_nir
   real :: orient_tree
   real :: orient_grass
   real :: clump_tree
   real :: clump_grass
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Fraction of diffuse solar radiation in the PAR band.  Used when you don't know    !
   ! the direct/diffuse breakdown. (parameters)                                            !
   !---------------------------------------------------------------------------------------!
   real :: fvis_beam_def
   real :: fvis_diff_def
   real :: fnir_beam_def
   real :: fnir_diff_def
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Structure with scratch variables for radiation (Thanks RGK!).                      !
   !---------------------------------------------------------------------------------------!
   type radscrtype
      integer         , pointer, dimension(:)   :: pft_array
      real(kind=8)    , pointer, dimension(:)   :: leaf_temp_array
      real(kind=8)    , pointer, dimension(:)   :: wood_temp_array
      real(kind=8)    , pointer, dimension(:)   :: lai_array
      real(kind=8)    , pointer, dimension(:)   :: wai_array 
      real(kind=8)    , pointer, dimension(:)   :: CA_array
      real(kind=8)    , pointer, dimension(:)   :: htop_array
      real(kind=8)    , pointer, dimension(:)   :: hbot_array
      real(kind=8)    , pointer, dimension(:)   :: par_level_beam
      real(kind=8)    , pointer, dimension(:)   :: par_level_diffd
      real(kind=8)    , pointer, dimension(:)   :: par_level_diffu
      real(kind=8)    , pointer, dimension(:)   :: light_level_array
      real(kind=8)    , pointer, dimension(:)   :: light_beam_level_array
      real(kind=8)    , pointer, dimension(:)   :: light_diff_level_array
      real            , pointer, dimension(:)   :: par_v_beam_array
      real            , pointer, dimension(:)   :: rshort_v_beam_array
      real            , pointer, dimension(:)   :: par_v_diffuse_array
      real            , pointer, dimension(:)   :: rshort_v_diffuse_array
      real            , pointer, dimension(:)   :: lw_v_array
      real            , pointer, dimension(:,:) :: radprof_array
   end type radscrtype
   type(radscrtype)   , pointer,dimension(:)    :: radscr(:)
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Factors that define the orientation and clumping of leaves.                       !
   ! CLUMPING FACTOR - factor indicating the degree of clumpiness of leaves.               !
   ! ORIENT_FACTOR   - mean leaf orientation.                                              !
   !                     0 -- leaves are randomly oriented                                 !
   !                     1 -- all leaves are perfectly horizontal                          !
   !                    -1 -- all leaves are perfectly vertical.                           !
   ! PHI1            - The phi1 term from the CLM technical manual                         !
   ! PHI2            - The phi2 term from the CLM technical manual                         !
   ! MU_BAR          - average cosine of incidence angle for hemispheric (diffuse)         !
   !                   radiation (for both short wave and long wave)                       !
   !---------------------------------------------------------------------------------------!
   real(kind=8), dimension(n_pft) :: clumping_factor
   real(kind=8), dimension(n_pft) :: orient_factor
   real(kind=8), dimension(n_pft) :: phi1
   real(kind=8), dimension(n_pft) :: phi2
   real(kind=8), dimension(n_pft) :: mu_bar
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Reflectance coefficients.                                                         !
   !---------------------------------------------------------------------------------------!
   !----- Visible (PAR). ------------------------------------------------------------------!
   real(kind=8), dimension(n_pft) :: leaf_reflect_vis
   real(kind=8), dimension(n_pft) :: wood_reflect_vis
   !----- Near infrared. ------------------------------------------------------------------!
   real(kind=8), dimension(n_pft) :: leaf_reflect_nir
   real(kind=8), dimension(n_pft) :: wood_reflect_nir
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Transmittance coefficients.                                                       !
   !---------------------------------------------------------------------------------------!
   !----- Visible (PAR). ------------------------------------------------------------------!
   real(kind=8), dimension(n_pft) :: leaf_trans_vis
   real(kind=8), dimension(n_pft) :: wood_trans_vis
   !----- Near infrared. ------------------------------------------------------------------!
   real(kind=8), dimension(n_pft) :: leaf_trans_nir
   real(kind=8), dimension(n_pft) :: wood_trans_nir
   !---------------------------------------------------------------------------------------!




   !----- Emissivity of the vegetation (TIR). ---------------------------------------------!
   real(kind=8), dimension(n_pft) :: leaf_emiss_tir
   real(kind=8), dimension(n_pft) :: wood_emiss_tir
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Scattering coefficients.                                                          !
   !---------------------------------------------------------------------------------------!
   !----- Visible (PAR). ------------------------------------------------------------------!
   real(kind=8), dimension(n_pft) :: leaf_scatter_vis
   real(kind=8), dimension(n_pft) :: wood_scatter_vis
   !----- Near infrared. ------------------------------------------------------------------!
   real(kind=8), dimension(n_pft) :: leaf_scatter_nir
   real(kind=8), dimension(n_pft) :: wood_scatter_nir
   !----- Thermal infrared. ---------------------------------------------------------------!
   ! real(kind=8), dimension(n_pft) :: leaf_scatter_tir
   ! real(kind=8), dimension(n_pft) :: wood_scatter_tir
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Fraction of diffuse radiation that is upscattered.                                !
   !---------------------------------------------------------------------------------------!
   !----- Visible (PAR). ------------------------------------------------------------------!
   real(kind=8), dimension(n_pft) :: leaf_backscatter_vis
   real(kind=8), dimension(n_pft) :: wood_backscatter_vis
   !----- Near infrared. ------------------------------------------------------------------!
   real(kind=8), dimension(n_pft) :: leaf_backscatter_nir
   real(kind=8), dimension(n_pft) :: wood_backscatter_nir
   !----- Backscattering of thermal infrared. ---------------------------------------------!
   real(kind=8), dimension(n_pft) :: leaf_backscatter_tir
   real(kind=8), dimension(n_pft) :: wood_backscatter_tir
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Snow pack properties.                                                              !
   !---------------------------------------------------------------------------------------!
   real(kind=4) :: snow_albedo_vis
   real(kind=4) :: snow_albedo_nir
   real(kind=4) :: snow_emiss_tir
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     The following variables control whether to call things that should be called      !
   ! when there is still some light.                                                       !
   !---------------------------------------------------------------------------------------!
   real(kind=4)    :: rshort_twilight_min
   real(kind=4)    :: cosz_min
   real(kind=8)    :: cosz_min8
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     The following variables control the method that allow light redistribution based  !
   ! on patch neighbourhood.  These are initialised in ed_xml_config.f90 or ed_params.f90. !
   !---------------------------------------------------------------------------------------!
   real(kind=4)    :: cci_radius   ! Maximum radius to calculate CCI               [     m]
   real(kind=4)    :: cci_pixres   ! Pixel resolution for TCH and CCI              [     m]
   real(kind=4)    :: cci_gapsize  ! Gap size                                      [     m]
   real(kind=4)    :: cci_gapmin   ! # of gaps associated with the smallest area   [   ---]
   integer         :: cci_nretn    ! "Return density" to generate the TCH map      [  1/m2]
   real(kind=4)    :: cci_hmax     ! Maximum height allowed in the CCI scheme      [     m]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     These variables are derived from the properties above, they will be allocated     !
   ! during the initialisation step.                                                       !
   !---------------------------------------------------------------------------------------!
   !----- Total area of each single gap. --------------------------------------------------!
   real(kind=4)                              :: cci_gaparea
   !----- Number of grid points in x and y direction (pseudo-landscape). ------------------!
   integer                                   :: rls_nxy
   !----- Number of pixels in the pseudo-landscape. ---------------------------------------!
   integer                                   :: rls_npixel
   !----- Number of gaps in the pseudo-landscape. -----------------------------------------!
   integer                                   :: rls_ngap
   !----- 'raster' length along x/y axes. -------------------------------------------------!
   real(kind=4)                              :: rls_length
   !----- Number of pixels in each gap. ---------------------------------------------------!
   integer                                   :: gap_npixel
   !----- Total 'raster' landscape area. --------------------------------------------------!
   real(kind=4)                              :: rls_area
   !----- Use fixed thresholds to split patches by illumination classes? ------------------!
   logical                                   :: fixed_hrz_classes
   !----- Default thresholds in case fixed classes are to be used. ------------------------!
   real(kind=4)                              :: at_bright_def
   real(kind=4)                              :: at_dark_def
   !----- x of the 'raster' landscape. ----------------------------------------------------!
   real(kind=4), dimension(:,:), allocatable :: rls_x
   !----- y of the 'raster' landscape. ----------------------------------------------------!
   real(kind=4), dimension(:,:), allocatable :: rls_y
   !----- Top-of-canopy height. -----------------------------------------------------------!
   real(kind=4), dimension(:,:), allocatable :: rls_ztch
   !----- Crown closure index. ------------------------------------------------------------!
   real(kind=4), dimension(:,:), allocatable :: rls_cci
   !----- Absorption correction for incident beam radiation. ------------------------------!
   real(kind=4), dimension(:,:), allocatable :: rls_fbeam
   !----- Gap indices (zero is the default). ----------------------------------------------!
   integer     , dimension(:,:), allocatable :: rls_igp0
   integer     , dimension(:,:), allocatable :: rls_igp
   integer     , dimension(:,:), allocatable :: rls_ipa
   !----- Mask to decide which gaps can be used for any patch. ----------------------------!
   logical     , dimension(:,:), allocatable :: rls_mask
   !----- Gap origin. ---------------------------------------------------------------------!
   real(kind=4), dimension(:)  , allocatable :: gap_x0
   real(kind=4), dimension(:)  , allocatable :: gap_y0
   !----- Mean absorption correction for incident beam radiation. -------------------------!
   real(kind=4), dimension(:)  , allocatable :: gap_fbeam
   integer     , dimension(:)  , allocatable :: gap_nuse
   !----- Patch associated with the gap. --------------------------------------------------!
   integer     , dimension(:)  , allocatable :: gap_ipa
   !----- Auxiliary variable, patch index before shuffling, gap index after shuffling. ----!
   integer     , dimension(:)  , allocatable :: gap_idx
   !----- Mask to decide which gaps can be used for any patch. ----------------------------!
   logical     , dimension(:)  , allocatable :: gap_mask
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Hold these parameters as constants, the functional form may change soon.          !
   !---------------------------------------------------------------------------------------!
   real(kind=4) :: at0
   real(kind=4) :: at1
   real(kind=8) :: at08
   real(kind=8) :: at18
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Define variables for computing the modified Chapman function.                     !
   !---------------------------------------------------------------------------------------!
   !----- Zenith angle resolution of the look-up table for the Chapman function. ----------!
   real(kind=8)                            :: dzen_ref
   !----- Dimension of the look-up table bins. --------------------------------------------!
   integer                                 :: nzen_ref
   !----- Reference zenith angles for the Chapman function. -------------------------------!
   real(kind=8), dimension(:), allocatable :: zend_ref
   !----- Reference values for the modified Chapman function (Huestis et al. 2001). -------!
   real(kind=8), dimension(:), allocatable :: huestis_ref
   !---------------------------------------------------------------------------------------!


   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine allocates the scratch variables after all pointers have been     !
   ! nullified.                                                                            !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_radscratch(cradscr,maxcohort)
      
      use ed_max_dims          , only : n_radprof            ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(radscrtype), target     :: cradscr
      integer         , intent(in) :: maxcohort
      !------------------------------------------------------------------------------------!

      call nullify_radscratch(cradscr)

      allocate (cradscr%pft_array                 (          maxcohort))
      allocate (cradscr%leaf_temp_array           (          maxcohort))
      allocate (cradscr%wood_temp_array           (          maxcohort))
      allocate (cradscr%lai_array                 (          maxcohort))
      allocate (cradscr%wai_array                 (          maxcohort))
      allocate (cradscr%CA_array                  (          maxcohort))
      allocate (cradscr%htop_array                (          maxcohort))
      allocate (cradscr%hbot_array                (          maxcohort))
      allocate (cradscr%par_level_beam            (          maxcohort))
      allocate (cradscr%par_level_diffu           (          maxcohort))
      allocate (cradscr%par_level_diffd           (          maxcohort))
      allocate (cradscr%light_level_array         (          maxcohort))
      allocate (cradscr%light_beam_level_array    (          maxcohort))
      allocate (cradscr%light_diff_level_array    (          maxcohort))
      allocate (cradscr%par_v_beam_array          (          maxcohort))
      allocate (cradscr%rshort_v_beam_array       (          maxcohort))
      allocate (cradscr%par_v_diffuse_array       (          maxcohort))
      allocate (cradscr%rshort_v_diffuse_array    (          maxcohort))
      allocate (cradscr%lw_v_array                (          maxcohort))
      allocate (cradscr%radprof_array             (n_radprof,maxcohort))
      
      return
   end subroutine alloc_radscratch
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine de-allocates the scratch variables.                              !
   !---------------------------------------------------------------------------------------!
   subroutine dealloc_radscratch(cradscr)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(radscrtype), target :: cradscr
      !------------------------------------------------------------------------------------!

      if (associated(cradscr%pft_array             )) deallocate(cradscr%pft_array        )
      if (associated(cradscr%leaf_temp_array       )) deallocate(cradscr%leaf_temp_array  )
      if (associated(cradscr%wood_temp_array       )) deallocate(cradscr%wood_temp_array  )
      if (associated(cradscr%lai_array             )) deallocate(cradscr%lai_array        )
      if (associated(cradscr%wai_array             )) deallocate(cradscr%wai_array        )
      if (associated(cradscr%CA_array              )) deallocate(cradscr%CA_array         )
      if (associated(cradscr%htop_array            )) deallocate(cradscr%htop_array       )
      if (associated(cradscr%hbot_array            )) deallocate(cradscr%hbot_array       )
      if (associated(cradscr%par_level_beam        )) deallocate(cradscr%par_level_beam   )
      if (associated(cradscr%par_level_diffu       )) deallocate(cradscr%par_level_diffu  )
      if (associated(cradscr%par_level_diffd       )) deallocate(cradscr%par_level_diffd  )
      if (associated(cradscr%light_level_array     )) deallocate(cradscr%light_level_array)
      if (associated(cradscr%light_beam_level_array))                                      &
                                                 deallocate(cradscr%light_beam_level_array)
      if (associated(cradscr%light_diff_level_array))                                      &
                                                 deallocate(cradscr%light_diff_level_array)
      if (associated(cradscr%par_v_beam_array      )) deallocate(cradscr%par_v_beam_array )
      if (associated(cradscr%rshort_v_beam_array   ))                                      &
                                                 deallocate(cradscr%rshort_v_beam_array   )
      if (associated(cradscr%par_v_diffuse_array   ))                                      &
                                                 deallocate(cradscr%par_v_diffuse_array   )
      if (associated(cradscr%rshort_v_diffuse_array))                                      &
                                                 deallocate(cradscr%rshort_v_diffuse_array)
      if (associated(cradscr%lw_v_array            )) deallocate(cradscr%lw_v_array       )
      if (associated(cradscr%radprof_array         )) deallocate(cradscr%radprof_array    )
      return
   end subroutine dealloc_radscratch
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine nullifies the pointers, for a safe allocation.                   !
   !---------------------------------------------------------------------------------------!
   subroutine nullify_radscratch(cradscr)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(radscrtype), target :: cradscr
      !------------------------------------------------------------------------------------!

      nullify(cradscr%pft_array                 )
      nullify(cradscr%leaf_temp_array           )
      nullify(cradscr%wood_temp_array           )
      nullify(cradscr%lai_array                 )
      nullify(cradscr%wai_array                 )
      nullify(cradscr%CA_array                  )
      nullify(cradscr%htop_array                )
      nullify(cradscr%hbot_array                )
      nullify(cradscr%par_level_beam            )
      nullify(cradscr%par_level_diffu           )
      nullify(cradscr%par_level_diffd           )
      nullify(cradscr%light_level_array         )
      nullify(cradscr%light_beam_level_array    )
      nullify(cradscr%light_diff_level_array    )
      nullify(cradscr%par_v_beam_array          )
      nullify(cradscr%rshort_v_beam_array       )
      nullify(cradscr%par_v_diffuse_array       )
      nullify(cradscr%rshort_v_diffuse_array    )
      nullify(cradscr%lw_v_array                )
      nullify(cradscr%radprof_array             )
      return
   end subroutine nullify_radscratch
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function calculates the effect of sun angle increasing the optical depth in a !
   ! sphere.  This allows accounting for the diffuse irradiance at twilight, using the     !
   ! modified Chapman function following H01.  The computation of the Chapman function is  !
   ! computationally demanding, therefore, we pre-calculate the values for a broad range   !
   ! of angles.  We do not account for the terrain altitude as this effect is minimal and  !
   ! would require substantial changes in the code.                                        !
   !                                                                                       !
   !    Reference:                                                                         !
   !                                                                                       !
   !    Huestis DL. 2001. Accurate evaluation of the chapman function for atmospheric      !
   !       attenuation. J. Quant. Spectrosc. Radiat. Transf., 69: 709-721.                 !
   !       doi:10.1016/S0022-4073(00)00107-2                                               !
   !---------------------------------------------------------------------------------------!
   subroutine set_huestis_lut(nzen,dzen,zend,huestis)
      use consts_coms, only : erad       & ! intent(in)
                            , ehgt       & ! intent(in)
                            , pio1808    & ! intent(in)
                            , tiny_num8  & ! intent(in)
                            , lnexp_min8 & ! intent(in)
                            , lnexp_max8 ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer                      , intent(in)  :: nzen    ! Number of bins
      real(kind=8)                 , intent(in)  :: dzen    ! Bin width
      real(kind=8), dimension(nzen), intent(out) :: zend    ! Reference zenith angle
      real(kind=8), dimension(nzen), intent(out) :: huestis ! Modified Chapman function
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8) :: xcurve       ! Curvature ratio
      real(kind=8) :: lambda       ! Integrating element for zenith angles
      real(kind=8) :: dlambda      ! Integrating width for zenith angles
      real(kind=8) :: sin_zend     ! Sine of zend
      real(kind=8) :: sin_lambda   ! Sine of lambda
      real(kind=8) :: cos_lambda   ! Cosine of lambda
      real(kind=8) :: ln_kernel    ! Natural logarithm of the kernel
      real(kind=8) :: kernel       ! Kernel
      real(kind=8) :: integ_kernel ! Integral of the kernel from 0 to the current angle.
      integer      :: i            ! Row counter
      integer      :: j            ! Column counter
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Dimensionless curvature ratio.                                                !
      !------------------------------------------------------------------------------------!
      xcurve = dble( erad / ehgt )
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Define the reference zenith angles.                                           !
      !------------------------------------------------------------------------------------!
      do i=1,nzen
         zend    (i) = min(1.80d2, dzen * dble(i-1) )
      end do
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Initialise integrator and the modified Chapman function for the first         !
      ! element, which should be always 1.                                                 !
      !------------------------------------------------------------------------------------!
      integ_kern = 0.d0
      huestis(1) = 1.d0
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Integrate kernel for the first bin. We use a staggered approach for computing !
      ! the kernels and deriving the modified Chapman function evaluation.                 !
      !------------------------------------------------------------------------------------!
      do i=2,nzen
         !----- Find the mid points for integrand, and the trigonometric functions. -------!
         lambda     = 5.d-1 * ( zend(i-1) + zend(i) )
         dlambda    = zend(i) - zend(i-1)
         sin_zend   = sin( zend(i) * pio1808 )
         sin_lambda = sin( lambda  * pio1808 )
         cos_lambda = cos( lambda  * pio1808 )
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find kernel. When the sine of lambda approaches zero, the kernel function   !
         ! becomes undefined, so we use the limit values instead, as our goal is to        !
         ! integrate the function.  Also, we cap the natural logarithm of the kernel to    !
         ! avoid floating point exceptions (though it should be safe with double           !
         ! precision).                                                                     !
         !---------------------------------------------------------------------------------!
         if ( abs(sin_lambda) < tiny_num8 ) then
            kernel    = 0.d0
         elseif (abs(1.d0+cos_lambda) < tiny_num8) then
            kernel    = exp(lnexp_max8)
         else
            ln_kernel = xcurve * ( 1.d0 - sin_zend / sin_lambda )
            ln_kernel = max(lnexp_min8,min(lnexp_max8,ln_kernel))
            kernel    = exp(ln_kernel) / ( 1.d0 + cos_lambda )
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Update the kernel integral and find the modified Chapman function.          !
         !---------------------------------------------------------------------------------!
         integ_kern = integ_kern + kernel * dlambda_ref * pio1808
         huestis(i) = min( exp(lnexp_max8), 1.d0 + xcurve*sin(zend*pio1808) * integ_kern)
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!

      return
   end subroutine set_huestis_lut
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This sub-routine interpolates the modified Chapman function (after H01) by inter-  !
   ! polating the pre-calculated values from the lookup table, then finds the effective    !
   ! cosine of the zenith angle.                                                           !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function find_eff_cosz(cosz)
      use consts_coms, only : pio1808     & ! intent(in)
                            , tiny_offset ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in ) :: cosz
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)              :: zen
      integer                   :: iprev
      integer                   :: inext
      real(kind=8)              :: pwr_prev
      real(kind=8)              :: pwr_next
      real(kind=8)              :: chapman
      !---- External functions. -----------------------------------------------------------!
      real(kind=4), external    :: sngloff
      !------------------------------------------------------------------------------------!



      !----- Find the zenith angle and the nearest look-up table points. ------------------!
      zen   = acos( dble(cosz) ) / pio1808
      iprev = max(1       ,floor  ( zen / dzen_ref ))
      inext = min(nzen_ref,ceiling( zen / dzen_ref ))
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Retrieve either the exact value from the look-up table (when the zenith        !
      ! angle matches a value in the look-up table or when both values of the modified     !
      ! Chapman function are identical) or log-linearly interpolate the values.            !
      !------------------------------------------------------------------------------------!
      if (huestis(iprev) == huestis(inext)) then
         chapman  = huestis_ref(iprev)
      else
         pwr_next = ( zen - zen_ref(iprev) ) / ( zen_ref(inext) - zen_ref(iprev))
         pwr_prev = 1.d0 - pwr_next
         chapman  = huestis_ref(iprev) ** pwr_prev * huestis_ref(inext) ** pwr_next
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     The effective cosine of the zenith angle is the inverse of the Chapman         !
      ! function.                                                                          !
      !------------------------------------------------------------------------------------!
      find_eff_cosz = sngloff( 1.d0 / chapman, tiny_num)
      !------------------------------------------------------------------------------------!

      return
   end function find_eff_cosz
   !=======================================================================================!
   !=======================================================================================!
end module canopy_radiation_coms
!==========================================================================================!
!==========================================================================================!
