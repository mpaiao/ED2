!==========================================================================================!
!==========================================================================================!
module mortality
   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the total PFT-dependent mortality rate:                   !
   !---------------------------------------------------------------------------------------!
   subroutine mortality_rates(cpatch,ico,avg_daily_temp, patch_age,dist_type)
      use ed_state_vars , only : patchtype                  ! ! Structure
      use pft_coms      , only : mort0                      & ! intent(in)
                               , mort1                      & ! intent(in)
                               , mort2                      & ! intent(in)
                               , mort3                      & ! intent(in)
                               , hydro_mort0                & ! intent(in)
                               , hydro_mort1                & ! intent(in)
                               , plant_min_temp             & ! intent(in)
                               , frost_mort                 & ! intent(in)
                               , cbr_severe_stress          ! ! intent(in)
      use disturb_coms  , only : treefall_disturbance_rate  & ! intent(in)
                               , treefall_hite_threshold    & ! intent(in)
                               , time2canopy                ! ! intent(in)
      use ed_max_dims   , only : n_pft                      ! ! intent(in)
      use physiology_coms,only : carbon_mortality_scheme
      use consts_coms   , only : lnexp_min                  & ! intent(in)
                               , lnexp_max                  ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(patchtype), target     :: cpatch          ! Current patch
      integer        , intent(in) :: ico             ! Current cohort ID
      real           , intent(in) :: avg_daily_temp  ! Mean temperature yesterday
      real           , intent(in) :: patch_age       ! Patch age
      integer        , intent(in) :: dist_type       ! Disturbance type.
      !----- Local variables --------------------------------------------------------------!
      integer                     :: ipft            ! PFT 
      real                        :: temp_dep        ! Temp. function  (frost mortality)
      real                        :: expmort         ! Carbon-balance term
      real                        :: cbr_use         ! Bounded carbon balance.
      real                        :: growth_past_year! DBH growth rates in the past year
      !------------------------------------------------------------------------------------!


      !----- Assume happy end, all plants survive... --------------------------------------!
      cpatch%mort_rate(1:5,ico) = 0.0
      ipft = cpatch%pft(ico)

      !------------------------------------------------------------------------------------!
      ! 1.  Ageing, PFT-dependent but otherwise constant.                                  !
      !------------------------------------------------------------------------------------!
      cpatch%mort_rate(1,ico) = mort3(ipft)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      ! 2.  Mortality rates due to negative carbon balance.   Note that the functional     !
      !     form is slightly different for economics_scheme = 1.                           !
      !------------------------------------------------------------------------------------!
      cbr_use = max(cpatch%cbr_bar(ico),cbr_severe_stress(ipft))
      expmort = max( lnexp_min, min( lnexp_max,mort2(ipft) * ( cbr_use - mort0(ipft) ) ) )
      select case (carbon_mortality_scheme)
      case (2)
         !----- Camac et al (2017).  But use absolute growth rates ------------------------!
         growth_past_year = sum(cpatch%ddbh_monthly(1:12,ico)) / 12.0
         expmort = max( lnexp_min, min( lnexp_max, mort2(ipft) * growth_past_year))
         cpatch%mort_rate(2,ico) = mort1(ipft) * exp(-expmort)

      case (1)
         !----- Camac et al (2017).  Mind the minus sign. ---------------------------------!
         cpatch%mort_rate(2,ico) = mort1(ipft) * exp(-expmort)
         !---------------------------------------------------------------------------------!
      case default
         !----- Moorcroft et al (2001).  Exponential should not have minus sign. ----------!
         cpatch%mort_rate(2,ico) = mort1(ipft) / (1. + exp(expmort))
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      ! 3.  Background mortality, defined by treefall_disturbance_rate.  For small trees,  !
      !     this mortality is applied here because they aren't big enough to generate a    !
      !     new gap.  Likewise, this rate is applied to all trees in a forest plantation,  !
      !     because otherwise treefall would create a new patch and the land would         !
      !     transition from managed to unmanaged.                                          !
      !------------------------------------------------------------------------------------!
      if (dist_type == 2) then
         cpatch%mort_rate(3,ico) = treefall_disturbance_rate
      elseif ( cpatch%hite(ico) <= treefall_hite_threshold .and.                           &
               patch_age        >  time2canopy             ) then
         cpatch%mort_rate(3,ico) = treefall_disturbance_rate
      else
         cpatch%mort_rate(3,ico) = 0.
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! 4.   Mortality due to cold, after:                                                 !
      !      Albani, M.; D. Medvigy; G. C. Hurtt; P. R. Moorcroft, 2006: The contributions !
      !           of land-use change, CO2 fertilization, and climate variability to the    !
      !           Eastern US carbon sink.  Glob. Change Biol., 12, 2370-2390,              !
      !           doi: 10.1111/j.1365-2486.2006.01254.x                                    !
      !------------------------------------------------------------------------------------!
      temp_dep = max( 0.0, min( 1.0, 1.0 - (avg_daily_temp - plant_min_temp(ipft)) / 5.0) )
      cpatch%mort_rate(4,ico) = frost_mort(ipft) * temp_dep
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      ! 5. Hydraulic failure moratlity. Exponential increases of mortality rates with PLC  !
      !------------------------------------------------------------------------------------!
      cpatch%mort_rate(5,ico) = sum( hydro_mort0(ipft) *                                   &
                                     cpatch%plc_monthly(1:12,ico) ** hydro_mort1(ipft)     &
                                   ) / 12.0

      !------------------------------------------------------------------------------------!
      ! 6. Disturbance rate mortality.  This is not used by the cohort dynamics, instead   !
      !    this is just to account for the lost density due to the patch creation.  This   !
      !    mortality will be determined by the disturbance_mortality subroutine, not here. !
      !------------------------------------------------------------------------------------!
      !cpatch%mort_rate(6,ico) = TBD
      !------------------------------------------------------------------------------------!

      return
   end subroutine mortality_rates
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine determines the mortality rates associated with the current        !
   ! disturbance.                                                                          !
   !---------------------------------------------------------------------------------------!
   subroutine disturbance_mortality(csite,ipa,area_loss,mindbh_harvest                     &
                                   ,avg_fire_intensity,avg_fire_tlethal)
      use ed_state_vars, only : sitetype      & ! structure
                              , patchtype     ! ! structure
      use ed_max_dims  , only : n_pft         & ! intent(in)
                              , n_dist_types  ! ! intent(in)
      use consts_coms  , only : lnexp_max     & ! intent(in)
                              , tiny_num      & ! intent(in)
                              , almost_one    ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)                         , target      :: csite
      integer                                , intent(in)  :: ipa
      real          , dimension(n_dist_types), intent(in)  :: area_loss
      real          , dimension(n_pft)       , intent(in)  :: mindbh_harvest
      real          , dimension(12)          , intent(in)  :: avg_fire_intensity
      real          , dimension(12)          , intent(in)  :: avg_fire_tlethal
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype)                        , pointer     :: cpatch
      integer                                              :: ico
      integer                                              :: new_lu
      real                                                 :: f_survival
      real           , dimension(:)          , allocatable :: a_factor
      !------------------------------------------------------------------------------------!


      !----- Current patch, in case it is empty, return. ----------------------------------!
      cpatch => csite%patch(ipa)
      if (cpatch%ncohorts == 0) return
      !------------------------------------------------------------------------------------!


      !----- Allocate the "a_factor", which will integrate all disturbances. --------------!
      allocate(a_factor(cpatch%ncohorts))
      a_factor(:) = 0.0
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Loop over new disturbance types, add survivors from each disturbance type.     !
      !------------------------------------------------------------------------------------!
      do new_lu=1,n_dist_types
         if (area_loss(new_lu) > tiny_num) then
            do ico=1,cpatch%ncohorts
              f_survival    = survivorship(new_lu,csite%dist_type(ipa),mindbh_harvest      &
                                          ,avg_fire_intensity,avg_fire_tlethal,cpatch,ico)
              a_factor(ico) = a_factor(ico)                                                &
                            + ( 1.0 - f_survival ) * area_loss(new_lu) / csite%area(ipa)
            end do
         end if
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Loop over cohorts, and find mortality.                                         !
      !------------------------------------------------------------------------------------!
      do ico=1,cpatch%ncohorts
         if ( a_factor(ico) < almost_one ) then
            cpatch%mort_rate(6,ico) = log( 1.0 / (1.0 - a_factor(ico)) )
         else
            cpatch%mort_rate(6,ico) = lnexp_max
         end if
      end do
      !------------------------------------------------------------------------------------!

      deallocate(a_factor)
      return
   end subroutine disturbance_mortality
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the survivorship rate associated with a disturbance.       !
   !  Input variables:                                                                     !
   !  -- new_lu/old_lu: the disturbance/land use type after/before disturbance             !
   !     1. Pasture.                                                                       !
   !     2. Forest plantation.                                                             !
   !     3. Tree fall.                                                                     !
   !     4. Fire.                                                                          !
   !     5. Forest regrowth.                                                               !
   !     6. Logging (tree felling).                                                        !
   !     7. Logging (mechanical damage).                                                   !
   !     8. Cropland.                                                                      !
   !  -- mindbh_harvest: minimum DBH for selective logging.  All trees above threshold     !
   !                     will be logged in the tree felling patch.                         !
   !  -- cpatch: current patch.                                                            !
   !  -- ico: index for current cohort.                                                    !
   !---------------------------------------------------------------------------------------!
   real function survivorship(new_lu,old_lu,mindbh_harvest                                 &
                             ,avg_fire_intensity,avg_fire_tlethal,cpatch,ico)
      use ed_state_vars, only : patchtype                ! ! structure
      use disturb_coms , only : include_fire             & ! intent(in)
                              , treefall_hite_threshold  & ! intent(in)
                              , fx_tlc_slope             & ! intent(in)
                              , fx_pmtau_di              & ! intent(in)
                              , fx_pmtau_ds              ! ! intent(in)
      use pft_coms     , only : treefall_s_ltht          & ! intent(in)
                              , treefall_s_gtht          & ! intent(in)
                              , fire_s_min               & ! intent(in)
                              , fire_s_max               & ! intent(in)
                              , fire_s_inter             & ! intent(in)
                              , fire_s_slope             & ! intent(in)
                              , felling_s_gtharv         & ! intent(in)
                              , felling_s_ltharv         & ! intent(in)
                              , skid_s_ltharv            & ! intent(in)
                              , skid_s_gtharv            & ! intent(in)
                              , fscorch                  & ! intent(in)
                              , fx_rck_pft               & ! intent(in)
                              , fx_pck_pft               ! ! intent(in)
      use ed_max_dims  , only : n_pft                    ! ! intent(in)
      use consts_coms  , only : lnexp_min                & ! intent(in)
                              , lnexp_max                & ! intent(in)
                              , tiny_num                 & ! intent(in)
                              , onetwelfth               ! ! intent(in)
      use allometry    , only : h2crownbh                ! ! function
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype)                 , target     :: cpatch
      real          , dimension(n_pft), intent(in) :: mindbh_harvest
      real          , dimension(12)   , intent(in) :: avg_fire_intensity
      real          , dimension(12)   , intent(in) :: avg_fire_tlethal
      integer                         , intent(in) :: ico
      integer                         , intent(in) :: new_lu
      integer                         , intent(in) :: old_lu
      !----- Local variables. -------------------------------------------------------------!
      integer                                      :: ipft
      real                                         :: lnexp
      real          , dimension(12)                :: avg_scorch_height
      real          , dimension(12)                :: avg_crown_damage
      real          , dimension(12)                :: avg_tl_ratio
      real          , dimension(12)                :: avg_pmtau
      real          , dimension(12)                :: avg_pmck
      real          , dimension(12)                :: avg_mortality
      real                                         :: tlethal_crit
      real                                         :: chbase
      real                                         :: clength
      !------------------------------------------------------------------------------------!


      !----- Alias for current PFT. -------------------------------------------------------!
      ipft = cpatch%pft(ico)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      Survivorship depends on the disturbance type (new_lu), and within each        !
      ! disturbance type, it may depend on the disturbance pathway (dist_path), the PFT,   !
      ! and size.                                                                          !
      !------------------------------------------------------------------------------------!
      select case(new_lu)
      case (1,2,8)
         !---------------------------------------------------------------------------------!
         !     Clear cut (cropland/pasture/forest plantation).  Nothing survives.          !
         !---------------------------------------------------------------------------------!
         survivorship = 0.0
         !---------------------------------------------------------------------------------!
      case (3)
         !---------------------------------------------------------------------------------!
         !     Tree fall.  Mortality depends on the cohort height and PFT.                 !
         !---------------------------------------------------------------------------------!
         if (cpatch%hite(ico) < treefall_hite_threshold) then
            survivorship = treefall_s_ltht(ipft)
         else
            survivorship = treefall_s_gtht(ipft)
         end if
         !---------------------------------------------------------------------------------!

      case (4)
         !---------------------------------------------------------------------------------!
         !     Fire.  Survivorship will depend on the fire model.                          !
         !---------------------------------------------------------------------------------!
         select case (include_fire)
         case (4)
            !------------------------------------------------------------------------------!
            !      SPITFIRE-based approach.  We calculate the average mortality of each    !
            ! month, based on the average fire intensity and scorch height.                !
            !------------------------------------------------------------------------------!

            !------ Average scorch height [m]. --------------------------------------------!
            avg_scorch_height(:) = fscorch(ipft) * avg_fire_intensity(:)
            !------------------------------------------------------------------------------!

            !------ Average fire damage. --------------------------------------------------!
            chbase               = h2crownbh(cpatch%hite(ico),ipft)
            clength              = cpatch%hite(ico) - chbase
            avg_crown_damage (:) = ( avg_scorch_height(:) - chbase ) / clength
            avg_crown_damage (:) = max(0.,avg_crown_damage(:))
            !------------------------------------------------------------------------------!


            !------ Find critical duration of lethal heating. -----------------------------!
            tlethal_crit = fx_tlc_slope * cpatch%thbark(ico) * cpatch%thbark(ico)
            !------------------------------------------------------------------------------!


            !------ Find how close to critical fire duration the system has been. ---------!
            if ( tlethal_crit > tiny_num ) then
               !------ Find the tl:tc ratio. ----------------------------------------------!
               avg_tl_ratio (:) =  avg_fire_tlethal(:) / tlethal_crit
               !---------------------------------------------------------------------------!
            else
               !---------------------------------------------------------------------------!
               !     Critical duration is zero, assume tl:tc for maximum mortality in case !
               ! there is any fire, or zero otherwise.                                     !
               !---------------------------------------------------------------------------!
               avg_tl_ratio (:) = merge( (1.0-fx_pmtau_di)/fx_pmtau_ds                     &
                                       , 0.0                                               &
                                       , avg_fire_intensity(:) > 0.    )
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!


            !------ Find the mortality probability due to cambial damage. -----------------!
            avg_pmtau(:) = max(0.,min(1.,fx_pmtau_di + fx_pmtau_ds * avg_tl_ratio(:)))
            !------------------------------------------------------------------------------!


            !------ Find the mortality probability due to crown damage. -------------------!
            avg_pmck(:) = fx_rck_pft(ipft) * avg_crown_damage(:) ** fx_pck_pft(ipft)
            avg_pmck(:) = max(0.,min(1.,avg_pmck(:)))
            !------------------------------------------------------------------------------!


            !------ Find the probability of mortality. ------------------------------------!
            avg_mortality(:) = avg_pmtau(:) + avg_pmck(:) - avg_pmtau(:) * avg_pmck(:)
            survivorship     = sum(1. - avg_mortality(:)) * onetwelfth
            !------------------------------------------------------------------------------!

         case default
            !------------------------------------------------------------------------------!
            !     Old fire approach.  Fire survival rates are not dependent upon fire      !
            ! intensity or flame height.  Survival rates are a function of bark thickness  !
            ! (and size as BT depends on DBH and height).   The original scheme kills all  !
            ! individuals and this is maintained by setting both fire_s_min and fire_s_max !
            ! to 0.                                                                        !
            !------------------------------------------------------------------------------!
            lnexp        = fire_s_inter(ipft) + fire_s_slope(ipft) * cpatch%thbark(ico)
            lnexp        = max(lnexp_min,min(lnexp_max,lnexp))
            survivorship = fire_s_min(ipft)                                                &
                         + (fire_s_max(ipft) - fire_s_min(ipft)) / (1. + exp(lnexp))
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!

      case (5)
         !---------------------------------------------------------------------------------!
         !     Abandonment (secondary regrowth).  Two paths are possible: abandonment      !
         ! occurs after one last harvest, or the field/plantation is left as is.           !
         !---------------------------------------------------------------------------------!
         select case (old_lu)
         case (8)
            !----- Cropland: final harvest. -----------------------------------------------!
            survivorship = 0.0
            !------------------------------------------------------------------------------!
         case (2)
            !------------------------------------------------------------------------------!
            !     Forest plantation.  Assume fire causes abandonment.   See fire           !
            ! explanation above.                                                           !
            !------------------------------------------------------------------------------!
            lnexp        = fire_s_inter(ipft) + fire_s_slope(ipft) * cpatch%thbark(ico)
            lnexp        = max(lnexp_min,min(lnexp_max,lnexp))
            survivorship = fire_s_min(ipft)                                                &
                         + (fire_s_max(ipft) - fire_s_min(ipft)) / (1. + exp(lnexp))
            !------------------------------------------------------------------------------!
         case default
            !------------------------------------------------------------------------------!
            !    The only other option is pasture.  Leaving it as a default: everything    !
            ! survives.                                                                    !
            !------------------------------------------------------------------------------!
            survivorship = 1.0
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!
      case (6)
         !---------------------------------------------------------------------------------!
         !     Tree felling.                                                               !
         !---------------------------------------------------------------------------------!
         if (cpatch%dbh(ico) >= mindbh_harvest(ipft)) then
            survivorship = felling_s_gtharv(ipft)
         else
            survivorship = felling_s_ltharv(ipft)
         end if
         !---------------------------------------------------------------------------------!
      case (7)
         !---------------------------------------------------------------------------------!
         !     Collateral damage from logging (skid trails, roads).  The damage currently  !
         ! takes into account the minimum harvest size, because presumably loggers avoid   !
         ! damaging trees that may be potentially harvested.                               !
         !---------------------------------------------------------------------------------!
         if (cpatch%dbh(ico) >= mindbh_harvest(ipft)) then
            survivorship = skid_s_gtharv(ipft)
         else
            survivorship = skid_s_ltharv(ipft)
         end if
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!

      return
   end function survivorship
   !=======================================================================================!
   !=======================================================================================!
end module mortality
!==========================================================================================!
!==========================================================================================!
