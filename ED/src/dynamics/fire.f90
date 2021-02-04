!==========================================================================================!
!==========================================================================================!
! MODULE FIRE
!
!> \brief This module contains routines to obtain fire disturbance rates
!> \details These subroutines are intended to calculate fire intensity and burned area
!!          which are or can be used to obtain fire disturbance rate and survivorship
!> \author  Paul Moorcroft, converted to fortran by David Medvigy
!> \author  10 Jan 2018.  MLO converted it into module so the code compiles with ifort 17.
!!          Also implementing process-based model, step by step.
!------------------------------------------------------------------------------------------!
module fire

   contains 

   !=======================================================================================!
   !=======================================================================================!
   ! SUB-ROUTINE FIRE_FREQUENCY
   !> This subroutine will evaluate whether fire conditions exist, and if that is the
   !! case, it will calculate the disturbance rate due to fire.
   !---------------------------------------------------------------------------------------!
   subroutine fire_frequency(cgrid)
      use ed_state_vars , only : edtype                 & ! structure
                               , polygontype            & ! structure
                               , sitetype               & ! structure
                               , patchtype              ! ! structure
      use ed_misc_coms  , only : simtime                & ! intent(in)
                               , current_time           & ! intent(in)
                               , dtlsm                  ! ! intent(in)
      use pft_coms      , only : is_grass               ! ! intent(in)
      use disturb_coms  , only : include_fire           & ! intent(in)
                               , fire_parameter         & ! intent(in)
                               , fuel_height_max        & ! intent(in)
                               , f_combusted_fast_c     & ! intent(in)
                               , f_combusted_struct_c   ! ! intent(in)
      use consts_coms   , only : wdns                   & ! intent(in)
                               , wdnsi                  & ! intent(in)
                               , day_sec                & ! intent(in)
                               , almost_one             & ! intent(in)
                               , lnexp_max              ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)      , target     :: cgrid
      !----- Local variables --------------------------------------------------------------!
      type(polygontype) , pointer    :: cpoly
      type(sitetype)    , pointer    :: csite
      type(patchtype)   , pointer    :: cpatch
      type(simtime)                  :: lastmonth
      integer                        :: ipy
      integer                        :: isi
      integer                        :: ipa
      integer                        :: ico
      integer                        :: imo
      integer                        :: ipft
      real                           :: ndaysi
      real                           :: normfac
      real                           :: fire_wmass_threshold
      real                           :: fire_intensity
      real                           :: fuel
      real                           :: ignition_rate
      real                           :: mean_fire_intensity
      real                           :: sum_accp
      logical                        :: people_around
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the number of days of last month so we can normalise the integrated       !
      ! ground water.                                                                      !
      !------------------------------------------------------------------------------------!
      call lastmonthdate(current_time,lastmonth,ndaysi)
      normfac = dtlsm * ndaysi / (day_sec)
      !------------------------------------------------------------------------------------!


      !----- Current month. ---------------------------------------------------------------!
      imo = current_time%month
      !------------------------------------------------------------------------------------!


      !----- Loop over polygons and sites. ------------------------------------------------!
      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)


         !---------------------------------------------------------------------------------!
         !     Loop over all sites.                                                        !
         !---------------------------------------------------------------------------------!
         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            !------------------------------------------------------------------------------!
            !      Decide how to compute fire disturbance, based on the method.            !
            !------------------------------------------------------------------------------!
            select case (include_fire)
            case (4)
               !---------------------------------------------------------------------------!
               !     HESFIRE/SPITFIRE:  Fires have already been integrated over the month, !
               ! calculate disturbance area from burnt area.                               !
               !---------------------------------------------------------------------------!


               !----- Loop over patches. --------------------------------------------------!
               resetloop: do ipa=1,csite%npatches
                  !----- Reset the ground water for next month. ---------------------------!
                  csite%avg_monthly_gndwater(ipa) = 0.
                  !------------------------------------------------------------------------!
               end do resetloop
               !---------------------------------------------------------------------------!


               !----- Use burnt area to find disturbance rate. ----------------------------!
               if (cpoly%burnt_area(isi) >= almost_one) then
                  cpoly%lambda_fire(imo,isi) = lnexp_max
               else
                  cpoly%lambda_fire(imo,isi) = log( 1.0 / ( 1.0 - cpoly%burnt_area(isi) ) )
               end if
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
            case default
               !---------------------------------------------------------------------------!
               !     Other approaches, use monthly time step.                              !
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Find the total rainfall of the past year and reset the counter for    !
               ! this month.                                                               !
               !---------------------------------------------------------------------------!
               sum_accp                        = sum(cpoly%avg_monthly_accp(:,isi))
               cpoly%avg_monthly_accp(imo,isi) = 0.
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Initialize ignition rate and mean fire intensity (site variables).    !
               !---------------------------------------------------------------------------!
               ignition_rate       = 0.0
               mean_fire_intensity = 0.0
               !---------------------------------------------------------------------------!



               !----- Temporary check for human activities. -------------------------------!
               people_around = .false.
               humanloop: do ipa=1,csite%npatches
                  select case (csite%dist_type(ipa))
                  case (3)
                     continue
                  case default
                     people_around = .true.
                     exit humanloop
                  end select
               end do humanloop
               !----- Allow fires to ignite in intact forests. ----------------------------!
               people_around = .true.
               !---------------------------------------------------------------------------!



               !----- Loop over patches. --------------------------------------------------!
               patchloop: do ipa=1,csite%npatches
                  cpatch => csite%patch(ipa)

                  !----- Normalise the monthly mean ground water. -------------------------!
                  csite%avg_monthly_gndwater(ipa) = csite%avg_monthly_gndwater(ipa)        &
                                                  * normfac
                  !------------------------------------------------------------------------!


                  !----- Normalise the monthly mean ground water. -------------------------!
                  csite%avg_monthly_waterdef(ipa) = max( 0.0                               &
                                                       , csite%avg_monthly_waterdef(ipa) )
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Obtain fuel stocks.  The original fire model would consider all    !
                  ! above-ground biomass and no litter.  When include_fire is set to 3,    !
                  ! fuel is the sum of all above-ground fast soil C, and biomass from      !
                  ! grasses and small individuals (up to fuel_max_height).                 !
                  !------------------------------------------------------------------------!
                  select case (include_fire)
                  case (3)
                     fuel = csite%fast_grnd_C(ipa) + csite%structural_grnd_C(ipa)
                     fuelcohloop_3: do ico = 1,cpatch%ncohorts
                        ipft = cpatch%pft(ico)
                        if (is_grass(ipft) .or. cpatch%hite(ico) <= fuel_height_max) then
                           fuel = fuel + cpatch%nplant(ico) * cpatch%agb(ico)
                        end if
                     end do fuelcohloop_3
                  case default
                     fuel = 0.0
                     fuelcohloop_d: do ico = 1,cpatch%ncohorts
                        fuel = fuel + cpatch%nplant(ico) * cpatch%agb(ico)
                     end do fuelcohloop_d
                  end select
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Determine the correct threshold to ignite fires according to the   !
                  ! fire method.                                                           !
                  !------------------------------------------------------------------------!
                  select case (include_fire)
                  case (0)
                     !------ Set fire intensity to zero (fires should not happen). --------!
                     fire_intensity       = 0.
                     !---------------------------------------------------------------------!

                  case (1,2)
                     !---------------------------------------------------------------------!
                     !    Fires occur when average soil water goes below the soil moisture !
                     ! threshold.                                                          !
                     !---------------------------------------------------------------------!
                     if (csite%avg_monthly_gndwater(ipa)<cpoly%fire_wmass_threshold(isi))  &
                     then
                        fire_intensity      = fire_parameter
                        mean_fire_intensity = mean_fire_intensity                          &
                                            + fire_intensity * csite%area(ipa)
                     else
                        fire_intensity      = 0.0
                     end if
                     !---------------------------------------------------------------------!

                  case (3)
                     !---------------------------------------------------------------------!
                     !     Set fire intensity the same as method 2, except that we prevent !
                     ! fires until the site has any anthropogenic disturbance.  This will  !
                     ! change in the future to allow natural fires.                        !
                     !---------------------------------------------------------------------!
                     if (people_around .and.                                               &
                         csite%avg_monthly_gndwater(ipa)<cpoly%fire_wmass_threshold(isi) ) &
                     then
                        fire_intensity      = fire_parameter
                        mean_fire_intensity = mean_fire_intensity                          &
                                            + fire_intensity * csite%area(ipa)
                     else
                        !------ Set fire intensity to zero, so fires won't happen. --------!
                        fire_wmass_threshold = huge(1.)
                        fire_intensity       = 0.
                        !------------------------------------------------------------------!
                     end if
                     !---------------------------------------------------------------------!
                  end select
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !    If the soil is dry, then calculate patch contribution to the        !
                  ! ignition rate.                                                         !
                  !------------------------------------------------------------------------!
                  ignition_rate = ignition_rate + fire_intensity * fuel * csite%area(ipa)
                  !------------------------------------------------------------------------!


                  !----- Reset the ground water for next month. ---------------------------!
                  csite%avg_monthly_gndwater(ipa) = 0.
                  !------------------------------------------------------------------------!

               end do patchloop
               !---------------------------------------------------------------------------!



               !----- Calculate fire disturbance rate [1/year]. ---------------------------!
               cpoly%lambda_fire  (imo,isi) = ignition_rate
               if (mean_fire_intensity > 0.) then
                  cpoly%ignition_rate (isi) = ignition_rate / mean_fire_intensity
               else
                  cpoly%ignition_rate (isi) = 0.0
               end if
               !---------------------------------------------------------------------------!


               !----- Set the combusted fraction based on default values. -----------------!
               cpoly%avg_fire_f_bherb  (imo,isi) = f_combusted_fast_c
               cpoly%avg_fire_f_bwoody (imo,isi) = f_combusted_struct_c
               cpoly%avg_fire_f_fgc    (imo,isi) = f_combusted_fast_c
               cpoly%avg_fire_f_stgc   (imo,isi) = f_combusted_struct_c
               !---------------------------------------------------------------------------!
            end select
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine fire_frequency
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !       Sub-routine that integrates the Nesterov index.  Although not used by the model !
   ! unless INCLUDE_FIRE = 4, we still define  the index in case it is useful as a         !
   ! diagnostic variable.                                                                  !
   !---------------------------------------------------------------------------------------!
   subroutine integ_nesterov(cgrid)
      use ed_state_vars , only : edtype                 & ! structure
                               , polygontype            & ! structure
                               , sitetype               ! ! structure
      use consts_coms   , only : t00                    ! ! intent(in)
      use disturb_coms  , only : fh_pcpg_ni0            ! ! intent(in)
      implicit none
      !----- -Arguments. ------------------------------------------------------------------!
      type(edtype)     , target     :: cgrid
      !------ Local variables. ------------------------------------------------------------!
      type(polygontype), pointer    :: cpoly
      type(sitetype)   , pointer    :: csite
      integer                       :: ipy
      integer                       :: isi
      integer                       :: ipa
      real                          :: today_nesterov
      !------------------------------------------------------------------------------------!


      !------ Loop over polygons. ---------------------------------------------------------!
      poly_loop: do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)


         !----- Loop over sites. ----------------------------------------------------------!
         site_loop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)


            !------------------------------------------------------------------------------!
            !      Decide whether or not to reset the index.                               !
            !------------------------------------------------------------------------------!
            if (cpoly%today_pcpg(isi) <= fh_pcpg_ni0) then
               !---------------------------------------------------------------------------!
               !     Update Nesterov index.  Find the site average by looping through      !
               ! patches.                                                                  !
               !---------------------------------------------------------------------------!
               today_nesterov = 0.0
               patch_loop: do ipa=1,csite%npatches

                  !----- Add patch contribution. ------------------------------------------!
                  today_nesterov = today_nesterov                                          &
                                 + (csite%tdmax_can_temp(ipa) - t00                      ) &
                                 * (csite%tdmax_can_temp(ipa) - csite%today_can_tdew(ipa)) &
                                 * csite%area(ipa)
                  !------------------------------------------------------------------------!
               end do patch_loop
               !---------------------------------------------------------------------------!


               !----- Update polygon-level Nesterov Index. --------------------------------!
               cpoly%nesterov_index(isi) = cpoly%nesterov_index(isi) + today_nesterov
               !---------------------------------------------------------------------------!
            else
               !----- Rainy day, reset the index. -----------------------------------------!
               cpoly%nesterov_index(isi) = 0.
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!

         end do site_loop
         !---------------------------------------------------------------------------------!
      end do poly_loop
      !------------------------------------------------------------------------------------!


      return
   end subroutine integ_nesterov
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !       Sub-routine that integrates the fire disturbance rate when using (mostly)       !
   ! HESFIRE, as described in LP15/LP17.  The calculation of maximum rate of spread is     !
   ! based on R72 model revised by A18 and using the very dry fuel conditions described    !
   ! in SB05.  The model is complemented by a few additional equations to describe fire    !
   ! intensity, scorch height, and the impact on survivorship, which are all based on      !
   ! SPITFIRE (T10).                                                                       !
   !                                                                                       !
   ! References:                                                                           !
   !                                                                                       !
   ! Andrews PL. 2018. The Rothermel surface fire spread model and associated develop-     !
   !    ments: A compre- hensive explanation. Gen. Tech. Rep. RMRS-GTR-371, U.S.           !
   !    Department of Agriculture, Forest Service, Rocky Mountain Research Station, Fort   !
   !    Collins, CO, U.S.A. https://www.fs.usda.gov/treesearch/pubs/55928 (A18).           !
   !                                                                                       !
   ! Le Page Y, Morton D, Bond-Lamberty B, Pereira JMC , Hurtt G. 2015. HESFIRE: a global  !
   !    fire model to explore the role of anthropogenic and weather drivers. Biogeo-       !
   !    sciences, 12: 887-903. doi:10.5194/bg-12-887-2015 (LP15).                          !
   !                                                                                       !
   ! Le Page Y, Morton D, Hartin C, Bond-Lamberty B, Pereira JMC, Hurtt G , Asrar G. 2017. !
   !    Synergy between land use and climate change increases future fire risk in Amazon   !
   !    forests. Earth Syst. Dynam., 8: 1237-1246. doi:10.5194/esd-8-1237-2017 (LP17).     !
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
   ! Thonicke K, Spessa A, Prentice IC, Harrison SP, Dong L, Carmona-Moreno C. 2010. The   !
   !    influence of vegetation, fire spread and fire behaviour on biomass burning and     !
   !    trace gas emissions: results from a process-based model. Biogeosciences, 7:        !
   !    1991-2011. doi:10.5194/bg-7-1991-2010 (T10).                                       !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   subroutine integ_hesfire(cgrid,dtfire)
      use ed_state_vars , only : edtype                 & ! structure
                               , polygontype            & ! structure
                               , sitetype               & ! structure
                               , patchtype              ! ! structure
      use ed_misc_coms  , only : simtime                & ! structure
                               , current_time           ! ! intent(in)
      use disturb_coms  , only : fh_f0001               & ! intent(in)
                               , fh_f0010               & ! intent(in)
                               , fh_f0100               & ! intent(in)
                               , fh_f1000               & ! intent(in)
                               , fh_grid                & ! intent(in)
                               , fi_cg_ignp             & ! intent(in)
                               , fi_hdi_exp             & ! intent(in)
                               , fi_hdi_upr             & ! intent(in)
                               , fi_lu_ignd             & ! intent(in)
                               , fi_lu_exp              & ! intent(in)
                               , fi_lu_off              & ! intent(in)
                               , fi_lu_upr              & ! intent(in)
                               , fi_sf_maxage           & ! intent(in)
                               , fr_h                   & ! intent(in)
                               , fs_ba_frag             & ! intent(in)
                               , fs_gw_infty            & ! intent(in)
                               , fs_lbr_exp             & ! intent(in)
                               , fs_lbr_slp             & ! intent(in)
                               , fs_rhv_dti             & ! intent(in)
                               , fs_rhv_exp             & ! intent(in)
                               , fs_rhv_lwr             & ! intent(in)
                               , fs_smpot_dti           & ! intent(in)
                               , fs_smpot_exp           & ! intent(in)
                               , fs_smpot_lwr           & ! intent(in)
                               , fs_temp_dti            & ! intent(in)
                               , fs_temp_exp            & ! intent(in)
                               , fs_temp_lwr            & ! intent(in)
                               , ft_fint_dti            & ! intent(in)
                               , ft_fint_exp            & ! intent(in)
                               , ft_fint_lwr            & ! intent(in)
                               , ft_frag_exp            & ! intent(in)
                               , ft_lu_exp              & ! intent(in)
                               , ft_lu_upr              & ! intent(in)
                               , ft_hdi_exp             & ! intent(in)
                               , ft_hdi_upr             & ! intent(in)
                               , fuel_height_max        & ! intent(in)
                               , fx_a0001               & ! intent(in)
                               , fx_a0010               & ! intent(in)
                               , fx_a0100               & ! intent(in)
                               , fx_rmfac               & ! intent(in)
                               , fx_tlh_slope           & ! intent(in)
                               , include_fire           & ! intent(in)
                               , n_fst                  ! ! strcuture
      use pft_coms      , only : agf_bs                 & ! intent(in)
                               , C2B                    & ! intent(in)
                               , f_labile_leaf          & ! intent(in)
                               , f_labile_stem          & ! intent(in)
                               , is_grass               ! ! intent(in)
      use consts_coms   , only : almost_one             & ! intent(in)
                               , almost_zero            & ! intent(in)
                               , pio4                   & ! intent(in)
                               , tiny_num               & ! intent(in)
                               , lnexp_min              & ! intent(in)
                               , lnexp_max              ! ! intent(in)
      implicit none
      !----- -Arguments. ------------------------------------------------------------------!
      type(edtype)  , target     :: cgrid             ! Current grid              [    ---]
      real          , intent(in) :: dtfire            ! Fire full time step       [      s]
      !------ Local variables. ------------------------------------------------------------!
      type(polygontype), pointer :: cpoly             ! Current polygon           [    ---]
      type(sitetype)   , pointer :: csite             ! Current site              [    ---]
      type(patchtype)  , pointer :: cpatch            ! Current patch             [    ---]
      type(simtime)              :: hier              ! Yesterdays' date info.    [    ---]
      integer                    :: ico               ! Cohort index              [    ---]
      integer                    :: iflash            ! Idx to use from lightning [    ---]
      integer                    :: imonth            ! Month matrix index        [    ---]
      integer                    :: ipa               ! Patch index               [    ---]
      integer                    :: ipft              ! PFT index                 [    ---]
      integer                    :: ipy               ! Polygon index             [    ---]
      integer                    :: isei              ! Idx to use from SEI       [    ---]
      integer                    :: isi               ! Site index                [    ---]
      integer                    :: iwhen             ! Time loop index           [    ---]
      integer                    :: iyear             ! Year matrix index         [    ---]
      integer                    :: ndays             ! # days in month           [    ---]
      logical                    :: night_1st         ! Assume night for 1st step [    T|F]
      logical                    :: night             ! Nighttime step?           [    T|F]
      real, dimension(n_fst)     :: Mx_i              ! Fuel moist. extinct.      [  kg/kg]
      real                       :: anth_ign_rate     ! Anthropogenic ignt. rate  [ 1/m2/s]
      real                       :: apy_area          ! Grid area                 [     m2]
      real                       :: bherb             ! Cohort Herbaceous fuels   [ kgC/pl]
      real                       :: bherb_pat         ! Patch Herbaceous fuels    [ kgC/m2]
      real                       :: bherb_tot         ! Site Herbaceous fuels     [ kgC/m2]
      real                       :: bfuel_d0001_pat   ! Patch 1-hr dead fuels     [ kgC/m2]
      real                       :: bfuel_d0001_tot   ! Site 1-hr dead fuels      [ kgC/m2]
      real                       :: bfuel_d0010_pat   ! Patch 10-hr dead fuels    [ kgC/m2]
      real                       :: bfuel_d0010_tot   ! Site 10-hr dead fuels     [ kgC/m2]
      real                       :: bfuel_d0100_pat   ! Patch 100-hr dead fuels   [ kgC/m2]
      real                       :: bfuel_d0100_tot   ! Site 100-hr dead fuels    [ kgC/m2]
      real                       :: bfuel_d1000_pat   ! Patch 1000-hr dead fuels  [ kgC/m2]
      real                       :: bfuel_d1000_tot   ! Site 1000-hr dead fuels   [ kgC/m2]
      real                       :: bfuel_d0111_tot_i ! 1./(1+10+100-hr fuels)    [ m2/kgC]
      real                       :: burnt_area_step   ! Burnt area (time step)    [  m2/m2]
      real                       :: bwoody            ! Cohort living woody fuels [ kgC/pl]
      real                       :: bwoody_pat        ! Patch living woody fuels  [ kgC/m2]
      real                       :: bwoody_tot        ! Site living woody fuels   [ kgC/m2]
      real                       :: can_rhvn          ! Norm. CAS relative hum.   [    ---]
      real                       :: can_tempn         ! Normalised CAS temp.      [    ---]
      real                       :: fintn             ! Norm. fire intensity      [    ---]
      real                       :: fire_density_mid  ! Fire density (midstep)    [   1/m2]
      real                       :: fragn             ! Norm. fragmentation       [    ---]
      real                       :: fs_iarea          ! Individual fire area      [     m2]
      real                       :: fs_length         ! Length of fire ellipse    [    ---]
      real                       :: fs_rhv_fun        ! Rel. Hum. control funct.  [    ---]
      real                       :: fs_smpot_fun      ! Soil Potl. ctrl. funct.   [    ---]
      real                       :: fs_temp_fun       ! Temp. control function    [    ---]
      real                       :: fs_wind_fun       ! Wind control function     [    ---]
      real                       :: ft_frag_fun       ! Fragmentation limit. fac. [    ---]
      real                       :: ft_fuel_fun       ! Fuel limitation factor    [    ---]
      real                       :: ft_rhv_fun        ! Rel. Hum. control funct.  [    ---]
      real                       :: ft_smpot_fun      ! Soil Potl. ctrl. funct.   [    ---]
      real                       :: ft_supp_fun       ! Suppression limit. factor [    ---]
      real                       :: ft_suppress       ! Fire suppression factor   [    ---]
      real                       :: ft_temp_fun       ! Temp. control function    [    ---]
      real                       :: ft_decay_fun      ! Fire termination factor   [    ---]
      real                       :: ft_wind_fun       ! Wind control function     [    ---]
      real                       :: ft_wthr_fun       ! Weather limitation factor [    ---]
      real                       :: fx_b0001          ! Fuel consumption 1-h      [ kgC/m2]
      real                       :: fx_b0010          ! Fuel consumption 10-h     [ kgC/m2]
      real                       :: fx_b0100          ! Fuel consumption 100-h    [ kgC/m2]
      real                       :: fx_b1000          ! Fuel consumption 1000-h   [ kgC/m2]
      real                       :: fx_bherb          ! Fuel consumption herb     [ kgC/m2]
      real                       :: fx_bwoody         ! Fuel consumpton woody     [ kgC/m2]
      real                       :: fx_f_b0001        ! Rel. fuel consumpt. 1-h   [    ---]
      real                       :: fx_f_b0010        ! Rel. fuel consumpt. 10-h  [    ---]
      real                       :: fx_f_b0100        ! Rel. fuel consumpt. 100-h [    ---]
      real                       :: fx_f_b1000        ! Rel. fuel cons. 1000-h    [    ---]
      real                       :: fx_f_bherb        ! Rel. fuel consumpt. herb  [    ---]
      real                       :: fx_f_bwoody       ! Rel. fuel consumpt. woody [    ---]
      real                       :: fx_f_fgc          ! Rel. fuel cons. fast C    [    ---]
      real                       :: fx_f_stgc         ! Rel. fuel cons. struct C  [    ---]
      real                       :: fx_f_wn1000       ! Rel. f. cons. woody-1000h [    ---]
      real                       :: fx_intensity      ! Step fire intensity       [    W/m]
      real                       :: fx_tlethal        ! Step lethal heat duration [      s]
      real                       :: fx_wn1000         ! F. consumpt. woody-1000h  [ kgC/m2]
      real                       :: hb_ratio          ! Head:back ratio           [    ---]
      real                       :: hdin              ! Norm. human develop. idx  [    ---]
      real                       :: lb_ratio          ! Length:breadth ratio      [    ---]
      real                       :: lnexp             ! Aux. var. for safe exp    [    ---]
      real                       :: lu_area           ! LU area                   [    ---]
      real                       :: lu_effect         ! LU effect on ignition     [    ---]
      real                       :: lu_norm           ! Norm. LU area             [    ---]
      real                       :: moist_bfuel       ! Dead fuel moisture        [    ---]
      real                       :: moist_bherb       ! Herbaceous fuel moisture  [    ---]
      real                       :: moist_bwoody      ! Living woody fuel moist.  [    ---]
      real                       :: nat_ign_rate      ! Natural ignition rate     [ 1/m2/s]
      real                       :: ndaysi            ! 1/# days in a month       [  1/day]
      real                       :: rmoist_b0001      ! 1-hr dead rel. moisture   [    ---]
      real                       :: rmoist_b0010      ! 1-hr dead rel. moisture   [    ---]
      real                       :: rmoist_b0100      ! 1-hr dead rel. moisture   [    ---]
      real                       :: rmoist_b1000      ! 1-hr dead rel. moisture   [    ---]
      real                       :: rmoist_bherb      ! Herbaceous rel. moisture  [    ---]
      real                       :: rmoist_bwoody     ! Living woody rel. moist.  [    ---]
      real                       :: rosmax            ! Maximum rate of spread    [    m/s]
      real                       :: rosnow            ! Actual rate of spread     [    m/s]
      real                       :: sfc_smpotn        ! Norm. soil matrix potl.   [    ---]
      real                       :: total_ignition    ! Number of ignitions       [   1/m2]
      !------ External functions. ---------------------------------------------------------!
      real, external                :: solid_area     ! Solid-angle area          [     m2]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Run this subroutine only if using HESFIRE/SPITFIRE.                           !
      !------------------------------------------------------------------------------------!
      select case (include_fire)
      case (4)
         !------ Run the routine. ---------------------------------------------------------!
         continue
         !---------------------------------------------------------------------------------!
      case default
         !------ Skip the routine. --------------------------------------------------------!
         return
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the number of days of last day so we can normalise the integrated fire    !
      ! intensity and retrieve lightning and HDI information.                              !
      !------------------------------------------------------------------------------------!
      call yesterday_info(current_time,hier,ndays,ndaysi)
      imonth = hier%month
      iyear  = hier%year
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Loop over all polygons.                                                       !
      !------------------------------------------------------------------------------------!
      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)


         !---------------------------------------------------------------------------------!
         !     Find the absolute area of this polygon.                                     !
         !---------------------------------------------------------------------------------!
         apy_area = solid_area( cgrid%lon(ipy) - 0.5 * fh_grid                             &
                              , cgrid%lat(ipy) - 0.5 * fh_grid                             &
                              , cgrid%lon(ipy) + 0.5 * fh_grid                             &
                              , cgrid%lat(ipy) + 0.5 * fh_grid )
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Flag to decide whether to have night step followed by day step, or the      !
         ! other way round.  On average minimum temperature occurs around 6am LT, and      !
         ! maximum temperature occurs at around 2-3pm LT.  We use these as a rough         !
         ! guidance to decide the order.                                                   !
         !---------------------------------------------------------------------------------!
         night_1st =      ( cgrid%lon(ipy) >= -90. .and. cgrid%lon(ipy) < 127.5 )          &
                     .or.   cgrid%lon(ipy) >= 210.
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Loop over sites.                                                           !
         !---------------------------------------------------------------------------------!
         siteloop: do isi=1,cpoly%nsites
            csite => cpoly%site(isi)




            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
            !                                   IGNITION                                   !
            !------------------------------------------------------------------------------!
            !                                                                              !
            !      We find the ignition rates at site level, but acknowledging that the    !
            ! landscape may be fragmented (to be consistent with HESFIRE).  However, ED2   !
            ! does not have geographic information at the site level, therefore, for       !
            ! simplicity, we assume that fragmentation is homogeneous across sites. Within !
            ! sites, the following patches are considered areas that cannot sustain a      !
            ! fire:                                                                        !
            !                                                                              !
            !  1 -- Croplands                                                              !
            !  2 -- Areas that recently burnt                                              !
            !  3 -- Deserts (patches with bare soil)                                       !
            !  4 -- Fractions of patches covered in water (flooded) or snow.               !
            !                                                                              !
            !     The polygon area also accounts for the fraction of the polygon that      !
            ! cannot sustain a fire (oceans, glaciers, inland water, urban/built-up).      !
            ! These areas are not technically part of the site, so we only account for     !
            ! them in the empirical equations.                                             !
            !------------------------------------------------------------------------------!

            !------ Initialise local variables. -------------------------------------------!
            fragn    = 0.
            lu_area  = 0.
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !      Patch loop.                                                             !
            !------------------------------------------------------------------------------!
            patch_ignt_loop: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)

               !---------------------------------------------------------------------------!
               !     Check fragmentation cases.                                            !
               !---------------------------------------------------------------------------!
               if (csite%dist_type(ipa) == 8) then
                  !----- 1. Croplands, exclude the area. ----------------------------------!
                  fragn = fragn + csite%area(ipa)
                  !------------------------------------------------------------------------!
               elseif (csite%dist_type(ipa) == 4 .and. csite%age(ipa) <= fs_ba_frag) then
                  !----- 2. Recently burnt area. ------------------------------------------!
                  fragn = fragn + csite%area(ipa)
                  !------------------------------------------------------------------------!
               elseif (cpatch%ncohorts == 0) then
                  !----- 3. Deserts, exclude the area. ------------------------------------!
                  fragn = fragn + csite%area(ipa)
                  !------------------------------------------------------------------------!
               else
                  !----- 4. The patch may burn, but we exclude flooded/snowpack fraction. -!
                  fragn = fragn + csite%snowfac(ipa) * csite%area(ipa)
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Check land use area.                                                 !
               !---------------------------------------------------------------------------!
               select case (csite%dist_type(ipa))
               case (1,8)
                  !------ Pastures and  croplands. ----------------------------------------!
                  lu_area = lu_area + csite%area(ipa)
                  !------------------------------------------------------------------------!
               case (2,5,6,7)
                  !------------------------------------------------------------------------!
                  !     "Secondary" forests (forest plantations, abandoned lands, and      !
                  ! logged forests).  Add only when they are recently disturbed.           !
                  !------------------------------------------------------------------------!
                  if (csite%age(ipa) <= fi_sf_maxage) then
                     lu_area = lu_area + csite%area(ipa)
                  end if
                  !------------------------------------------------------------------------!
               end select
               !---------------------------------------------------------------------------!
            end do patch_ignt_loop
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !       Match the SEI time (currently SEI data are yearly).                    !
            !------------------------------------------------------------------------------!
            isei = cpoly%num_sei_times  (isi)
            find_sei_time: do iwhen=1,cpoly%num_sei_times(isi)
               if (iyear == cpoly%seitimes(iwhen,isi)%sei_year) then
                  isei = iwhen
                  exit find_sei_time
               end if
            end do find_sei_time
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Match the lightning time (current lightning data are climatological     !
            ! average by month, but this could be updated to account for IAV).             !
            !--------------------------------------------------------------------- --------!
            iflash     = cpoly%num_flash_times(isi)
            !----- Match year with lightning data. ----------------------------------------!
            find_flash_time: do iwhen=1,cpoly%num_sei_times(isi)
               if (imonth == cpoly%flashtimes(iwhen,isi)%flash_month) then
                  iflash = iwhen
                  exit find_flash_time
               end if
            end do find_flash_time
            !------------------------------------------------------------------------------!




            !------ Find natural ignition rate. -------------------------------------------!
            nat_ign_rate = cpoly%flashtimes(iflash,isi)%c2g * fi_cg_ignp * (1. - fragn)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !       Find the land use effect.  This is the analytical integral of LP15'    !
            ! Equation 3 (note that there is a "dLU" missing).  We account for the         !
            ! fraction of the polygon that cannot sustain a fire as LP15 also did.  We cap !
            ! the land use area to the maximum land use area that contributes to           !
            ! anthropogenic ignitions.                                                     !
            !------------------------------------------------------------------------------!
            lu_area       = min(fi_lu_upr,lu_area * cgrid%landfrac(ipy))
            lu_norm       = max(almost_zero,min(almost_one,1. - lu_area / fi_lu_upr))
            lu_effect     = fi_lu_off                                                      &
                          + ( fi_lu_upr + ( lu_area - fi_lu_upr ) * lu_norm ** fi_lu_exp ) &
                          / ( 1. + fi_lu_exp )
            lu_effect     = max( almost_zero, lu_effect )
            hdin          = max( almost_zero                                               &
                               , min(almost_one,cpoly%seitimes(isei,isi)%hdi/fi_hdi_upr) )
            anth_ign_rate = (1. - hdin ** fi_hdi_exp) * fi_lu_ignd * lu_effect
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !      Find total number of ignitions (1/m2).  This is applied twice inside    !
            ! the model time step, because the ignitions do not change during the same     !
            ! day.                                                                         !
            !------------------------------------------------------------------------------!
            cpoly%ignition_rate(isi) = nat_ign_rate + anth_ign_rate
            total_ignition           = cpoly%ignition_rate(isi) * 0.5 * dtfire
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Reset fire variables that will be average during the day.                !
            !------------------------------------------------------------------------------!
            cpoly%fire_spread   (isi) = 0.0
            cpoly%fire_intensity(isi) = 0.0
            cpoly%fire_tlethal  (isi) = 0.0
            cpoly%fire_f_bherb  (isi) = 0.0
            cpoly%fire_f_bwoody (isi) = 0.0
            cpoly%fire_f_fgc    (isi) = 0.0
            cpoly%fire_f_stgc   (isi) = 0.0
            !------------------------------------------------------------------------------!




            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
            !                                  FIRE DYNAMICS                               !
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !       Run the model twice, using either nighttime or daytime values to       !
            ! ensure to capture the very large differences between the time of the day and !
            ! the  chance of continuous fire spread or termination.                        !
            !------------------------------------------------------------------------------!
            timestep_loop: do iwhen=1,2
               !---------------------------------------------------------------------------!
               !       Select time step order.  This is just so most of the night/day time !
               ! steps are chronologically consistent with the longitude (as ED2 calls     !
               ! this routine at midnight UTC).                                            !
               !---------------------------------------------------------------------------!
               night = (iwhen == 1) .eqv. night_1st
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Initialise the functions that control fire spread and termination.   !
               !---------------------------------------------------------------------------!
               !----- Functions shared by fire spread. ------------------------------------!
               fs_temp_fun  = 0.
               fs_rhv_fun   = 0.
               fs_smpot_fun = 0.
               fs_wind_fun  = 0.
               fs_iarea     = 0.
               !----- Functions used for fire termination. --------------------------------!
               ft_supp_fun  = 0.
               ft_wthr_fun  = 0.
               !---------------------------------------------------------------------------!


               !------ Initialise fuel stocks. --------------------------------------------!
               bfuel_d0001_tot = 0.
               bfuel_d0010_tot = 0.
               bfuel_d0100_tot = 0.
               bfuel_d1000_tot = 0.
               bherb_tot       = 0.
               bwoody_tot      = 0.
               !---------------------------------------------------------------------------!


               !------ Initialise live woody moisture. ------------------------------------!
               moist_bherb    = 0.
               moist_bwoody   = 0.
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Ignite fires.                                                         !
               !---------------------------------------------------------------------------!
               cpoly%fire_density(isi) = cpoly%fire_density(isi) + total_ignition
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Loop through patches.                                                 !
               !---------------------------------------------------------------------------!
               fst_patch_loop: do ipa=1,csite%npatches
                  !------------------------------------------------------------------------!
                  !    For temperature and humidity, decide between minimum (night) or     !
                  ! maximum (day) temperature, and associated relative humidity (maximum   !
                  ! during the night, minimum during the day).                             !
                  !------------------------------------------------------------------------!
                  if (night) then
                     !------ Night time, use Tmin and RHmax. ------------------------------!
                     can_tempn = ( csite%tdmin_can_temp(ipa) - fs_temp_lwr ) * fs_temp_dti
                     can_rhvn  = ( csite%tdmax_can_rhv (ipa) - fs_rhv_lwr  ) * fs_rhv_dti
                     can_tempn = max( almost_zero, min( almost_one, can_tempn ) )
                     can_rhvn  = max( almost_zero, min( almost_one, can_rhvn  ) )
                     !---------------------------------------------------------------------!
                  else
                     !------ Day time, use Tmax and RHmin. --------------------------------!
                     can_tempn = ( csite%tdmax_can_temp(ipa) - fs_temp_lwr ) * fs_temp_dti
                     can_rhvn  = ( csite%tdmin_can_rhv (ipa) - fs_rhv_lwr  ) * fs_rhv_dti
                     can_tempn = max( almost_zero, min( almost_one, can_tempn ) )
                     can_rhvn  = max( almost_zero, min( almost_one, can_rhvn  ) )
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!


                  !----- Normalised soil potential. ---------------------------------------!
                  sfc_smpotn = ( csite%today_sfc_mstpot(ipa) - fs_smpot_lwr) * fs_smpot_dti
                  sfc_smpotn = max( almost_zero, min( almost_one, sfc_smpotn) )
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Integrate functions of temperature, relative humidity, and soil   !
                  ! potential.                                                             !
                  !------------------------------------------------------------------------!
                  fs_temp_fun  = fs_temp_fun                                               &
                               + ( 1. - can_tempn  ** fs_temp_exp ) * csite%area(ipa)
                  fs_rhv_fun   = fs_rhv_fun                                                &
                               + ( 1. - can_rhvn   ** fs_rhv_exp  ) * csite%area(ipa)
                  fs_smpot_fun = fs_smpot_fun                                              &
                               + ( 1. - sfc_smpotn ** fs_smpot_exp) * csite%area(ipa)
                  !------------------------------------------------------------------------!


                  !----- Find the wind influence function. --------------------------------!
                  lnexp       = max( lnexp_min                                             &
                                   , min( lnexp_max                                        &
                                        , fs_lbr_exp * csite%today_can_vels(ipa)) )
                  lb_ratio    = 1. + fs_lbr_slp * (1. - exp(lnexp))
                  hb_ratio    = ( lb_ratio + sqrt( lb_ratio * lb_ratio - 1. ) )            &
                              / ( lb_ratio - sqrt( lb_ratio * lb_ratio - 1. ) )
                  fs_wind_fun = fs_wind_fun                                                &
                              + 2. * lb_ratio / ( 1. + 1. / hb_ratio ) * fs_gw_infty       &
                              * csite%area(ipa)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !       Allocate fuels.                                                  !
                  !------------------------------------------------------------------------!
                  bfuel_d0001_pat = csite%fast_grnd_C(ipa)                                 &
                                  + fh_f0001 * csite%structural_grnd_C(ipa)
                  bfuel_d0010_pat = fh_f0010 * csite%structural_grnd_C(ipa)
                  bfuel_d0100_pat = fh_f0100 * csite%structural_grnd_C(ipa)
                  bfuel_d1000_pat = fh_f1000 * csite%structural_grnd_C(ipa)
                  bherb_pat       = 0.
                  bwoody_pat      = 0.
                  spread_cohort_loop: do ico=1,cpatch%ncohorts
                     ipft = cpatch%pft(ico)
                     if (is_grass(ipft) .or. cpatch%hite(ico) <= fuel_height_max) then
                        !------ Herbaceous fuel.  AG labile biomass + AG storage. ---------!
                        bherb  = f_labile_leaf(ipft) * cpatch%bleaf(ico)                   &
                               + f_labile_stem(ipft)                                       &
                               * ( cpatch%bsapwooda(ipft)                                  &
                                 + cpatch%bbarka   (ipft) + cpatch%bdeada(ico) )           &
                               + agf_bs(ipft) * cpatch%bstorage(ico)
                        !------------------------------------------------------------------!



                        !------ Woody living fuel.  AG lignified biomass. -----------------!
                        bwoody = (1. - f_labile_leaf(ipft)) * cpatch%bleaf(ico)            &
                               + (1. - f_labile_stem(ipft))                                &
                               * ( cpatch%bsapwooda(ipft)                                  &
                                 + cpatch%bbarka   (ipft) + cpatch%bdeada(ico) )
                        !------------------------------------------------------------------!


                        !------ Accumulate fuels to the patch level. ----------------------!
                        bherb_pat  = bherb_pat  + cpatch%nplant(ico) * bherb
                        bwoody_pat = bwoody_pat + cpatch%nplant(ico) * bwoody
                        !------------------------------------------------------------------!
                     end if
                  end do spread_cohort_loop
                  bfuel_d0001_tot = bfuel_d0001_tot + bfuel_d0001_pat * csite%area(ipa)
                  bfuel_d0010_tot = bfuel_d0010_tot + bfuel_d0010_pat * csite%area(ipa)
                  bfuel_d0100_tot = bfuel_d0100_tot + bfuel_d0100_pat * csite%area(ipa)
                  bfuel_d1000_tot = bfuel_d1000_tot + bfuel_d1000_pat * csite%area(ipa)
                  bherb_tot       = bherb_tot       + bherb_pat       * csite%area(ipa)
                  bwoody_tot      = bwoody_tot      + bwoody_pat      * csite%area(ipa)
                  !------------------------------------------------------------------------!


                  !------ Find the maximum rate of spread and the actual rate of spread. --!
                  rosmax = find_rosmax(bfuel_d0001_pat,bfuel_d0010_pat,bfuel_d0100_pat     &
                                      ,bfuel_d1000_pat,bherb_pat,bwoody_pat)
                  rosnow = rosmax * fs_rhv_fun * fs_temp_fun * fs_smpot_fun * fs_wind_fun
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !      Find the length of the major axis, following T10 and LP17.  This  !
                  ! can be applied to both new fires and existing fires.                   !
                  !------------------------------------------------------------------------!
                  fs_length = rosnow * dtfire
                  fs_iarea  = fs_iarea                                                     &
                            + pio4 * fs_length * fs_length / lb_ratio * csite%area(ipa)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !       Normalised weather control on termination.                       !
                  !------------------------------------------------------------------------!
                  if (    (can_tempn  <= almost_zero) .or. (can_rhvn >= almost_one)        &
                     .or. (sfc_smpotn >= almost_one )                               ) then
                     ft_wthr_fun = ft_wthr_fun + csite%area(ipa)
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Find herb and woody fuel wetness, using the wetness of the top    !
                  ! soil.                                                                  !
                  !------------------------------------------------------------------------!
                  moist_bherb  = moist_bherb  + csite%today_sfc_wetness(ipa)               &
                                              * bherb  * csite%area(ipa)
                  moist_bwoody = moist_bwoody + csite%today_sfc_wetness(ipa)               &
                                              * bwoody * csite%area(ipa)
                  !------------------------------------------------------------------------!


               end do fst_patch_loop
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !       Normalise moisture for herbs and living woody materials.            !
               !---------------------------------------------------------------------------!
               if (bherb_tot  > tiny_num) then
                  moist_bherb  = moist_bherb  / bherb_tot
               else
                  moist_bherb  = 1.0
               end if
               if (bwoody_tot > tiny_num) then
                  moist_bwoody = moist_bwoody / bwoody_tot
               else
                  moist_bwoody = 1.0
               end if
               !------ Apply correction factor for fuel moisture. -------------------------!
               moist_bherb  = ( (1.+fx_rmfac) * moist_bherb  - 1. ) / fx_rmfac
               moist_bwoody = ( (1.+fx_rmfac) * moist_bwoody - 1. ) / fx_rmfac
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !       Compute fuel moisture, based on SPITFIRE (T10).                     !
               !---------------------------------------------------------------------------!
               bfuel_d0111_tot_i = 1.                                                      &
                                 / (bfuel_d0100_tot + bfuel_d0010_tot + bfuel_d0001_tot)
               lnexp             = - ( fx_a0001 * bfuel_d0001_tot                          &
                                     + fx_a0010 * bfuel_d0010_tot                          &
                                     + fx_a0100 * bfuel_d0100_tot ) * bfuel_d0111_tot_i    &
                                 * cpoly%nesterov_index(isi)
               moist_bfuel       = exp(max(lnexp_min,min(lnexp_max,lnexp)))
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !    Find moisture of extinction for fuels and the relative moisture.       !
               !---------------------------------------------------------------------------!
               call find_mextinct(bfuel_d0001_tot,bfuel_d0010_tot,bfuel_d0100_tot          &
                                 ,bfuel_d1000_tot,bherb_tot,bwoody_tot                     &
                                 ,moist_bfuel,moist_bfuel,moist_bfuel,moist_bfuel          &
                                 ,moist_bherb,moist_bwoody,Mx_i)
               rmoist_b0001  = moist_bfuel  / Mx_i(1)
               rmoist_b0010  = moist_bfuel  / Mx_i(1)
               rmoist_b0100  = moist_bfuel  / Mx_i(1)
               rmoist_b1000  = moist_bfuel  / Mx_i(1)
               rmoist_bherb  = moist_bherb  / Mx_i(2)
               rmoist_bwoody = moist_bwoody / Mx_i(2)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !      Find fuel consumption.                                               !
               !---------------------------------------------------------------------------!
               !----- Find the fuel consumption factors for all fuel classes. -------------!
               call find_fx_factors(rmoist_bherb,rmoist_bwoody,rmoist_b0001,rmoist_b0010   &
                                   ,rmoist_b0100,rmoist_b1000,fx_f_bherb,fx_f_bwoody       &
                                   ,fx_f_wn1000,fx_f_b0001,fx_f_b0010,fx_f_b0100           &
                                   ,fx_f_b1000)
               !----- Find fuel consumption. ----------------------------------------------!
               fx_b0001  = fx_f_b0001  * bfuel_d0001_tot              * burnt_area_step
               fx_b0010  = fx_f_b0010  * bfuel_d0010_tot              * burnt_area_step
               fx_b0100  = fx_f_b0100  * bfuel_d0100_tot              * burnt_area_step
               fx_b1000  = fx_f_b1000  * bfuel_d1000_tot              * burnt_area_step
               fx_bherb  = fx_f_bherb  * bherb_tot                    * burnt_area_step
               fx_bwoody = fx_f_bwoody * bwoody_tot                   * burnt_area_step
               fx_wn1000 = fx_f_wn1000 * bwoody_tot * (1. - fh_f1000) * burnt_area_step
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Find combusted fraction for fast and structural pools.  This is a     !
               ! simplification that ought to be revisited at some point.  Ideally the     !
               ! above-ground structural pool should be split into the fuel classes, so    !
               ! different fractions can be burnt for each class, independently.           !
               !---------------------------------------------------------------------------!
               fx_f_fgc  = fx_f_b0001
               fx_f_stgc = fh_f0001 * fx_f_b0001 + fh_f0010 * fx_f_b0010                   &
                         + fh_f0100 * fx_f_b0100 + fh_f1000 * fx_f_b1000
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !       Find fire intensity of this step.                                   !
               !---------------------------------------------------------------------------!
               if (burnt_area_step > tiny_num) then
                  !----- Find fire intensity. ---------------------------------------------!
                  fx_intensity = fr_h * rosnow * C2B                                       &
                               * ( fx_b0001 + fx_b0010 + fx_b0100 + fx_bherb + fx_wn1000)  &
                               / burnt_area_step
                  if (fx_intensity < ft_fint_lwr) fx_intensity = 0.0
                  !------------------------------------------------------------------------!
               else
                  !----- No burnt area, set it to zero. -----------------------------------!
                  fx_intensity = 0.0
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !       Find duration of lethal bole heating.                               !
               !---------------------------------------------------------------------------!
               if (fx_intensity == 0.0) then
                  !----- No burning, set it to zero. --------------------------------------!
                  fx_tlethal = 0.0
                  !------------------------------------------------------------------------!
               else
                  !----- Find lethal duration. --------------------------------------------!
                  fx_tlethal = fx_tlh_slope * C2B                                          &
                             * ( bfuel_d0001_tot * (1. - fx_f_b0001 ) * (1. - fx_f_b0001 ) &
                               + bfuel_d0010_tot * (1. - fx_f_b0010 ) * (1. - fx_f_b0010 ) &
                               + bfuel_d0100_tot * (1. - fx_f_b0100 ) * (1. - fx_f_b0100 ) &
                               + bherb_tot       * (1. - fx_f_bherb ) * (1. - fx_f_bherb ) &
                               + bwoody_tot      * (1. - fh_f1000)                         &
                                                 * (1. - fx_f_wn1000) * (1. - fx_f_wn1000) )
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Compute the fuel limitation control on termination.  This replaces   !
               ! the precipitation term in LP15 with a fire intensity term that accounts   !
               ! for fuel loads and fuel moisture.  The first guess parameters are based   !
               ! on typical scorch heights for tropical broadleaf evergreen forests, using !
               ! T10 parameters.                                                           !
               !---------------------------------------------------------------------------!
               fintn          = ( fx_intensity - ft_fint_lwr ) * ft_fint_dti
               fintn          = max(almost_zero,min(almost_one,fintn))
               ft_fuel_fun    = 1. - fintn ** ft_fint_exp
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !      Compute fragmentation control on termination.  This is an empirical  !
               ! function and thus should include the fragmentation due to areas that      !
               ! cannot sustain vegetation. 
               !---------------------------------------------------------------------------!
               ft_frag_fun    = ( 1. - (1. - fragn) * cgrid%landfrac(ipy)) ** ft_frag_exp
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !       The termination functions of temperature, relative humidity, soil   !
               ! matric potential, and wind are taken as 1 - their spread counterparts.    !
               ! The rationale is that conditions that make it easy for fire to spread are !
               ! the same conditions that make it hard to suppress the fire.               !
               !---------------------------------------------------------------------------!
               ft_temp_fun  = 1. - fs_temp_fun
               ft_rhv_fun   = 1. - fs_rhv_fun
               ft_smpot_fun = 1. - fs_smpot_fun
               ft_wind_fun  = 1. - fs_wind_fun
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !       Compute the fire suppression function.                              !
               !---------------------------------------------------------------------------!
               !----- Fire suppressibility. -----------------------------------------------!
               ft_suppress = ft_rhv_fun * ft_smpot_fun * ft_temp_fun * ft_wind_fun         &
                           * ft_fuel_fun
               !----- Normalised land use . -----------------------------------------------!
               lu_norm       = max(almost_zero,min(almost_one,1. - lu_area / ft_lu_upr))
               !----- Normalised HDI . ----------------------------------------------------!
               hdin          = max( almost_zero                                            &
                                  , min( almost_one                                        &
                                       , cpoly%seitimes(isei,isi)%hdi / ft_hdi_upr ) )
               !----- Fire suppression function. ------------------------------------------!
               ft_supp_fun   = lu_norm ** ft_lu_exp * hdin ** ft_hdi_exp * ft_suppress
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !       Fire decay due to termination. \                                    !
               !---------------------------------------------------------------------------!
               ft_decay_fun = ( 1. - ft_fuel_fun ) * ( 1. - ft_frag_fun )                  &
                            * ( 1. - ft_supp_fun ) * ( 1. - ft_wthr_fun )
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Update the burnt area using an intermediate value to account for      !
               ! failed ignitions and early termination.                                   !
               !---------------------------------------------------------------------------!
               if (cgrid%landfrac(ipy) > tiny_num) then
                  !------ Increment burnt area until it is saturated. ---------------------!
                  fire_density_mid      = 0.5 * ( 1. + ft_decay_fun )                      &
                                        * cpoly%fire_density(isi)
                  burnt_area_step       = min( fragn - cpoly%burnt_area(isi)               &
                                             , fire_density_mid * fs_iarea                 &
                                             / ( apy_area * cgrid%landfrac(ipy) ) )
                  cpoly%burnt_area(isi) = cpoly%burnt_area(isi) + burnt_area_step
                  !------------------------------------------------------------------------!
               else
                  !------ No land to burn, set burnt area to zero... ----------------------!
                  burnt_area_step       = 0.0
                  cpoly%burnt_area(isi) = 0.0
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Terminate fires.                                                      !
               !---------------------------------------------------------------------------!
               cpoly%fire_density(isi) = cpoly%fire_density(isi) * ft_decay_fun
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Integrate the daily average fire spread.                              !
               !---------------------------------------------------------------------------!
               cpoly%fire_spread   (isi) = cpoly%fire_spread   (isi) + 0.5 * rosnow
               cpoly%fire_intensity(isi) = cpoly%fire_intensity(isi) + 0.5 * fx_intensity
               cpoly%fire_tlethal  (isi) = cpoly%fire_tlethal  (isi) + 0.5 * fx_tlethal
               cpoly%fire_f_bherb  (isi) = cpoly%fire_f_bherb  (isi) + 0.5 * fx_f_bherb
               cpoly%fire_f_bwoody (isi) = cpoly%fire_f_bwoody (isi) + 0.5 * fx_f_bwoody
               cpoly%fire_f_fgc    (isi) = cpoly%fire_f_fgc    (isi) + 0.5 * fx_f_fgc
               cpoly%fire_f_stgc   (isi) = cpoly%fire_f_stgc   (isi) + 0.5 * fx_f_fgc
               !---------------------------------------------------------------------------!
            end do timestep_loop
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !       Update the average fire intensity and lethal heating duration.         !
            !------------------------------------------------------------------------------!
            cpoly%avg_fire_intensity(imonth,isi) = cpoly%avg_fire_intensity(imonth,isi)    &
                                                 + cpoly%fire_intensity           (isi)    &
                                                 * ndaysi
            cpoly%avg_fire_tlethal  (imonth,isi) = cpoly%avg_fire_tlethal  (imonth,isi)    &
                                                 + cpoly%fire_tlethal             (isi)    &
                                                 * ndaysi
            cpoly%avg_fire_f_bherb  (imonth,isi) = cpoly%avg_fire_f_bherb  (imonth,isi)    &
                                                 + cpoly%fire_f_bherb             (isi)    &
                                                 * ndaysi
            cpoly%avg_fire_f_bwoody (imonth,isi) = cpoly%avg_fire_f_bwoody (imonth,isi)    &
                                                 + cpoly%fire_f_bwoody            (isi)    &
                                                 * ndaysi
            cpoly%avg_fire_f_fgc    (imonth,isi) = cpoly%avg_fire_f_fgc    (imonth,isi)    &
                                                 + cpoly%fire_f_fgc               (isi)    &
                                                 * ndaysi
            cpoly%avg_fire_f_stgc   (imonth,isi) = cpoly%avg_fire_f_stgc   (imonth,isi)    &
                                                 + cpoly%fire_f_stgc              (isi)    &
                                                 * ndaysi
            !------------------------------------------------------------------------------!



            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!



      return
   end subroutine integ_hesfire
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !       Sub-routine that resets the fire-related variables when using HESFIRE.  This    !
   ! subroutine is needed because HESFIRE timestep is daily (actually 12 hr).               
   !---------------------------------------------------------------------------------------!
   subroutine reset_hesfire(cgrid)
      use ed_state_vars , only : edtype                 & ! structure
                               , polygontype            & ! structure
                               , sitetype               ! ! structure
      use ed_misc_coms  , only : current_time           ! ! intent(in)
      use disturb_coms  , only : include_fire           ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)      , target     :: cgrid
      !----- Local variables --------------------------------------------------------------!
      type(polygontype) , pointer    :: cpoly
      type(sitetype)    , pointer    :: csite
      integer                        :: ipy
      integer                        :: isi
      integer                        :: ipa
      integer                        :: imo
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Run this subroutine only if using HESFIRE/SPITFIRE.                           !
      !------------------------------------------------------------------------------------!
      select case (include_fire)
      case (4);
         !------ Run the routine. ---------------------------------------------------------!
         continue
         !---------------------------------------------------------------------------------!
      case default
         !------ Skip the routine. --------------------------------------------------------!
         return
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!



      !----- Current month. ---------------------------------------------------------------!
      imo = current_time%month
      !------------------------------------------------------------------------------------!


      !----- Loop over polygons and sites. ------------------------------------------------!
      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         !---- Loop over all sites. -------------------------------------------------------!
         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)


            !---- Loop over all sites. ----------------------------------------------------!
            patchloop: do ipa=1,csite%npatches
               !----- Reset the ground water for next month. ------------------------------!
               csite%avg_monthly_gndwater(ipa) = 0.
               !---------------------------------------------------------------------------!
            end do patchloop
            !------------------------------------------------------------------------------!

            !----- Reset fire disturbance rates. ------------------------------------------!
            cpoly%lambda_fire       (imo,isi) = 0.0
            cpoly%avg_fire_intensity(imo,isi) = 0.0
            cpoly%avg_fire_tlethal  (imo,isi) = 0.0
            cpoly%avg_fire_f_bherb  (imo,isi) = 0.0
            cpoly%avg_fire_f_bwoody (imo,isi) = 0.0
            cpoly%avg_fire_f_fgc    (imo,isi) = 0.0
            cpoly%avg_fire_f_stgc   (imo,isi) = 0.0
            !------------------------------------------------------------------------------!


            !------ Resetting these variables just to play it safe. -----------------------!
            cpoly%fire_intensity         (isi) = 0.0
            cpoly%fire_spread            (isi) = 0.0
            cpoly%burnt_area             (isi) = 0.0
            cpoly%ignition_rate          (isi) = 0.0
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!


      return
   end subroutine reset_hesfire
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This function calculates the maximum spread rate, using the A18's revision of    !
   ! the R72 fire spread model, using "very dry" moisture conditions as defined by SB05.   !
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
   real function find_rosmax(bfuel_d0001,bfuel_d0010,bfuel_d0100,bfuel_d1000,bherb,bwoody)
      use disturb_coms, only : n_fst         & ! intent(in)
                             , n_fcl         & ! intent(in)
                             , n_sbmax       & ! intent(in)
                             , n_sbins       & ! intent(in)
                             , fr_sigma_ij   & ! intent(in)
                             , fr_moist_ij   & ! intent(in)
                             , fr_sgclss_ij  & ! intent(in)
                             , fr_ST         & ! intent(in)
                             , fr_dead_j     & ! intent(in)
                             , fr_rhop       & ! intent(in)
                             , fr_sigma_00   & ! intent(in)
                             , fr_hh_ij      & ! intent(in)
                             , fr_epsil_ee   & ! intent(in)
                             , fr_g_W_00     & ! intent(in)
                             , fr_Mxdead     & ! intent(in)
                             , fr_mxl_aa     & ! intent(in)
                             , fr_eta_m_aa   & ! intent(in)
                             , fr_Se_i       & ! intent(in)
                             , fr_eta_s_uu   & ! intent(in)
                             , fr_depth      & ! intent(in)
                             , fr_beta_op_uu & ! intent(in)
                             , fr_gamma_xx   & ! intent(in)
                             , fr_AA_uu      & ! intent(in)
                             , fr_xi_xx      & ! intent(in)
                             , fr_Umax_uu    & ! intent(in)
                             , fr_BB_uu      & ! intent(in)
                             , fr_CC_xx      & ! intent(in)
                             , fr_EE_ee      & ! intent(in)
                             , fr_Qig_aa     ! ! intent(in)
      use pft_coms    , only : C2B           ! ! intent(in)
      use consts_coms , only : tiny_num      & ! intent(in)
                             , lnexp_min     & ! intent(in)
                             , lnexp_max     ! ! intent(in)
      implicit none

      !------ Arguments. ------------------------------------------------------------------!
      real, intent(in)               :: bfuel_d0001  !  1-hr dead fuel load     [   kgC/m2]
      real, intent(in)               :: bfuel_d0010  !  10-hr dead fuel load    [   kgC/m2]
      real, intent(in)               :: bfuel_d0100  !  100-hr dead fuel load   [   kgC/m2]
      real, intent(in)               :: bfuel_d1000  !  1000-hr dead fuel load  [   kgC/m2]
      real, intent(in)               :: bherb        !  Herbaceous fuel load    [   kgC/m2]
      real, intent(in)               :: bwoody       !  Living Woody fuel load  [   kgC/m2]
      !----- Local variables (by fuel class and status). ----------------------------------!
      real, dimension(n_fst,n_fcl)   :: wood_ij      ! Fuel load                [   kgB/m2]
      real, dimension(n_fst,n_fcl)   :: fai_ij       ! Fuel area index          [  m2_f/m2]
      real, dimension(n_fst,n_fcl)   :: fwgt_ij      ! Weighting factor         [       --]
      real, dimension(n_fst,n_fcl)   :: wnod_ij      ! Net fuel load            [   kgB/m2]
      real, dimension(n_fst,n_fcl)   :: epsil_ij     ! Effective heating number [       --]
      real, dimension(n_fst,n_fcl)   :: rmoist_ij    ! Relative moisture        [       --]
      real, dimension(n_fst,n_fcl)   :: Qig_ij       ! Heat of pre-ignition     [     J/kg]
      !----- Local variables (by fuel size and status). -----------------------------------!
      real, dimension(n_fst,n_sbmax) :: gwgt_ik      ! Weighting factor (net)   [       --]
      real, dimension(n_fst,n_sbmax) :: gwnod_ik     ! Net fuel load (size bin) [       --]
      !----- Local variables (by fuel status, aggregated across fuel classes). ------------!
      real, dimension(n_fst)         :: wood_i       ! Fuel load                [   kgB/m2]
      real, dimension(n_fst)         :: fai_i        ! Fuel area index          [  m2_f/m2]
      real, dimension(n_fst)         :: fwgt_i       ! Weighting factor         [       --]
      real, dimension(n_fst)         :: wnod_i       ! Net fuel load            [   kgB/m2]
      real, dimension(n_fst)         :: sigma_i      ! Effective SAV            [      1/m]
      real, dimension(n_fst)         :: hh_i         ! Heat content             [      1/m]
      real, dimension(n_fst)         :: Mx_i         ! Fuel moist. of extinct.  [    kg/kg]
      real, dimension(n_fst)         :: moist_i      ! Fuel moisture            [    kg/kg]
      real, dimension(n_fst)         :: rmoist_i     ! Relative moisture        [       --]
      real, dimension(n_fst)         :: eta_m_i      ! Moisture dampening coef. [       --]
      real, dimension(n_fst)         :: eta_s_i      ! Mineral dampening coef.  [       --]
      real, dimension(n_fst)         :: epsilQig_i   ! Effective heat pre-ign.  [     J/kg]
      !----- Local variables (aggregated across fuel classes and statuses). ---------------!
      real                           :: g_wood       ! Fuel load                [   kgB/m2]
      real                           :: g_fai        ! Fuel area index          [  m2_f/m2]
      real                           :: g_wnod       ! Net fuel load            [   kgB/m2]
      real                           :: g_sigma      ! Effective SAV            [      1/m]
      real                           :: g_hh         ! Heat content             [      1/m]
      real                           :: g_W          ! Dead-to-live load ratio  [    kg/kg]
      real                           :: g_rhob       ! Effective bulk density   [    kg/m3]
      real                           :: g_beta       ! Mean packing ratio       [       --]
      real                           :: g_beta_op    ! Optimum packing ratio    [       --]
      real                           :: g_r_beta     ! Relative packing ratio   [       --]
      real                           :: g_gamma_max  ! Maximum reaction veloc.  [      1/s]
      real                           :: g_gamma      ! Optimum reaction veloc.  [      1/s]
      real                           :: g_AA         ! Aux. variable            [       --]
      real                           :: g_BB         ! Aux. variable            [       --]
      real                           :: g_CC         ! Aux. variable            [       --]
      real                           :: g_EE         ! Aux. variable            [       --]
      real                           :: g_xi         ! Propagating flux ratio   [       --]
      real                           :: g_Ir         ! Reaction intensity       [     W/m2]
      real                           :: g_Umax       ! Max. wind (corrected)    [      m/s]
      real                           :: g_HeatSink   ! Heat sink                [     J/m3]
      real                           :: g_psiw       ! R72's wind factor        [       --]
      real                           :: g_ROSmax_w0  ! Max. ROS (no wind)       [      m/s]
      !----- Additional local variables. --------------------------------------------------!
      integer                        :: i            ! Counter                  [       --]
      integer                        :: j            ! Counter                  [       --]
      integer                        :: k            ! Counter                  [       --]
      real                           :: mf_dead      ! "Fine" dead fuel moist.  [    kg/kg]
      real                           :: lnexp        ! Temp var to avoid FPE    [       --]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !       Define fuel load.                                                            !
      !------------------------------------------------------------------------------------!
      wood_ij(1,:) = (/ bfuel_d0001, bfuel_d0010, bfuel_d0100, bfuel_d1000, bherb, bwoody /)
      wood_ij(1,:) = wood_ij(1,:) * C2B
      wood_ij(2,:) = (1. - fr_dead_j(:)) * wood_ij(1,:)
      wood_ij(1,:) = fr_dead_j(:) * wood_ij(1,:)
      do i=1,n_fst
         wood_i(i) = sum(wood_ij(i,:))
      end do
      g_wood = sum(wood_i)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find area and weighting factors.                                               !
      !------------------------------------------------------------------------------------!
      !----- Fuel area index. -------------------------------------------------------------!
      fai_ij(:,:) = fr_sigma_ij(:,:) * wood_ij(:,:) / fr_rhop
      do i=1,n_fst
         fai_i(i) = sum(fai_ij(i,:))
      end do
      g_fai = sum(fai_i)
      !----- Element-wise weighting factor. -----------------------------------------------!
      do i=1,n_fst
         !------ Make sure data are bounded. ----------------------------------------------!
         if (fai_i(i) > tiny_num) then
            fwgt_ij(i,:) = fai_ij(i,:) / fai_i(i)
         elseif (i == 1) then
            fwgt_ij(i,:) = fr_dead_j(:) / sum(fr_dead_j(:))
         else
            fwgt_ij(i,:) = (1. - fr_dead_j(:)) / sum(1. - fr_dead_j(:))
         end if
         !---------------------------------------------------------------------------------!
      end do
      !----- Weighting factor for fuel status. --------------------------------------------!
      if (g_fai == 0.) then
         fwgt_i(:) = 1. / real(n_fst)
      else
         fwgt_i(:) = fai_i(:) / g_fai
      end if
      !------------------------------------------------------------------------------------!



      !------ Calculate the net fuel load. ------------------------------------------------!
      wnod_ij(:,:) = wood_ij(:,:) * (1. - fr_ST)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the weighting factors for the net fuel load.                              !
      !------------------------------------------------------------------------------------!
      gwgt_ik (:,:) = 0.
      gwnod_ik(:,:) = 0.
      do i=1,n_fst
         !---------------------------------------------------------------------------------!
         !      Loop through bins, populate the bins according to their sizes.  Following  !
         ! A18, we ignore fuels at the lowest SAV size class.                              !
         !---------------------------------------------------------------------------------!
         do k=2,n_sbins
            gwgt_ik(i,k) = sum(fwgt_ij(i,:),mask=fr_sgclss_ij(i,:) == k)
         end do
         !----- Net fuel loads, add all elements for completeness. ------------------------!
         do k=1,n_sbins
            gwnod_ik(i,k) = sum(wnod_ij(i,:),mask=fr_sgclss_ij(i,:) == k)
         end do
         !---------------------------------------------------------------------------------!

         !------ Standardise weights so they add up to 1. ---------------------------------!
         gwgt_ik(i,:) = gwgt_ik(i,:) / sum(gwgt_ik(i,:))
         !---------------------------------------------------------------------------------!
      end do
      !------ Aggregate net fuel load by size class. --------------------------------------!
      do i=1,n_fst
         wnod_i(i) = sum(gwnod_ik(i,:)*gwgt_ik(i,:))
      end do
      g_wnod = sum(wnod_i)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the effective SAV.                                                        !
      !------------------------------------------------------------------------------------!
      !------ Effective SAV by fuel status. -----------------------------------------------!
      do i=1,n_fst
         !------ Ensure that fuel area is not zero (in case it is, assume dummy values). --!
         if (fai_i(i) > tiny_num) then
            !----- Weighted mean. ---------------------------------------------------------!
            sigma_i(i) = sum(fr_sigma_ij(i,:)*fwgt_ij(i,:))
            !------------------------------------------------------------------------------!
         else
            !----- Dummy value. -----------------------------------------------------------!
            sigma_i(i) = fr_sigma_00
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end do
      !------ Effective SAV. --------------------------------------------------------------!
      g_sigma = sum(sigma_i*fwgt_i) / sum(fwgt_i)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the heat content.                                                         !
      !------------------------------------------------------------------------------------!
      !------ Heat content by fuel status. ------------------------------------------------!
      do i=1,n_fst
         hh_i(i) = sum(fr_hh_ij(i,:)*fwgt_ij(i,:))
      end do
      !------ Heat content. ---------------------------------------------------------------!
      g_hh = sum(hh_i*fwgt_i) / sum(fwgt_i)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the relative moisture.                                                    !
      !------------------------------------------------------------------------------------!
      !----- Find the effective heating number. -------------------------------------------!
      do i=1,n_fst
         do j=1,n_fcl
            lnexp         = max(lnexp_min,min(lnexp_max,fr_epsil_ee(i)/fr_sigma_ij(i,j)))
            epsil_ij(i,j) = exp(lnexp)
         end do
      end do
      !----- Fine dead fuel moisture. -----------------------------------------------------!
      mf_dead = sum(fr_moist_ij(1,:)*wnod_ij(1,:)*epsil_ij(1,:))                           &
              / sum(wnod_ij(1,:)*epsil_ij(1,:))
      !----- Find the dead-to-live load ratio (or use dummy in case live load is 0). ------!
      if (wnod_i(2) > tiny_num) then
         g_W = sum(wnod_ij(1,:)*epsil_ij(1,:)) / sum(wnod_ij(2,:)*epsil_ij(2,:))
      else
         g_W = fr_g_W_00
      end if
      !----- Live fuel moisture of extinction. --------------------------------------------!
      Mx_i(1) = fr_Mxdead
      Mx_i(2) = min( fr_Mxdead                                                             &
                   , fr_mxl_aa(1) + fr_mxl_aa(2) * g_W * (1. - mf_dead/fr_Mxdead) )
      !----- Fuel moisture. ---------------------------------------------------------------!
      do i=1,n_fst
         moist_i(i) = sum(fr_moist_ij(i,:)*fwgt_ij(i,:))
      end do
      !----- Relative moisture. -----------------------------------------------------------!
      do i=1,n_fst
         !----- Moisture by status and class. ---------------------------------------------!
         do j=1,n_fcl
            rmoist_ij(i,j) = min(1.,fr_moist_ij(i,j)/Mx_i(i))
         end do
         !---------------------------------------------------------------------------------!

         !----- Aggregated moisture by status. --------------------------------------------!
         rmoist_i(i) = min(1.,moist_i(i)/Mx_i(i))
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!


      !----- Moisture dampening coefficient. ----------------------------------------------!
      eta_m_i(:) = fr_eta_m_aa(1) + rmoist_i(:)                                            &
                                  * ( fr_eta_m_aa(2) + rmoist_i(:)                         &
                                                     * ( fr_eta_m_aa(3) + fr_eta_m_aa(4)   &
                                                                        * rmoist_i(:)    ) )
      !------------------------------------------------------------------------------------!


      !----- Mineral dampening coefficient. -----------------------------------------------!
      do i=1,n_fst
         eta_s_i(i) = min(1.,fr_eta_s_uu(1) * fr_Se_i(i) ** fr_eta_s_uu(2))
      end do
      !------------------------------------------------------------------------------------!


      !----- Effective bulk density. ------------------------------------------------------!
      g_rhob = g_wood / fr_depth
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find packing ratios.                                                           !
      !------------------------------------------------------------------------------------!
      !----- Mean packing ratio. ----------------------------------------------------------!
      g_beta    = g_rhob / fr_rhop
      !----- Optimum packing ratio. -------------------------------------------------------!
      g_beta_op = fr_beta_op_uu(1) * g_sigma ** fr_beta_op_uu(2)
      !----- Relative packing ratio. ------------------------------------------------------!
      g_r_beta  = g_beta / g_beta_op
      !------------------------------------------------------------------------------------!


      !---- Maximum reaction velocity  [1/s]. ---------------------------------------------!
      g_gamma_max = fr_gamma_xx(1)                                                         &
                  / ( fr_gamma_xx(2) + fr_gamma_xx(3) * g_sigma ** fr_gamma_xx(4) )
      !------------------------------------------------------------------------------------!


      !---- Optimum reaction velocity  [1/s]. ---------------------------------------------!
      g_AA    = fr_AA_uu(1) * g_sigma ** fr_AA_uu(2)
      g_gamma = g_gamma_max * g_r_beta ** g_AA * exp(g_AA * (1. - g_r_beta))
      !------------------------------------------------------------------------------------!


      !---- Propagating flux ratio. -------------------------------------------------------!
      g_xi = exp( ( fr_xi_xx(1) + fr_xi_xx(2) * g_sigma **fr_xi_xx(3) )                    &
                * ( g_beta + fr_xi_xx(4) )                              )                  &
           / ( fr_xi_xx(5) + fr_xi_xx(6) * g_sigma )
      !------------------------------------------------------------------------------------!


      !---- Reaction intensity [W/m2]. ----------------------------------------------------!
      g_Ir = g_gamma * sum(wnod_i(:) * hh_i(:) * eta_m_i(:) * eta_s_i(:))
      !------------------------------------------------------------------------------------!


      !---- Maximum wind (corrected for saturation) [ m/s]. -------------------------------!
      g_Umax = fr_Umax_uu(1) * g_Ir ** fr_Umax_uu(2)
      !------------------------------------------------------------------------------------!


      !---- Ancillary parameters that depend on surface-area-to-volume-ratio. -------------!
      g_BB  = fr_BB_uu(1) * g_sigma ** fr_BB_uu(2)
      lnexp = max( lnexp_min, min( lnexp_max, fr_CC_xx(3) * g_sigma ** fr_CC_xx(4) ) )
      g_CC  = fr_CC_xx(1) * fr_CC_xx(2) ** g_BB * exp(lnexp)
      lnexp = max( lnexp_min, min( lnexp_max, fr_EE_ee(2) * g_sigma ) )
      g_EE  = fr_EE_ee(1) * exp(lnexp)
      !------------------------------------------------------------------------------------!


      !----- Heat of pre-ignition [ J/kg]. ------------------------------------------------!
      Qig_ij(:,:) = fr_Qig_aa(1) + fr_Qig_aa(2) * fr_moist_ij(:,:)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the heat sink.                                                            !
      !------------------------------------------------------------------------------------!
      !----- Effective heat of pre-ignition [J/kg]. ---------------------------------------!
      do i=1,n_fst
         epsilQig_i(i) = sum(epsil_ij(i,:)*Qig_ij(i,:)*fwgt_ij(i,:))
      end do
      !----- Heat sink [J/m3]. ------------------------------------------------------------!
      g_HeatSink = g_rhob * sum(epsilQig_i(:)*fwgt_i(:))
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the wind factor based on R72, assuming the maximum possible wind.         !
      !------------------------------------------------------------------------------------!
      g_psiw = g_CC * g_Umax ** g_BB / g_beta ** g_EE
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the maximum rate of spread.                                               !
      !------------------------------------------------------------------------------------!
      !----- Maximum rate of spread in the absence of wind  [m/s]. ------------------------!
      g_ROSmax_w0 = g_Ir * g_xi / g_HeatSink
      !----- Maximum rate of spread                         [m/s]. ------------------------!
      find_rosmax = g_ROSmax_w0 * (1. + g_psiw)
      !------------------------------------------------------------------------------------!

      return
   end function find_rosmax
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This function calculates the moisture of extinction, using the A18's revision of !
   ! the R72 fire spread model, for any given fuel content and moisture.                   !
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
   subroutine find_mextinct(bfuel_d0001,bfuel_d0010,bfuel_d0100,bfuel_d1000,bherb,bwoody   &
                           ,moist_b0001,moist_b0010,moist_b0100,moist_b1000,moist_bherb    &
                           ,moist_bwoody,Mx_i)
      use disturb_coms, only : n_fst         & ! intent(in)
                             , n_fcl         & ! intent(in)
                             , n_sbmax       & ! intent(in)
                             , n_sbins       & ! intent(in)
                             , fr_sigma_ij   & ! intent(in)
                             , fr_sgclss_ij  & ! intent(in)
                             , fr_ST         & ! intent(in)
                             , fr_dead_j     & ! intent(in)
                             , fr_rhop       & ! intent(in)
                             , fr_epsil_ee   & ! intent(in)
                             , fr_g_W_00     & ! intent(in)
                             , fr_Mxdead     & ! intent(in)
                             , fr_mxl_aa     ! ! intent(in)
      use pft_coms    , only : C2B           ! ! intent(in)
      use consts_coms , only : tiny_num      & ! intent(in)
                             , lnexp_min     & ! intent(in)
                             , lnexp_max     ! ! intent(in)
      implicit none

      !------ Arguments. ------------------------------------------------------------------!
      real, intent(in)                    :: bfuel_d0001  ! 1-hr dead fuel       [  kgC/m2]
      real, intent(in)                    :: bfuel_d0010  ! 10-hr dead fuel      [  kgC/m2]
      real, intent(in)                    :: bfuel_d0100  ! 100-hr dead fuel     [  kgC/m2]
      real, intent(in)                    :: bfuel_d1000  ! 1000-hr dead fuel    [  kgC/m2]
      real, intent(in)                    :: bherb        ! Herb. fuel           [  kgC/m2]
      real, intent(in)                    :: bwoody       ! Woody fuel           [  kgC/m2]
      real, intent(in)                    :: moist_b0001  ! 1-hr fuel moist.     [     ---]
      real, intent(in)                    :: moist_b0010  ! 10-hr fuel moist.    [     ---]
      real, intent(in)                    :: moist_b0100  ! 100-hr fuel moist.   [     ---]
      real, intent(in)                    :: moist_b1000  ! 1000-hr fuel moist.  [     ---]
      real, intent(in)                    :: moist_bherb  ! Herb. fuel moist.    [     ---]
      real, intent(in)                    :: moist_bwoody ! Woody fuel moisture  [     ---]
      real, dimension(n_fst), intent(out) :: Mx_i         ! Fuel moist. extinct. [   kg/kg]
      !----- Local variables (by fuel class and status). ----------------------------------!
      real, dimension(n_fst,n_fcl)   :: wood_ij      ! Fuel load                [   kgB/m2]
      real, dimension(n_fst,n_fcl)   :: moist_ij     ! Fuel moisture            [      ---]
      real, dimension(n_fst,n_fcl)   :: fai_ij       ! Fuel area index          [  m2_f/m2]
      real, dimension(n_fst,n_fcl)   :: fwgt_ij      ! Weighting factor         [       --]
      real, dimension(n_fst,n_fcl)   :: wnod_ij      ! Net fuel load            [   kgB/m2]
      real, dimension(n_fst,n_fcl)   :: epsil_ij     ! Effective heating number [       --]
      !----- Local variables (by fuel size and status). -----------------------------------!
      real, dimension(n_fst,n_sbmax) :: gwgt_ik      ! Weighting factor (net)   [       --]
      real, dimension(n_fst,n_sbmax) :: gwnod_ik     ! Net fuel load (size bin) [       --]
      !----- Local variables (by fuel status, aggregated across fuel classes). ------------!
      real, dimension(n_fst)         :: fai_i        ! Fuel area index          [  m2_f/m2]
      real, dimension(n_fst)         :: fwgt_i       ! Weighting factor         [       --]
      real, dimension(n_fst)         :: wnod_i       ! Net fuel load            [   kgB/m2]
      !----- Local variables (aggregated across fuel classes and statuses). ---------------!
      real                           :: g_fai        ! Fuel area index          [  m2_f/m2]
      real                           :: g_W          ! Dead-to-live load ratio  [    kg/kg]
      !----- Additional local variables. --------------------------------------------------!
      integer                        :: i            ! Counter                  [       --]
      integer                        :: j            ! Counter                  [       --]
      integer                        :: k            ! Counter                  [       --]
      real                           :: mf_dead      ! "Fine" dead fuel moist.  [    kg/kg]
      real                           :: lnexp        ! Temp var to avoid FPE    [       --]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !       Define fuel load.                                                            !
      !------------------------------------------------------------------------------------!
      wood_ij(1,:) = (/ bfuel_d0001, bfuel_d0010, bfuel_d0100, bfuel_d1000, bherb, bwoody /)
      wood_ij(1,:) = wood_ij(1,:) * C2B
      wood_ij(2,:) = (1. - fr_dead_j(:)) * wood_ij(1,:)
      wood_ij(1,:) = fr_dead_j(:) * wood_ij(1,:)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !       Define fuel moisture.                                                        !
      !------------------------------------------------------------------------------------!
      moist_ij(1,:) = (/  moist_b0001,  moist_b0010,  moist_b0100,  moist_b1000            &
                       ,  moist_bherb, moist_bwoody /)
      moist_ij(2,:) = moist_ij(1,:)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find area and weighting factors.                                               !
      !------------------------------------------------------------------------------------!
      !----- Fuel area index. -------------------------------------------------------------!
      fai_ij(:,:) = fr_sigma_ij(:,:) * wood_ij(:,:) / fr_rhop
      do i=1,n_fst
         fai_i(i) = sum(fai_ij(i,:))
      end do
      g_fai = sum(fai_i)
      !----- Element-wise weighting factor. -----------------------------------------------!
      do i=1,n_fst
         !------ Make sure data are bounded. ----------------------------------------------!
         if (fai_i(i) > tiny_num) then
            fwgt_ij(i,:) = fai_ij(i,:) / fai_i(i)
         elseif (i == 1) then
            fwgt_ij(i,:) = fr_dead_j(:) / sum(fr_dead_j(:))
         else
            fwgt_ij(i,:) = (1. - fr_dead_j(:)) / sum(1. - fr_dead_j(:))
         end if
         !---------------------------------------------------------------------------------!
      end do
      !----- Weighting factor for fuel status. --------------------------------------------!
      if (g_fai == 0.) then
         fwgt_i(:) = 1. / real(n_fst)
      else
         fwgt_i(:) = fai_i(:) / g_fai
      end if
      !------------------------------------------------------------------------------------!



      !------ Calculate the net fuel load. ------------------------------------------------!
      wnod_ij(:,:) = wood_ij(:,:) * (1. - fr_ST)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the weighting factors for the net fuel load.                              !
      !------------------------------------------------------------------------------------!
      gwgt_ik (:,:) = 0.
      gwnod_ik(:,:) = 0.
      do i=1,n_fst
         !---------------------------------------------------------------------------------!
         !      Loop through bins, populate the bins according to their sizes.  Following  !
         ! A18, we ignore fuels at the lowest SAV size class.                              !
         !---------------------------------------------------------------------------------!
         do k=2,n_sbins
            gwgt_ik(i,k) = sum(fwgt_ij(i,:),mask=fr_sgclss_ij(i,:) == k)
         end do
         !----- Net fuel loads, add all elements for completeness. ------------------------!
         do k=1,n_sbins
            gwnod_ik(i,k) = sum(wnod_ij(i,:),mask=fr_sgclss_ij(i,:) == k)
         end do
         !---------------------------------------------------------------------------------!

         !------ Standardise weights so they add up to 1. ---------------------------------!
         gwgt_ik(i,:) = gwgt_ik(i,:) / sum(gwgt_ik(i,:))
         !---------------------------------------------------------------------------------!
      end do
      !------ Aggregate net fuel load by size class. --------------------------------------!
      do i=1,n_fst
         wnod_i(i) = sum(gwnod_ik(i,:)*gwgt_ik(i,:))
      end do
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the relative moisture.                                                    !
      !------------------------------------------------------------------------------------!
      !----- Find the effective heating number. -------------------------------------------!
      do i=1,n_fst
         do j=1,n_fcl
            lnexp         = max(lnexp_min,min(lnexp_max,fr_epsil_ee(i)/fr_sigma_ij(i,j)))
            epsil_ij(i,j) = exp(lnexp)
         end do
      end do
      !----- Fine dead fuel moisture. -----------------------------------------------------!
      mf_dead = sum(moist_ij(1,:)*wnod_ij(1,:)*epsil_ij(1,:))                              &
              / sum(wnod_ij(1,:)*epsil_ij(1,:))
      !----- Find the dead-to-live load ratio (or use dummy in case live load is 0). ------!
      if (wnod_i(2) > tiny_num) then
         g_W = sum(wnod_ij(1,:)*epsil_ij(1,:)) / sum(wnod_ij(2,:)*epsil_ij(2,:))
      else
         g_W = fr_g_W_00
      end if
      !----- Live fuel moisture of extinction. --------------------------------------------!
      Mx_i(1) = fr_Mxdead
      Mx_i(2) = min( fr_Mxdead                                                             &
                   , fr_mxl_aa(1) + fr_mxl_aa(2) * g_W * (1. - mf_dead/fr_Mxdead) )
      !------------------------------------------------------------------------------------!

      return
   end subroutine find_mextinct
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This function calculates the fuel consumption factor from the relative moisture  !
   ! of each flux class.  This is based on appendix B of T10, with the following           !
   ! modifications:                                                                        !
   !                                                                                       !
   ! 1. We keep dead 1-hr fuels separated from herbaceous fuels (i.e., same equations,     !
   !    but different moistures), and we solve the woody fuel moisture based on the        !
   !    fractions attributed to 1-hr, 10-hr, 100-hr, and 1000-hr.                          !
   ! 2. We changed the thresholds for the different linear relationships, so the           !
   !    functions are continuous.                                                          !
   ! 3. We ignore the "very dry" case, and instead ensure that the "dry" function would    !
   !    be capped at 1.  We also make sure the consumption factor cannot be negative.      !
   !                                                                                       !
   ! References:                                                                           !
   !                                                                                       !
   ! Thonicke K, Spessa A, Prentice IC, Harrison SP, Dong L, Carmona-Moreno C. 2010. The   !
   !    influence of vegetation, fire spread and fire behaviour on biomass burning and     !
   !    trace gas emissions: results from a process-based model. Biogeosciences, 7:        !
   !    1991-2011. doi:10.5194/bg-7-1991-2010 (T10).                                       !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   subroutine find_fx_factors(rmoist_bherb,rmoist_bwoody,rmoist_b0001,rmoist_b0010         &
                             ,rmoist_b0100,rmoist_b1000,fx_f_bherb,fx_f_bwoody,fx_f_wn1000 &
                             ,fx_f_b0001,fx_f_b0010,fx_f_b0100,fx_f_b1000)
      use disturb_coms, only : fh_f0001      & ! intent(in)
                             , fh_f0010      & ! intent(in)
                             , fh_f0100      & ! intent(in)
                             , fh_f1000      & ! intent(in)
                             , fx_c0001_di   & ! intent(in)
                             , fx_c0001_ds   & ! intent(in)
                             , fx_c0001_mi   & ! intent(in)
                             , fx_c0001_ms   & ! intent(in)
                             , fx_c0010_di   & ! intent(in)
                             , fx_c0010_ds   & ! intent(in)
                             , fx_c0010_mi   & ! intent(in)
                             , fx_c0010_ms   & ! intent(in)
                             , fx_c0100_di   & ! intent(in)
                             , fx_c0100_ds   & ! intent(in)
                             , fx_c0100_mi   & ! intent(in)
                             , fx_c0100_ms   & ! intent(in)
                             , fx_c1000_di   & ! intent(in)
                             , fx_c1000_ds   & ! intent(in)
                             , fx_c1000_mi   & ! intent(in)
                             , fx_c1000_ms   ! ! intent(in)
      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      real, intent(in)  :: rmoist_bherb  ! Rel. moisture - Herb. fuels               [ ---]
      real, intent(in)  :: rmoist_bwoody ! Rel. moisture - Woody fuels               [ ---]
      real, intent(in)  :: rmoist_b0001  ! Rel. moisture - 1-hr fuels                [ ---]
      real, intent(in)  :: rmoist_b0010  ! Rel. moisture - 10-hr fuels               [ ---]
      real, intent(in)  :: rmoist_b0100  ! Rel. moisture - 100-hr fuels              [ ---]
      real, intent(in)  :: rmoist_b1000  ! Rel. moisture - 100-hr fuels              [ ---]
      real, intent(out) :: fx_f_bherb    ! Fuel consumpt. factor - Herb. fuels       [ ---]
      real, intent(out) :: fx_f_bwoody   ! Fuel consumpt. factor - Woody fuels       [ ---]
      real, intent(out) :: fx_f_wn1000   ! Fuel consumpt. factor - Woody (1-100hr)   [ ---]
      real, intent(out) :: fx_f_b0001    ! Fuel consumpt. factor - 1-hr fuels        [ ---]
      real, intent(out) :: fx_f_b0010    ! Fuel consumpt. factor - 10-hr fuels       [ ---]
      real, intent(out) :: fx_f_b0100    ! Fuel consumpt. factor - 100-hr fuels      [ ---]
      real, intent(out) :: fx_f_b1000    ! Fuel consumpt. factor - 100-hr fuels      [ ---]
      !------ Local variables. ------------------------------------------------------------!
      real              :: fx_f_w0001    ! Fuel consumpt. fac. - 1-hr woody fuels    [ ---]
      real              :: fx_f_w0010    ! Fuel consumpt. fac. - 10-hr woody fuels   [ ---]
      real              :: fx_f_w0100    ! Fuel consumpt. fac. - 100-hr woody fuels  [ ---]
      real              :: fx_f_w1000    ! Fuel consumpt. fac. - 1000-hr woody fuels [ ---]
      !------------------------------------------------------------------------------------!


      !----- Dead 1-hr fuels. -------------------------------------------------------------!
      fx_f_b0001  = max(0., min(1., min( fx_c0001_di + fx_c0001_ds * rmoist_b0001          &
                                       , fx_c0001_mi + fx_c0001_ms * rmoist_b0001  ) ) )
      !------------------------------------------------------------------------------------!


      !----- Dead 10-hr fuels. ------------------------------------------------------------!
      fx_f_b0010  = max(0., min(1., min( fx_c0010_di + fx_c0010_ds * rmoist_b0010          &
                                       , fx_c0010_mi + fx_c0010_ms * rmoist_b0010  ) ) )
      !------------------------------------------------------------------------------------!


      !----- Dead 100-hr fuels. -----------------------------------------------------------!
      fx_f_b0100  = max(0., min(1., min( fx_c0100_di + fx_c0100_ds * rmoist_b0100          &
                                       , fx_c0100_mi + fx_c0100_ms * rmoist_b0100  ) ) )
      !------------------------------------------------------------------------------------!


      !----- Dead 1000-hr fuels. ----------------------------------------------------------!
      fx_f_b1000  = max(0., min(1., min( fx_c1000_di + fx_c1000_ds * rmoist_b1000          &
                                       , fx_c1000_mi + fx_c1000_ms * rmoist_b1000  ) ) )
      !------------------------------------------------------------------------------------!


      !----- Herbaceous fuels.  Use the same function as 1-hr fuels. ----------------------!
      fx_f_bherb  = max(0., min(1., min( fx_c0001_di + fx_c0001_ds * rmoist_bherb          &
                                       , fx_c0001_mi + fx_c0001_ms * rmoist_bherb  ) ) )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Woody fuels.  Partition the fuels across all categories, and weight based on  !
      ! the fraction allocated for each category.                                          !
      !------------------------------------------------------------------------------------!
      !----- Woody 1-hr fuels. ------------------------------------------------------------!
      fx_f_w0001  = max(0., min(1., min( fx_c0001_di + fx_c0001_ds * rmoist_bwoody         &
                                       , fx_c0001_mi + fx_c0001_ms * rmoist_bwoody ) ) )
      !----- Woody 10-hr fuels. -----------------------------------------------------------!
      fx_f_w0010  = max(0., min(1., min( fx_c0010_di + fx_c0010_ds * rmoist_bwoody         &
                                       , fx_c0010_mi + fx_c0010_ms * rmoist_bwoody ) ) )
      !----- Woody 100-hr fuels. ----------------------------------------------------------!
      fx_f_w0100  = max(0., min(1., min( fx_c0100_di + fx_c0100_ds * rmoist_bwoody         &
                                       , fx_c0100_mi + fx_c0100_ms * rmoist_bwoody ) ) )
      !----- Woody 1000-hr fuels. ---------------------------------------------------------!
      fx_f_w1000  = max(0., min(1., min( fx_c1000_di + fx_c1000_ds * rmoist_bwoody         &
                                       , fx_c1000_mi + fx_c1000_ms * rmoist_bwoody ) ) )
      !----- Weight fuel consumption factor (without 1000-hr fuels). ----------------------!
      fx_f_wn1000 = ( fh_f0001 * fx_f_w0001 + fh_f0010 * fx_f_w0010                        &
                    + fh_f0100 * fx_f_w0100 )                                              &
                  / (1. - fh_f1000 )
      !----- Weight fuel consumption factor. ----------------------------------------------!
      fx_f_bwoody = fh_f0001 * fx_f_w0001 + fh_f0010 * fx_f_w0010                          &
                  + fh_f0100 * fx_f_w0100 + fh_f1000 * fx_f_w1000
      !------------------------------------------------------------------------------------!


      return
   end subroutine find_fx_factors
   !=======================================================================================!
   !=======================================================================================!
end module fire
!==========================================================================================!
!==========================================================================================!

