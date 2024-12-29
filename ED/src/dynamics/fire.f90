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
      use disturb_coms  , only : include_fire           & ! intent(in)
                               , fire_parameter         & ! intent(in)
                               , fe_combusted_fast_c    & ! intent(in)
                               , fe_combusted_struct_c  ! ! intent(in)
      use pft_coms      , only : fire_s_max             & ! intent(in)
                               , fire_s_efac            ! ! intent(in)
      use consts_coms   , only : wdns                   & ! intent(in)
                               , wdnsi                  & ! intent(in)
                               , day_sec                & ! intent(in)
                               , almost_one             & ! intent(in)
                               , onetwelfth             & ! intent(in)
                               , tiny_num               & ! intent(in)
                               , lnexp_min              & ! intent(in)
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
      real                           :: fire_intensity
      real                           :: fuel
      real                           :: ignition_rate
      real                           :: mean_fire_intensity
      real                           :: sum_accp
      real                           :: mean_gndwater_si
      real                           :: mean_fuel_si
      real                           :: lnexp
      real                           :: prev_not_burnt
      real                           :: curr_not_burnt
      real                           :: fire_lethal_now
      !----- Local parameters. ------------------------------------------------------------!
      character(len=18) , parameter  :: firefile = 'edfire_details.txt'
      logical           , parameter  :: printout = .false.
      !----- Locally saved variables. -----------------------------------------------------!
      logical           , save       :: first_time = .true.
      !------------------------------------------------------------------------------------!


      !----- First time, and the user wants to print the output.  Make a header. ----------!
      if (first_time) then

         !----- Make the header. ----------------------------------------------------------!
         if (printout .and. (include_fire > 0 .and. include_fire < 4)) then
            open (unit=35,file=firefile,status='replace',action='write')
            write (unit=35,fmt='(9(a,1x))')                                                &
                     '  YEAR',      ' MONTH',      '   DAY',      '   ISI','   INTENSITY'  &
                             ,'    IGNITION','        FUEL','  SOIL_WATER','SW_THRESHOLD'
            close (unit=35,status='keep')
         end if
         !---------------------------------------------------------------------------------!

         first_time = .false.
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the number of days of last month so we can normalise the integrated       !
      ! ground water.                                                                      !
      !------------------------------------------------------------------------------------!
      call lastmonthdate(current_time,lastmonth,ndaysi)
      normfac = dtlsm * ndaysi / (day_sec)
      !------------------------------------------------------------------------------------!


      !----- Current month. ---------------------------------------------------------------!
      imo = lastmonth%month
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
            case (0)
               !---------------------------------------------------------------------------!
               !    No fire.  Reset all variables and ensure disturbance remains zero.     !
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Reset the precipitation counter for this month.                       !
               !---------------------------------------------------------------------------!
               cpoly%avg_monthly_accp(imo,isi) = 0.
               !---------------------------------------------------------------------------!


               !----- Loop over patches. --------------------------------------------------!
               resetloop_0: do ipa=1,csite%npatches
                  !----- Reset the ground water for next month. ---------------------------!
                  csite%avg_monthly_gndwater(ipa) = 0.
                  !------------------------------------------------------------------------!
               end do resetloop_0
               !---------------------------------------------------------------------------!


               !----- Reset disturbance rates. --------------------------------------------!
               cpoly%lambda_fire  (imo,isi) = 0.
               cpoly%ignition_rate    (isi) = 0.
               cpoly%burnt_area       (isi) = 0.
               !---------------------------------------------------------------------------!


               !----- Set the combusted fraction based on default values. -----------------!
               cpoly%avg_fire_f_bherb  (imo,isi) = 0.
               cpoly%avg_fire_f_bwoody (imo,isi) = 0.
               cpoly%avg_fire_f_fgc    (imo,isi) = 0.
               cpoly%avg_fire_f_stgc   (imo,isi) = 0.
               !---------------------------------------------------------------------------!

            case (1,2)
               !---------------------------------------------------------------------------!
               !     ED-1/ED-2 approaches, use monthly time step.                          !
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Find the total rainfall of the past year and reset the counter for    !
               ! this month.                                                               !
               !---------------------------------------------------------------------------!
               sum_accp                        = sum(cpoly%avg_monthly_accp(:,isi))
               cpoly%avg_monthly_accp(imo,isi) = 0.
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Initialise variables that must be integrated across patches.          !
               !---------------------------------------------------------------------------!
               ignition_rate       = 0.0
               mean_fire_intensity = 0.0
               mean_gndwater_si    = 0.0
               mean_fuel_si        = 0.0
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


                  !----- Integrate site-level average monthly_gndwater. -------------------!
                  mean_gndwater_si = mean_gndwater_si                                      &
                                   + csite%avg_monthly_gndwater(ipa) * csite%area(ipa)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Obtain fuel stocks.  The original fire model would consider all    !
                  ! above-ground biomass and no litter.                                    !
                  !------------------------------------------------------------------------!
                  fuel = 0.0
                  fuel_cohort_loop: do ico = 1,cpatch%ncohorts
                     fuel = fuel + cpatch%nplant(ico) * cpatch%agb(ico)
                  end do fuel_cohort_loop
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Determine the correct threshold to ignite fires according to the   !
                  ! fire method.  Fires occur when average soil water goes below the soil  !
                  ! moisture threshold.                                                    !
                  !------------------------------------------------------------------------!
                  if (csite%avg_monthly_gndwater(ipa) < cpoly%fire_wmass_threshold(isi))   &
                  then
                     fire_intensity      = fire_parameter
                     mean_fire_intensity = mean_fire_intensity                             &
                                         + fire_intensity * csite%area(ipa)
                  else
                     fire_intensity      = 0.0
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !    If the soil is dry, then calculate patch contribution to the        !
                  ! ignition rate.                                                         !
                  !------------------------------------------------------------------------!
                  ignition_rate = ignition_rate + fire_intensity * fuel * csite%area(ipa)
                  !------------------------------------------------------------------------!


                  !----- Integrate fuels. -------------------------------------------------!
                  mean_fuel_si = mean_fuel_si + fuel * csite%area(ipa)
                  !------------------------------------------------------------------------!


                  !----- Reset the ground water for next month. ---------------------------!
                  csite%avg_monthly_gndwater(ipa) = 0.
                  !------------------------------------------------------------------------!

               end do patchloop
               !---------------------------------------------------------------------------!



               !----- Calculate fire disturbance rate [1/month]. --------------------------!
               cpoly%lambda_fire  (imo,isi) = min(lnexp_max,onetwelfth * ignition_rate)
               if (mean_fire_intensity > 0.) then
                  cpoly%ignition_rate (isi)     = cpoly%lambda_fire  (imo,isi)             &
                                                / mean_fire_intensity
               else
                  cpoly%ignition_rate (isi)     = 0.0
               end if
               !---------------------------------------------------------------------------!



               !------ Fire intensity (reporting only). -----------------------------------!
               cpoly%fire_intensity        (isi) = mean_fire_intensity
               cpoly%avg_fire_intensity(imo,isi) = mean_fire_intensity
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Find burnt area.                                                     !
               !---------------------------------------------------------------------------!
               lnexp                         = max( lnexp_min                              &
                                                  , min( lnexp_max, - ignition_rate ) )
               cpoly%avg_burnt_area(imo,isi) = (1. - cpoly%burnt_area(isi))                &
                                             * (1. - exp(lnexp))
               cpoly%burnt_area        (isi) = min(1., cpoly%burnt_area        (isi)       &
                                                     + cpoly%avg_burnt_area(imo,isi) )
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Find fire-related rates for this month.  For most of them, we must    !
               ! take the average by burnt area (so they are average values given a fire), !
               ! and thus we must check that there was any fire in this month.             !
               !---------------------------------------------------------------------------!
               if (cpoly%avg_burnt_area(imo,isi) > tiny_num) then
                  !------------------------------------------------------------------------!
                  !     Fire occurred this month.                                          !
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Find fire lethality rate (fire mortality given that fire has       !
                  ! occurred).  Here we use the conversion between the discrete lethality  !
                  ! (akin to m in SM96) and the instantaneous, exponential lethality       !
                  ! (akin to lambda in SM96).                                              !
                  !                                                                        !
                  ! Reference:                                                             !
                  !                                                                        !
                  ! Sheil D , May RM. 1996. Mortality and recruitment rate evaluations in  !
                  !    heterogeneous tropical forests. J. Ecol., 84: 91-100.               !
                  !    doi:10.2307/2261703 (SM96).                                         !
                  !------------------------------------------------------------------------!
                  !------ Loop through patches. -------------------------------------------!
                  patchlethal_12: do ipa=1,csite%npatches
                     cpatch => csite%patch(ipa)

                     !------ Loop through cohorts. ----------------------------------------!
                     cohortlethal_12: do ico=1,cpatch%ncohorts
                        !------------------------------------------------------------------!
                        !      In ED1, fire kills all the cohorts.  Set lethality          !
                        ! probability and lethality rate accordingly.                      !
                        !------------------------------------------------------------------!
                        cpatch%fire_lethal_rate (13,ico) = cpoly%avg_burnt_area(imo,isi)
                        cpatch%fire_lethal_prob    (ico) = cpatch%fire_lethal_prob   (ico) &
                                                         + cpatch%fire_lethal_rate(13,ico)
                        cpatch%fire_lethal_rate(imo,ico) = lnexp_max
                        !------------------------------------------------------------------!
                     end do cohortlethal_12
                     !---------------------------------------------------------------------!
                  end do patchlethal_12
                  !------------------------------------------------------------------------!

               else
                  !------------------------------------------------------------------------!
                  !     No fire, set variables to zero.                                    !
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !      Fire lethality rate.                                              !
                  !------------------------------------------------------------------------!
                  !------ Loop through patches. -------------------------------------------!
                  patchnofire_12: do ipa=1,csite%npatches
                     cpatch => csite%patch(ipa)
                     !------ Loop through cohorts. ----------------------------------------!
                     cohortnofire_12: do ico=1,cpatch%ncohorts
                        !----- Lethality rate. --------------------------------------------!
                        cpatch%fire_lethal_rate(imo,ico) = 0.0
                        !------------------------------------------------------------------!
                     end do cohortnofire_12
                     !---------------------------------------------------------------------!
                  end do patchnofire_12
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !     Relative consumption rates.  Set them to zero as this is the default  !
               ! in the original ED1 fire model.                                           !
               !---------------------------------------------------------------------------!
               cpoly%avg_fire_f_bherb (imo,isi) = 0.0
               cpoly%avg_fire_f_bwoody(imo,isi) = 0.0
               cpoly%avg_fire_f_fgc   (imo,isi) = 0.0
               cpoly%avg_fire_f_stgc  (imo,isi) = 0.0
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Print the output if needed.                                           !
               !---------------------------------------------------------------------------!
               if (printout) then
                  open(unit=35,file=firefile,status='old',position='append',action='write')
                  write(unit=35,fmt='(4(i6,1x),5(f12.6,1x))')                              &
                             current_time%year,current_time%month,current_time%date,isi    &
                            ,mean_fire_intensity,ignition_rate,mean_fuel_si                &
                            ,mean_gndwater_si,cpoly%fire_wmass_threshold(isi)
                  close(unit=35,status='keep')
               end if
               !---------------------------------------------------------------------------!

            case (3)
               !---------------------------------------------------------------------------!
               !     EMBERFIRE:  Fires have already been integrated over the month,        !
               ! calculate disturbance area from average fire intensity.                   !
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Reset the precipitation counter for this month.                       !
               !---------------------------------------------------------------------------!
               cpoly%avg_monthly_accp(imo,isi) = 0.
               !---------------------------------------------------------------------------!


               !----- Loop over patches. --------------------------------------------------!
               resetloop_3: do ipa=1,csite%npatches
                  !----- Reset the ground water for next month. ---------------------------!
                  csite%avg_monthly_gndwater(ipa) = 0.
                  !------------------------------------------------------------------------!
               end do resetloop_3
               !---------------------------------------------------------------------------!


               !----- Use fire "intensity" to find disturbance rate. ----------------------!
               cpoly%lambda_fire  (imo,isi) = min( lnexp_max                               &
                                                 , onetwelfth                              &
                                                 * cpoly%avg_fire_intensity(imo,isi) )
               cpoly%ignition_rate    (isi) = cpoly%avg_fire_intensity(imo,isi)            &
                                            / fire_parameter
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Find burnt area.                                                     !
               !---------------------------------------------------------------------------!
               lnexp = max( lnexp_min, min( lnexp_max, - cpoly%lambda_fire(imo,isi) ) )
               cpoly%avg_burnt_area(imo,isi) = (1. - cpoly%burnt_area(isi))                &
                                             * (1. - exp(lnexp))
               cpoly%burnt_area        (isi) = min(1., cpoly%burnt_area        (isi)       &
                                                     + cpoly%avg_burnt_area(imo,isi) )
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Find fire-related rates for this month.  For most of them, we must    !
               ! take the average by burnt area (so they are average values given a fire), !
               ! and thus we must check that there was any fire in this month.             !
               !---------------------------------------------------------------------------!
               if (cpoly%avg_burnt_area(imo,isi) > tiny_num) then
                  !------------------------------------------------------------------------!
                  !     Fire occurred this month.                                          !
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Find fire lethality rate (fire mortality given that fire has       !
                  ! occurred).  Here we use the conversion between the discrete lethality  !
                  ! (akin to m in SM96) and the instantaneous, exponential lethality       !
                  ! (akin to lambda in SM96).                                              !
                  !                                                                        !
                  ! Reference:                                                             !
                  !                                                                        !
                  ! Sheil D , May RM. 1996. Mortality and recruitment rate evaluations in  !
                  !    heterogeneous tropical forests. J. Ecol., 84: 91-100.               !
                  !    doi:10.2307/2261703 (SM96).                                         !
                  !------------------------------------------------------------------------!
                  !------ Loop through patches. -------------------------------------------!
                  patchlethal_3: do ipa=1,csite%npatches
                     cpatch => csite%patch(ipa)

                     !------ Loop through cohorts. ----------------------------------------!
                     cohortlethal_3: do ico=1,cpatch%ncohorts
                        ipft = cpatch%pft(ico)

                        !------------------------------------------------------------------!
                        !      Normalise lethality by the area burnt this month.  This     !
                        ! gives the "discrete" lethality (assuming delta t = 1 month).     !
                        !------------------------------------------------------------------!
                        lnexp           = fire_s_efac(ipft) * cpatch%thbark(ico)
                        lnexp           = max(lnexp_min,min(lnexp_max,lnexp))
                        fire_lethal_now = ( 1. - fire_s_max(ipft) )                        &
                                        / ( 1. - fire_s_max(ipft) * exp(lnexp) ) 
                        fire_lethal_now = max(0.,min(1.,fire_lethal_now))
                        !------------------------------------------------------------------!



                        !------- Update fire lethality. -----------------------------------!
                        cpatch%fire_lethal_rate(13,ico) = fire_lethal_now                  &
                                                        * cpoly%avg_burnt_area(imo,isi)
                        cpatch%fire_lethal_prob   (ico) = cpatch%fire_lethal_prob   (ico)  &
                                                        + cpatch%fire_lethal_rate(13,ico)
                        !------------------------------------------------------------------!




                        !------------------------------------------------------------------!
                        !     Check to see if the lethality was excessive (which could     !
                        ! cause singularities in the general conversion equation).         !
                        !------------------------------------------------------------------!
                        if (fire_lethal_now > almost_one) then
                           !---------------------------------------------------------------!
                           !     Cataclysmic fire, all affected plants died. The lethality !
                           ! rate should be infinity, due to numeric precision we set it   !
                           ! to the maximum number we can find exponentials without        !
                           ! trigger FPE.                                                  !
                           !---------------------------------------------------------------!
                           cpatch%fire_lethal_rate(imo,ico) = lnexp_max
                           !---------------------------------------------------------------!
                        else
                           !------ Set monthly fire lethality rate. -----------------------!
                           cpatch%fire_lethal_rate(imo,ico) = log(1./(1.-fire_lethal_now))
                           !---------------------------------------------------------------!
                        end if
                        !------------------------------------------------------------------!
                     end do cohortlethal_3
                     !---------------------------------------------------------------------!
                  end do patchlethal_3
                  !------------------------------------------------------------------------!




                  !----- Set the combusted fraction based on default values. --------------!
                  cpoly%avg_fire_f_bherb  (imo,isi) = fe_combusted_fast_c
                  cpoly%avg_fire_f_bwoody (imo,isi) = fe_combusted_struct_c
                  cpoly%avg_fire_f_fgc    (imo,isi) = fe_combusted_fast_c
                  cpoly%avg_fire_f_stgc   (imo,isi) = fe_combusted_struct_c
                  !------------------------------------------------------------------------!

               else
                  !------------------------------------------------------------------------!
                  !     No fire, set variables to zero.                                    !
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !      Fire lethality rate.                                              !
                  !------------------------------------------------------------------------!
                  !------ Loop through patches. -------------------------------------------!
                  patchnofire_3: do ipa=1,csite%npatches
                     cpatch => csite%patch(ipa)
                     !------ Loop through cohorts. ----------------------------------------!
                     cohortnofire_3: do ico=1,cpatch%ncohorts
                        !----- Lethality rate. --------------------------------------------!
                        cpatch%fire_lethal_rate(imo,ico) = 0.0
                        !------------------------------------------------------------------!
                     end do cohortnofire_3
                     !---------------------------------------------------------------------!
                  end do patchnofire_3
                  !------------------------------------------------------------------------!


                  !----- Relative consumption rates. --------------------------------------!
                  cpoly%avg_fire_f_bherb (imo,isi) = 0.0
                  cpoly%avg_fire_f_bwoody(imo,isi) = 0.0
                  cpoly%avg_fire_f_fgc   (imo,isi) = 0.0
                  cpoly%avg_fire_f_stgc  (imo,isi) = 0.0
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!

            case (4)
               !---------------------------------------------------------------------------!
               !     FIRESTARTER:  Fires have already been integrated over the month,      !
               ! calculate disturbance area from burnt area.                               !
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Reset the precipitation counter for this month.                       !
               !---------------------------------------------------------------------------!
               cpoly%avg_monthly_accp(imo,isi) = 0.
               !---------------------------------------------------------------------------!
               !----- Loop over patches. --------------------------------------------------!
               patchloop_4: do ipa=1,csite%npatches
                  !----- Reset the ground water for next month. ---------------------------!
                  csite%avg_monthly_gndwater(ipa) = 0.
                  !------------------------------------------------------------------------!
               end do patchloop_4
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Find fire-related rates for this month.  For most of them, we must    !
               ! take the average by burnt area (so they are average values given a fire), !
               ! and thus we must check that there was any fire in this month.             !
               !---------------------------------------------------------------------------!
               if (cpoly%avg_burnt_area(imo,isi) > tiny_num) then
                  !------------------------------------------------------------------------!
                  !     Fire occurred this month.                                          !
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Find fire lethality rate (fire mortality given that fire has       !
                  ! occurred).  Here we use the conversion between the discrete lethality  !
                  ! (akin to m in SM96) and the instantaneous, exponential lethality       !
                  ! (akin to lambda in SM96).                                              !
                  !                                                                        !
                  ! Reference:                                                             !
                  !                                                                        !
                  ! Sheil D , May RM. 1996. Mortality and recruitment rate evaluations in  !
                  !    heterogeneous tropical forests. J. Ecol., 84: 91-100.               !
                  !    doi:10.2307/2261703 (SM96).                                         !
                  !------------------------------------------------------------------------!
                  !------ Loop through patches. -------------------------------------------!
                  patchlethal_4: do ipa=1,csite%npatches
                     cpatch => csite%patch(ipa)

                     !------ Loop through cohorts. ----------------------------------------!
                     cohortlethal_4: do ico=1,cpatch%ncohorts
                        !------------------------------------------------------------------!
                        !      Normalise lethality by the area burnt this month.  This     !
                        ! gives the "discrete" lethality (assuming delta t = 1 month).     !
                        !------------------------------------------------------------------!
                        fire_lethal_now = cpatch%fire_lethal_rate(13,ico)                  &
                                        / cpoly%avg_burnt_area(imo,isi)
                        fire_lethal_now = max(0.,min(1.,fire_lethal_now))
                        !------------------------------------------------------------------!


                        !------------------------------------------------------------------!
                        !     Check to see if the lethality was excessive (which could     !
                        ! cause singularities in the general conversion equation).         !
                        !------------------------------------------------------------------!
                        if (fire_lethal_now > almost_one) then
                           !---------------------------------------------------------------!
                           !     Cataclysmic fire, all affected plants died. The lethality !
                           ! rate should be infinity, due to numeric precision we set it   !
                           ! to the maximum number we can find exponentials without        !
                           ! trigger FPE.                                                  !
                           !---------------------------------------------------------------!
                           cpatch%fire_lethal_rate(imo,ico) = lnexp_max
                           !---------------------------------------------------------------!
                        else
                           !------ Set monthly fire lethality rate. -----------------------!
                           cpatch%fire_lethal_rate(imo,ico) = log(1./(1.-fire_lethal_now))
                           !---------------------------------------------------------------!
                        end if
                        !------------------------------------------------------------------!
                     end do cohortlethal_4
                     !---------------------------------------------------------------------!
                  end do patchlethal_4
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !      Normalise fractional fire fuel consumption.  Similarly to         !
                  ! lethality, we want to obtain the fractions given that a fire has       !
                  ! occurred, because this is applied to new patches (by definition, the   !
                  ! areas where disturbance occurred).                                     !
                  !------------------------------------------------------------------------!
                  !----- Divide fraction by area burnt this month. ------------------------!
                  cpoly%avg_fire_f_bherb (imo,isi) = cpoly%avg_fire_f_bherb (imo,isi)      &
                                                   / cpoly%avg_burnt_area   (imo,isi)
                  cpoly%avg_fire_f_bwoody(imo,isi) = cpoly%avg_fire_f_bwoody(imo,isi)      &
                                                   / cpoly%avg_burnt_area   (imo,isi)
                  cpoly%avg_fire_f_fgc   (imo,isi) = cpoly%avg_fire_f_fgc   (imo,isi)      &
                                                   / cpoly%avg_burnt_area   (imo,isi)
                  cpoly%avg_fire_f_stgc  (imo,isi) = cpoly%avg_fire_f_stgc  (imo,isi)      &
                                                   / cpoly%avg_burnt_area   (imo,isi)
                  !----- Ensure consumption terms are bounded. ----------------------------!
                  cpoly%avg_fire_f_bherb (imo,isi) =                                       &
                                           max(0.,min(1.,cpoly%avg_fire_f_bherb (imo,isi)))
                  cpoly%avg_fire_f_bwoody(imo,isi) =                                       &
                                           max(0.,min(1.,cpoly%avg_fire_f_bwoody(imo,isi)))
                  cpoly%avg_fire_f_fgc   (imo,isi) =                                       &
                                           max(0.,min(1.,cpoly%avg_fire_f_fgc   (imo,isi)))
                  cpoly%avg_fire_f_stgc  (imo,isi) =                                       &
                                           max(0.,min(1.,cpoly%avg_fire_f_stgc  (imo,isi)))
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !     Fire disturbance rate.  This also follows SM96, but we ought to    !
                  ! consider the burnt area relative to the area not previously burnt.     !
                  !------------------------------------------------------------------------!
                  if (cpoly%burnt_area(isi) > almost_one) then
                     !---------------------------------------------------------------------!
                     !       The entire grid cell burned this month, assume maximum        !
                     ! disturbance rate.                                                   !
                     !---------------------------------------------------------------------!
                     cpoly%lambda_fire(imo,isi) = lnexp_max
                     !---------------------------------------------------------------------!
                  else
                     !---------------------------------------------------------------------!
                     !      Compute the disturbance rate based on the new burnt area       !
                     ! relative to the area not previously burnt.                          !
                     !---------------------------------------------------------------------!
                     curr_not_burnt             = 1. - cpoly%burnt_area(isi)
                     prev_not_burnt             = curr_not_burnt                           &
                                                + cpoly%avg_burnt_area(imo,isi)
                     cpoly%lambda_fire(imo,isi) = log(prev_not_burnt/curr_not_burnt)
                     cpoly%lambda_fire(imo,isi) = min(lnexp_max,cpoly%lambda_fire(imo,isi))
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!

               else
                  !------------------------------------------------------------------------!
                  !     No fire, set variables to zero.                                    !
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !      Fire lethality rate.                                              !
                  !------------------------------------------------------------------------!
                  !------ Loop through patches. -------------------------------------------!
                  patchnofire_4: do ipa=1,csite%npatches
                     cpatch => csite%patch(ipa)
                     !------ Loop through cohorts. ----------------------------------------!
                     cohortnofire_4: do ico=1,cpatch%ncohorts
                        !----- Lethality rate. --------------------------------------------!
                        cpatch%fire_lethal_rate(imo,ico) = 0.0
                        !------------------------------------------------------------------!
                     end do cohortnofire_4
                     !---------------------------------------------------------------------!
                  end do patchnofire_4
                  !------------------------------------------------------------------------!

                  !----- Fire disturbance rate. -------------------------------------------!
                  cpoly%lambda_fire(imo,isi) = 0.0
                  !------------------------------------------------------------------------!

                  !----- Relative consumption rates. --------------------------------------!
                  cpoly%avg_fire_f_bherb (imo,isi) = 0.0
                  cpoly%avg_fire_f_bwoody(imo,isi) = 0.0
                  cpoly%avg_fire_f_fgc   (imo,isi) = 0.0
                  cpoly%avg_fire_f_stgc  (imo,isi) = 0.0
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Reset the precipitation counter for this month.                       !
               !---------------------------------------------------------------------------!
               cpoly%avg_monthly_accp(imo,isi) = 0.
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
   ! SUB-ROUTINE INTEG_FIRE_DANGER
   !\brief Subroutine that integrates variables used to assess fire risk.
   !\details This sub-routine that integrates the fire disturbance rate when using the 
   !!       Sub-routine that integrates the Nesterov index (T10) and the VPD-based fire
   !!        danger index (D19), but replacing their approximation for VPD with the 
   !!        calculated VPD scaled by the reference pressure at sea level (see rationale in
   !!        PS09).  The foliar projective cover is based on (S18), but using accounting
   !!        for leaf orientation factor. Although not used by every fire model, we still
   !!        define the indices in case they are useful as diagnostic variables.
   !!
   !! Reference:
   !! 
   !! Druke  M, Forkel M, von Bloh W, Sakschewski B, Cardoso M, Bustamante M, Kurths J,
   !!    Thonicke K. 2019. Improving the LPJmL4-SPITFIRE vegetation--fire model for South
   !!    America using satellite data. Geosci. Model Dev., 12: 5029-5054.
   !!    doi:10.5194/gmd-12-5029-2019 (D19).
   !!
   !! Schaphoff S, von Bloh W, Rammig A, Thonicke K, Biemans H, Forkel M, Gerten D, 
   !!    Heinke J, Jagermeyr J, Knauer J et al. 2018. LPJmL4 -- a dynamic global 
   !!    vegetation model with managed land -- part 1: Model description. 
   !!    Geosci. Model Dev., 11: 13431375. doi:10.5194/gmd-11-1343-2018 (S18).
   !!
   !! Pechony O, Shindell DT. 2009. Fire parameterization on a global scale. J. Geophys.
   !!    Res.-Atmos., 114(D16): D16115. doi:10.1029/2009JD011927 (PS09).
   !!
   !! Thonicke K, Spessa A, Prentice IC, Harrison SP, Dong L, Carmona-Moreno C. 2010. The
   !!    influence of vegetation, fire spread and fire behaviour on biomass burning and
   !!    trace gas emissions: results from a process-based model. Biogeosciences, 7:
   !!    1991-2011. doi:10.5194/bg-7-1991-2010 (T10).
   !---------------------------------------------------------------------------------------!
   subroutine integ_fire_danger(cgrid)
      use ed_state_vars        , only : edtype       & ! structure
                                      , polygontype  & ! structure
                                      , sitetype     & ! structure
                                      , patchtype    ! ! structure
      use ed_misc_coms         , only : current_time & ! intent(in)
                                      , ndfire       ! ! intent(in)
      use consts_coms          , only : t00          & ! intent(in)
                                      , day_sec      & ! intent(in)
                                      , tiny_num     & ! intent(in)
                                      , lnexp_min    & ! intent(in)
                                      , lnexp_max    & ! intent(in)
                                      , prefsea      ! ! intent(in)
      use disturb_coms         , only : fh_pcpg_ni0  & ! intent(in)
                                      , fh_pcpg_edi  & ! intent(in)
                                      , fh_pcpg_win  ! ! intent(in)
      use canopy_radiation_coms, only : eproj_light  ! ! intent(in)
      use pft_coms             , only : alpha_fdivpd ! ! intent(in)
      implicit none
      !----- -Arguments. ------------------------------------------------------------------!
      type(edtype)     , target     :: cgrid
      !------ Local variables. ------------------------------------------------------------!
      type(polygontype), pointer    :: cpoly
      type(sitetype)   , pointer    :: csite
      type(patchtype)  , pointer    :: cpatch
      integer                       :: ipy
      integer                       :: isi
      integer                       :: ipa
      integer                       :: ico
      integer                       :: ipft
      integer                       :: k
      logical                       :: dry_day
      real                          :: avgrun_accp
      real                          :: frain_fdivpd
      real                          :: today_nesterov
      real                          :: today_accp
      real                          :: tdmax_atm_temp
      real                          :: tdmin_atm_temp
      real                          :: today_atm_tdew
      real                          :: tdmax_atm_vpdef
      real                          :: tdmin_atm_vpdef
      real                          :: tdmax_can_temp
      real                          :: tdmin_can_temp
      real                          :: tdmax_can_vpdef
      real                          :: tdmin_can_vpdef
      real                          :: today_can_tdew
      real                          :: today_pcpg
      real                          :: fdi_vpd_max
      real                          :: fdi_vpd_min
      real                          :: lnexp
      real                          :: lai_ind
      real                          :: fpc_pat
      real                          :: fpc_coh
      real                          :: alpha_pat
      real                          :: light_above
      !----- Local parameters. ------------------------------------------------------------!
      character(len=22) , parameter :: firefile = 'firedanger_details.txt'
      logical           , parameter :: printout = .false.
      !----- Locally saved variables. -----------------------------------------------------!
      logical           , save      :: first_time = .true.
      real              , save      :: wgt_running
      real              , save      :: wgt_today
      real              , save      :: ndfirei
      !------------------------------------------------------------------------------------!


      !----- First time, and the user wants to print the output.  Make a header. ----------!
      if (first_time) then

         !----- Make the header. ----------------------------------------------------------!
         if (printout) then
            open (unit=35,file=firefile,status='replace',action='write')
            write (unit=35,fmt='(27(a,1x))')                                               &
                     '  YEAR',      ' MONTH',      '   DAY',      '   ISI',      '   IPA'  &
              ,'        AREA','   CAN_DEPTH','   ALPHA_VPD','     FPC_PAT','      PRECIP'  &
              ,' RUNAVG_PREC','   WGT_TODAY',' WGT_RUNNING','ATM_TEMP_MAX','CAN_TEMP_MAX'  &
              ,'ATM_TEMP_MIN','CAN_TEMP_MIN','    ATM_TDEW','    CAN_TDEW',' ATM_VPD_MAX'  &
              ,' CAN_VPD_MAX',' ATM_VPD_MIN',' CAN_VPD_MIN','NESTEROV_PAT','NESTEROV_INT'  &
              ,' FDI_VPD_MAX',' FDI_VPD_MIN'
            close (unit=35,status='keep')
         end if
         !---------------------------------------------------------------------------------!


         !----- Define weight for running average. ----------------------------------------!
         if (fh_pcpg_win > 1.) then
            !----- Weighting factor for today is the inverse of running average window. ---!
            wgt_today   = max(0.,min(1.,1.  / fh_pcpg_win))
            wgt_running = 1. - wgt_today
            !------------------------------------------------------------------------------!
         else
            !----- Running average window is too short, use today's value only. -----------!
            wgt_today   = 1.0
            wgt_running = 0.0
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!


         !----- Inverse of number of bins in a day. ---------------------------------------!
         ndfirei = 1. / real(ndfire)
         !---------------------------------------------------------------------------------!


         first_time = .false.
      end if
      !------------------------------------------------------------------------------------!


      !------ Loop over polygons. ---------------------------------------------------------!
      poly_loop: do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)



         !---------------------------------------------------------------------------------!
         !      Loop over sites.                                                           !
         !---------------------------------------------------------------------------------!
         site_loop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)


            !----- Update precipitation running average. ----------------------------------!
            today_pcpg                  = ndfirei * sum(cpoly%tdfire_pcpg(:,isi),dim=1)
            cpoly%avg_running_pcpg(isi) = wgt_running * cpoly%avg_running_pcpg(isi)        &
                                        + wgt_today   * today_pcpg
            avgrun_accp                 = cpoly%avg_running_pcpg(isi) * day_sec
            !------------------------------------------------------------------------------!



            !----- Find the rainfall down-regulation term for VPD-based FDI. --------------!
            lnexp        = max( lnexp_min                                                  &
                              , min( lnexp_max                                             &
                                   , fh_pcpg_edi * cpoly%avg_running_pcpg(isi) ) )
            frain_fdivpd = exp(lnexp)
            !------------------------------------------------------------------------------!



            !----- Convert air temperature to degC (useful for the report). ---------------!
            tdmax_atm_temp  = maxval(cpoly%tdfire_atm_temp (:,isi),dim=1)           - t00
            tdmin_atm_temp  = minval(cpoly%tdfire_atm_temp (:,isi),dim=1)           - t00
            today_atm_tdew  = sum   (cpoly%tdfire_atm_tdew (:,isi),dim=1) * ndfirei - t00
            tdmax_atm_vpdef = maxval(cpoly%tdfire_atm_vpdef(:,isi),dim=1)           * 0.01
            tdmin_atm_vpdef = minval(cpoly%tdfire_atm_vpdef(:,isi),dim=1)           * 0.01
            !------------------------------------------------------------------------------!



            !----- Check whether this is a dry day. ---------------------------------------!
            today_accp = today_pcpg * day_sec
            dry_day    = today_pcpg <= fh_pcpg_ni0
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Update Nesterov index.  Find the site average by looping through         !
            ! patches.                                                                     !
            !------------------------------------------------------------------------------!
            today_nesterov = 0.0
            patch_loop: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)


               !----- Convert temperature to degC (also useful for the report). -----------!
               tdmax_can_temp  = maxval(csite%tdfire_can_temp (:,ipa),dim=1) - t00
               tdmin_can_temp  = minval(csite%tdfire_can_temp (:,ipa),dim=1) - t00
               today_can_tdew  = sum   (csite%tdfire_can_tdew (:,ipa),dim=1) * ndfirei - t00
               tdmax_can_vpdef = maxval(csite%tdfire_can_vpdef(:,ipa),dim=1) * 0.01
               tdmin_can_vpdef = minval(csite%tdfire_can_vpdef(:,ipa),dim=1) * 0.01
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !      Loop through cohorts.                                                !
               !---------------------------------------------------------------------------!
               light_above = 1.
               alpha_pat   = 0.
               cohort_loop: do ico=1,cpatch%ncohorts
                  !------ Handy aliases. --------------------------------------------------!
                  ipft = cpatch%pft(ico)
                  !------------------------------------------------------------------------!


                  !------ Find the individual leaf cover (S18). ---------------------------!
                  if (cpatch%leaf_resolvable(ico)) then
                     lai_ind = cpatch%lai(ico) / cpatch%crown_area(ico)
                     lnexp   = max(lnexp_min,min(lnexp_max,-eproj_light(ipft)*lai_ind))
                     fpc_coh = light_above * cpatch%crown_area(ico) * (1. - exp(lnexp))
                  else
                     fpc_coh = 0.
                  end if
                  !------------------------------------------------------------------------!

                  !------ Integrate foliar projective cover, and the alpha term. ----------!
                  light_above = max(0.,light_above - fpc_coh)
                  alpha_pat   = alpha_pat + alpha_fdivpd(ipft) * fpc_coh
                  !------------------------------------------------------------------------!
               end do cohort_loop
               !---------------------------------------------------------------------------!



               !----- Normalise the patch-level alpha factor. -----------------------------!
               fpc_pat = max(0.,min(1.,1. - light_above))
               if (fpc_pat > tiny_num) then
                  alpha_pat = alpha_pat / fpc_pat
               else
                  alpha_pat = 0.
               end if
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !      Find the fire danger index, following D19, but using the VPD from    !
               ! ED2 equations (Murphy-Koop) instead of the Goff-Gratch equation.  Follow- !
               ! ing the rationale of PS09, we make VPD dimensionless by dividing it by    !
               ! the reference pressure at the sea level.  We also make sure that values   !
               ! are bounded between 0. and 1.                                             !
               !---------------------------------------------------------------------------!
               do k=1,ndfire
                  csite%tdfire_fdi_vpd(k,ipa) = alpha_pat * fpc_pat * frain_fdivpd         &
                                              * csite%tdfire_can_vpdef(k,ipa) / prefsea
                  csite%tdfire_fdi_vpd(k,ipa) = max( 0.                                    &
                                                   , min( 1., csite%tdfire_fdi_vpd(k,ipa)) )
               end do
               fdi_vpd_max                 = maxval(csite%tdfire_fdi_vpd(:,ipa),dim=1)
               fdi_vpd_min                 = minval(csite%tdfire_fdi_vpd(:,ipa),dim=1)
               !---------------------------------------------------------------------------!


               !---- Patch Nesterov index, and make sure it is never negative. ------------!
               if (dry_day) then
                  today_nesterov = max(0.,tdmax_can_temp*(tdmax_can_temp-today_can_tdew))
               else
                  today_nesterov = 0.
               end if
               !---------------------------------------------------------------------------!


               !----- Add patch contribution to today's Nesterov index. -------------------!
               csite%nesterov_index(ipa) = csite%nesterov_index(ipa)                       &
                                         + today_nesterov
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Print the output if needed.                                           !
               !---------------------------------------------------------------------------!
               if (printout) then
                  open(unit=35,file=firefile,status='old',position='append',action='write')
                  write(unit=35,fmt='(5(i6,1x),22(f12.4,1x))')                             &
                             current_time%year,current_time%month,current_time%date,isi    &
                            ,ipa,csite%area(ipa),csite%can_depth(ipa),alpha_pat            &
                            ,fpc_pat,today_accp,avgrun_accp,wgt_today,wgt_running          &
                            ,tdmax_atm_temp,tdmax_can_temp,tdmin_atm_temp,tdmin_can_temp   &
                            ,today_atm_tdew,today_can_tdew,tdmax_atm_vpdef,tdmax_can_vpdef &
                            ,tdmin_atm_vpdef,tdmin_can_vpdef,today_nesterov                &
                            ,csite%nesterov_index(ipa),fdi_vpd_max,fdi_vpd_min
                  close(unit=35,status='keep')
               end if
               !---------------------------------------------------------------------------!
           end do patch_loop
            !------------------------------------------------------------------------------!
         end do site_loop
         !---------------------------------------------------------------------------------!
      end do poly_loop
      !------------------------------------------------------------------------------------!


      return
   end subroutine integ_fire_danger
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   ! SUB-ROUTINE INTEG_EMBERFIRE
   !\brief Integration step of the EMBERFIRE model
   !\details This subroutine will integrate the fuel loads and fire intensity for the
   !!        "Empirical Model for Basic Ecosystem Response to FIRE" (EMBERFIRE model).
   !---------------------------------------------------------------------------------------!
   subroutine integ_emberfire(cgrid)
      use ed_state_vars , only : edtype            & ! structure
                               , polygontype       & ! structure
                               , sitetype          & ! structure
                               , patchtype         ! ! structure
      use ed_misc_coms  , only : simtime           & ! intent(in)
                               , current_time      & ! intent(in)
                               , ndfire            ! ! intent(in)
      use pft_coms      , only : agf_bs            & ! intent(in)
                               , is_grass          & ! intent(in)
                               , f_labile_leaf     & ! intent(in)
                               , f_labile_stem     ! ! intent(in)
      use disturb_coms  , only : fire_parameter    & ! intent(in)
                               , fuel_height_max   & ! intent(in)
                               , fe_anth_ignt_only & ! intent(in)
                               , fe_fdivpd_exp     & ! intent(in)
                               , fe_fdivpd_slp     & ! intent(in)
                               , fe_use_fdivpd     & ! intent(in)
                               , fh_f0001          & ! intent(in)
                               , fh_f0010          & ! intent(in)
                               , fh_f0100          & ! intent(in)
                               , fh_f1000          & ! intent(in)
                               , fx_a0001          & ! intent(in)
                               , fx_a0010          & ! intent(in)
                               , fx_a0100          & ! intent(in)
                               , fx_rmfac          & ! intent(in)
                               , n_fst             ! ! intent(in)
      use consts_coms   , only : day_sec           & ! intent(in)
                               , t00               & ! intent(in)
                               , tiny_num          & ! intent(in)
                               , lnexp_min         & ! intent(in)
                               , lnexp_max         ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)  , target     :: cgrid             ! Current grid              [    ---]
      !----- Local variables --------------------------------------------------------------!
      type(polygontype), pointer :: cpoly             ! Current polygon           [    ---]
      type(sitetype)   , pointer :: csite             ! Current site              [    ---]
      type(patchtype)  , pointer :: cpatch            ! Current patch             [    ---]
      type(simtime)              :: hier              ! Yesterdays' date info.    [    ---]
      integer                    :: ico               ! Cohort index              [    ---]
      integer                    :: imonth            ! Month matrix index        [    ---]
      integer                    :: ipa               ! Patch index               [    ---]
      integer                    :: ipft              ! PFT index                 [    ---]
      integer                    :: ipy               ! Polygon index             [    ---]
      integer                    :: isi               ! Site index                [    ---]
      integer                    :: ndays             ! # days in month           [    ---]
      logical                    :: people_around     ! Any anth. disturb. type   [    T|F]
      real, dimension(n_fst)     :: Mx_i              ! Fuel moist. extinct.      [  kg/kg]
      real                       :: bfuel_d0001_pat   ! Patch 1-hr dead fuels     [ kgC/m2]
      real                       :: bfuel_d0001_tot   ! Site 1-hr dead fuels      [ kgC/m2]
      real                       :: bfuel_d0010_pat   ! Patch 10-hr dead fuels    [ kgC/m2]
      real                       :: bfuel_d0010_tot   ! Site 10-hr dead fuels     [ kgC/m2]
      real                       :: bfuel_d0100_pat   ! Patch 100-hr dead fuels   [ kgC/m2]
      real                       :: bfuel_d0100_tot   ! Site 100-hr dead fuels    [ kgC/m2]
      real                       :: bfuel_d1000_pat   ! Patch 1000-hr dead fuels  [ kgC/m2]
      real                       :: bfuel_d1000_tot   ! Site 1000-hr dead fuels   [ kgC/m2]
      real                       :: bfuel_d0111_pat   ! (1+10+100-hr fuels)       [ kgC/m2]
      real                       :: bfuel_d0111_tot   ! (1+10+100-hr fuels)       [ kgC/m2]
      real                       :: bfuel_dead_tot    ! Total dead fuels          [ kgC/m2]
      real                       :: bfuel_live_tot    ! Total live fuels          [ kgC/m2]
      real                       :: bherb             ! Cohort Herbaceous fuels   [ kgC/pl]
      real                       :: bherb_pat         ! Patch Herbaceous fuels    [ kgC/m2]
      real                       :: bherb_tot         ! Site Herbaceous fuels     [ kgC/m2]
      real                       :: bwoody            ! Cohort living woody fuels [ kgC/pl]
      real                       :: bwoody_pat        ! Patch living woody fuels  [ kgC/m2]
      real                       :: bwoody_tot        ! Site living woody fuels   [ kgC/m2]
      real                       :: fdivpd_pat        ! Patch fire danger index   [     --]
      real                       :: fdivpd_avg        ! Average fire danger index [     --]
      real                       :: ignition_rate     ! Ignition probability rate [     --]
      real                       :: lnexp             ! Aux. var. for safe exp    [    ---]
      real                       :: moist_bfuel_avg   ! Dead fuel moisture        [    ---]
      real                       :: moist_bfuel_pat   ! Dead fuel moisture        [    ---]
      real                       :: moist_bherb_avg   ! Herbaceous fuel moisture  [    ---]
      real                       :: moist_bherb_pat   ! Herbaceous fuel moisture  [    ---]
      real                       :: moist_bwoody_avg  ! Living woody fuel moist.  [    ---]
      real                       :: moist_bwoody_pat  ! Living woody fuel moist.  [    ---]
      real                       :: ndaysi            ! 1/# days in a month       [  1/day]
      real                       :: nesterov_avg      ! Average Nesterov index    [  degC2]
      real                       :: rmoist_avg        ! Relative fuel moisture    [    ---]
      real                       :: tdmax_can_temp    ! Maximum temperature       [   degC]
      real                       :: today_can_tdew    ! Dew point temperature     [   degC]
      real                       :: tdmax_can_vpdef   ! Max. Vapour press. def.   [    hPa]
      real                       :: tdmin_can_vpdef   ! Max. Vapour press. def.   [    hPa]
      real                       :: today_pcpg        ! Daily rainfall            [     mm]
      !------ External functions. ---------------------------------------------------------!
      real              , external  :: bpow01         ! Power funct. for [0-1]    [    ---]
      !----- Local parameters. ------------------------------------------------------------!
      character(len=21) , parameter :: firefile = 'emberfire_details.txt'
      logical           , parameter :: printout = .false.
      !----- Locally saved variables. -----------------------------------------------------!
      logical           , save      :: first_time = .true. ! First time calling routine
      real              , save      :: ndfirei             ! 1. / ndfire
      !------------------------------------------------------------------------------------!


      !----- First time, and the user wants to print the output.  Make a header. ----------!
      if (first_time) then

         !----- Make the header. ----------------------------------------------------------!
         if (printout) then
            open (unit=35,file=firefile,status='replace',action='write')
            write (unit=35,fmt='(23(a,1x))')                                               &
                     '  YEAR',      ' MONTH',      '   DAY',      '   ISI',      'PEOPLE'  &
              ,'      PRECIP','CAN_TEMP_MAX','    CAN_TDEW',' CAN_VPD_MAX',' CAN_VPD_MIN'  &
              ,'    NESTEROV','     FDI_VPD','  BFUEL_HERB','  BFUEL_WOOD','  BFUEL_DEAD'  &
              ,' MOIST_BHERB','MOIST_BWOODY',' MOIST_BFUEL',' MSTEXT_DEAD',' MSTEXT_LIVE'  &
              ,' RMOIST_FUEL','    IGNITION','   INTENSITY'
            close (unit=35,status='keep')
         end if
         !---------------------------------------------------------------------------------!

         !------ Inverse of sub-daily bins. -----------------------------------------------!
         ndfirei = 1. / real(ndfire)
         !---------------------------------------------------------------------------------!

         first_time = .false.
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the number of days of last day so we can normalise the integrated fire    !
      ! intensity and retrieve lightning and HDI information.                              !
      !------------------------------------------------------------------------------------!
      call yesterday_info(current_time,hier,ndays,ndaysi)
      imonth = hier%month
      !------------------------------------------------------------------------------------!



      !----- Loop over polygons and sites. ------------------------------------------------!
      main_polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         !---------------------------------------------------------------------------------!
         !     Loop through all sites and patches and check whether or not there are signs !
         ! of anthropogenic activities (In case the user wants to check this).             !
         !---------------------------------------------------------------------------------!
         if (fe_anth_ignt_only) then
            !------------------------------------------------------------------------------!
            !     Assume that fires can occur if at least one patch has anthropogenic      !
            ! disturbance (any disturbance type other than tree fall as of now).           !
            !------------------------------------------------------------------------------!
            !------ Assume no people. -----------------------------------------------------!
            people_around = .false.
            !----- Loop through sites, leave if we find anthropogenic disturbance. --------!
            anth_siteloop: do isi = 1,cpoly%nsites
               csite => cpoly%site(isi)
               !----- Loop through patches. -----------------------------------------------!
               anth_patchloop: do ipa=1,csite%npatches
                  select case (csite%dist_type(ipa))
                  case (3)
                     !------ Natural disturbance. -----------------------------------------!
                     continue
                     !---------------------------------------------------------------------!
                  case default
                     !----- Anthropogenic disturbance, update flag and leave outer loop. --!
                     people_around = .true.
                     exit anth_siteloop
                     !---------------------------------------------------------------------!
                  end select
               end do anth_patchloop
               !---------------------------------------------------------------------------!
            end do anth_siteloop
            !------------------------------------------------------------------------------!
         else
            !----- Do not check for anthropogenic disturbance, allow fires everywhere. ----!
            people_around = .true.
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!






         !---------------------------------------------------------------------------------!
         !     Loop over all sites.                                                        !
         !---------------------------------------------------------------------------------!
         main_siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            !------ Initialise fuel stocks. -----------------------------------------------!
            bfuel_d0001_tot = 0.
            bfuel_d0010_tot = 0.
            bfuel_d0100_tot = 0.
            bfuel_d1000_tot = 0.
            bfuel_d0111_tot = 0.
            bherb_tot       = 0.
            bwoody_tot      = 0.
            !------------------------------------------------------------------------------!


            !------ Initialise dead and live fuel moisture. -------------------------------!
            moist_bfuel_avg    = 0.
            moist_bherb_avg    = 0.
            moist_bwoody_avg   = 0.
            !------------------------------------------------------------------------------!


            !------ Initialise average indices. -------------------------------------------!
            nesterov_avg    = 0.
            fdivpd_avg      = 0.
            !------------------------------------------------------------------------------!



            !----- Initialise site-average temperatures. ----------------------------------!
            tdmax_can_temp  = 0.
            today_can_tdew  = 0.
            tdmax_can_vpdef = 0.
            tdmin_can_vpdef = 0.
            !------------------------------------------------------------------------------!


            !----- Loop over patches. -----------------------------------------------------!
            patchloop: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)

               !---------------------------------------------------------------------------!
               !       Allocate fuels.                                                     !
               !---------------------------------------------------------------------------!
               bfuel_d0001_pat = csite%fast_grnd_C(ipa)                                    &
                               + fh_f0001 * csite%structural_grnd_C(ipa)
               bfuel_d0010_pat = fh_f0010 * csite%structural_grnd_C(ipa)
               bfuel_d0100_pat = fh_f0100 * csite%structural_grnd_C(ipa)
               bfuel_d1000_pat = fh_f1000 * csite%structural_grnd_C(ipa)
               bfuel_d0111_pat = bfuel_d0100_pat + bfuel_d0010_pat + bfuel_d0001_pat
               bherb_pat       = 0.
               bwoody_pat      = 0.
               cohort_loop: do ico=1,cpatch%ncohorts
                  ipft = cpatch%pft(ico)
                  if (is_grass(ipft) .or. cpatch%hite(ico) <= fuel_height_max) then
                     !------ Herbaceous fuel.  AG labile biomass + AG storage. ------------!
                     bherb  = f_labile_leaf(ipft) * cpatch%bleaf(ico)                      &
                            + f_labile_stem(ipft)                                          &
                            * ( cpatch%bsapwooda(ico)                                      &
                              + cpatch%bbarka   (ico) + cpatch%bdeada(ico) )               &
                            + agf_bs(ipft) * cpatch%bstorage(ico)
                     !---------------------------------------------------------------------!



                     !------ Woody living fuel.  AG lignified biomass. --------------------!
                     bwoody = (1. - f_labile_leaf(ipft)) * cpatch%bleaf(ico)               &
                            + (1. - f_labile_stem(ipft))                                   &
                            * ( cpatch%bsapwooda(ico)                                      &
                              + cpatch%bbarka   (ico) + cpatch%bdeada(ico) )
                     !---------------------------------------------------------------------!


                     !------ Accumulate fuels to the patch level. -------------------------!
                     bherb_pat  = bherb_pat  + cpatch%nplant(ico) * bherb
                     bwoody_pat = bwoody_pat + cpatch%nplant(ico) * bwoody
                     !---------------------------------------------------------------------!
                  end if
               end do cohort_loop
               !----- Integrate fuel. -----------------------------------------------------!
               bfuel_d0001_tot = bfuel_d0001_tot + bfuel_d0001_pat * csite%area(ipa)
               bfuel_d0010_tot = bfuel_d0010_tot + bfuel_d0010_pat * csite%area(ipa)
               bfuel_d0100_tot = bfuel_d0100_tot + bfuel_d0100_pat * csite%area(ipa)
               bfuel_d1000_tot = bfuel_d1000_tot + bfuel_d1000_pat * csite%area(ipa)
               bfuel_d0111_tot = bfuel_d0111_tot + bfuel_d0111_pat * csite%area(ipa)
               bherb_tot       = bherb_tot       + bherb_pat       * csite%area(ipa)
               bwoody_tot      = bwoody_tot      + bwoody_pat      * csite%area(ipa)
               !---------------------------------------------------------------------------!



               !----- Integrate site-average temperatures and VPD. ------------------------!
               tdmax_can_temp  = tdmax_can_temp                                            &
                               + maxval(csite%tdfire_can_temp (:,ipa),dim=1)               &
                               * csite%area(ipa)
               today_can_tdew  = today_can_tdew                                            &
                               + sum   (csite%tdfire_can_tdew (:,ipa),dim=1) * ndfirei     &
                               * csite%area(ipa)
               tdmax_can_vpdef = tdmax_can_vpdef                                           &
                               + maxval(csite%tdfire_can_vpdef(:,ipa),dim=1)               &
                               * csite%area(ipa)
               tdmin_can_vpdef = tdmin_can_vpdef                                           &
                               + minval(csite%tdfire_can_vpdef(:,ipa),dim=1)               &
                               * csite%area(ipa)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !      For the fire danger index, we use the average between maximum and    !
               ! minimum.                                                                  !
               !---------------------------------------------------------------------------!
               fdivpd_pat = ndfirei * sum(csite%tdfire_fdi_vpd(:,ipa),dim=1)
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Find herb and woody fuel wetness, using the wetness of the top       !
               ! soil.                                                                     !
               !---------------------------------------------------------------------------!
               moist_bherb_pat  = ndfirei * sum(csite%tdfire_sfc_wetness(:,ipa),dim=1)
               moist_bwoody_pat = ndfirei * sum(csite%tdfire_sfc_wetness(:,ipa),dim=1)
               moist_bherb_pat  = max(0.,min(1.,moist_bherb_pat ))
               moist_bwoody_pat = max(0.,min(1.,moist_bwoody_pat))
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !    Compute fuel moisture.  Decide whether to use the original SPITFIRE    !
               ! approach (based on Nesterov index) or the VPD-based fire danger index.    !
               !---------------------------------------------------------------------------!
               if (fe_use_fdivpd) then
                  !------------------------------------------------------------------------!
                  !       Use the fire danger index to estimate fuel moisture.             !
                  !------------------------------------------------------------------------!
                  moist_bfuel_pat = fe_fdivpd_slp * bpow01(1. - fdivpd_pat,fe_fdivpd_exp)
                  moist_bfuel_pat = max(0.,min(1.,moist_bfuel_pat))
                  !------------------------------------------------------------------------!
               else
                  !------------------------------------------------------------------------!
                  !       Compute fuel moisture, based on SPITFIRE (T10).                  !
                  !------------------------------------------------------------------------!
                  if (bfuel_d0111_pat > tiny_num) then
                     lnexp             = - ( fx_a0001 * bfuel_d0001_pat                    &
                                           + fx_a0010 * bfuel_d0010_pat                    &
                                           + fx_a0100 * bfuel_d0100_pat )                  &
                                           / bfuel_d0111_pat * csite%nesterov_index(ipa)
                     moist_bfuel_pat   = exp(max(lnexp_min,min(lnexp_max,lnexp)))
                  else
                     moist_bfuel_pat   = 1.0
                  end if
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Integrate fuel moisture, scaling by area and fuel stocks in this     !
               ! patch.                                                                    !
               !---------------------------------------------------------------------------!
               moist_bfuel_avg  = moist_bfuel_avg                                          &
                                + moist_bfuel_pat  * bfuel_d0111_pat * csite%area(ipa)
               moist_bherb_avg  = moist_bherb_avg                                          &
                                + moist_bherb_pat  * bherb_pat       * csite%area(ipa)
               moist_bwoody_avg = moist_bwoody_avg                                         &
                                + moist_bwoody_pat * bwoody_pat      * csite%area(ipa)
               !---------------------------------------------------------------------------!


               !------ Integrate average indices. -----------------------------------------!
               nesterov_avg    = nesterov_avg + csite%nesterov_index(ipa) * csite%area(ipa)
               fdivpd_avg      = fdivpd_avg   + fdivpd_pat                * csite%area(ipa)
               !---------------------------------------------------------------------------!

            end do patchloop
            !------------------------------------------------------------------------------!


            !----- Total live/dead fuel loads. --------------------------------------------!
            bfuel_live_tot = bherb_tot + bwoody_tot
            bfuel_dead_tot = bfuel_d0001_tot + bfuel_d0010_tot                             &
                           + bfuel_d0100_tot + bfuel_d1000_tot
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !       Normalise fuel moisture for dead and live components.                  !
            !------------------------------------------------------------------------------!
            if (bfuel_d0111_tot  > tiny_num) then
               moist_bfuel_avg  = moist_bfuel_avg  / bfuel_d0111_tot
            else
               moist_bfuel_avg  = 1.0
            end if
            if (bherb_tot  > tiny_num) then
               moist_bherb_avg  = moist_bherb_avg  / bherb_tot
            else
               moist_bherb_avg  = 1.0
            end if
            if (bwoody_tot > tiny_num) then
               moist_bwoody_avg = moist_bwoody_avg / bwoody_tot
            else
               moist_bwoody_avg = 1.0
            end if
            !------ Apply correction factor for fuel moisture. ----------------------------!
            moist_bherb_avg  = max(0., (1.+fx_rmfac) * moist_bherb_avg  - 1. ) / fx_rmfac
            moist_bwoody_avg = max(0., (1.+fx_rmfac) * moist_bwoody_avg - 1. ) / fx_rmfac
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Find moisture of extinction for fuels and the relative moisture.          !
            !------------------------------------------------------------------------------!
            call find_mextinct(bfuel_d0001_tot,bfuel_d0010_tot,bfuel_d0100_tot             &
                              ,bfuel_d1000_tot,bherb_tot,bwoody_tot                        &
                              ,moist_bfuel_avg,moist_bfuel_avg,moist_bfuel_avg             &
                              ,moist_bfuel_avg,moist_bherb_avg,moist_bwoody_avg,Mx_i)
            if ( (bfuel_d0111_tot + bfuel_live_tot) > tiny_num ) then
               rmoist_avg = ( moist_bfuel_avg  * bfuel_d0111_tot                           &
                            + moist_bherb_avg  * bherb_tot                                 &
                            + moist_bwoody_avg * bwoody_tot           )                    &
                          / ( Mx_i(1) * bfuel_d0111_tot + Mx_i(2) * bfuel_live_tot )
            else
               rmoist_avg = 1.0
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !       Use T10's spread probability as a proxy for the probability that       !
            ! fire ignitions could turn into spreading fire.  Depending on the user        !
            ! settings,we assume no ignitions in case the polygon does not have any        !
            ! existing signs of land use.                                                  !
            !------------------------------------------------------------------------------!
            if (people_around) then
               ignition_rate = max(0.,1. - rmoist_avg )
            else
               ignition_rate = 0.
            end if
            !------------------------------------------------------------------------------!


            !----- Compute fire intensity. ------------------------------------------------!
            cpoly%fire_intensity(isi) = fire_parameter * ignition_rate                     &
                                      * ( bfuel_live_tot  + bfuel_dead_tot )
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !       Update the average fire intensity.                                     !
            !------------------------------------------------------------------------------!
            cpoly%avg_fire_intensity(imonth,isi) = cpoly%avg_fire_intensity(imonth,isi)    &
                                                 + cpoly%fire_intensity           (isi)    &
                                                 * ndaysi
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Print the output if needed.                                              !
            !------------------------------------------------------------------------------!
            if (printout) then
               !------ Convert units for output. ------------------------------------------!
               today_pcpg      = today_pcpg      * day_sec
               tdmax_can_temp  = tdmax_can_temp  - t00
               today_can_tdew  = today_can_tdew  - t00
               tdmax_can_vpdef = tdmax_can_vpdef * 100.
               tdmin_can_vpdef = tdmin_can_vpdef * 100.
               !---------------------------------------------------------------------------!


               open(unit=35,file=firefile,status='old',position='append',action='write')
               write(unit=35,fmt='(4(i6,1x),1(5x,l1,1x),18(f12.4,1x))')                    &
                          current_time%year,current_time%month,current_time%date,isi       &
                         ,people_around,today_pcpg,tdmax_can_temp,today_can_tdew           &
                         ,tdmax_can_vpdef,tdmin_can_vpdef,nesterov_avg,fdivpd_avg          &
                         ,bherb_tot,bwoody_tot,bfuel_dead_tot,moist_bherb_avg              &
                         ,moist_bwoody_avg,moist_bfuel_avg,Mx_i(1),Mx_i(2),rmoist_avg      &
                         ,ignition_rate,cpoly%fire_intensity(isi)
               close(unit=35,status='keep')
            end if
            !------------------------------------------------------------------------------!
         end do main_siteloop
         !---------------------------------------------------------------------------------!
      end do main_polyloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine integ_emberfire
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   ! SUB-ROUTINE INTEG_FIRESTARTER
   !\brief Main integrator of the FIRESTARTER model
   !\details This sub-routine that integrates the fire disturbance rate when using the 
   !!        "Fire Ignition, Rate of Elliptical Spread, and Termination Approaches to
   !!        Represent the Terrestrial Ecosystem Responses (to fires)" (FIRESTARTER) model.
   !!        The FIRESTARTER model builds on the HESFIRE model (LP15/LP17) for ignitions 
   !!        and termination, and the SPITFIRE (T10) for ecosystem response to fires.  The
   !!        calculation of maximum rate of spread is based on R72 model revised by A18
   !!        and using the very dry fuel conditions described in SB05.
   !!
   !> References:
   !!
   !!  Andrews PL. 2018. The Rothermel surface fire spread model and associated
   !!     developments: A compre- hensive explanation. Gen. Tech. Rep. RMRS-GTR-371, U.S.
   !!     Department of Agriculture, Forest Service, Rocky Mountain Research Station, Fort
   !!    Collins, CO, U.S.A. https://www.fs.usda.gov/treesearch/pubs/55928 (A18).
   !!
   !! Le Page Y, Morton D, Bond-Lamberty B, Pereira JMC , Hurtt G. 2015. HESFIRE: a global
   !!    fire model to explore the role of anthropogenic and weather drivers.
   !!    Biogeosciences, 12: 887-903. doi:10.5194/bg-12-887-2015 (LP15).
   !!
   !! Le Page Y, Morton D, Hartin C, Bond-Lamberty B, Pereira JMC, Hurtt G , Asrar G. 2017.
   !!    Synergy between land use and climate change increases future fire risk in Amazon
   !!    forests. Earth Syst. Dynam., 8: 1237-1246. doi:10.5194/esd-8-1237-2017 (LP17).
   !!
   !! Peterson, D. L., and K. C. Ryan, 1986: Modeling postfire conifer mortality for 
   !!    long-range planning. Environ. Manage., 10 (6), 797-808, doi:10.1007/BF01867732
   !!    (PR86).
   !!
   !! Rothermel RC. 1972. A mathematical model for predicting fire spread in wildland
   !!    fuels. Res. Pap. INT- 115, U.S. Department of Agriculture, Intermountain Forest
   !!    and Range Experiment Station, Ogden, UT, U. S. A.,
   !!    https://www.fs.usda.gov/treesearch/pubs/32533 (R72).
   !!
   !! Schaphoff S, von Bloh W, Rammig A, Thonicke K, Biemans H, Forkel M, Gerten D, 
   !!    Heinke J, Jagermeyr J, Knauer J et al. 2018. LPJmL4 -- a dynamic global 
   !!    vegetation model with managed land -- part 1: Model description. 
   !!    Geosci. Model Dev., 11: 13431375. doi:10.5194/gmd-11-1343-2018 (S18).
   !!
   !! Scott JH , Burgan RE. 2005. Standard fire behavior fuel models: a comprehensive set
   !!    for use with Rothermel's surface fire spread model. Gen. Tech. Rep. RMRS-GTR-153,
   !!    U.S. Department of Agriculture, Forest Service, Rocky Mountain Research Station,
   !!    Fort Collins, CO, U.S.A. doi:10.2737/RMRS-GTR-153 (SB05).
   !!
   !! Thonicke K, Spessa A, Prentice IC, Harrison SP, Dong L, Carmona-Moreno C. 2010. The
   !!    influence of vegetation, fire spread and fire behaviour on biomass burning and
   !!    trace gas emissions: results from a process-based model. Biogeosciences, 7:
   !!    1991-2011. doi:10.5194/bg-7-1991-2010 (T10).
   !!
   !---------------------------------------------------------------------------------------!
   subroutine integ_firestarter(cgrid,dtfull)
      use ed_state_vars , only : edtype                 & ! structure
                               , polygontype            & ! structure
                               , sitetype               & ! structure
                               , patchtype              ! ! structure
      use ed_misc_coms  , only : simtime                & ! structure
                               , current_time           & ! intent(in)
                               , ndfire                 ! ! intent(in)
      use disturb_coms  , only : fe_fdivpd_exp          & ! intent(in)
                               , fe_fdivpd_slp          & ! intent(in)
                               , fe_use_fdivpd          & ! intent(in)
                               , fh_f0001               & ! intent(in)
                               , fh_f0010               & ! intent(in)
                               , fh_f0100               & ! intent(in)
                               , fh_f1000               & ! intent(in)
                               , fh_grid                & ! intent(in)
                               , fi_cg_ignp             & ! intent(in)
                               , fi_hdi_exp             & ! intent(in)
                               , fi_hdi_dti             & ! intent(in)
                               , fi_hdi_lwr             & ! intent(in)
                               , fi_lu_ignd             & ! intent(in)
                               , fi_lu_exp              & ! intent(in)
                               , fi_lu_upr              & ! intent(in)
                               , fi_sf_maxage           & ! intent(in)
                               , fr_h                   & ! intent(in)
                               , fr_Mxdead              & ! intent(in)
                               , fs_bck_exp             & ! intent(in)
                               , fs_gw_infty            & ! intent(in)
                               , fs_gw_upr              & ! intent(in)
                               , fs_lbr_exp             & ! intent(in)
                               , fs_lbr_slp             & ! intent(in)
                               , ft_fint_dti            & ! intent(in)
                               , ft_fint_exp            & ! intent(in)
                               , ft_fint_lwr            & ! intent(in)
                               , ft_frag_exp            & ! intent(in)
                               , ft_lu_exp              & ! intent(in)
                               , ft_lu_upr              & ! intent(in)
                               , ft_hdi_exp             & ! intent(in)
                               , ft_hdi_dti             & ! intent(in)
                               , ft_hdi_lwr             & ! intent(in)
                               , ft_fdi_upr             & ! intent(in)
                               , ft_fdi_exp             & ! intent(in)
                               , fuel_height_max        & ! intent(in)
                               , fx_a0001               & ! intent(in)
                               , fx_a0010               & ! intent(in)
                               , fx_a0100               & ! intent(in)
                               , fx_tlh_slope           & ! intent(in)
                               , fx_rmfac               & ! intent(in)
                               , n_fst                  ! ! intent(in)
      use pft_coms      , only : agf_bs                 & ! intent(in)
                               , C2B                    & ! intent(in)
                               , f_labile_leaf          & ! intent(in)
                               , f_labile_stem          & ! intent(in)
                               , is_grass               ! ! intent(in)
      use consts_coms   , only : pio4                   & ! intent(in)
                               , tiny_num               & ! intent(in)
                               , almost_zero            & ! intent(in)
                               , almost_one             & ! intent(in)
                               , lnexp_min              & ! intent(in)
                               , lnexp_max              ! ! intent(in)
      implicit none
      !----- -Arguments. ------------------------------------------------------------------!
      type(edtype)  , target     :: cgrid             ! Current grid              [    ---]
      real          , intent(in) :: dtfull            ! Fire full time step       [      s]
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
      real, dimension(n_fst)     :: Mx_i              ! Fuel moist. extinct.      [  kg/kg]
      real                       :: anth_ign_rate     ! Anthropogenic ignt. rate  [ 1/m2/s]
      real                       :: apy_area          ! Grid area                 [     m2]
      real                       :: bherb             ! Cohort Herbaceous fuels   [ kgC/pl]
      real                       :: bherb_pat         ! Patch Herbaceous fuels    [ kgC/m2]
      real                       :: bherb_tot         ! Site Herbaceous fuels     [ kgC/m2]
      real                       :: bfuel_all_pat     ! Patch Total fuels         [ kgC/m2]
      real                       :: bfuel_all_tot     ! Site Total fuels          [ kgC/m2]
      real                       :: bfuel_d0001_pat   ! Patch 1-hr dead fuels     [ kgC/m2]
      real                       :: bfuel_d0001_tot   ! Site 1-hr dead fuels      [ kgC/m2]
      real                       :: bfuel_d0010_pat   ! Patch 10-hr dead fuels    [ kgC/m2]
      real                       :: bfuel_d0010_tot   ! Site 10-hr dead fuels     [ kgC/m2]
      real                       :: bfuel_d0100_pat   ! Patch 100-hr dead fuels   [ kgC/m2]
      real                       :: bfuel_d0100_tot   ! Site 100-hr dead fuels    [ kgC/m2]
      real                       :: bfuel_d1000_pat   ! Patch 1000-hr dead fuels  [ kgC/m2]
      real                       :: bfuel_d1000_tot   ! Site 1000-hr dead fuels   [ kgC/m2]
      real                       :: bfuel_d0111_pat   ! (1+10+100-hr fuels)       [ kgC/m2]
      real                       :: bfuel_d0111_tot   ! (1+10+100-hr fuels)       [ kgC/m2]
      real                       :: burnt_area_deja   ! Burnt area so far         [  m2/m2]
      real                       :: burnt_area_potl   ! Potl. Burnt area (step)   [  m2/m2]
      real                       :: burnt_area_step   ! Burnt area (time step)    [  m2/m2]
      real                       :: burnt_area_max    ! Maximum burnt area        [  m2/m2]
      real                       :: bwoody            ! Cohort living woody fuels [ kgC/pl]
      real                       :: bwoody_pat        ! Patch living woody fuels  [ kgC/m2]
      real                       :: bwoody_tot        ! Site living woody fuels   [ kgC/m2]
      real                       :: bwn1000           ! Cohort 1-100hr woody fuel [ kgC/pl]
      real                       :: bwn1000_pat       ! Patch 1-100hr woody fuel  [ kgC/m2]
      real                       :: bwn1000_tot       ! Site 1-100hr woody fuels  [ kgC/m2]
      real                       :: can_vels          ! Canopy air velocity       [    m/s]
      real                       :: ell_length        ! Length of main ell. axis  [      m]
      real                       :: fdi_pat           ! Patch fire danger index   [     --]
      real                       :: fdin              ! Norm. fire danger index   [     --]
      real                       :: fdivpd_pat        ! Patch fire danger index   [     --]
      real                       :: fdivpd_avg        ! Average fire danger index [     --]
      real                       :: fintn             ! Norm. fire intensity      [    ---]
      real                       :: fragn             ! Norm. fragmentation       [    ---]
      real                       :: flamn             ! Flammable area            [    ---]
      real                       :: fp_anth_fun       ! LU/HDI-mediated persist.  [     --]
      real                       :: fp_cntg_fun       ! Contiguity factor         [    ---]
      real                       :: fp_fuel_fun       ! Fuel availability factor  [    ---]
      real                       :: fp_fdi_fun        ! FDI control function      [    ---]
      real                       :: fp_fdi_loc        ! FDI control function      [    ---]
      real                       :: fp_wild_fun       ! Wildfire risk function    [    ---]
      real                       :: fp_wind_fun       ! Wind control function     [    ---]
      real                       :: fp_wind_loc       ! Wind control function     [    ---]
      real                       :: fs_iarea_pat      ! Individual fire area      [     m2]
      real                       :: fs_iarea_avg      ! Individual fire area      [     m2]
      real                       :: fx_b0001          ! Fuel consumption 1-h      [ kgC/m2]
      real                       :: fx_b0001_potl     ! Potl. fuel cons. 1-h      [ kgC/m2]
      real                       :: fx_b0010          ! Fuel consumption 10-h     [ kgC/m2]
      real                       :: fx_b0010_potl     ! Potl. fuel cons. 10-h     [ kgC/m2]
      real                       :: fx_b0100          ! Fuel consumption 100-h    [ kgC/m2]
      real                       :: fx_b0100_potl     ! Potl. fuel cons. 100-h    [ kgC/m2]
      real                       :: fx_b1000          ! Fuel consumption 1000-h   [ kgC/m2]
      real                       :: fx_b1000_potl     ! Potl. fuel cons. 1000-h   [ kgC/m2]
      real                       :: fx_bherb          ! Fuel consumption herb     [ kgC/m2]
      real                       :: fx_bherb_potl     ! Potl. fuel cons. herb     [ kgC/m2]
      real                       :: fx_bwoody         ! Fuel consumption woody    [ kgC/m2]
      real                       :: fx_bwoody_potl    ! Potl. fuel cons. woody    [ kgC/m2]
      real                       :: fx_duration       ! Fire duration correction  [    ---]
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
      real                       :: fx_wn1000_potl    ! Potl. F. C. woody-1000h   [ kgC/m2]
      real                       :: g_Umax            ! Maximum wind for ROS      [    m/s]
      real                       :: gw_factor         ! Wind speed effect on ROS  [    m/s]
      real                       :: hb_ratio          ! Head:back ratio           [    ---]
      real                       :: hdin              ! Norm. human develop. idx  [    ---]
      real                       :: lb_ratio          ! Length:breadth ratio      [    ---]
      real                       :: lnexp             ! Aux. var. for safe exp    [    ---]
      real                       :: lu_area           ! LU area                   [    ---]
      real                       :: lu_effect         ! LU effect on ignition     [    ---]
      real                       :: lu_norm           ! Norm. LU area             [    ---]
      real                       :: moist_bfuel_avg   ! Dead fuel moisture        [    ---]
      real                       :: moist_bfuel_pat   ! Dead fuel moisture        [    ---]
      real                       :: moist_bherb_avg   ! Herbaceous fuel moisture  [    ---]
      real                       :: moist_bherb_pat   ! Herbaceous fuel moisture  [    ---]
      real                       :: moist_bwoody_avg  ! Living woody fuel moist.  [    ---]
      real                       :: moist_bwoody_pat  ! Living woody fuel moist.  [    ---]
      real                       :: nat_ign_rate      ! Natural ignition rate     [ 1/m2/s]
      real                       :: ndaysi            ! 1/# days in a month       [  1/day]
      real                       :: nesterov_avg      ! Average Nesterov index    [  degC2]
      real                       :: prob_persist      ! Persistence probability   [    ---]
      real                       :: rmoist_b0001      ! 1-hr dead rel. moisture   [    ---]
      real                       :: rmoist_b0010      ! 10-hr dead rel. moisture  [    ---]
      real                       :: rmoist_b0100      ! 100-hr dead rel. moisture [    ---]
      real                       :: rmoist_b1000      ! 1000-hr dead rel. mst.    [    ---]
      real                       :: rmoist_bherb      ! Herbaceous rel. moisture  [    ---]
      real                       :: rmoist_bwoody     ! Living woody rel. moist.  [    ---]
      real                       :: rosbwd            ! Backward rate of spread   [    m/s]
      real                       :: rosfwd            ! Forward rate of spread    [    m/s]
      real                       :: rosbwd_avg        ! Site-avg bwd. spread rate [    m/s]
      real                       :: rosfwd_avg        ! Site-avg fwd. spread rate [    m/s]
      real                       :: sfc_wetness       ! Sfc. soil wetness         [    ---]
      real                       :: suppressibility   ! Fire suppressibility      [    ---]
      real                       :: total_ignition    ! Number of ignitions       [   1/m2]
      !------ External functions. ---------------------------------------------------------!
      real              , external  :: solid_area     ! Solid-angle area          [     m2]
      real              , external  :: bpow01         ! Power funct. for [0-1]    [    ---]
      real              , external  :: cbrt           ! Cube root                 [    ---]
      !----- Local parameters. ------------------------------------------------------------!
      character(len=23) , parameter :: firefile = 'firestarter_details.txt'
      logical           , parameter :: printout = .false.
      !----- Locally saved variables. -----------------------------------------------------!
      logical           , save      :: first_time = .true. ! First time calling   [    T|F]
      real              , save      :: ndfirei             ! 1./ndfire            [     --]
      real              , save      :: dtfire              ! Fire time step       [      s]
      !------------------------------------------------------------------------------------!


      !----- First time, and the user wants to print the output.  Make a header. ----------!
      if (first_time) then

         !----- Make the header. ----------------------------------------------------------!
         if (printout) then
            open (unit=35,file=firefile,status='replace',action='write')
            write (unit=35,fmt='(37(a,1x))')                                               &
                     '  YEAR',      ' MONTH',      '   DAY',      '  STEP',      '   ISI'  &
              ,'    APY_AREA','    LANDFRAC','       FRAGN','     LU_AREA','         HDI'  &
              ,'   C2G_FLASH','TOT_IGNITION','    NESTEROV','     FDI_VPD',' MOIST_BHERB'  &
              ,'MOIST_BWOODY',' MOIST_BFUEL','      ROSFWD','  FCOMB_FAST','FCOMB_STRUCT'  &
              ,' FCOMB_BHERB','FCOMB_BWOODY',' FX_DURATION','FX_INTENSITY','  FX_TLETHAL'  &
              ,'  FP_FDI_FUN',' FP_WIND_FUN','  SUPPRESSIB',' FP_ANTH_FUN',' FP_WILD_FUN'  &
              ,' FP_FUEL_FUN',' FP_CNTG_FUN','PROB_PERSIST','  BAREA_STEP','  BURNT_AREA'  &
              ,'FIRE_DENSITY','FIRE_EXTINCT'
            close (unit=35,status='keep')
         end if
         !---------------------------------------------------------------------------------!


         !------ Set inverse of number of sub-daily bins/steps. ---------------------------!
         ndfirei    = 1. / real(ndfire)
         !---------------------------------------------------------------------------------!


         !----- Set fire time step. -------------------------------------------------------!
         dtfire = ndfirei * dtfull
         !---------------------------------------------------------------------------------!

         first_time = .false.
      end if
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
               elseif (cpatch%ncohorts == 0) then
                  !----- 2. Deserts, exclude the area. ------------------------------------!
                  fragn = fragn + csite%area(ipa)
                  !------------------------------------------------------------------------!
               else
                  !----- 3. The patch may burn, but we exclude flooded/snowpack fraction. -!
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
            !     The maximum are that can burn is scaled by the landscape fragmentation   !
            ! (including the permanent fragmentation such as lakes, oceans, and glaciers). !
            !------------------------------------------------------------------------------!
            burnt_area_max = cgrid%landfrac(ipy) * (1. - fragn)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Land use area relative to the total grid area (including permanent       !
            ! areas not solved by ED2 such as lakes, oceans, and glaciers).                !
            !------------------------------------------------------------------------------!
            lu_area       = lu_area * cgrid%landfrac(ipy)
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
            ! Equation 3, without multiplying by fi_lu_upr / (1. + fi_lu_exp), so the      !
            ! land use effects reaches 1 at the maximum land use effect.                   !
            ! We account for the fraction of the polygon that cannot sustain a fire as     !
            ! LP15 also did.  We cap the land use area to the maximum land use area that   !
            ! contributes to anthropogenic ignitions.                                      !
            !------------------------------------------------------------------------------!
            lu_norm       = max(0., min(1.,lu_area / fi_lu_upr))
            lu_effect     = max(0., (1. - bpow01(1.- lu_norm, fi_lu_exp + 1. ) ) )
            hdin          = ( cpoly%seitimes(isei,isi)%hdi - fi_hdi_lwr ) * fi_hdi_dti
            hdin          = max(0., min(1.,hdin) )
            anth_ign_rate = (1. - bpow01(hdin,fi_hdi_exp)) * fi_lu_ignd * lu_effect
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !      Find total number of ignitions (1/m2).  This is applied twice inside    !
            ! the model time step, because the ignitions do not change during the same     !
            ! day.                                                                         !
            !------------------------------------------------------------------------------!
            cpoly%ignition_rate(isi) = nat_ign_rate + anth_ign_rate
            total_ignition           = cpoly%ignition_rate(isi) * dtfire
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
            !       Run the model for every sub-daily step.                                !
            !------------------------------------------------------------------------------!
            timestep_loop: do iwhen=1,ndfire
 
 

               !---------------------------------------------------------------------------!
               !      Initialise the functions that control fire spread and termination.   !
               !---------------------------------------------------------------------------!
               !----- Functions shared by fire spread. ------------------------------------!
               fs_iarea_avg = 0.
               !----- Functions used for fire persistence. --------------------------------!
               fp_fdi_fun   = 0.
               fp_wind_fun  = 0.
               !---------------------------------------------------------------------------!


               !------ Initialise fuel stocks. --------------------------------------------!
               bfuel_d0001_tot = 0.
               bfuel_d0010_tot = 0.
               bfuel_d0100_tot = 0.
               bfuel_d1000_tot = 0.
               bfuel_d0111_tot = 0.
               bherb_tot       = 0.
               bwoody_tot      = 0.
               bwn1000_tot     = 0.
               bfuel_all_tot   = 0.
               !---------------------------------------------------------------------------!


               !------ Initialise dead and live fuel moisture. ----------------------------!
               moist_bfuel_avg    = 0.
               moist_bherb_avg    = 0.
               moist_bwoody_avg   = 0.
               !---------------------------------------------------------------------------!


               !------ Initialise the average rate of spread (forward and backward). ------!
               rosfwd_avg     = 0.
               rosbwd_avg     = 0.
               !---------------------------------------------------------------------------!


               !------ Scale burnt area up to now by removing permanent fragmentation. ----!
               burnt_area_deja = cpoly%burnt_area  (isi) * cgrid%landfrac(ipy)
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Update fire density.  Add fire ignitions, and apply last step's       !
               ! extinction rate.                                                          !
               !---------------------------------------------------------------------------!
               prob_persist = exp(max( lnexp_min                                           &
                                     ,min(lnexp_max,-cpoly%fire_extinction(isi)*dtfire) ) )
               if (prob_persist < almost_zero) then
                  cpoly%fire_density(isi) = total_ignition
               else
                  cpoly%fire_density(isi) = cpoly%fire_density(isi) * prob_persist         &
                                          + total_ignition
               end if
               !---------------------------------------------------------------------------!


               !------ Initialise average indices. ----------------------------------------!
               nesterov_avg    = 0.
               fdivpd_avg      = 0.
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Loop through patches.                                                 !
               !---------------------------------------------------------------------------!
               fst_patch_loop: do ipa=1,csite%npatches
                  cpatch => csite%patch(ipa)

                  !------------------------------------------------------------------------!
                  !    Select patch-level variables for the current step.                  !
                  !------------------------------------------------------------------------!
                  fdivpd_pat  = csite%tdfire_fdi_vpd    (iwhen,ipa)
                  can_vels    = csite%tdfire_can_vels   (iwhen,ipa)
                  sfc_wetness = csite%tdfire_sfc_wetness(iwhen,ipa)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !       Allocate fuels.                                                  !
                  !------------------------------------------------------------------------!
                  bfuel_d0001_pat = csite%fast_grnd_C(ipa)                                 &
                                  + fh_f0001 * csite%structural_grnd_C(ipa)
                  bfuel_d0010_pat = fh_f0010 * csite%structural_grnd_C(ipa)
                  bfuel_d0100_pat = fh_f0100 * csite%structural_grnd_C(ipa)
                  bfuel_d1000_pat = fh_f1000 * csite%structural_grnd_C(ipa)
                  bfuel_d0111_pat = bfuel_d0100_pat + bfuel_d0010_pat + bfuel_d0001_pat
                  bherb_pat       = 0.
                  bwoody_pat      = 0.
                  bwn1000_pat     = 0.
                  spread_cohort_loop: do ico=1,cpatch%ncohorts
                     ipft = cpatch%pft(ico)
                     if (is_grass(ipft) .or. cpatch%hite(ico) <= fuel_height_max) then
                        !------ Herbaceous fuel.  AG labile biomass + AG storage. ---------!
                        bherb  = f_labile_leaf(ipft) * cpatch%bleaf(ico)                   &
                               + f_labile_stem(ipft)                                       &
                               * ( cpatch%bsapwooda(ico)                                   &
                                 + cpatch%bbarka   (ico) + cpatch%bdeada(ico) )            &
                               + agf_bs(ipft) * cpatch%bstorage(ico)
                        !------------------------------------------------------------------!



                        !------ Woody living fuel.  AG lignified biomass. -----------------!
                        bwoody = (1. - f_labile_leaf(ipft)) * cpatch%bleaf(ico)            &
                               + (1. - f_labile_stem(ipft))                                &
                               * ( cpatch%bsapwooda(ico)                                   &
                                 + cpatch%bbarka   (ico) + cpatch%bdeada(ico) )
                        !------------------------------------------------------------------!


                        !----- Woody living fuel (excluding 1000-hr fuels). ---------------!
                        bwn1000 = bwoody * (1. - fh_f1000)
                        !------------------------------------------------------------------!


                        !------ Accumulate fuels to the patch level. ----------------------!
                        bherb_pat   = bherb_pat   + cpatch%nplant(ico) * bherb
                        bwoody_pat  = bwoody_pat  + cpatch%nplant(ico) * bwoody
                        bwn1000_pat = bwn1000_pat + cpatch%nplant(ico) * bwn1000
                        !------------------------------------------------------------------!
                     end if
                  end do spread_cohort_loop
                  !------------------------------------------------------------------------!


                  !------ Find total fuel loads. ------------------------------------------!
                  bfuel_all_pat = bfuel_d0111_pat + bherb_pat + bwn1000_pat
                  !------------------------------------------------------------------------!



                  !------ Integrate site-level fuel loads. --------------------------------!
                  bfuel_d0001_tot = bfuel_d0001_tot + bfuel_d0001_pat * csite%area(ipa)
                  bfuel_d0010_tot = bfuel_d0010_tot + bfuel_d0010_pat * csite%area(ipa)
                  bfuel_d0100_tot = bfuel_d0100_tot + bfuel_d0100_pat * csite%area(ipa)
                  bfuel_d1000_tot = bfuel_d1000_tot + bfuel_d1000_pat * csite%area(ipa)
                  bfuel_d0111_tot = bfuel_d0111_tot + bfuel_d0111_pat * csite%area(ipa)
                  bherb_tot       = bherb_tot       + bherb_pat       * csite%area(ipa)
                  bwoody_tot      = bwoody_tot      + bwoody_pat      * csite%area(ipa)
                  bwn1000_tot     = bwn1000_tot     + bwn1000_pat     * csite%area(ipa)
                  bfuel_all_tot   = bfuel_all_tot   + bfuel_all_pat   * csite%area(ipa)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !    Compute fuel moisture.  Decide whether to use the original SPITFIRE !
                  ! approach (based on Nesterov index) or the VPD-based fire danger index. !
                  !------------------------------------------------------------------------!
                  if (fe_use_fdivpd) then
                     !---------------------------------------------------------------------!
                     !       Use the fire danger index to estimate fuel moisture.          !
                     !---------------------------------------------------------------------!
                     moist_bfuel_pat = fe_fdivpd_slp                                       &
                                     * bpow01(1. - fdivpd_pat,fe_fdivpd_exp)
                     moist_bfuel_pat = max(0.,min(1.,moist_bfuel_pat))
                     fdi_pat         = fdivpd_pat
                     !---------------------------------------------------------------------!
                  else
                     !---------------------------------------------------------------------!
                     !       Compute fuel moisture, based on SPITFIRE (T10).               !
                     !---------------------------------------------------------------------!
                     if (bfuel_d0111_pat > tiny_num) then
                        lnexp             = - ( fx_a0001 * bfuel_d0001_pat                 &
                                              + fx_a0010 * bfuel_d0010_pat                 &
                                              + fx_a0100 * bfuel_d0100_pat )               &
                                              / bfuel_d0111_pat * csite%nesterov_index(ipa)
                        moist_bfuel_pat   = exp(max(lnexp_min,min(lnexp_max,lnexp)))
                     else
                        moist_bfuel_pat   = 1.0
                     end if
                     fdi_pat              = max(0.,1.-moist_bfuel_pat/fr_Mxdead)
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Integrate persistence functions for fire danger index (using an   !
                  ! approach similar to HESFIRE (LP15).   This substitutes the more        !
                  ! convoluted function of soil matric potential, relative humidity and    !
                  ! temperature with a single, normalised metric.                          !
                  !------------------------------------------------------------------------!
                  fdin       = min(1.,fdi_pat / ft_fdi_upr)
                  fp_fdi_loc = bpow01(fdin   ,ft_fdi_exp )
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !      Find the wind influence function on termination.                  !
                  !------------------------------------------------------------------------!
                  lnexp       = max( lnexp_min, min( lnexp_max, fs_lbr_exp * can_vels) )
                  lb_ratio    = 1. + fs_lbr_slp * (1. - exp(lnexp))
                  hb_ratio    = ( lb_ratio + sqrt( lb_ratio * lb_ratio - 1. ) )            &
                              / ( lb_ratio - sqrt( lb_ratio * lb_ratio - 1. ) )
                  gw_factor   = 2. * lb_ratio / ( 1. + 1. / hb_ratio ) * fs_gw_infty
                  fp_wind_loc = max( 0., min(1., gw_factor / fs_gw_upr ) )
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !      Integrate persistence functions.                                  !
                  !------------------------------------------------------------------------!
                  fp_fdi_fun   = fp_fdi_fun   + fp_fdi_loc   * csite%area(ipa)
                  fp_wind_fun  = fp_wind_fun  + fp_wind_loc  * csite%area(ipa)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Find herb and woody fuel wetness, using the wetness of the top    !
                  ! soil and T10 correction factor.                                        !
                  !------------------------------------------------------------------------!
                  moist_bherb_pat  = max(0., (1.+fx_rmfac) * sfc_wetness - 1.) / fx_rmfac
                  moist_bwoody_pat = max(0., (1.+fx_rmfac) * sfc_wetness - 1.) / fx_rmfac
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Integrate fuel moisture, scaling by area and fuel stocks in this  !
                  ! patch.                                                                 !
                  !------------------------------------------------------------------------!
                  moist_bfuel_avg  = moist_bfuel_avg                                       &
                                   + moist_bfuel_pat  * bfuel_d0111_pat * csite%area(ipa)
                  moist_bherb_avg  = moist_bherb_avg                                       &
                                   + moist_bherb_pat  * bherb_pat       * csite%area(ipa)
                  moist_bwoody_avg = moist_bwoody_avg                                      &
                                   + moist_bwoody_pat * bwn1000_pat     * csite%area(ipa)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !       Find the fire rates of spread using the R72 model modified by    !
                  ! A18 (forward) and the SPITFIRE parametrisation (backward).             !
                  !------------------------------------------------------------------------!
                  !------ Forward. --------------------------------------------------------!
                  call rate_of_spread(isi,bfuel_d0001_pat,bfuel_d0010_pat,bfuel_d0100_pat  &
                                     ,0.*bfuel_d1000_pat,bherb_pat,bwn1000_pat             &
                                     ,moist_bfuel_pat,moist_bfuel_pat,moist_bfuel_pat      &
                                     ,moist_bfuel_pat,moist_bherb_pat,moist_bwoody_pat     &
                                     ,can_vels,.false.,g_Umax,rosfwd     )
                  !------ Backward. -------------------------------------------------------!
                  lnexp  = max( lnexp_min, min( lnexp_max, fs_bck_exp * can_vels ) )
                  rosbwd = rosfwd * exp(lnexp)
                  !------------------------------------------------------------------------!


                  !------ Integrate the rate of spread (use kinetic energy for average). --!
                  rosfwd_avg = rosfwd_avg + rosfwd * rosfwd * csite%area(ipa)
                  rosbwd_avg = rosbwd_avg + rosbwd * rosbwd * csite%area(ipa)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !      Find the length of the major axis, following S18.                 !
                  !------------------------------------------------------------------------!
                  ell_length   = ( rosfwd + rosbwd ) * dtfire
                  fs_iarea_pat = pio4 * rosfwd * rosfwd * dtfire * dtfire / lb_ratio       &
                               * ( 1. + 1./hb_ratio ) * ( 1. + 1./hb_ratio )
                  fs_iarea_avg = fs_iarea_avg + fs_iarea_pat * csite%area(ipa)
                  !------------------------------------------------------------------------!


                  !------ Integrate average indices. --------------------------------------!
                  nesterov_avg    = nesterov_avg                                           &
                                  + csite%nesterov_index(ipa) * csite%area(ipa)
                  fdivpd_avg      = fdivpd_avg                                             &
                                  + fdivpd_pat                * csite%area(ipa)
                  !------------------------------------------------------------------------!


               end do fst_patch_loop
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !       Normalise fuel moisture for dead and live components.               !
               !---------------------------------------------------------------------------!
               if (bfuel_d0111_tot  > tiny_num) then
                  moist_bfuel_avg  = moist_bfuel_avg  / bfuel_d0111_tot
               else
                  moist_bfuel_avg  = 1.0
               end if
               if (bherb_tot  > tiny_num) then
                  moist_bherb_avg  = moist_bherb_avg  / bherb_tot
               else
                  moist_bherb_avg  = 1.0
               end if
               if (bwn1000_tot > tiny_num) then
                  moist_bwoody_avg = moist_bwoody_avg / bwn1000_tot
               else
                  moist_bwoody_avg = 1.0
               end if
               !------ Ensure fuel moisture is bounded. -----------------------------------!
               moist_bfuel_avg  = max(0., min(1., moist_bfuel_avg ))
               moist_bherb_avg  = max(0., min(1., moist_bherb_avg ))
               moist_bwoody_avg = max(0., min(1., moist_bwoody_avg))
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !      Make sure the sum is bounded (it could go off the 0-1 interval due   !
               ! to truncation errors).                                                    !
               !---------------------------------------------------------------------------!
               fp_fdi_fun  = max(0.,min(1.,fp_fdi_fun ))
               fp_wind_fun = max(0.,min(1.,fp_wind_fun))
               !---------------------------------------------------------------------------!


               !------ Normalise rate of spread (square root is needed). ------------------!
               rosfwd_avg = sqrt(rosfwd_avg)
               rosbwd_avg = sqrt(rosbwd_avg)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Find the potential burnt area for this step (potential means if fires !
               ! persist the entire time step).                                            !
               !---------------------------------------------------------------------------!
               if (cgrid%landfrac(ipy) > tiny_num) then
                  !------ Increment burnt area until it is saturated. ---------------------!
                  burnt_area_potl = cpoly%fire_density(isi) * fs_iarea_avg
                  burnt_area_potl = max( 0., min( burnt_area_max - burnt_area_deja         &
                                                , burnt_area_potl                    ) )
                  !------------------------------------------------------------------------!
               else
                  !------ No land to burn, set burnt area to zero... ----------------------!
                  burnt_area_potl = 0.0
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !    Find moisture of extinction for fuels and the relative moisture.       !
               !---------------------------------------------------------------------------!
               call find_mextinct(bfuel_d0001_tot,bfuel_d0010_tot,bfuel_d0100_tot          &
                                 ,bfuel_d1000_tot,bherb_tot,bwn1000_tot                    &
                                 ,moist_bfuel_avg,moist_bfuel_avg,moist_bfuel_avg          &
                                 ,moist_bfuel_avg,moist_bherb_avg,moist_bwoody_avg,Mx_i)
               rmoist_b0001  = moist_bfuel_avg  / Mx_i(1)
               rmoist_b0010  = moist_bfuel_avg  / Mx_i(1)
               rmoist_b0100  = moist_bfuel_avg  / Mx_i(1)
               rmoist_b1000  = moist_bfuel_avg  / Mx_i(1)
               rmoist_bherb  = moist_bherb_avg  / Mx_i(2)
               rmoist_bwoody = moist_bwoody_avg / Mx_i(2)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !      Find potential fuel consumption.                                     !
               !---------------------------------------------------------------------------!
               !----- Find the fuel consumption factors for all fuel classes. -------------!
               call find_fx_factors(rmoist_bherb,rmoist_bwoody,rmoist_b0001,rmoist_b0010   &
                                   ,rmoist_b0100,rmoist_b1000,fx_f_bherb,fx_f_bwoody       &
                                   ,fx_f_wn1000,fx_f_b0001,fx_f_b0010,fx_f_b0100           &
                                   ,fx_f_b1000)
               !----- Find fuel consumption. ----------------------------------------------!
               fx_b0001_potl  = fx_f_b0001  * bfuel_d0001_tot * burnt_area_potl
               fx_b0010_potl  = fx_f_b0010  * bfuel_d0010_tot * burnt_area_potl
               fx_b0100_potl  = fx_f_b0100  * bfuel_d0100_tot * burnt_area_potl
               fx_b1000_potl  = fx_f_b1000  * bfuel_d1000_tot * burnt_area_potl
               fx_bherb_potl  = fx_f_bherb  * bherb_tot       * burnt_area_potl
               fx_bwoody_potl = fx_f_bwoody * bwoody_tot      * burnt_area_potl
               fx_wn1000_potl = fx_f_wn1000 * bwn1000_tot     * burnt_area_potl
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !       Find fire intensity of this step.                                   !
               !---------------------------------------------------------------------------!
               if (burnt_area_potl > tiny_num) then
                  !----- Find fire intensity. ---------------------------------------------!
                  fx_intensity = fr_h * rosfwd_avg * C2B                                   &
                               * ( fx_b0001_potl + fx_b0010_potl + fx_b0100_potl           &
                                 + fx_bherb_potl + fx_wn1000_potl                )         &
                               / burnt_area_potl
                  !------------------------------------------------------------------------!


                  !------ If potential intensity is zero, set burnt area to zero as well. -!
                  if (fx_intensity <= tiny_num) then
                     fx_intensity    = 0.0
                     burnt_area_potl = 0.0
                  end if
                  !------------------------------------------------------------------------!
               else
                  !------------------------------------------------------------------------!
                  !      No burnt area, set it to zero, and zero fire intensity and        !
                  ! combustion.                                                            !
                  !------------------------------------------------------------------------!
                  burnt_area_potl   = 0.0
                  fx_intensity      = 0.0
                  fx_b0001_potl     = 0.0
                  fx_b0010_potl     = 0.0
                  fx_b0100_potl     = 0.0
                  fx_b1000_potl     = 0.0
                  fx_bherb_potl     = 0.0
                  fx_bwoody_potl    = 0.0
                  fx_wn1000_potl    = 0.0
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Compute the fuel limitation control on fire persistence.  This       !
               ! replaces the precipitation term in LP15 with a potential fire intensity   !
               ! term that accounts for fuel loads and fuel moisture.  The first guess     !
               ! parameters are based on typical scorch heights for tropical broadleaf     !
               ! evergreen forests, using T10 parameters.                                  !
               !---------------------------------------------------------------------------!
               fintn       = ( fx_intensity - ft_fint_lwr ) * ft_fint_dti
               fintn       = max(0.,min(1.,fintn))
               fp_fuel_fun = bpow01(fintn,ft_fint_exp)
               if (fp_fuel_fun < almost_zero) fp_fuel_fun = 0.0
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !      Compute fragmentation control on termination.  This is an empirical  !
               ! function and thus should include the fragmentation due to areas that      !
               ! cannot sustain vegetation.                                                !
               !---------------------------------------------------------------------------!
               flamn       = max(0.,min(1.,(1. - fragn) * cgrid%landfrac(ipy)))
               fp_cntg_fun = 1. - bpow01(1. - flamn,ft_frag_exp)
               if (fp_cntg_fun < almost_zero) fp_cntg_fun = 0.0
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !       Compute the fire suppression function.                              !
               !---------------------------------------------------------------------------!
               !----- Fire suppressibility. -----------------------------------------------!
               suppressibility = sqrt((1. - fp_fdi_fun) * (1. - fp_wind_fun))
               suppressibility = max(0.,min(1.,suppressibility))
               !----- Normalised land use . -----------------------------------------------!
               lu_norm     = max(0.,min(1.,lu_area / ft_lu_upr))
               !----- Normalised HDI . ----------------------------------------------------!
               hdin        = ( cpoly%seitimes(isei,isi)%hdi - ft_hdi_lwr ) * ft_hdi_dti
               hdin        = max(0., min(1., hdin) )
               !----- Anthropogenic-controlled persistence driven by HDI and land use. ----!
               fp_anth_fun = 1. - bpow01(lu_norm,ft_lu_exp) * bpow01(hdin,ft_hdi_exp)
               !----- Fire wildfire persistence function. ---------------------------------!
               fp_wild_fun = fp_anth_fun * (1. - suppressibility)
               if (fp_wild_fun < almost_zero) fp_wild_fun = 0.0
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !       Probability of fire persistence (i.e., probability that  fires will !
               ! not extinguish).  Because the fire time step is flexible and we don't     !
               ! want parameters to be strongly dependent upon the fire time step, we      !
               ! assume this probability refers to the probability of fires to persist for !
               ! 24 hours.  This is different from HESFIRE (which uses 12-h steps) but     !
               ! 24 hours is more convenient for ED2.                                      !
               !---------------------------------------------------------------------------!
               prob_persist = cbrt( fp_fuel_fun * fp_cntg_fun * fp_wild_fun )
               if (prob_persist < almost_zero) then
                  cpoly%fire_extinction(isi) = - lnexp_min / dtfull
               else
                  cpoly%fire_extinction(isi) = - log(prob_persist) / dtfull
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Find the scaling factor average time duration.  The average time      !
               ! duration is given by 1/(extinction rate).  In case the average time       !
               ! exceeds the time step, we do not amplify fires as we will continue to     !
               ! integrate them over the next step.                                        !
               !---------------------------------------------------------------------------!
               if (cpoly%fire_extinction(isi) < (1. / dtfire)) then
                  fx_duration = 1.
               else
                  fx_duration = 1. / ( cpoly%fire_extinction(isi) * dtfire )
               end if
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !      Scale burnt area and combustion factors by the square of fire        !
               ! duration (as burnt area is proportional to the square of time).  Fire     !
               ! intensity does not need to be rescaled because the burnt area appears in  !
               ! the numerator and denominator, so it is effectively independent on the    !
               ! fire duration.                                                            !
               !---------------------------------------------------------------------------!
               burnt_area_step = burnt_area_potl * fx_duration * fx_duration
               fx_b0001        = fx_b0001_potl   * fx_duration * fx_duration
               fx_b0010        = fx_b0010_potl   * fx_duration * fx_duration
               fx_b0100        = fx_b0100_potl   * fx_duration * fx_duration
               fx_b1000        = fx_b1000_potl   * fx_duration * fx_duration
               fx_bherb        = fx_bherb_potl   * fx_duration * fx_duration
               fx_bwoody       = fx_bwoody_potl  * fx_duration * fx_duration
               fx_wn1000       = fx_wn1000_potl  * fx_duration * fx_duration
               !---------------------------------------------------------------------------!







               !---------------------------------------------------------------------------!
               !       Find duration of lethal bole heating, following PR86.               !
               !---------------------------------------------------------------------------!
               if (fx_intensity > tiny_num) then
                  !----- Find lethal duration. --------------------------------------------!
                  fx_tlethal = fx_tlh_slope * C2B                                          &
                             * ( bfuel_d0001_tot      * (1. - sqrt(1. - fx_f_b0001 ) )     &
                               + bfuel_d0010_tot      * (1. - sqrt(1. - fx_f_b0010 ) )     &
                               + bfuel_d0100_tot      * (1. - sqrt(1. - fx_f_b0100 ) )     &
                               + bherb_tot            * (1. - sqrt(1. - fx_f_bherb ) )     &
                               + bwn1000_tot          * (1. - sqrt(1. - fx_f_wn1000) ) )
                  !------------------------------------------------------------------------!
               else
                  !----- No burning, set it to zero. --------------------------------------!
                  fx_tlethal = 0.0
                  !------------------------------------------------------------------------!
               end if
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
               !      Add burnt area from this step, then calculate the lethality and      !
               ! combustion fraction associated with this time step.                       !
               !---------------------------------------------------------------------------!
               if (cgrid%landfrac(ipy) > tiny_num) then
                  !----- Integrate area. --------------------------------------------------!
                  cpoly%burnt_area(isi) = cpoly%burnt_area(isi)                            &
                                        + burnt_area_step / cgrid%landfrac(ipy)
                  !------------------------------------------------------------------------!



                  !----- Integrate monthly burnt area. ------------------------------------!
                  cpoly%avg_burnt_area(imonth,isi) = cpoly%avg_burnt_area(imonth,isi)      &
                                                   + burnt_area_step / cgrid%landfrac(ipy)
                  !------------------------------------------------------------------------!



                  !----- Integrate lethality. ---------------------------------------------!
                  call integ_fire_lethality( cpoly,isi,iwhen,fx_intensity,fx_tlethal       &
                                           , burnt_area_step / cgrid%landfrac(ipy) )
                  !------------------------------------------------------------------------!



                  !----- Integrate combustion fraction. -----------------------------------!
                  cpoly%fire_f_bherb (isi) = cpoly%fire_f_bherb (isi)                      &
                                           + fx_f_bherb                                    &
                                           * burnt_area_step / cgrid%landfrac(ipy)
                  cpoly%fire_f_bwoody(isi) = cpoly%fire_f_bwoody(isi)                      &
                                           + fx_f_bwoody                                   &
                                           * burnt_area_step / cgrid%landfrac(ipy)
                  cpoly%fire_f_fgc   (isi) = cpoly%fire_f_fgc   (isi)                      &
                                           + fx_f_fgc                                      &
                                           * burnt_area_step / cgrid%landfrac(ipy)
                  cpoly%fire_f_stgc  (isi) = cpoly%fire_f_stgc  (isi)                      &
                                           + fx_f_stgc                                     &
                                           * burnt_area_step / cgrid%landfrac(ipy)
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Integrate the daily average fire spread.                              !
               !---------------------------------------------------------------------------!
               cpoly%today_fire_density   (isi) = cpoly%today_fire_density       (isi)     &
                                                + ndfirei * cpoly%fire_density   (isi)
               cpoly%today_fire_extinction(isi) = cpoly%today_fire_extinction    (isi)     &
                                                + ndfirei * cpoly%fire_extinction(isi)
               cpoly%fire_spread          (isi) = cpoly%fire_spread              (isi)     &
                                                + ndfirei * rosfwd_avg
               cpoly%fire_intensity       (isi) = cpoly%fire_intensity           (isi)     &
                                                + ndfirei * fx_intensity
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !       Lethal time is averaged using the inverse time, to account for      !
               ! periods with no fires.                                                    !
               !---------------------------------------------------------------------------!
               if (fx_tlethal > tiny_num) then
                  cpoly%fire_tlethal(isi) = cpoly%fire_tlethal(isi) + ndfirei / fx_tlethal
               end if
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Print the output if needed.                                           !
               !---------------------------------------------------------------------------!
               if (printout) then
                  open(unit=35,file=firefile,status='old',position='append',action='write')
                  write(unit=35,fmt='(5(i6,1x),33(es12.3,1x))')                            &
                             current_time%year,current_time%month,current_time%date,iwhen  &
                            ,isi,apy_area,cgrid%landfrac(ipy),fragn,lu_area                &
                            ,cpoly%seitimes(isei,isi)%hdi,cpoly%flashtimes(iflash,isi)%c2g &
                            ,total_ignition,nesterov_avg,fdivpd_avg,moist_bherb_avg        &
                            ,moist_bwoody_avg,moist_bfuel_avg,rosfwd_avg,fx_f_fgc          &
                            ,fx_f_stgc,fx_f_bherb,fx_f_bwoody,fx_duration,fx_intensity     &
                            ,fx_tlethal,fp_fdi_fun,fp_wind_fun,suppressibility,fp_anth_fun &
                            ,fp_wild_fun,fp_fuel_fun,fp_cntg_fun,prob_persist              &
                            ,burnt_area_step,cpoly%burnt_area(isi),cpoly%fire_density(isi) &
                            ,cpoly%fire_extinction(isi)
                  close(unit=35,status='keep')
               end if
               !---------------------------------------------------------------------------!


            end do timestep_loop
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !       Update the average fire intensity.                                     !
            !------------------------------------------------------------------------------!
            cpoly%avg_fire_intensity(imonth,isi) = cpoly%avg_fire_intensity(imonth,isi)    &
                                                 + cpoly%fire_intensity           (isi)    &
                                                 * ndaysi
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      For fire consumption, we do not weight by time step; instead, we use    !
            ! the burnt area as weights.  We normalise the consumption rate at the monthly !
            ! time step.                                                                   !
            !------------------------------------------------------------------------------!
            cpoly%avg_fire_f_bherb  (imonth,isi) = cpoly%avg_fire_f_bherb  (imonth,isi)    &
                                                 + cpoly%fire_f_bherb             (isi)
            cpoly%avg_fire_f_bwoody (imonth,isi) = cpoly%avg_fire_f_bwoody (imonth,isi)    &
                                                 + cpoly%fire_f_bwoody            (isi)
            cpoly%avg_fire_f_fgc    (imonth,isi) = cpoly%avg_fire_f_fgc    (imonth,isi)    &
                                                 + cpoly%fire_f_fgc               (isi)
            cpoly%avg_fire_f_stgc   (imonth,isi) = cpoly%avg_fire_f_stgc   (imonth,isi)    &
                                                 + cpoly%fire_f_stgc              (isi)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Integrate lethal heating duration. We use the inverse so we account for  !
            ! times with no fire (when time should go to infinity).  We correct the        !
            ! averaged values at subroutine fire_frequency.  We fix the instantaneous      !
            ! value here though.                                                           !
            !------------------------------------------------------------------------------!
            if (cpoly%fire_tlethal(isi) > tiny_num) then
               !------ Invert instantaneous lethal heating duration. ----------------------!
               cpoly%fire_tlethal(isi) = 1. / cpoly%fire_tlethal(isi)
               !---------------------------------------------------------------------------!
            else
               !------ No fire. -----------------------------------------------------------!
               cpoly%fire_tlethal(isi) = 0.
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!



            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!



      return
   end subroutine integ_firestarter
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine integrates fire lethality (i.e., fire mortality given that there  !
   ! was fire) when running FIRESTARTER.                                                   !
   !---------------------------------------------------------------------------------------!
   subroutine integ_fire_lethality(cpoly,isi,iwhen,fx_intensity,fx_tlethal,burnt_area_step)
      use disturb_coms , only : fx_tlc_slope  & ! intent(in)
                              , fx_pmtau_di   & ! intent(in)
                              , fx_pmtau_ds   ! ! intent(in)
      use ed_misc_coms , only : current_time  ! ! intent(in)
      use ed_state_vars, only : polygontype   & ! structure
                              , sitetype      & ! structure
                              , patchtype     ! ! structure
      use ed_max_dims  , only : n_pft         & ! intent(in)
                              , n_dist_types  ! ! intent(in)
      use consts_coms  , only : lnexp_max     & ! intent(in)
                              , tiny_num      & ! intent(in)
                              , almost_one    ! ! intent(in)
      use pft_coms     , only : escorch       & ! intent(in)
                              , fscorch       & ! intent(in)
                              , fx_rck_pft    & ! intent(in)
                              , fx_pck_pft    ! ! intent(in)
      use allometry    , only : h2crownbh     ! ! function
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(polygontype), target      :: cpoly
      integer          , intent(in)  :: isi
      integer          , intent(in)  :: iwhen
      real             , intent(in)  :: fx_intensity
      real             , intent(in)  :: fx_tlethal
      real             , intent(in)  :: burnt_area_step
      !----- Local variables. -------------------------------------------------------------!
      type(sitetype)   , pointer     :: csite
      type(patchtype)  , pointer     :: cpatch
      logical                        :: has_tlethal
      logical                        :: has_tlcrit
      logical                        :: has_intensity
      integer                        :: ipa
      integer                        :: ico
      integer                        :: ipft
      real                           :: scorch_height
      real                           :: crown_damage
      real                           :: pmtau
      real                           :: pmck
      real                           :: p_mort
      real                           :: tlethal_crit
      real                           :: chbase
      real                           :: clength
      !----- Local parameters. ------------------------------------------------------------!
      character(len=23) , parameter  :: firefile = 'lethalfire_details.txt'
      logical           , parameter  :: printout = .false.
      !----- Locally saved variables. -----------------------------------------------------!
      logical           , save       :: first_time = .true. ! First time calling  [    T|F]
      !------------------------------------------------------------------------------------!


      !----- First time, and the user wants to print the output.  Make a header. ----------!
      if (first_time) then

         !----- Make the header. ----------------------------------------------------------!
         if (printout) then
            open (unit=35,file=firefile,status='replace',action='write')
            write (unit=35,fmt='(19(a,1x))')                                               &
                     '  YEAR',      ' MONTH',      '   DAY',      '  STEP',      '   ISI'  &
              ,      '   IPA',      '   ICO',      '   PFT','         DBH','      HEIGHT'  &
              ,'      THBARK','FX_INTENSITY','     HSCORCH','CROWN_DAMAGE','  FX_TLETHAL'  &
              ,'TLETHAL_CRIT','       PMTAU','        PMCK','      P_MORT'
            close (unit=35,status='keep')
         end if
         !---------------------------------------------------------------------------------!


         first_time = .false.
      end if
      !------------------------------------------------------------------------------------!


      !------ Alias for current site. -----------------------------------------------------!
      csite => cpoly%site(isi)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Loop through patches.                                                          !
      !------------------------------------------------------------------------------------!
      patch_loop: do ipa=1,csite%npatches
         cpatch => csite%patch(ipa)

         !---------------------------------------------------------------------------------!
         !     Loop through cohorts.                                                       !
         !---------------------------------------------------------------------------------!
         cohort_loop: do ico=1,cpatch%ncohorts
            !------ Handy aliases. --------------------------------------------------------!
            ipft = cpatch%pft(ico)
            !------------------------------------------------------------------------------!


            !------ Scorch height [m]. ----------------------------------------------------!
            scorch_height = fscorch(ipft) * fx_intensity ** escorch(ipft)
            !------------------------------------------------------------------------------!


            !------ Compute fire damage from scorch height relative to cohort height. -----!
            chbase       = h2crownbh(cpatch%hite(ico),ipft)
            clength      = cpatch%hite(ico) - chbase
            crown_damage = ( scorch_height - chbase ) / clength
            crown_damage = max( 0., crown_damage )
            !------------------------------------------------------------------------------!


            !------ Find critical duration of lethal heating. -----------------------------!
            tlethal_crit = fx_tlc_slope * cpatch%thbark(ico) * cpatch%thbark(ico)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Flags for critical fire duration and sufficient fire intensity.         !
            !------------------------------------------------------------------------------!
            has_tlethal   = fx_tlethal   > tiny_num
            has_tlcrit    = tlethal_crit > tiny_num
            has_intensity = fx_intensity > tiny_num
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !      Find the mortality probability due to cambial damage.  Check the        !
            ! critical fire duration to avoid singularities.                               !
            !------------------------------------------------------------------------------!
            if ( has_tlcrit ) then
               !------ Find the tl:tc ratio. ----------------------------------------------!
               pmtau = fx_pmtau_di + fx_pmtau_ds * fx_tlethal / tlethal_crit
               pmtau = max(0.,min(1., pmtau))
               !---------------------------------------------------------------------------!
            else if ( has_intensity ) then
               !------ Critical duration is zero, assume maximum mortality. ---------------!
               pmtau = 1.0
               !---------------------------------------------------------------------------!
            else
               !------ No fire, no fire mortality... --------------------------------------!
               pmtau = 0.0
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!


            !------ Find the mortality probability due to crown damage. -------------------!
            pmck = fx_rck_pft(ipft) * crown_damage ** fx_pck_pft(ipft)
            pmck = max(0.,min(1.,pmck))
            !------------------------------------------------------------------------------!


            !------ Find the probability of mortality. ------------------------------------!
            p_mort = max(0.,min(1.,pmtau + pmck - pmtau * pmck))
            !------------------------------------------------------------------------------!


            !------ Integrate lethality probability due to fire. --------------------------!
            cpatch%fire_lethal_rate(13,ico) = cpatch%fire_lethal_rate(13,ico)              &
                                            + p_mort * burnt_area_step
            cpatch%fire_lethal_prob   (ico) = cpatch%fire_lethal_prob   (ico)              &
                                            + p_mort * burnt_area_step
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !     Print the output if needed.                                              !
            !------------------------------------------------------------------------------!
            if (printout .and. ( has_tlethal .or. has_intensity ) ) then
               open(unit=35,file=firefile,status='old',position='append',action='write')
               write(unit=35,fmt='(8(i6,1x),3(f12.4,1x),f12.3,1x,7(f12.4,1x))')            &
                          current_time%year,current_time%month,current_time%date,iwhen     &
                         ,isi,ipa,ico,ipft,cpatch%dbh(ico),cpatch%hite(ico)                &
                         ,cpatch%thbark(ico),fx_intensity,scorch_height,crown_damage       &
                         ,fx_tlethal,tlethal_crit,pmtau,pmck,p_mort
               close(unit=35,status='keep')
            end if
            !------------------------------------------------------------------------------!


         end do cohort_loop
         !---------------------------------------------------------------------------------!
      end do patch_loop
      !------------------------------------------------------------------------------------!

      return
   end subroutine integ_fire_lethality
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !       Sub-routine that resets the monthly fire-related variables when using fire      !
   ! models that require daily or sub-daily integration.  Currently this applies to        !
   ! INCLUDE_FIRE = 3 (EMBERFIRE) or INCLUDE_FIRE = 4 (FIRESTARTER).                       !
   !---------------------------------------------------------------------------------------!
   subroutine reset_monthly_fire(cgrid)
      use ed_state_vars , only : edtype                 & ! structure
                               , polygontype            & ! structure
                               , sitetype               & ! structure
                               , patchtype              ! ! structure
      use ed_misc_coms  , only : current_time           ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)      , target     :: cgrid
      !----- Local variables --------------------------------------------------------------!
      type(polygontype) , pointer    :: cpoly
      type(sitetype)    , pointer    :: csite
      type(patchtype)   , pointer    :: cpatch
      integer                        :: ipy
      integer                        :: isi
      integer                        :: ipa
      integer                        :: ico
      integer                        :: imo
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
               cpatch => csite%patch(ipa)


               !----- Reset the ground water for next month. ------------------------------!
               csite%avg_monthly_gndwater(ipa) = 0.
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Reset fire lethality by month.                                        !
               !---------------------------------------------------------------------------!
               cohortloop: do ico=1,cpatch%ncohorts
                  cpatch%fire_lethal_rate(13,ico) = 0.0
               end do cohortloop
               !---------------------------------------------------------------------------!
            end do patchloop
            !------------------------------------------------------------------------------!


            !------ Resetting these variables just to play it safe. -----------------------!
            cpoly%fire_intensity         (isi) = 0.0
            cpoly%fire_spread            (isi) = 0.0
            cpoly%ignition_rate          (isi) = 0.0
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!


      return
   end subroutine reset_monthly_fire
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !       Sub-routine that resets the monthly fire-related variables when using fire      !
   ! models that require daily or sub-daily integration.  Currently this applies to        !
   ! INCLUDE_FIRE = 3 (EMBERFIRE) or INCLUDE_FIRE = 4 (FIRESTARTER).                       !
   !---------------------------------------------------------------------------------------!
   subroutine reset_yearly_fire(cgrid)
      use ed_state_vars , only : edtype                 & ! structure
                               , polygontype            & ! structure
                               , sitetype               & ! structure
                               , patchtype              ! ! structure
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)      , target     :: cgrid
      !----- Local variables --------------------------------------------------------------!
      type(polygontype) , pointer    :: cpoly
      type(sitetype)    , pointer    :: csite
      type(patchtype)   , pointer    :: cpatch
      integer                        :: ipy
      integer                        :: isi
      integer                        :: ipa
      integer                        :: ico
      !------------------------------------------------------------------------------------!


      !----- Loop over polygons and sites. ------------------------------------------------!
      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         !---- Loop over all sites. -------------------------------------------------------!
         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            !----- Reset burnt area. ------------------------------------------------------!
            cpoly%burnt_area(isi) = 0.0
            !------------------------------------------------------------------------------!



            !----- Reset fire disturbance rates. ------------------------------------------!
            cpoly%lambda_fire       (:,isi) = 0.0
            cpoly%avg_burnt_area    (:,isi) = 0.0
            cpoly%avg_fire_intensity(:,isi) = 0.0
            cpoly%avg_fire_f_bherb  (:,isi) = 0.0
            cpoly%avg_fire_f_bwoody (:,isi) = 0.0
            cpoly%avg_fire_f_fgc    (:,isi) = 0.0
            cpoly%avg_fire_f_stgc   (:,isi) = 0.0
            !------------------------------------------------------------------------------!


            !---- Loop over all patches. --------------------------------------------------!
            patchloop: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)
               !----- Loop over all cohorts. ----------------------------------------------!
               cohortloop: do ico=1,cpatch%ncohorts
                  cpatch%fire_lethal_prob  (ico) = 0.0
                  cpatch%fire_lethal_rate(:,ico) = 0.0
               end do cohortloop
               !---------------------------------------------------------------------------!
            end do patchloop
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!


      return
   end subroutine reset_yearly_fire
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !      This function calculates the forward spread rate, using the A18's revision of    !
   ! the R72 fire spread model.  If needed, it is possible to obtain the maximum rate of   !
   ! spread  using maximum winds and "very dry" moisture conditions as defined by SB05.    !
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
   subroutine rate_of_spread(isi                                                           &
                            ,bfuel_d0001,bfuel_d0010,bfuel_d0100,bfuel_d1000,bherb,bwoody  &
                            ,moist_b0001,moist_b0010,moist_b0100,moist_b1000,moist_bherb   &
                            ,moist_bwoody,can_wind,use_max,g_Umax,rosfwd)
      use disturb_coms, only : n_fst         & ! intent(in)
                             , n_fcl         & ! intent(in)
                             , n_sbmax       & ! intent(in)
                             , n_sbins       & ! intent(in)
                             , fr_moist_ij   & ! intent(in)
                             , fr_sigma_ij   & ! intent(in)
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
      use ed_misc_coms, only : current_time  ! ! intent(in)
      use pft_coms    , only : C2B           ! ! intent(in)
      use consts_coms , only : tiny_num      & ! intent(in)
                             , lnexp_min     & ! intent(in)
                             , lnexp_max     & ! intent(in)
                             , almost_zero   ! ! intent(in)
      implicit none

      !------ Arguments. ------------------------------------------------------------------!
      integer, intent(in)            :: isi          !  Site number             [      ---]
      real   , intent(in)            :: bfuel_d0001  !  1-hr dead fuel load     [   kgC/m2]
      real   , intent(in)            :: bfuel_d0010  !  10-hr dead fuel load    [   kgC/m2]
      real   , intent(in)            :: bfuel_d0100  !  100-hr dead fuel load   [   kgC/m2]
      real   , intent(in)            :: bfuel_d1000  !  1000-hr dead fuel load  [   kgC/m2]
      real   , intent(in)            :: bherb        !  Herbaceous fuel load    [   kgC/m2]
      real   , intent(in)            :: bwoody       !  Living Woody fuel load  [   kgC/m2]
      real   , intent(in)            :: moist_b0001  ! 1-hr fuel moisture       [      ---]
      real   , intent(in)            :: moist_b0010  ! 10-hr fuel moisture      [      ---]
      real   , intent(in)            :: moist_b0100  ! 100-hr fuel moisture     [      ---]
      real   , intent(in)            :: moist_b1000  ! 1000-hr fuel moisture    [      ---]
      real   , intent(in)            :: moist_bherb  ! Herb. fuel moisture      [      ---]
      real   , intent(in)            :: moist_bwoody ! Woody fuel moisture      [      ---]
      real   , intent(in)            :: can_wind     ! Canopy air space wind    [      m/s]
      logical, intent(in)            :: use_max      ! Find maximum spread      [      T|F]
      real   , intent(out)           :: g_Umax       ! Max. wind (corrected)    [      m/s]
      real   , intent(out)           :: rosfwd       ! Forward rate of spread   [      m/s]
      !----- Local variables (by fuel class and status). ----------------------------------!
      real, dimension(n_fst,n_fcl)   :: wood_ij      ! Fuel load                [   kgB/m2]
      real, dimension(n_fst,n_fcl)   :: moist_ij     ! Fuel moisture            [      ---]
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
      real                           :: g_HeatSink   ! Heat sink                [     J/m3]
      real                           :: g_psiw       ! R72's wind factor        [       --]
      real                           :: g_ros_w0     ! rate of spread (no wind) [      m/s]
      !----- Additional local variables. --------------------------------------------------!
      integer                        :: i            ! Counter                  [       --]
      integer                        :: j            ! Counter                  [       --]
      integer                        :: k            ! Counter                  [       --]
      real                           :: mf_dead      ! "Fine" dead fuel moist.  [    kg/kg]
      real                           :: lnexp        ! Temp var to avoid FPE    [       --]
      !----- Additional parameters. -------------------------------------------------------!
      character(len= 3), parameter   :: fmth = '(a)'
      character(len=21), parameter   :: fmtt = '(a,1x,i4.4,2(a,i2.2))'
      character(len= 9), parameter   :: fmti = '(a,1x,i5)'
      character(len= 9), parameter   :: fmtl = '(a,1x,l1)'
      character(len=13), parameter   :: fmte = '(a,1x,es10.3)'
      !----- External functions. ----------------------------------------------------------!
      logical, external              :: isnan_real   ! Number is NaN            [      T|F]
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
      !       Define fuel moisture.                                                        !
      !------------------------------------------------------------------------------------!
      if (use_max) then
         moist_ij(:,:) = fr_moist_ij(:,:)
      else
         moist_ij(1,:) = (/  moist_b0001,  moist_b0010,  moist_b0100,  moist_b1000         &
                          ,  moist_bherb, moist_bwoody /)
         moist_ij(2,:) = moist_ij(1,:)
      end if
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
      Mx_i(2) = max( fr_Mxdead                                                             &
                   , fr_mxl_aa(1) + fr_mxl_aa(2) * g_W * (1. - mf_dead/fr_Mxdead) )
      !----- Fuel moisture. ---------------------------------------------------------------!
      do i=1,n_fst
         moist_i(i) = sum(moist_ij(i,:)*fwgt_ij(i,:))
      end do
      !----- Relative moisture. -----------------------------------------------------------!
      do i=1,n_fst
         !----- Moisture by status and class. ---------------------------------------------!
         do j=1,n_fcl
            rmoist_ij(i,j) = min(1.,moist_ij(i,j)/Mx_i(i))
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
      Qig_ij(:,:) = fr_Qig_aa(1) + fr_Qig_aa(2) * moist_ij(:,:)
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
      !     Find the wind factor based on R72, but capping at the maximum wind.            !
      !------------------------------------------------------------------------------------!
      if (use_max) then
         g_psiw = g_CC * g_Umax ** g_BB / g_beta ** g_EE
      else
         g_psiw = g_CC * min(can_wind,g_Umax) ** g_BB / g_beta ** g_EE
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the forward rate of spread.                                               !
      !------------------------------------------------------------------------------------!
      !----- Forward rate of spread in the absence of wind  [m/s]. ------------------------!
      g_ros_w0 = g_Ir * g_xi / g_HeatSink
      !----- Forward rate of spread                         [m/s]. ------------------------!
      rosfwd   = g_ros_w0 * (1. + g_psiw)
      !----- Make rate of spread bounded. -------------------------------------------------!
      if (isnan_real(rosfwd) .or. (rosfwd < - almost_zero)) then
         write(unit=*,fmt=fmth) '---------------------------------------------------------'
         write(unit=*,fmt=fmth) '  Incorrect rate of spread detected in FIRESTARTER.'
         write(unit=*,fmt=fmth) '---------------------------------------------------------'
         write(unit=*,fmt=fmtt) 'Time:',current_time%year,'-',current_time%month,'-'       &
                               ,current_time%date
         write(unit=*,fmt=fmti) 'Site:',isi
         write(unit=*,fmt=fmth) ''
         write(unit=*,fmt=fmte) ' BFUEL_D0001     =',bfuel_d0001
         write(unit=*,fmt=fmte) ' BFUEL_D0010     =',bfuel_d0010
         write(unit=*,fmt=fmte) ' BFUEL_D0100     =',bfuel_d0100
         write(unit=*,fmt=fmte) ' BFUEL_D1000     =',bfuel_d1000
         write(unit=*,fmt=fmte) ' BHERB           =',bherb
         write(unit=*,fmt=fmte) ' BWOODY          =',bwoody
         write(unit=*,fmt=fmte) ' MOIST_B0001     =',moist_b0001
         write(unit=*,fmt=fmte) ' MOIST_B0010     =',moist_b0010
         write(unit=*,fmt=fmte) ' MOIST_B0100     =',moist_b0100
         write(unit=*,fmt=fmte) ' MOIST_B1000     =',moist_b1000
         write(unit=*,fmt=fmte) ' MOIST_BHERB     =',moist_bherb
         write(unit=*,fmt=fmte) ' MOIST_BWOODY    =',moist_bwoody
         write(unit=*,fmt=fmte) ' CAN_WIND        =',can_wind
         write(unit=*,fmt=fmtl) ' USE_MAX         =',use_max
         write(unit=*,fmt=fmte) ' G_FAI           =',g_fai
         write(unit=*,fmt=fmte) ' G_WOOD          =',g_wood
         write(unit=*,fmt=fmte) ' G_WNOD          =',g_wnod
         write(unit=*,fmt=fmte) ' G_SIGMA         =',g_sigma
         write(unit=*,fmt=fmte) ' G_HH            =',g_hh
         write(unit=*,fmt=fmte) ' MF_DEAD         =',mf_dead
         write(unit=*,fmt=fmte) ' G_W             =',g_W
         write(unit=*,fmt=fmte) ' MX_DEAD         =',Mx_i(1)
         write(unit=*,fmt=fmte) ' MX_LIVE         =',Mx_i(2)
         write(unit=*,fmt=fmte) ' G_RHOB          =',g_rhob
         write(unit=*,fmt=fmte) ' G_BETA          =',g_beta
         write(unit=*,fmt=fmte) ' G_BETA_OP       =',g_beta_op
         write(unit=*,fmt=fmte) ' G_R_BETA        =',g_r_beta
         write(unit=*,fmt=fmte) ' G_GAMMA_MAX     =',g_gamma_max
         write(unit=*,fmt=fmte) ' G_GAMMA         =',g_gamma
         write(unit=*,fmt=fmte) ' G_XI            =',g_xi
         write(unit=*,fmt=fmte) ' G_IR            =',g_Ir
         write(unit=*,fmt=fmte) ' G_UMAX          =',g_Umax
         write(unit=*,fmt=fmte) ' G_AA            =',g_AA
         write(unit=*,fmt=fmte) ' G_BB            =',g_BB
         write(unit=*,fmt=fmte) ' G_CC            =',g_CC
         write(unit=*,fmt=fmte) ' G_EE            =',g_EE
         write(unit=*,fmt=fmte) ' G_HEATSINK      =',g_HeatSink
         write(unit=*,fmt=fmte) ' G_PSIW          =',g_psiw
         write(unit=*,fmt=fmte) ' G_ROS_W0        =',g_ros_w0
         write(unit=*,fmt=fmte) ' ROSFWD          =',rosfwd
         write(unit=*,fmt=fmth) '---------------------------------------------------------'
         call fatal_error('Invalid rate of spread in FIRESTARTER.'                         &
                         ,'rate_of_spread','fire.f90')
      else if (rosfwd < almost_zero) then
         rosfwd = 0.0
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine rate_of_spread
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
      real                  , intent(in)  :: bfuel_d0001  ! 1-hr dead fuel       [  kgC/m2]
      real                  , intent(in)  :: bfuel_d0010  ! 10-hr dead fuel      [  kgC/m2]
      real                  , intent(in)  :: bfuel_d0100  ! 100-hr dead fuel     [  kgC/m2]
      real                  , intent(in)  :: bfuel_d1000  ! 1000-hr dead fuel    [  kgC/m2]
      real                  , intent(in)  :: bherb        ! Herb. fuel           [  kgC/m2]
      real                  , intent(in)  :: bwoody       ! Woody fuel           [  kgC/m2]
      real                  , intent(in)  :: moist_b0001  ! 1-hr fuel moist.     [     ---]
      real                  , intent(in)  :: moist_b0010  ! 10-hr fuel moist.    [     ---]
      real                  , intent(in)  :: moist_b0100  ! 100-hr fuel moist.   [     ---]
      real                  , intent(in)  :: moist_b1000  ! 1000-hr fuel moist.  [     ---]
      real                  , intent(in)  :: moist_bherb  ! Herb. fuel moist.    [     ---]
      real                  , intent(in)  :: moist_bwoody ! Woody fuel moisture  [     ---]
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
      Mx_i(2) = max( fr_Mxdead                                                             &
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

