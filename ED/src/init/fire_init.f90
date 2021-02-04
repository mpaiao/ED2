!==========================================================================================!
!==========================================================================================!
! MODULE: FIRE_INIT
!> \brief This module contains sub-routines used to initialise data used to drive fire 
!>        ignition and fire suppression.
!------------------------------------------------------------------------------------------!
module fire_init
   contains

   !=======================================================================================!
   !=======================================================================================!
   !     This file will read the landuse files and assign the anthropogenic disturbance    !
   ! matrices.                                                                             !
   !---------------------------------------------------------------------------------------!
   subroutine read_fire_ignition

      use ed_state_vars , only : edtype            & ! structure
                               , polygontype       & ! structure
                               , sitetype          & ! structure
                               , edgrid_g          ! ! structure
      use consts_coms   , only : erad              & ! intent(in)
                               , pio180            & ! intent(in)
                               , huge_num          ! ! intent(in)
      use disturb_coms  , only : seitime           & ! structure
                               , flashtime         & ! structure
                               , include_fire      & ! structure
                               , max_lu_years      & ! intent(in)
                               , sei_database      & ! intent(in)
                               , flash_database    ! ! intent(in)
      use ed_misc_coms  , only : iyeara            & ! intent(in)
                               , iyearz            ! ! intent(in)
      use grid_coms     , only : ngrids            ! ! intent(in)
      use ed_max_dims   , only : str_len           & ! intent(in)
                               , maxlist           ! ! intent(in)

      implicit none
      !----- Local variables --------------------------------------------------------------!
      type(edtype)          , pointer                 :: cgrid
      type(polygontype)     , pointer                 :: cpoly
      type(sitetype)        , pointer                 :: csite
      type(seitime)         , pointer                 :: cseitime
      type(seitime)         , pointer                 :: rseitime
      type(flashtime)       , pointer                 :: cflashtime
      type(flashtime)       , pointer                 :: rflashtime
      character(len=str_len), dimension(maxlist)      :: full_list
      character(len=str_len), dimension(maxlist)      :: fire_list
      real                  , dimension(maxlist)      :: llon_list
      real                  , dimension(maxlist)      :: llat_list
      real                  , dimension(maxlist)      :: file_ldist
      character(len=6)                                :: hform
      character(len=str_len)                          :: fire_name
      character(len=str_len)                          :: cdum
      integer                                         :: nf
      integer                                         :: nflist
      integer                                         :: nflfire
      integer                                         :: ncl
      integer                                         :: iya
      integer                                         :: iyr
      integer                                         :: iyz
      integer                                         :: imo
      integer                                         :: igr
      integer                                         :: ipy
      integer                                         :: isi
      integer                                         :: sim_years
      integer                                         :: yfirst
      integer                                         :: yd_1st
      integer                                         :: yd_this
      integer                                         :: ylast
      integer                                         :: yd_last
      logical                                         :: inside
      real                                            :: fire_area
      real                                            :: wlon
      real                                            :: elon
      real                                            :: slat
      real                                            :: nlat
      !----- Local constants. -------------------------------------------------------------!
      character(len=12)     , parameter :: fffmt    = '(a,1x,f12.5)'
      character(len=13)     , parameter :: esfmt    = '(a,1x,es12.5)'
      integer               , parameter :: hoff     = 18
      !----- External function. -----------------------------------------------------------!
      real                  , external  :: dist_gc
      real                  , external  :: solid_area
      !------------------------------------------------------------------------------------!



      !----- Find number of simulation years ----------------------------------------------!
      sim_years = iyearz-iyeara+1
      !------------------------------------------------------------------------------------!


      !----- Crashing the run if the user set up a very long run... -----------------------!
      if (include_fire == 4 .and. sim_years > max_lu_years) then
         write (unit=*,fmt='(a,1x,i5)') 'IYEARA       (From namelist)        :',iyeara
         write (unit=*,fmt='(a,1x,i5)') 'IYEARZ       (From namelist)        : ',iyearz
         write (unit=*,fmt='(a,1x,i5)') 'MAX_LU_YEARS (From disturb_coms.f90): '           &
                                                                              ,max_lu_years
         write (unit=*,fmt='(a)') ' Your run is too long.  Try increasing max_lu_years,'
         write (unit=*,fmt='(a)') ' so MAX_LU_YEARS >= IYEARZ-IYEARA+1, then recompile'
         write (unit=*,fmt='(a)') ' your model.'
         call fatal_error ('Simulation is too long for fire disturbance.'                  &
                          ,'read_fire_ignition','fire_init.f90')
      end if
      !------------------------------------------------------------------------------------!



      !----- Define the format for the header. --------------------------------------------!
      write(hform,fmt='(a,i3.3,a)') '(a',str_len,')'
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !     Find the list of files containing socio-economic indices.                      !
      !------------------------------------------------------------------------------------!
      select case (include_fire)
      case (4)
         !---------------------------------------------------------------------------------!
         !     Loop through all grids.                                                     !
         !---------------------------------------------------------------------------------!
         fire_gridloop: do igr = 1,ngrids
            cgrid=>edgrid_g(igr)


            !..............................................................................!
            !..............................................................................!
            !..............................................................................!
            !..............................................................................!
            !     Socio-economic indices.                                                  !
            !------------------------------------------------------------------------------!



            !------ Extract SEI file list information. ------------------------------------!
            call ed_filelist(full_list,sei_database(igr),nflist)
            call ed1_fileinfo('.sei',nflist,full_list,nflfire,fire_list,llon_list,llat_list)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Loop through polygons.                                                   !
            !------------------------------------------------------------------------------!
            seipoly_loop: do ipy=1,cgrid%npolygons
               cpoly => cgrid%polygon(ipy)



               !---------------------------------------------------------------------------!
               !     Compute the distance between the current polygon and all the files.   !
               !---------------------------------------------------------------------------!
               do nf=1,nflfire
                  file_ldist(nf) = dist_gc(cgrid%lon(ipy),llon_list(nf)                    &
                                          ,cgrid%lat(ipy),llat_list(nf) )
               end do
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !    Pick the closest file.  This is not a guarantee that it will be used   !
               ! because the closest polygon area must contain the point associated to the !
               ! current polygon.                                                          !
               !---------------------------------------------------------------------------!
               ncl       = minloc(file_ldist(1:nflfire),dim=1)
               fire_name = fire_list(ncl)
               write (unit=*,fmt='(2a)') 'Using socio-economic file: ',trim(fire_name)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    Open the socio-economic file and read in all years.                    !
               !---------------------------------------------------------------------------!
               open(unit=12,file=trim(fire_name),form='formatted',status='old'             &
                   ,action='read')
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    Read the header.                                                       !
               !---------------------------------------------------------------------------!
               !----- Westernmost longitude. ----------------------------------------------!
               read (unit=12,fmt=hform) cdum
               cdum = cdum(hoff:)
               read (cdum, fmt=*) wlon
               !----- Easternmost longitude. ----------------------------------------------!
               read (unit=12,fmt=hform) cdum
               cdum = cdum(hoff:)
               read (cdum, fmt=*) elon
               !----- Southernmost latitude. ----------------------------------------------!
               read (unit=12,fmt=hform) cdum
               cdum = cdum(hoff:)
               read (cdum, fmt=*) slat
               !----- Northernmost latitude. ----------------------------------------------!
               read (unit=12,fmt=hform) cdum
               cdum = cdum(hoff:)
               read (cdum, fmt=*) nlat
               !----- Grid cell area (currently not used). --------------------------------!
               read (unit=12,fmt=hform) cdum
               cdum = cdum(hoff:)
               read (cdum, fmt=*) fire_area
               !----- First year with socio-economic data. --------------------------------!
               read (unit=12,fmt=hform) cdum
               cdum = cdum(hoff:)
               read (cdum, fmt=*) yd_1st
               !----- First year with socio-economic data. --------------------------------!
               read (unit=12,fmt=hform) cdum
               cdum = cdum(hoff:)
               read (cdum, fmt=*) yd_last
               !----- Table header, skip it. ----------------------------------------------!
               read (unit=12,fmt=*) 
               !---------------------------------------------------------------------------!


               !----- Determine whether this block contains the current polygon. ----------!
               inside = cgrid%lon(ipy) >= wlon .and. cgrid%lon(ipy) <= elon .and.          &
                        cgrid%lat(ipy) >= slat .and. cgrid%lat(ipy) <= nlat
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Here we will only use the land information if the polygon centre is   !
               ! inside the block of land use disturbance we are about to read.            !
               !---------------------------------------------------------------------------!
               if (inside) then

                  !----- Define the first and last year to be read. -----------------------!
                  yfirst = min(yd_1st ,iyeara)
                  ylast  = max(yd_last,iyearz)
                  iya    = yd_1st  - yfirst + 1
                  iyz    = yd_last - yfirst + 1
                  !------------------------------------------------------------------------!


                  !----- Determine the number of disturbance years. -----------------------!
                  cpoly%num_sei_times(isi) = ylast - yfirst + 1
                  !------------------------------------------------------------------------!



                  !----- Allocate the number of simulation years. -------------------------!
                  allocate(cpoly%seitimes(cpoly%num_sei_times(isi),cpoly%nsites))
                  !------------------------------------------------------------------------!


                  !----- Copy the file information to the first site. ---------------------!
                  isi = 1
                  csite => cpoly%site(isi)
                  !------------------------------------------------------------------------!




                  !------------------------------------------------------------------------!
                  !    First, populate the years with available data.                      !
                  !------------------------------------------------------------------------!
                  do yd_this = yd_1st,yd_last
                     iyr      =  yd_this - yfirst + 1
                     cseitime => cpoly%seitimes(iyr,isi)

                     !----- Read data. ----------------------------------------------------!
                     read(unit=12,fmt=*)  cseitime%sei_year,cseitime%hdi,cseitime%gdpc
                     !---------------------------------------------------------------------!
                  end do
                  !------------------------------------------------------------------------!




                  !------------------------------------------------------------------------!
                  !     For years before the first one with data, copy information from    !
                  ! the first year.                                                        !
                  !------------------------------------------------------------------------!
                  rseitime => cpoly%seitimes(iya,isi)
                  do yd_this = iyeara,(yd_1st-1)
                     iyr      =  yd_this - yfirst + 1
                     cseitime => cpoly%seitimes(iyr,isi)

                     !----- Copy data from first year with data. --------------------------!
                     cseitime%sei_year = yd_this
                     cseitime%hdi      = rseitime%hdi
                     cseitime%gdpc     = rseitime%gdpc
                     !---------------------------------------------------------------------!
                  end do
                  !------------------------------------------------------------------------!





                  !------------------------------------------------------------------------!
                  !     For years after the last one with data, copy information from      !
                  ! the last year.                                                         !
                  !------------------------------------------------------------------------!
                  rseitime => cpoly%seitimes(iyz,isi)
                  do yd_this = (yd_last+1),iyearz
                     iyr      =  yd_this - yfirst + 1
                     cseitime => cpoly%seitimes(iyr,isi)

                     !----- Copy data from first year with data. --------------------------!
                     cseitime%sei_year = yd_this
                     cseitime%hdi      = rseitime%hdi
                     cseitime%gdpc     = rseitime%gdpc
                     !---------------------------------------------------------------------!
                  end do
                  !------------------------------------------------------------------------!


                  !----- Close the land use file, outside the if statement. ---------------!
                  close(unit=12,status='keep')
                  !------------------------------------------------------------------------!





                  !------------------------------------------------------------------------!
                  !      Copy the information from the first site to the other ones.       !
                  !------------------------------------------------------------------------!
                  seiland_siteloop: do isi=2,cpoly%nsites
                     !----- Assume maximum HDI and infinite GDP. --------------------------!
                     cpoly%num_sei_times(isi) = cpoly%num_sei_times(1)
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !     For years after the last one with data, copy information from   !
                     ! the last year.                                                      !
                     !---------------------------------------------------------------------!
                     do iyr = 1,cpoly%num_sei_times(isi)
                        cseitime => cpoly%seitimes(iyr,isi)
                        rseitime => cpoly%seitimes(iyr,1)

                        !----- Copy data from first year with data. -----------------------!
                        cseitime%sei_year = rseitime%sei_year
                        cseitime%hdi      = rseitime%hdi
                        cseitime%gdpc     = rseitime%gdpc
                        !------------------------------------------------------------------!
                     end do
                     !---------------------------------------------------------------------!
                  end do seiland_siteloop
                  !------------------------------------------------------------------------!


               else
                  !------------------------------------------------------------------------!
                  !      No SEI data for this site.  Probably water.                       !
                  !------------------------------------------------------------------------!
                  write (unit=*,fmt='(a)') '----------------------------------------------'
                  write (unit=*,fmt='(a)') ' The closest SEI point is too far away.'
                  write (unit=*,fmt='(a)') ' - File:                     ',trim(fire_name)
                  write (unit=*,fmt=fffmt) ' - Polygon longitude:        ',cgrid%lon(ipy)
                  write (unit=*,fmt=fffmt) ' - Polygon latitude:         ',cgrid%lat(ipy)
                  write (unit=*,fmt=fffmt) ' - Closest central long.:    ',llon_list(ncl)
                  write (unit=*,fmt=fffmt) ' - Closest central lat.:     ',llat_list(ncl)
                  write (unit=*,fmt=fffmt) ' - Closest western edge:     ',wlon
                  write (unit=*,fmt=fffmt) ' - Closest eastern edge:     ',elon
                  write (unit=*,fmt=fffmt) ' - Closest southern edge:    ',slat
                  write (unit=*,fmt=fffmt) ' - Closest northern edge:    ',nlat
                  write (unit=*,fmt=esfmt) ' - Distance:                 ',file_ldist(ncl)
                  write (unit=*,fmt=*)     ' '
                  write (unit=*,fmt='(a)') ' We will assign no SEI-based ignition rate.'
                  write (unit=*,fmt='(a)') '----------------------------------------------'
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Create dummy data set.                                            !
                  !------------------------------------------------------------------------!
                  allocate(cpoly%seitimes  (sim_years,cpoly%nsites))
                  !------------------------------------------------------------------------!



                  !----- Set the parameters in a way that fires will never happen. --------!
                  seiwater_siteloop: do isi = 1,cpoly%nsites
                     !----- Create structure with as many years as the simulation. --------!
                     cpoly%num_sei_times(isi)       = sim_years
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !      Loop through years, and set dummy indices.                     !
                     !---------------------------------------------------------------------!
                     do yd_this = iyeara,iyearz
                        iyr      = yd_this - iyeara + 1
                        cseitime => cpoly%seitimes(iyr,isi)

                        !----- Copy data from first year with data. -----------------------!
                        cseitime%sei_year = yd_this
                        cseitime%hdi      = 1.0
                        cseitime%gdpc     = huge_num
                        !------------------------------------------------------------------!
                     end do
                     !---------------------------------------------------------------------!
                  end do seiwater_siteloop
                  !---------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            end do seipoly_loop
            !------------------------------------------------------------------------------!

            !..............................................................................!
            !..............................................................................!
            !..............................................................................!
            !..............................................................................!






            !..............................................................................!
            !..............................................................................!
            !..............................................................................!
            !..............................................................................!
            !     Lightning data.                                                          !
            !------------------------------------------------------------------------------!




            !------ Extract lightning file list information. ------------------------------!
            call ed_filelist(full_list,flash_database(igr),nflist)
            call ed1_fileinfo('.sei',nflist,full_list,nflfire,fire_list,llon_list,llat_list)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Loop through polygons.                                                   !
            !------------------------------------------------------------------------------!
            flashpoly_loop: do ipy=1,cgrid%npolygons
               cpoly => cgrid%polygon(ipy)



               !---------------------------------------------------------------------------!
               !     Compute the distance between the current polygon and all the files.   !
               !---------------------------------------------------------------------------!
               do nf=1,nflfire
                  file_ldist(nf) = dist_gc(cgrid%lon(ipy),llon_list(nf)                    &
                                          ,cgrid%lat(ipy),llat_list(nf) )
               end do
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !    Pick the closest file.  This is not a guarantee that it will be used   !
               ! because the closest polygon area must contain the point associated to the !
               ! current polygon.                                                          !
               !---------------------------------------------------------------------------!
               ncl       = minloc(file_ldist(1:nflfire),dim=1)
               fire_name = fire_list(ncl)
               write (unit=*,fmt='(2a)') 'Using lightning file: ',trim(fire_name)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    Open the socio-economic file and read in all years.                    !
               !---------------------------------------------------------------------------!
               open(unit=12,file=trim(fire_name),form='formatted',status='old'             &
                   ,action='read')
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    Read the header.                                                       !
               !---------------------------------------------------------------------------!
               !----- Westernmost longitude. ----------------------------------------------!
               read (unit=12,fmt=hform) cdum
               cdum = cdum(hoff:)
               read (cdum, fmt=*) wlon
               !----- Easternmost longitude. ----------------------------------------------!
               read (unit=12,fmt=hform) cdum
               cdum = cdum(hoff:)
               read (cdum, fmt=*) elon
               !----- Southernmost latitude. ----------------------------------------------!
               read (unit=12,fmt=hform) cdum
               cdum = cdum(hoff:)
               read (cdum, fmt=*) slat
               !----- Northernmost latitude. ----------------------------------------------!
               read (unit=12,fmt=hform) cdum
               cdum = cdum(hoff:)
               read (cdum, fmt=*) nlat
               !----- Grid cell area (currently not used). --------------------------------!
               read (unit=12,fmt=hform) cdum
               cdum = cdum(hoff:)
               read (cdum, fmt=*) fire_area
               !----- First year with socio-economic data. --------------------------------!
               read (unit=12,fmt=hform) cdum
               cdum = cdum(hoff:)
               read (cdum, fmt=*) yd_1st
               !----- First year with socio-economic data. --------------------------------!
               read (unit=12,fmt=hform) cdum
               cdum = cdum(hoff:)
               read (cdum, fmt=*) yd_last
               !----- Table header, skip it. ----------------------------------------------!
               read (unit=12,fmt=*) 
               !---------------------------------------------------------------------------!


               !----- Determine whether this block contains the current polygon. ----------!
               inside = cgrid%lon(ipy) >= wlon .and. cgrid%lon(ipy) <= elon .and.          &
                        cgrid%lat(ipy) >= slat .and. cgrid%lat(ipy) <= nlat
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Here we will only use the land information if the polygon centre is   !
               ! inside the block of land use disturbance we are about to read.            !
               !---------------------------------------------------------------------------!
               if (inside) then

                  !----- Define the first and last year to be read. -----------------------!
                  yfirst = min(yd_1st ,iyeara)
                  ylast  = max(yd_last,iyearz)
                  iya    = yd_1st  - yfirst + 1
                  iyz    = yd_last - yfirst + 1
                  !------------------------------------------------------------------------!


                  !----- Currently, flash times are always 12. ----------------------------!
                  cpoly%num_flash_times(isi) = 12
                  !------------------------------------------------------------------------!



                  !----- Allocate the number of simulation years. -------------------------!
                  allocate(cpoly%flashtimes(cpoly%num_flash_times(isi),cpoly%nsites))
                  !------------------------------------------------------------------------!


                  !----- Copy the file information to the first site. ---------------------!
                  isi = 1
                  csite => cpoly%site(isi)
                  !------------------------------------------------------------------------!




                  !------------------------------------------------------------------------!
                  !    Currently, this should always have 12 entries, one for each month.  !
                  !------------------------------------------------------------------------!
                  do imo=1,12
                     cflashtime => cpoly%flashtimes(imo,isi)

                     !----- Read data. ----------------------------------------------------!
                     read(unit=12,fmt=*)  cflashtime%flash_month,cflashtime%frd            &
                                         ,cflashtime%c2g
                     !---------------------------------------------------------------------!
                  end do
                  !------------------------------------------------------------------------!


                  !----- Close the land use file, outside the if statement. ---------------!
                  close(unit=12,status='keep')
                  !------------------------------------------------------------------------!





                  !------------------------------------------------------------------------!
                  !      Copy the information from the first site to the other ones.       !
                  !------------------------------------------------------------------------!
                  flashland_siteloop: do isi=2,cpoly%nsites
                     !----- Assume maximum HDI and infinite GDP. --------------------------!
                     cpoly%num_flash_times(isi)       = cpoly%num_flash_times(1)
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !     For years after the last one with data, copy information from   !
                     ! the last year.                                                      !
                     !---------------------------------------------------------------------!
                     do imo = 1,12
                        cflashtime => cpoly%flashtimes(imo,isi)
                        rflashtime => cpoly%flashtimes(imo,1)

                        !----- Copy data from first year with data. -----------------------!
                        cflashtime%flash_month = rflashtime%flash_month
                        cflashtime%frd         = rflashtime%frd
                        cflashtime%c2g         = rflashtime%c2g
                        !------------------------------------------------------------------!
                     end do
                     !---------------------------------------------------------------------!
                  end do flashland_siteloop
                  !------------------------------------------------------------------------!


               else
                  !------------------------------------------------------------------------!
                  !      No SEI data for this site.  Probably water.                       !
                  !------------------------------------------------------------------------!
                  write (unit=*,fmt='(a)') '----------------------------------------------'
                  write (unit=*,fmt='(a)') ' The closest lightning point is too far away.'
                  write (unit=*,fmt='(a)') ' - File:                     ',trim(fire_name)
                  write (unit=*,fmt=fffmt) ' - Polygon longitude:        ',cgrid%lon(ipy)
                  write (unit=*,fmt=fffmt) ' - Polygon latitude:         ',cgrid%lat(ipy)
                  write (unit=*,fmt=fffmt) ' - Closest central long.:    ',llon_list(ncl)
                  write (unit=*,fmt=fffmt) ' - Closest central lat.:     ',llat_list(ncl)
                  write (unit=*,fmt=fffmt) ' - Closest western edge:     ',wlon
                  write (unit=*,fmt=fffmt) ' - Closest eastern edge:     ',elon
                  write (unit=*,fmt=fffmt) ' - Closest southern edge:    ',slat
                  write (unit=*,fmt=fffmt) ' - Closest northern edge:    ',nlat
                  write (unit=*,fmt=esfmt) ' - Distance:                 ',file_ldist(ncl)
                  write (unit=*,fmt=*)     ' '
                  write (unit=*,fmt='(a)') ' We will assign no lightning-based ignition.'
                  write (unit=*,fmt='(a)') '----------------------------------------------'
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Create dummy data set.                                            !
                  !------------------------------------------------------------------------!
                  allocate(cpoly%flashtimes  (12,cpoly%nsites))
                  !------------------------------------------------------------------------!



                  !----- Set the parameters in a way that fires will never happen. --------!
                  flashwater_siteloop: do isi = 1,cpoly%nsites
                     !----- Create structure with as many years as the simulation. --------!
                     cpoly%num_flash_times(isi)       = 12
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !      Loop through years, and set dummy indices.                     !
                     !---------------------------------------------------------------------!
                     do imo = 1,12
                        cflashtime => cpoly%flashtimes(imo,isi)

                        !----- Copy data from first year with data. -----------------------!
                        cflashtime%flash_month = imo
                        cflashtime%frd         = 0.0
                        cflashtime%c2g         = 0.0
                        !------------------------------------------------------------------!
                     end do
                  end do flashwater_siteloop
                  !---------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            end do flashpoly_loop
            !------------------------------------------------------------------------------!

            !..............................................................................!
            !..............................................................................!
            !..............................................................................!
            !..............................................................................!

         end do fire_gridloop
         !---------------------------------------------------------------------------------!

      case default
         !---------------------------------------------------------------------------------!
         !     Loop through all grids.                                                     !
         !---------------------------------------------------------------------------------!
         nofire_gridloop: do igr = 1,ngrids
            cgrid=>edgrid_g(igr)

            !------------------------------------------------------------------------------!
            !     Socio-economic indices are not needed. Set a single year with no data.   !
            !------------------------------------------------------------------------------!
            nofire_polyloop: do ipy = 1,cgrid%npolygons
               cpoly => cgrid%polygon(ipy)

               !---------------------------------------------------------------------------!
               !      Fire ignition data not needed, allocate only a single SEI year.  We  !
               ! allocate 12 months for flash rates, but make the climatology zero.        !
               !---------------------------------------------------------------------------!
               allocate(cpoly%seitimes  ( 1,cpoly%nsites))
               allocate(cpoly%flashtimes(12,cpoly%nsites))
               !---------------------------------------------------------------------------!

               !----- Set the parameters in a way that fires will never happen. -----------!
               nofire_siteloop: do isi = 1,cpoly%nsites
                  !----- Assume maximum HDI and infinite GDP. -----------------------------!
                  cpoly%num_sei_times(isi)       = 1
                  cpoly%seitimes(1,isi)%sei_year = iyeara
                  cpoly%seitimes(1,isi)%hdi      = 1.0
                  cpoly%seitimes(1,isi)%gdpc     = huge_num
                  !------------------------------------------------------------------------!


                  !----- Assume lightning flashes do not occur. ---------------------------!
                  cpoly%num_flash_times(isi)     = 12
                  do imo=1,12
                     cpoly%flashtimes(imo,isi)%flash_month = imo
                     cpoly%flashtimes(imo,isi)%frd         = 0.0
                     cpoly%flashtimes(imo,isi)%c2g         = 0.0
                  end do
                  !------------------------------------------------------------------------!
               end do nofire_siteloop
               !---------------------------------------------------------------------------!
            end do nofire_polyloop
            !------------------------------------------------------------------------------!
         end do nofire_gridloop
         !---------------------------------------------------------------------------------!
      end select
      !---------------------------------------------------------------------------------!

      return
   end subroutine read_fire_ignition
   !=======================================================================================!
   !=======================================================================================!
end module fire_init
!==========================================================================================!
!==========================================================================================!
