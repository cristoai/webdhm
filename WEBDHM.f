      module webdhm_core
         use metdata
         use output_sptl
         use damops ! Moiz: Damops
         private

         public webdhm


!#######################################################################
!
!     Topographic variable declarations.
!
!#######################################################################
!

! Static
         real, allocatable :: x  (:)            ! The x-coord. of the computational grid (m).    [nx]
         real, allocatable :: y  (:)            ! The y-coord. of the computational grid (m).    [ny]
! These are now defined in metdata (PRL)
!         real, allocatable :: xlon(:,:)         ! longtude at (x,y)                              [nx,ny]
!         real, allocatable :: ylat(:,:)         ! latitude at (x,y)                              [nx,ny]
         integer, allocatable :: inbasin(:,:)   ! flag to show the simulation boundary           [nx,ny]
         real, allocatable :: ele2d(:,:)        ! elevation                                      [nx,ny]
         real, allocatable :: sta_alt(:,:)      ! altitude of meteorological gauge station       [nx,ny]

         real, allocatable :: slp2d(:,:)        ! The slope of ground surface (ND)               [nx,ny]
         real, allocatable :: len2d(:,:)        ! hillslope length of each GRID (m)              [nx,ny]
         real, allocatable :: area(:,:)         ! area of each GRID (m2)                         [nx,ny]

         real, allocatable :: Ds2d(:,:)         ! depth of UZ soil of each GRID (m)              [nx,ny]
!> \todo Is it true that the depth to aquifer is fixed/static?
         real, allocatable :: Dg2d(:,:)         ! depth of unconfined aquifer of each GRID (m)   [nx,ny]

         real, allocatable :: Dr2d(:,:)         ! river depth in each GRID (m)                   [nx,ny]


!
!#######################################################################
!
!     Atmospheric forcing data
!
!#######################################################################
!
!
         real, allocatable :: zref(:,:)         ! The reference level of atmospheric forcing
! data from ground surface (m).
!         real    uv   (:,:)      ! wind velocity at reference level (m/s)
!         real    qv   (:,:)      ! specific humidity at reference level (kg/kg)
!         real    tair (:,:)      ! Air temperature at refernce level (K)
         real, allocatable :: diag_prc  (:,:)      ! convective precipitation rates (kg/(m**2*s))
         real, allocatable :: diag_prl  (:,:)      ! large-scale precipitation rates (kg/(m**2*s))
!         real    radsw(:,:)      ! Downward solar shortwave irradiance at bottom (W/m**2)
!         real    radlw(:,:)      ! Downward longwave irradiance at bottom (W/m**2)
!         real    cld  (:,:)      ! Cloud fraction
!         real    psfc (:,:)      ! surface air pressure
!         real    rh   (:,:)      ! relative humidity
!
!     secondary: rho
!
      real, allocatable :: rho  (:,:)          ! Surface air density (kg/m**3)

!
!#######################################################################
!
!     Vegetation parameters
!
!#######################################################################
!
!     primary
!
      integer, allocatable :: soiltyp  (:,:)    ! Soil type at each point
!      type(soil_params) :: soil(ns) ! info for each soil type
      real, allocatable ::    vtypfrct (:,:,:)  ! Fraction of vegetation types
      integer, allocatable :: vegtyp   (:,:,:)  ! Vegetation type at each point
      real, allocatable ::    ndvi     (:,:,:,:)! ndvi
!
!     Derived variables
!


      real, allocatable ::    roufns  (:,:,:) !< Surface roughness
      real, allocatable ::    green2d (:,:,:) !< Green leaf fraction
      real, allocatable ::    gmudmu2d(:,:,:) !< Time-mean radiation-weighted leaf projection (g(mu)/mu)

!
!#######################################################################
!
!     Prognostic variables
!
!#######################################################################
!

         real, allocatable :: Sstmax0(:,:)



!#######################################################################
!
!     State variables that are shared by multiple models
!
!#######################################################################

         real, allocatable :: Sst(:,:)
         real, allocatable :: Drw2d(:,:)        ! water depth in the river of each GRID (m)      [nx,ny]


!
!#######################################################################
!
!     LSM result variables
!
!#######################################################################
!

         real, allocatable :: runoff_inter(:,:)
         real, allocatable :: runoff_g(:,:)
         real, allocatable :: waterflx(:,:)


         type(meteorological_data), allocatable :: diagf(:,:)


      contains




!
!     ##################################################################
!     ##################################################################
!     ######                                                      ######
!     ######             Water & Energy Budget based              ######
!     ######            Distributed Hydrological Model            ######
!     ######                    (WEB-DHM)                         ######
!     ######                                                      ######
!     ######               Author:  WANG Lei                      ######
!     ######                                                      ######
!     ######                  Sep. 2007                           ######
!     ######                                                      ######
!     ######     River and Environmental Engineering Laboratory   ######
!     ######                University of Tokyo                   ######
!     ######                                                      ######
!     ##################################################################
!     ##################################################################

      subroutine WEBDHM
         use init_params
         use soilinit
         use hydro
         use forcing, only: init_read_meteo
         

         implicit none

         include 'dims.inc'
         include 'globcst.inc'
         include 'hydro.inc' ! inisub, finsub, numsub
         include 'incl_params.f'

         

! Separate logging file disabled (PRL)
!         OPEN(6,file='output.txt',status = 'unknown')

         CALL initpara

         if (inicon == -1) then
            inicon = INICON_CFG
         end if

         if (nx < 1) then
            call init_extents(trim(para_dir)//gridarea_map,nx,ny)
         end if

         ! Variables should be allocated before initializing data I/O
         call allocate_variables(nx,ny,nvtyps)

         allocate( diagf(nx,ny) )

         call init_projection(nx,ny,x,y)                    ! Necessary prior to setting up JRA forcing input
         call init_read_meteo(nx,ny,METEO_INPUT_TYPE,diagf) ! should be called after initpara

! Pre-initialization of subcatchments for GBHM
         call define_subcatchments
         inisub = fixed_inisub
         finsub = fixed_finsub
         if (finsub > num_sub) finsub = num_sub

         ! Read Dam and Hydropower Data
         !if (use_damops==0) then
         call init_damops(use_damops,dam_snapshot) ! Moiz: Damops
         !endif

         call WEBDHM_post_infile
!(nx,ny,nvtyps)
      end subroutine WEBDHM





      subroutine allocate_variables(nx,ny,nvtyps)
         use mod_flux2d, only:
     &                  lsm_allocate_variables => allocate_variables
         implicit none

! ===== PARAMETERS

         integer, intent(in) :: nx
         integer, intent(in) :: ny
         integer, intent(in) :: nvtyps


! ===== CODE

! Static
         allocate( x  (nx) )
         allocate( y  (ny) )
         allocate( xlon(nx,ny) )
         allocate( ylat(nx,ny) )
         allocate( slp2d(nx,ny) )
         allocate( len2d(nx,ny) )
         allocate( area(nx,ny) )
         allocate( Ds2d(nx,ny) )
!> \todo Is it true that the depth to aquifer is fixed/static?
         allocate( Dg2d(nx,ny) )
         allocate( Dr2d(nx,ny) )
         allocate( inbasin(nx,ny) )
         allocate( ele2d(nx,ny) )
         allocate( sta_alt(nx,ny) )

         allocate( Sstmax0(nx,ny) )

! Temporal
         allocate( Drw2d(nx,ny) )

!> This is a state variable shared by multiple models
         allocate( Sst(nx,ny) )

!#######################################################################
!     Atmospheric forcing data

         allocate( zref(nx,ny) )
! data from ground surface (m).
         allocate( diag_prc  (nx,ny) )
         allocate( diag_prl  (nx,ny) )
!
!     secondary: rho
!
         allocate( rho  (nx,ny) )


!#######################################################################
!     Vegetation parameters

!#######################################################################
!     primary

         allocate( soiltyp(nx,ny) )
         allocate( vtypfrct(nx,ny,nvtyps) )
         allocate( vegtyp (nx,ny,nvtyps) )
         allocate( ndvi   (nx,ny,nvtyps,12) )


!#######################################################################
!     Derived variables

         allocate( roufns (nx,ny,nvtyps) )
         allocate( green2d (nx,ny,nvtyps) )
         allocate( gmudmu2d(nx,ny,nvtyps) )


!#######################################################################
!     Prognostic variables

         allocate( runoff_inter(nx,ny) )
         allocate( runoff_g(nx,ny) )
         allocate( waterflx(nx,ny) )

         call lsm_allocate_variables(nx,ny,nvtyps)

      end subroutine allocate_variables




!
!#######################################################################
!
!     PURPOSE:
!
!     A coupling model of a revised simple
!     biosphere model and grid-based hydrological model.
!
!#######################################################################

      subroutine WEBDHM_post_infile
!(nx,ny,nvtyps)




!#######################################################################
!
!     Variable Declarations.
!
!#######################################################################
         use init_params
         use snapshot
         use snapshot_recv
!> \todo: check whether opc should be visible here
         use mod_flux2d, only: diagnostic_fluxes,
     &                  flux2d, init_flux2d, initsoil, vegtable,
     &                  opc, dgl2d,lai,fpar
         use sib2_public_vars
         use sib2_soil
         use topoinit
         use river_rtg, only: init_river_routing
         use hydro
         use metdata
         use output
         use output_bin_mod

         


         implicit none


! ===== INCLUDES

         include 'dims.inc'
         include 'globcst.inc'
         include 'hydro.inc'
         include 'incl_params.f'



! ===== LOCAL VARIABLES


!#######################################################################
! Diagnostic variables for land-surface model

      type(diagnostic_fluxes), allocatable :: diag2(:,:)


!#######################################################################
! Diagnostic variables for riverflow model


!#######################################################################
!     Output fluxes

      integer     ii, i, j, loop, ix,iy,k !,m,flag
      real        maxflow, totflow, tflow
      integer     maxi,maxj

      integer     flsoil

      character*4 chy
      character*2 chm,chd
      integer     outflag


!#######################################################################
!     Transient parameters for GBHM

      real q2(nc,300)


!#######################################################################
!     Misc. local variables:

      real    eps
      data eps/1.0e-6/
      integer lmet

      real         qin(300)

      character   recv_comm*256


 
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!     Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!************************************************************************


      write(6,'(/ 17(/5x,a)//)')

     :'###############################################################',
     :'###############################################################',
     :'#####                                                     #####',
     :'#####                      Welcome to                     #####',
     :'#####                                                     #####',
     :'#####       SiB2-GBHM coupled Model '//__SWVERSION__//
     &                                               '           #####',
     :'#####                      (SiB2-GBHM)                    #####',
     :'#####                                                     #####',
     :'#####                     Developed by                    #####',
     :'#####    River and Environmental Engineering Laboratory   #####',
     :'#####               University of Tokyo                   #####',
     :'#####                                                     #####',
     :'###############################################################',
     :'###############################################################'


!#######################################################################
! Allocate additional variables

      allocate( diag2(nx,ny) )


!#######################################################################
!
!     Initialize the model parameters. Most of them are global parameters
!     passed through common blocks.
!
!#######################################################################
!
      JustStart = 1
!      CALL initpara

!
!#######################################################################
!
!     Initialize look-up table related to soil type and vegetation type
!
!     These parameters are declared in the include files: lsmpar.inc
!
!#######################################################################
!
      CALL VEGTABLE(opc)
      CALL SOILTABLE
!
!#######################################################################
!
!     Initialize topography-related arrays
!
!#######################################################################
!
      CALL inittopo(nx,ny,ele2d,
     :     slp2d,len2d,area,Ds2d,Dg2d,inbasin)


!     open(1, file = 'xgrid.txt')
!     do i = 1, nx
!     write(1, "(i10, f15.2)") i, x(i)
!     end do
!     close(1)

!     open(1, file = 'ygrid.txt')
!     do j = 1, ny
!     write(1, "(i10, f15.2)") j, y(j)
!     end do
!     close(1)
!#######################################################################
!
!     Initialize hydrological parameter: 03/29/2006  Wang Lei
!
!#######################################################################

      call init_hydro

!#######################################################################
!
!     Initialize land and soil parameters and arrays for 2D flux
!     calculation
!
!#######################################################################
!

      CALL initsoil(
     &      nx,ny,nz,nvtyps, ! xlon,ylat,
     &      Sst, Sstmax0,
     &      ele2d,sta_alt,
     &      Ds2d,Dg2d,
     &      inbasin,
     &      zref,
     &      soiltyp,vtypfrct,vegtyp,ndvi,
     &      soil,diagf)

!#######################################################################
!
!     river initialization
!
!#######################################################################
      call init_river(Ds2d, Dr2d,Drw2d)
      do ii=inisub,finsub
         call surface_routing(ii,nx,ny,nvtyps,
     :        inbasin,area,slp2d,len2d,
     :        runoff_inter,runoff_g,
     :        Sst,Sstmax0,waterflx, qin)

         call init_river_routing(ii)
      end do
      JustStart = 0


! Spatial output should be initialized before flux2d so that
! variables are found correctly
         call init_output_sptl()
         call init_flux2d(inbasin)
         call init_grid_output(diag2,diagf)
         call init_output()

         call register_fast(Sst,'sst',2)
         call register_fast(waterflx,'waterflx',2)
         call register_fast(runoff_inter,'roff_inter',2)
         call register_fast(runoff_g,'roff_gw',2)
         call register_fast(drw2d,'dpth_rw',2)


!#######################################################################
!
!     Make some loops to initialize model input
!
!#######################################################################

      do loop = 1, 1            !3/18/2007: Calibrate

!#######################################################################
!
!     The current model time (initial time):
!
!#######################################################################
!
         curtim = tstart

         write(6,'(1x,a,i15,a)')
     :        'The initial model time is at ', curtim,' seconds.'

         nstep = curtim/dt_couple


! surface roughness changes as canopy cover changes

! Snapshot will only be loaded if the path is properly
! initialized.  An unintialized path will disable snapshot loading.
          call load_snapshot(
     &         nc,nx,ny,nz,nvtyps,  ! Not part of the record (only used to define size)
     &         Sst,waterflx
     &         )

         if( recv_mode.eq.1 ) then                                                                         
            call load_snapshot_recv(                                                                       
     &           nc,nx,ny,nz,nvtyps,  ! Not part of the record (only used to define size)                  
     &           Sst,waterflx                                                                              
     &           )                                                                                         
           endif 

         call distribute_river_water_depth(drw2d)

!#######################################################################
!
!     Time integration loop begins  ----------------------------------->
!
!#######################################################################
!

         do while(curtim+eps*dt_couple < tstop)

!
!#######################################################################
!
!     Define the time step counter nstep reference to time zero
!     (curtim=0).
!
!#######################################################################
!
            nstep   = nstep + 1
            abstcur = abstini + curtim

            CALL ABSS2CTIM(
     &               abstcur, year, month, day, hour, minute,second)
!     Convert the absolute time(sec) starting from 00:00:00 UTC, Jan.
!     1, 1960, to calendar day time.
            CALL JULDAY( year, month, day, jday )
!     Compute Julian day from year, month, and day
!     Start from 1 (Jan. 1) to 365, or 366 for leap year (Dec. 31)

!#######################################################################
!
!     Perform one time step integration for all equations:
!     On exit of this routine, all time dependent fields are
!     advanced by one time step.
!
!#######################################################################
!
            CALL flux2d(
     &         x,y, ! xlon,ylat,
     &         zref,slp2d,len2d,
     &         ele2d, sta_alt,
     &         Dg2d,Ds2d,Dr2d,Drw2d,
     &         inbasin,
     &         diagf,
     &         diag_prc,diag_prl,rho,
     &         soiltyp,vtypfrct,vegtyp,
     &         diag2,
     &         Sst,Sstmax0,
     &         waterflx, runoff_inter, runoff_g    ! OUTPUT
     &         )

! Parameters passed from flux2d to lateral_routing:
! runoff_inter, runoff_g, Sst, waterflx
! Sstmax0 is constant

            if (POINTSCALE_X == -1 .and. POINTSCALE_Y == -1) then
               call write_fast_sptl(2)
            end if

            CALL lateral_routing(
     &         inbasin,area,slp2d,len2d,Dr2d,               ! static topographic input
     &         runoff_inter,runoff_g,Sst,Sstmax0,waterflx,  ! INPUT
     &         q2,Drw2d                                     ! OUTPUT
     &     )




!#######################################################################
!     Update physical time and generate output files
          ! Moiz Damops
          ! This flag's role need to be clarified later
          if (use_damops>0) then
          !call damops_output
      
          ! Set values for 'previous' timestep 
          ! for calculation in next timestep
      
          ! Dam State Variables
          do idam=1,size(dam)
             call prep_dam(idam)
          enddo
      
          do ihpp=1,size(hpp)
             call prep_hpp(ihpp)
          enddo
      
          endif
          ! Moiz: Damops

            curtim = curtim + int(dt_couple)

!      call write_tprlavg_grid_data(pg3,diag2)
!      call write_tprlavg_flowint_data(q2)
!      call write_sptlavg_upstream_data(pg3,diag2)

            call write_sptlavg_upstream_data(inbasin,
     &                  diag2, diagf,lai,fpar,Dgl2d,diag_prc,diag_prl)
            call write_outflow_data(q2)
            call write_tprlavg_grid_data

            if (save_recv) then ! Moiz: Added
            if( hour.eq. 0 ) then    ! RECV
               write(recv_comm,"(a,4i8)")
     &            'cd '//trim(recv_dir)//'; sh rotate.sh',
     &            year, month, day, hour
               call system( recv_comm )
               call save_snapshot_recv(nc,nx,ny,nz,nvtyps,
     &               Sst,waterflx
     &               )
            endif
            endif





         end do ! while(curtim+eps*dt_couple < tstop)

         call save_snapshot(nc,nx,ny,nz,nvtyps,
     &         Sst,waterflx
     &         )

!
!#######################################################################
!
!     End of entire model time integration. The program stops.
!     Close opened files.
!
!#######################################################################
!

         write(6,'(//1x,a,/1x,a,i15,a,/1x,a)')
     &        'WEB-DHM stopped normally in the main program. ',
     &        'The ending time was ', curtim,' seconds.',
     &        'Thanks for using WEB-DHM.'

!#######################################################################
         end do                    ! loop
!#######################################################################

         close(6)

!     STOP
         return
      END subroutine WEBDHM_post_infile


      end module webdhm_core
