      module mod_flux2d
         use constants
         use time_and_fileio
         use lsmpar
         use vertical_profile,
     &            only: vertical_profile_in,vertical_profile_inout
         use sib2_public_vars
         use sib2_soil

! RASMY ADD WEB-S
         !use bats_vars, only: bats_state_t => state_t
         use bats, only:
     &            bats_initialize => initialize,   ! Moiz: Snapshot
     &            bats_state => state

         !use snow_vars, only: snow_state_t => state_t
         use snow, only:
     &            snow_initialize => initialize,
     &            snow_state => state



         private

! ======================================
! Initialization routines
         public :: init_flux2d
         public :: allocate_variables
         public :: initsoil
         public :: vegtable

         public :: flux2d
         public :: diagnostic_fluxes
         public :: diagnostic3d

!> These variables are made public for diagnostic purposes
         public :: lai
         public :: fpar
         public :: dgl2d


         public :: get_snapshot_size
         public :: get_snapshot_size_swe
         public :: save_snapshot
         public :: save_snapshot_swe
         public :: load_snapshot
         public :: load_snapshot_swe


         public :: in_analysis


         type(optical), save :: opc(nv)
         public :: opc

         logical, allocatable :: in_analysis(:,:)  !< flag to show the analyzed region

! These variables are used to support loop flattening for OpenMP speedup
         integer, allocatable :: omp_subbasin(:)
         integer, allocatable :: omp_flowinterval(:)
         integer :: isubflow_cnt

! For each sub-grid land types
         type diagnostic3d
            real wsfc      !< surface zone soil moisture content
            real wrt       !< root zone soil moisture content
            real wdp       !< deep zone soil moisture content
         end type diagnostic3d

! Avg for all land types in grid point
         type diagnostic_fluxes
            real wsfc      !< surface zone soil moisture content
            real wrt       !< root zone soil moisture content
            real wdp       !< deep zone soil moisture content

            real ptsflx    !< surface flux of heat (K*kg/(m**2*s))
            real qvsflx    !< surface flux of moisture (kg/(m**2*s))
            real Rnflx     !< net radiation
            real Ect_flx   !< canopy transpiration                      (mm)
            real Eci_flx   !< evap from canopy interception     (mm)
            real Egs_flx   !< evap from surface soil layer    (mm)
            real Egi_flx   !< evap from ground interception   (mm)

            real tsfc

            real tsoil
            real snowc
            real snowg
            real wetcanp
            real wetg
            real roff_i
            real roff_g
            real dgl       !< depth to groundwater level at current gridsquare (m)
            real drw       !< depth of river water at current gridsquare (m)

            real snowdens  ! Moiz: Added
            real snowdep   ! Moiz: Added

         end type diagnostic_fluxes



!#######################################################################
! Declaration of state variables

      type(remainder3dgrid), allocatable :: r3(:,:,:)

!> Vegetation variables
      real, allocatable ::    vegcover_mo (:,:,:)  !< monthly vegetation coverage
      real, allocatable ::    vcover     (:,:,:) !< (Instantaneous) Vegetation coverage
      real, allocatable ::    lai     (:,:,:) !< leaf area index
      real, allocatable ::    fpar    (:,:,:) !< the Fraction of Phyotosynthetically Active Radiation absorbed by the green vegetation canopy

!> Soil variables
      integer, allocatable :: nroot   (:,:,:)   !< number of layers in root zone
      integer, allocatable :: ndeep   (:,:,:)   !< number of layers in deep soil
      real, allocatable :: Dgl2d(:,:,:)         !< depth to groundwater level in each GRID (m)         [nx,ny]
      real, allocatable :: GWstorage(:,:,:)     !< volume of water in acquifer in each GRID     [nx,ny,nvtyp]

!> Prognostic variables
      type(prognostic3dgrid), allocatable :: pg3(:,:,:)

!> Water variables
      real, allocatable :: Sstmax0(:,:) !< TODO: How can tstmax change over time?
!RASMY WEB-S
      !type(snow_state_t), allocatable :: snow_state(:,:)
      !type(bats_state_t), allocatable :: bats_state(:,:)


!#######################################################################
! Precipitation data for current timestep

         real, allocatable, target :: rainh(:,:)          ! houry rainfall (mm/hour);
         real, pointer :: rainptr(:,:)

         ! Moiz: Added for READ_SEPARATE_RAIN_SNOW
         real, allocatable, target :: rainh_liquid(:,:)  ! hourly liquid precipitation (mm/hour)
         real, allocatable, target :: rainh_solid(:,:)   ! hourly solid precipitation (mm/hour)

         real, pointer :: rainptr_liquid(:,:)
         real, pointer :: rainptr_solid(:,:)


!#######################################################################
! Internal variables for diagnostics

         type(diagnostic3d), allocatable :: diag3(:,:,:) ! [NX,NY,NV]

         integer :: wsfc_hdl,wrt_hdl,wdp_hdl,dgl_hdl,swe_hdl
         integer :: eci_hdl,ect_hdl,egi_hdl,egs_hdl
         integer :: snowdep_hdl, snowdens_hdl   ! Moiz: Snow Vars
         integer :: rainh_hdl ! Moiz: precipitation
      contains



      subroutine init_flux2d(inbasin)
         use output_sptl
         use forcing, only: init_read_precip
         implicit none

! ===== INCLUDES
         include 'dims.inc'   !< nx,ny
         include 'hydro.inc'  !< num_flow,ngrid,grid_col,grid_row
         include 'incl_params.f' !< READ_SEPARATE_RAIN_SNOW

! ===== PARAMETERS
         integer, intent(in) :: inbasin(nx,ny)     !< flag to show the basin region

! ===== LOCAL VARIABLES
         integer i,j
         integer isub,iflow,ig
         integer tot_flowints

! ===== CODE
!      open(1001,file ='output/other/kurobe_pt.hourly',status='unknown')
!      open(1002,file='output/other/kurobe_pt_problem.hourly',
!     &                status='unknown')
!      open(1003,file='output/other/kurobe_pt_problem1.hourly',
!     &                status='unknown')
         open(123,   file='checkpoint.csv',
     &               status='unknown', action="write")
      write(123,'(" YY MM DD HH isat '//
     & 'RT DP SFCCAP SFCSNOW Tc Tg Td VLC1 VLC2 VLC3 '//
     & 'Dgl Hflx ETmass ECT ECI EGS EGI ROFFi ROFFg '//
     & 'rech qsub gwst s01 s02 s03 s04 s05 s06 s07 s08 s09 s10 '//
     & 's11 s12 s13 s14 s15 s16 s17 s18 s19 s20 s21 '//
     & 'ppl um Tm psur rad11 rad12 rad21 rad22 rad31 rad32 Drw")')

!#######################################################################
! Allocations for diagnostic and internal tracking purposes
         allocate( diag3(nx,ny,nv) )
         allocate(rainh(nx,ny))
         rainptr => rainh
         call init_read_precip(rainh)

         ! Moiz: Added for READ_SEPARATE_RAIN_SNOW
         if (READ_SEPARATE_RAIN_SNOW == .true.) then
            allocate(rainh_liquid(nx,ny))
            allocate(rainh_solid(nx,ny))

            rainptr_liquid => rainh_liquid
            rainptr_solid => rainh_solid
         endif

!#######################################################################
! Set up flag for whether gridpoint is within analysis region
! (analysis region may be a subset of entire basin)

         allocate( in_analysis(nx,ny) )

! RASMY WEB-S
         allocate (snow_state(nx,ny))
         allocate (bats_state(nx,ny))

         DO i = 1, nx
            DO j = 1, ny
               in_analysis(i,j) = .false.
! RASMY WEB-S
               call bats_initialize(bats_state(i,j))
               call snow_initialize(snow_state(i,j))
            end do
         end do


 



         tot_flowints = 0
         do isub = inisub, finsub
            do iflow = 1,num_flow(isub)
               tot_flowints = tot_flowints + 1
               do ig=1,ngrid(isub,iflow)
                  i=grid_col(isub,iflow,ig)
                  j=grid_row(isub,iflow,ig)
                  in_analysis(i,j) = .true.
               end do
            end do
         end do
         isubflow_cnt = tot_flowints

!#######################################################################
! Set up OpenMP [sub-basin,flowint] -> [INDEX] 2d -> 1d loop mapping

         allocate(omp_subbasin(tot_flowints))
         allocate(omp_flowinterval(tot_flowints))
         tot_flowints = 0
         do isub = inisub, finsub
            do iflow = 1,num_flow(isub)
               tot_flowints = tot_flowints + 1
               omp_subbasin(tot_flowints) = isub
               omp_flowinterval(tot_flowints) = iflow
            end do
         end do

         wsfc_hdl = register_slow('wsfc')
         wrt_hdl  = register_slow('wrt')
         wdp_hdl  = register_slow('wdp')

         dgl_hdl  = register_slow('dgl')
         swe_hdl  = register_slow('swe')

         eci_hdl  = register_slow('eci')
         ect_hdl  = register_slow('ect')
         egi_hdl  = register_slow('egi')
         egs_hdl  = register_slow('egs')

         ! Moiz: Snow vars
         snowdep_hdl = register_slow('snowdep')
         snowdens_hdl = register_slow('snowdens')

         ! Moiz: Precipitation
         rainh_hdl = register_slow('precipitation')


! RASMy NEED TO ADD register datA

         return
      end subroutine init_flux2d


      subroutine allocate_variables(nx,ny,nvtyps)
         implicit none

         integer, intent(in) :: nx
         integer, intent(in) :: ny
         integer, intent(in) :: nvtyps

         allocate( r3   (nx,ny,nvtyps) )

! Vegetation
         allocate( vegcover_mo(nx,ny,12) )
         allocate( vcover    (nx,ny,nvtyps) )
         allocate( lai    (nx,ny,nvtyps) )
         allocate( fpar   (nx,ny,nvtyps) )

! Soil
         allocate( nroot(nx,ny,nvtyps) )
         allocate( ndeep(nx,ny,nvtyps) )
         allocate( Dgl2d(nx,ny,0:nvtyps) )
         allocate( GWstorage(nx,ny,0:nvtyps) )

         allocate( pg3(nx,ny,0:nvtyps) )

         return
      end subroutine allocate_variables


      subroutine initsoil(
     &      nx,ny,nz,nvtyps, !xlon,ylat,
     &      Sst, Sstmax0,
     &      ele2d,sta_alt,
     &      Ds2d,Dg2d,
     &      inbasin,
     &      zref,
     &      soiltyp,vtypfrct,vegtyp,ndvi,
     &      soil,frc)
         use soilinit, only: sib2_initsoil => INITSOIL
         use metdata
         implicit none

      type(meteorological_data), target :: frc(:,:)

      integer, intent(in) :: nx           !< number of grids in x-dir.
      integer, intent(in) :: ny           !< number of grids in y-dir.
      integer, intent(in) :: nz           !< number of grids in z-dir of soil model.
      integer, intent(in) :: nvtyps       !< number of vegetation types
! These are now defined in metdata (PRL)
!      real, intent(in) :: xlon(nx,ny)     !< longitude at (x,y)
!      real, intent(in) :: ylat(nx,ny)     !< latitude at (x,y)

      real, intent(out) :: Sst(nx,ny)
      real, intent(out) :: Sstmax0(nx,ny)

      real, intent(in) :: ele2d(nx,ny)    !< Per-point land elevation <b>(m)</b>
      real, intent(out) :: sta_alt(nx,ny)

      integer, intent(in) :: inbasin(nx,ny) !< flag to show the simulation boundary


!#######################################################################
!     Soil parameters

      real, intent(in) :: Ds2d(nx,ny)              !< depth of UZ soil (m)
      real, intent(in) :: Dg2d(nx,ny)              !< depth of UZ soil (m)


!#######################################################################
!     Vegetation parameters

!>     primary
      integer, intent(out) :: soiltyp(nx,ny)       !< Soil type at each point
      real, intent(out) :: vtypfrct(nx,ny,nvtyps)  !< Fraction of vegetation types
      integer, intent(out) :: vegtyp(nx,ny,nvtyps) !< Vegetation type at each point
      real, intent(in) :: ndvi(nx,ny,nvtyps,12)    !< ndvi(Normalized Difference Vegetation Index)

!>     Derived variables
      real, intent(out) :: zref(nx,ny)       ! The reference level of atmospheric forcing


!#######################################################################
!     Prognostic variables

!     Soil / ground parameters / variables
      type(soil_params), intent(in) :: soil(:)


         call sib2_INITSOIL(
     &      nx,ny,nz,nvtyps, ! xlon,ylat, !2008.03.31
     &      Sst, Sstmax0,
     &      ele2d, sta_alt,
     &      Ds2d,Dg2d,Dgl2d,GWstorage,
     &      inbasin,nroot,ndeep,
     &      zref,
     &      soiltyp,vtypfrct,vegtyp,ndvi,lai,fpar,vegcover_mo,vcover,
     &      pg3,
     &      soil,mrf2,frc)
         return
      end subroutine initsoil


      SUBROUTINE vegtable(opc)
         use veginit, only: sib2_vegtable => VEGTABLE
         use lsmpar, only: mrf2
         implicit none

! ===== PARAMETERS
         type(optical), intent(out) :: opc(:)

         call sib2_VEGTABLE(opc,phy2)
         if (mrf2(2)%laimax <= 0) stop "LAIMAX <= 0"
         if (phy2(2)%anik <= 0) stop "ANIK1 <= 0"
         return
      end subroutine vegtable


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
! Begin snapshotting routines
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

      function get_snapshot_size()
         implicit none

         integer get_snapshot_size

! ===== CODE
         inquire(iolength=get_snapshot_size)
     &         Dgl2d,GWstorage,
     &         nroot,ndeep,
     &         lai,fpar,vcover,r3,pg3,
     &         vegcover_mo
         return
      end function get_snapshot_size

      ! Moiz: To override swe initial condition (for experiments)
      function get_snapshot_size_swe()
         implicit none

         integer get_snapshot_size_swe

! ===== CODE
         inquire(iolength=get_snapshot_size_swe)
     &         pg3%snowc, pg3%snowg
         return
      end function get_snapshot_size_swe

      subroutine save_snapshot(fileid,recid)
         implicit none

         include 'dims.inc' ! nx,ny,nvtyps

! ===== PARAMETERS
         integer, intent(in) :: fileid
         integer, intent(in) :: recid

         integer ii

! ===== CODE
         write(fileid, rec=recid)
     &         Dgl2d,GWstorage,
     &         nroot,ndeep,
     &         lai,fpar,vcover,r3,pg3,
     &         vegcover_mo

         return
      end subroutine save_snapshot

      subroutine save_snapshot_swe(fileid,recid)
         implicit none

         include 'dims.inc' ! nx,ny,nvtyps

! ===== PARAMETERS
         integer, intent(in) :: fileid
         integer, intent(in) :: recid

         integer ii

! ===== CODE
         write(fileid, rec=recid)
     &         pg3%snowc, pg3%snowg

         return
      end subroutine save_snapshot_swe


      subroutine load_snapshot(fileid,recid)
         implicit none

         include 'dims.inc' ! nx,ny,nvtyps

! ===== PARAMETERS
         integer, intent(in) :: fileid
         integer, intent(in) :: recid

! ===== LOCAL VARIABLES
         integer ii

! ===== CODE
         read(fileid, rec=recid)
     &         Dgl2d,GWstorage,
     &         nroot,ndeep,
     &         lai,fpar,vcover,r3,pg3,
     &         vegcover_mo

         return
      end subroutine load_snapshot

      subroutine load_snapshot_swe(fileid,recid)
         implicit none

         include 'dims.inc' ! nx,ny,nvtyps

! ===== PARAMETERS
         integer, intent(in) :: fileid
         integer, intent(in) :: recid

! ===== LOCAL VARIABLES
         integer ii

! ===== CODE
         read(fileid, rec=recid)
     &         pg3%snowc, pg3%snowg

         return
      end subroutine load_snapshot_swe

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! End snapshotting routines
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


!
!     ##################################################################
!     ##################################################################
!     ######                                                      ######
!     ######                SUBROUTINE FLUX2D                     ######
!     ######                                                      ######
!     ######                     Developed by                     ######
!     ######     River and Environmental Engineering Laboratory   ######
!     ######                University of Tokyo                   ######
!     ######                                                      ######
!     ##################################################################
!     ##################################################################
!

!> \details
!>     The subroutine is the driver of sib2-gbhm. It calculates
!>     surface momentum, heat, transevaporation and runoff. The runoff flux
!>     will be passed into kinematic hydraulic processes.

      SUBROUTINE flux2d (
     &      x,y, ! xlon,ylat,
     &      zref,slp2d,len2d,
     &      ele2d, sta_alt,
     &      Dg2d,Ds2d,Dr2d,Drw2d,
     &      inbasin,
     &      frc,
     &      prc,prl,rho,  ! rh
     &      soiltyp,vtypfrct,vegtyp,
     &      diag2,
     &      Sst, Sstmax0,
     &      waterflx, runoff_inter, runoff_g
     & )

!#######################################################################
!
!     02/21/2006 (Wang Lei)
!     Change the method to
!     derive aeodynamic and physiologic paramters once a day
!
!#######################################################################
         
         use metdata
         use forcing
         use vertical_profile
         use point_scale
         use soilinit, only: dailyveg
         use output_sptl

         implicit none
!
!#######################################################################


! ===== INCLUDES
         include 'globcst.inc'
         include 'sfcphycst.inc'
         include 'phycst.inc' ! cp,p0,rd
         include 'dims.inc' ! nx,ny,nvtyps
         include 'hydro.inc'
         include 'incl_params.f' ! MONTHLY_LAIFPAR,POINTSCALE_X,POINTSCALE_Y


! ===== PARAMETERS (BASIC)
      real, intent(in)     :: x  (nx)           !< The x-coord. of the computational grid (m).
      real, intent(in)     :: y  (ny)           !< The y-coord. of the computational grid (m).
! These are now defined in metdata (PRL)
!      real, intent(in)     :: xlon(nx,ny)       !< longtude at (x,y)
!      real, intent(in)     :: ylat(nx,ny)       !< latitude at (x,y)
      real, intent(in)     :: zref(nx,ny)       !< The reference level from ground surface (m).
      real, intent(in)     :: slp2d(nx,ny)      !< The slope of ground surface (ND)
      real, intent(in)     :: len2d(nx,ny)      !< hillslope length of each GRID (m)
      real, intent(in)     :: ele2d(nx,ny)
      real, intent(in)     :: sta_alt(nx,ny)
      real, intent(in)     :: Ds2d(nx,ny)       !< depth of UZ soil of each GRID (m)
      real, intent(in)     :: Dg2d(nx,ny)       !< depth of unconfined aquifer of each GRID (m)
      real, intent(in)     :: Dr2d(nx,ny)       !< river depth in each GRID (m)

      real, intent(inout)  :: Drw2d(nx,ny)      !< water depth in the river of each GRID (m)
!      real, intent(inout)  :: Dgl2d(nx,ny,0:nvtyps)       !< depth to groundwater level in each GRID (m)
!      real, intent(inout)  :: GWstorage(nx,ny,0:nvtyps)   !< depth to groundwater level in each GRID (m)

      integer, intent(in) :: inbasin(nx,ny)     !< flag to show the basin region

! ===== PARAMETERS (Atmospheric forcing data)
      type(meteorological_data), intent(inout) :: frc(:,:)
      real, intent(out) :: prc  (nx,ny)     !< convective precipitation rates (kg/(m**2*s))
      real, intent(out) :: prl  (nx,ny)     !< large-scale condensation rates (kg/(m**2*s))
!
!     secondary: rho
!
      real, intent(out) :: rho  (nx,ny)     !< Surface air density (kg/m**3)
!
!#######################################################################
!
!     Vegetation parameters
!
!#######################################################################
!
!     primary
!
      integer, intent(in) :: soiltyp(nx,ny)        !< Soil type at each point
!      type(soil_params), intent(in) :: soil(:)
      real, intent(in) :: vtypfrct(nx,ny,nvtyps)   !< Fraction of vegetation types
      integer, intent(in) :: vegtyp (nx,ny,nvtyps) !< Vegetation type at each point

!> \todo Move these elsewhere
!      type(diagnostic3d), intent(out) :: diag3(:,:,:)
      type(diagnostic_fluxes), intent(out) :: diag2(:,:)

!
!     Derived variables
!
!      integer, intent(inout) :: nroot  (nx,ny,nvtyps) !< number of layers in root zone
!      integer, intent(inout) :: ndeep  (nx,ny,nvtyps) !< number of layers in deep soil
!      real, intent(inout) :: lai    (nx,ny,nvtyps) !< leaf area index
!      real, intent(inout) :: fpar   (nx,ny,nvtyps) !< the Fraction of Phyotosynthetically Active Radiation absorbed by the green vegetation canopy
!a      real, intent(inout) :: roufns (nx,ny,nvtyps)  !< Surface roughness

!      type(remainder3dgrid), intent(inout) :: r3(nx,ny,nvtyps)

!
!#######################################################################
!
!     Prognostic variables
!
!#######################################################################
!
!      type(prognostic3dgrid), intent(inout) :: pg3(nx,ny,0:nvtyps)


! 'TEMPORARY' variable
! TODO: check if this is really temporary (PRL)
      real wetg_tmp(nx,ny)   !< for saving status


      real, intent(inout) :: Sst(nx,ny)
      real, intent(in) :: Sstmax0(nx,ny)

! ===== PARAMETERS (Output water flow)
      real, intent(out) :: runoff_inter(nx,ny)
      real, intent(out) :: runoff_g(nx,ny)
      real, intent(out) :: waterflx(nx,ny)


! ===== LOCAL VARIABLES

!#######################################################################
!     Snow coverage and extent


      integer snowc_daily(nx,ny)    ! daily fractional snow cover  (0~ 100)             !2/27/2007
      integer snowe_8daily(nx,ny)   ! 8-daily Maximum_Snow_Extent  (value: 200) !3/02/2007

      integer i, j,k,iv,marki,markj
      real    store             ! capac(2)  of previous step
      real    Dgltmp
      real    p0inv,eps
      integer isub, iflow , ig  ! 03/29/2006
      integer isubflow


! rasmy weB-s
      REAL CAPAC_IN   ! Moiz: Why is this defined? This is not used anywhere



      integer ireturn

      integer iii,jjj  ! just for testing (PRL)
      real tmax, tsum ! testing

      integer isoil, ia1, ia2, ib1, ib2
      character*10 varname

      integer outflag
!ysmod----------------------------------------------------------------72
      character yrstr*4, mostr*2, dystr*2, hrstr*2, rtstr*256
!ysmod----------------------------------------------------------------72

! URGENT: get rid of these (PRL)
!      include 'incl_ground_2d_params.f'
!      include 'incl_lsm_param_decl.f'
      type(ground_2d_params) :: i2d
!      type(o3d_vars) :: o3d
!      type(veg_cst) :: cst

!      type(prognostic3dgrid) pg3(nx,ny,0:nvtyps)
!      type(remainder3dgrid) r3i(nx,ny,0:nvtyps)

! TODO: These should be obtained from a memory pool
! Memory pool: input-output from point-scale
      type(ground_3d_vars) :: io3d
      type(ground_2d_vars) :: io2d
      type(vertical_profile_in) :: vp
      type(remainder) :: rm

! Memory pool: output from point-scale
      type(balance_info) :: bal


      integer ntotal_t

      real swdown !< incoming downward solar radiation (w/m2)
!      integer vsib2 !< id of vegetation type of current grid point
      real wsat
      real rhoair,rcp

!!!!!!!!!!!!!!! For OPENMP:
      integer tid, nthreads
      integer OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS


! ASSERTION:
      if (mrf2(2)%laimax <= 0) stop "FLUX2D1:LAIMAX=0"



! ===== CODE

!#######################################################################
!
!     Input hourly rainfall
!
!#######################################################################


         if (RAINFALL_INPUT_TYPE == 1) then
            stop "INTERNAL TEMPORAL DISAGGREGATION IS NOT SUPPORTED"
!            call read_grid_daily(nx,ny,year,jday,hour+1,rainh)
!
!            tmax = 0.
!            tsum = 0.
!            do i=1,nx
!               do j=1,ny
!                  tmax = max(rainh(i,j),tmax)
!                  tsum = tsum + rainh(i,j)
!               end do
!            end do
         endif

         if (RAINFALL_INPUT_TYPE == 2) then
               varname = 'Rain2'
            if (READ_SEPARATE_RAIN_SNOW == .false.) then
               rainptr => rainh
               call read_grid_old(varname, nx, ny, 0,
     &         FMT_PPT,'M', rainptr) !2008.04.11 ! Moiz: rainptr --> rainh
            else
               !Moiz: Added for READ_SEPARATE_RAIN_SNOW
               rainptr_liquid => rainh_liquid
               rainptr_solid => rainh_solid
               call read_grid_old(varname, nx, ny, 0,
     &         FMT_Rain,'M',rainptr_liquid)
               call read_grid_old(varname, nx, ny, 0,
     &         FMT_Snow,'M',rainptr_solid)
            endif
         end if

         if (RAINFALL_INPUT_TYPE == 3) then
! This is the case if precipitation is to be taken from JRA
         endif

         if (RAINFALL_INPUT_TYPE == 4) then
            varname = 'Rain4'
            if (READ_SEPARATE_RAIN_SNOW == .false.) then
               rainptr => rainh
               call read_grid_old(varname, nx, ny, 0,
     &         FMT_PPT,'H', rainptr) !2008.04.11 ! Moiz: rainptr --> rainh
            else
               !Moiz: Added for READ_SEPARATE_RAIN_SNOW
               rainptr_liquid => rainh_liquid
               rainptr_solid => rainh_solid
               call read_grid_old(varname, nx, ny, 0,
     &         FMT_Rain,'H',rainptr_liquid)
               call read_grid_old(varname, nx, ny, 0,
     &         FMT_Snow,'H',rainptr_solid)
            endif
         end if

         if (RAINFALL_INPUT_TYPE == 5) then
            varname = 'Rain5'
            if (READ_SEPARATE_RAIN_SNOW == .false.) then
               rainptr => rainh
               call read_grid_old(varname, nx, ny, 0,
     &         FMT_PPT,'H', rainptr) !2008.04.11 ! Moiz: rainptr --> rainh
            else
               !Moiz: Added for READ_SEPARATE_RAIN_SNOW
               rainptr_liquid => rainh_liquid
               rainptr_solid => rainh_solid
               call read_grid_old(varname, nx, ny, 0,
     &         FMT_Rain,'H',rainptr_liquid)
               call read_grid_old(varname, nx, ny, 0,
     &         FMT_Snow,'H',rainptr_solid)
            endif
         end if

!#######################################################################
!
!     Read surface meteorological data for each time step
!
!#######################################################################
!

         call read_meteo(
     &            nx,ny, ! xlon,ylat,
     &            year,month,jday,day,hour,minute,frc,
     &            inbasin, ele2d, sta_alt)

!
!#######################################################################
!
!     calculate air density
!
!#######################################################################
!
      DO i = 1, nx
         DO j = 1, ny
            if( inbasin(i,j) == 1) then
               rho (i,j) =
     &               frc(i,j)%psfc/
     &               (Rd*frc(i,j)%tair *(1+0.61*frc(i,j)%qv) )
            else
               rho(i,j) = 0
            end if
         END DO
      END DO
      p0inv=1./p0
      eps = 1.0e-6
!
!#######################################################################
!
!     Read 8-daily maximum snow extent.
!     Read 8-daily lai and fpar.
!     Lei Wang (3/09/3007)
!
!#######################################################################

!> \todo Should be based on month for monthly
      if(  (mod(jday,8)==1).and.hour==0  )then

!     call snowextent(nx, ny, year, jday, snowe_8daily)             !4/18/2007
      

      if (MONTHLY_LAIFPAR) then
         call LaiFpar_monthly(
     & nx, ny, year, month, lai(1,1,1), fpar(1,1,1))
      elseif (CONST_LAIFPAR) then
         call LaiFpar_const(
     & nx, ny, year, month, lai(1,1,1), fpar(1,1,1))
      else
         call LaiFpar_8day(nx, ny, year, jday, lai(1,1,1), fpar(1,1,1))
      end if

         do j = 1, ny
            do i = 1, nx
               if(  lai(i,j,1) < 0.1  )  lai(i,j,1) = 0.1
               if( fpar(i,j,1) < 0.01 ) fpar(i,j,1) = 0.01
            enddo
         enddo

      end if
!
!#######################################################################
!
!     Calculate vegetation coverage
!
!#######################################################################
!
#ifdef ALTERNATE_VEG
#else
      if(hour==0)  then
         call dailyveg(nx,ny,nvtyps,vegtyp,vegcover_mo,lai,vcover)
      end if
#endif

!#######################################################################
!     Set constants used in sib2

      CALL const2               ! set constants

!#######################################################################
!     Call sib2 to predict unknowns and output turbulent fluxes

!$OMP PARALLEL DEFAULT(none)
!$OMP&         PRIVATE(i2d,io3d,vp,bal,iv,iflow,ig,store,capac_in,
!$OMP&         rm,k,wsat,
!$OMP&         ntotal_t,p0inv,eps,marki,markj,Dgltmp,alfa,beta,swdown,
!$OMP&         rhoair,rcp,
!!!!$OMP&          i,j,    ! For loop type A
!$OMP&         i,j,isub,   ! For loop type A prime
!$OMP&         tid,nthreads)
!$OMP&         SHARED(pg3,inisub,finsub,num_flow,ngrid,
!$OMP&         grid_col,grid_row,inbasin,zref,
!$OMP&         bats_state, snow_state,
!$OMP&         ele2d, sta_alt,
!!!$OMP&       uv,tair,psfc,radsw,cld,radlw,qv,
!$OMP&         frc,
!$OMP&         rho,rainh,rainh_liquid, rainh_solid,
!$OMP&         prl,prc,year,month,day,hour,minute,second,xlon,ylat,
!$OMP&         sstmax0,Ds2d,Dg2d,Dr2d,len2d,slp2d,
!$OMP&         curtim,tstart,inicon,sst,dtlsm,
!$OMP&         runoff_inter,runoff_g,
!$OMP&         diag2,diag3,
!!!$OMP&         ptsflx,qvsflx,rnflx,ect_flx,eci_flx,egs_flx,egi_flx,
!!!$OMP&         wsfc_hour,wrt_hour,wdp_hour,
!$OMP&         waterflx,wetg_tmp,lai,vcover,vegtyp,soil,soiltyp,
!$OMP&         Dgl2d,GWstorage,Drw2d,r3,
!$OMP&         fpar,nroot,ndeep,vtypfrct,in_analysis
!$OMP&         ,prli,prlj
!$OMP&         ,isubflow_cnt,omp_subbasin,omp_flowinterval
!$OMP&         ,nx,ny   ! For loop type B
!$OMP&         ,pointscale_x,pointscale_y,
!$OMP&         rain_correct,snow_correct1,snow_correct2,thsnow,
!$OMP&         READ_SEPARATE_RAIN_SNOW)

!,csst,r3i,r3point,io3d,io2d,rm,vegt,vl,vp,bal)

      TID = OMP_GET_THREAD_NUM()
      IF (TID .EQ. 0) THEN
        NTHREADS = OMP_GET_NUM_THREADS()
!        PRINT *, 'Number of threads =', NTHREADS
      END IF

!!!! Parallelized loop type A
!!!!$OMP DO
!!!! collapse(3) not allowed because the loops aren't independent
!!!      do isub = inisub, finsub  !24, 40  for Agatsuma
!!!         do iflow = 1,num_flow(isub)
!!!            do ig=1,ngrid(isub,iflow)
!!!               i=grid_col(isub,iflow,ig)
!!!               j=grid_row(isub,iflow,ig)

!!!! Parallelized loop type A prime
!$OMP DO

      do isubflow = 1, isubflow_cnt
         isub = omp_subbasin(isubflow)
         iflow = omp_flowinterval(isubflow)
         do ig=1,ngrid(isub,iflow)
            i=grid_col(isub,iflow,ig)
            j=grid_row(isub,iflow,ig)


!!!! Parallelized loop type B
!!!!$OMP DO collapse(2)
!!!         do i=1,nx
!!!            do j=1,ny

!                  print *, "LOOP:",i,j
               IF(inbasin(i,j) == -9999)goto 200
               if (.not. in_analysis(i,j)) goto 200
!              if ((i==53.and.j==150)) then ! Moiz: pointsimulation
!               if ((i==26.and.j==208)) then ! Moiz: pointsimulation
!                if (i==2.and.j==204) then

               if ((POINTSCALE_X > 0) .and. (POINTSCALE_Y > 0) .and.
     &         ((i /= POINTSCALE_X) .or. (j /= POINTSCALE_Y)) ) goto 200

      prli = i
      prlj = j
      i2d%prli = i
      i2d%prlj = j

!#######################################################################
!
!     Begin preparing sib2 input ......
!
!#######################################################################

               rm%zwind=  zref(i,j) !02/08/2006
               rm%zwind= min(rm%zwind,z1limit)

               rm%um   = max( frc(i,j)%uv, vsfcmin ) ! total speed(m/s)
               rm%tmx   = frc(i,j)%tair ! temperature(K)

               rm%psur = frc(i,j)%psfc
               rhoair = rho(i,j)
               rcp    = rhoair * cp

!
!#######################################################################
!
!     For vegetation types, the following forcing data are required
!
!#######################################################################
!
               rm%psur = rm%psur /100 ! total pressure (hPa)

!  equation:   e = p*qv/(0.622+0.378*qv)   wanglei(12/23/2006)
               rm%em   = rm%psur*frc(i,j)%qv/
     &                     (epsfac+(1-epsfac)*frc(i,j)%qv) ! e (hPa)        !qv:wanglei 11/25/2006

              ! Added Moiz
              if (READ_SEPARATE_RAIN_SNOW == .false.) then
              !------------------------------------------------------------------------------
              ! Radar Precipitation Correction: Adjust JMA radar precipitation in some months based on ratio
              if (isnan(rainh(i,j))) then ! Sometimes there are missing values in the JMA radar product
                  rainh(i,j) = 0
                  write(*,*) 'NaN value in Radar Precipitation'
              endif
              !------------------------------------------------------------------------------

              !------------------------------------------------------------------------------
              ! Snow Precipitation Correction:
              !write(*,*) rain_correct, snow_correct1, snow_correct2
              rainh(i,j) = max(rainh(i,j),0.0)
              if (rm%tmx > (tf+thsnow)) then
                 rainh(i,j) =  rainh(i,j)*rain_correct    ! Radar-Amedas Rainfall Correction
              else
                 if (snow_correct1 > 0 .and. ele2d(i,j) <= 2000) then       ! Radar-Amedas Snowfall Correction [Shrestha et al., 2014]
                        rainh(i,j) = rainh(i,j)*snow_correct1
                        
                 elseif (snow_correct2 > 0 .and. ele2d(i,j) > 2000)then
                        rainh(i,j) = rainh(i,j)*snow_correct2
                 elseif (snow_correct1 < 0) then   ! Station Snowfall Correction [Shrestha et al., 2014]
            rainh(i,j) = (rainh(i,j))*(1.0+((ele2d(i,j) - sta_alt(i,j))
     &            *(-1)*snow_correct1))
                 endif
              endif
              
              else
               rainh(i,j) = rainh_liquid(i,j) + rainh_solid(i,j)
               endif
                
              
                 !------------------------------------------------------------------------------



!*****************************************************************
!     Climate change
!     Overall, land precipitation for the globe has increased by ~2% since 1900
!*****************************************************************
!     rainh(i,j) = rainh(i,j) * 1.02    !2% precipitation change for one century
!*****************************************************************

!     2/27/2007    Wang
!     change from
!     ppl =  0.3*rainh(i,j)     ! (mm)
!     ppc  = 0.7*rainh(i,j)     ! (mm)
!     to :                                 2007.12.22   SCE calibrate
!               rm%ppl =  rainh(i,j) * dtlsm / 3600. ! (mm) ! (scaled, PRL) ! Moiz: Commented out
!               rm%ppc  = 0.0       ! (mm)                                  ! Moiz: Commented out
                rm%ppl = 0.0
                rm%ppc = rainh(i,j) * dtlsm/3600

                ! Moiz: Added
                if (READ_SEPARATE_RAIN_SNOW == .true.) then
                rm%snowfall = rainh_solid(i,j) * dtlsm/3600
                rm%rainfall = rainh_liquid(i,j) * dtlsm/3600
                endif

!     end changing

               prl(i,j) = rm%ppl   ! (mm)    12/25/2006
               prc(i,j) = rm%ppc   ! (mm)    12/25/2006

               rm%po2m = 20900. * rm%psur / (p0/100) ! O2 pres. (Pa)

!*****************************************************************
!     CO2  change  study   280ppm (before 1800)  &  380ppm (now)
!     one standard air pressure (1.0e5 pa) :  28 pa  & 38 pa
!*****************************************************************
               rm%pco2m= 37. * rm%psur / (p0/100) ! CO2 pres.(Pa)
!     pco2m= 38. * psur / (p0/100)          ! CO2 pres.(Pa)
!     pco2m= 28. * psur / (p0/100)          ! CO2 pres.(Pa)
!     pco2m= 56. * psur / (p0/100)          ! CO2 pres.(Pa)
!     Climate change
!#######################################################################
!
!     partition solar radiation into diffuse and direct UV and IR
!     components
!     sunang = cosZ
!
               CALL zenith_loc(year, month, day, hour, minute, second,
     &              xlon(i,j),ylat(i,j),rm%sunang, alfa, beta) !02/22/2006: , alfa, beta

               swdown = amax1(frc(i,j)%radsw,0.1) ! Paritition solar radiation
               CALL radswpart(swdown,rm%sunang,frc(i,j)%cld,rm%radn)
               rm%radn(3,2) = frc(i,j)%radlw

               ! Moiz: Testing
               !swdown_pt = swdown
               !lwdown_pt = rm%radn(3,2)
               !cld_pt = frc(i,j)%cld

!
!#######################################################################
!
!     Soil parameters
!
!#######################################################################
!
               i2d%satcap_gnd = Sstmax0(i,j) !2007.12.28

! Fixed parameters (not time changing)
               vp%Ds     = Ds2d (i,j)
               vp%Dg     = Dg2d (i,j)
               vp%Dr     = Dr2d (i,j)
               vp%length = len2d(i,j)
               vp%slope  = slp2d(i,j)


!#######################################################################
!
!     we can change the hydraulic conductivity (ksat1, ksat2, kgs) for
!     different land use type   (e.g. ksat1->0 for urban area)
!
!     For urbanized area and barren (7)
!
!     Wang Lei
!     12/11/2006
!#######################################################################

! URGENT: Add the code for urban infiltration back in (PRL)
!               IF (vegtyp(i,j,1) == 7) THEN
!                  ksat1  =  0.1 * soil(soiltyp(i,j))%ksat1 ! 2008.04.09
!               end if

!#######################################################################
!
!     Initialize surface water storage for each grid
!
!     2007.10.30              2010.8.27
!
!********************************************************
!     change from

!     if(curtim == tstart .and. inicon==1) then
!     Sst(i,j) = amin1( wetg (i,j,0) , Sstmax0(i,j) )
!     end if
!     store      = Sst(i,j)        !(m)

!     to:
               if(curtim == tstart .and. inicon==0) then
                  store    =  amin1( pg3(i,j,0)%wetg, Sstmax0(i,j) )
               else
                  store    =  Sst(i,j)
               end if

!     end changing
!********************************************************

               diag2(i,j)%ptsflx = 0. ! surface flux of heat (K*kg/(m**2*s))
               diag2(i,j)%qvsflx = 0. ! surface flux of moisture (kg/(m**2*s))
               diag2(i,j)%Rnflx  = 0.

               runoff_inter(i,j) = 0.
               runoff_g  (i,j)   = 0.

               diag2(i,j)%Ect_flx = 0.
               diag2(i,j)%Eci_flx = 0.
               diag2(i,j)%Egs_flx = 0.
               diag2(i,j)%Egi_flx = 0.

               diag2(i,j)%tsfc = 0.

!     roff_thru (i,j)  = 0.        ! add: 02/15/2006   del: 07/02/2006

!RASMY
! WEB-S CHECK THIS peter swithched off
               Sst (i,j)   = 0. ! 02/08/2006: important
               waterflx  (i,j)   = 0.
               wetg_tmp(i,j)   = 0. ! 2007.10.30
!
!#######################################################################
!
!     Calculate prognostic variables and fluxes on each land use type
!
!#######################################################################
!
               Dgltmp = 0
               DO iv = 1, nvtyps

! ASSERTION:
	      if (pg3(i,j,iv)%tcanp <= 0) then
        	 print *,i,j,iv
	         stop "Tcanp <= 0"
	      end if

!RASMY
! WEB-S - WHY WE NEED TO CALCULATE CAPAC   ! Moiz: Is this needed or not?

!                 capac_in = pg3(i,j,iv)%wetg * store / !wetg(i,j,iv)                                                      
!     &                 max(1.0E-6,pg3(i,j,0)%wetg) !note 2010.8.27                                                         
!                 capac_in = amin1(capac_in,amax1(0.0,pg3(i,j,iv)%wetg))                                                   
!                 Sst(i,j) = Sst(i,j) - capac_in                                                                           
!                 Sst(i,j) = amax1(Sst(i,j), 0.0)  

      

!
!#######################################################################
!
!     For water body (10)
!
!#######################################################################
               IF (vegtyp(i,j,iv) == 10) THEN

                  call process_reservoir(rhoair,rcp,
     &                  rm,soil(soiltyp(i,j)),
     &                  io3d,bal,pg3(i,j,1),waterflx(i,j),
     &                  Dgl2d(i,j,iv),GWstorage(i,j,iv),
     &                  diag3(i,j,iv),
     &                  frc(i,j)%psfc,frc(i,j)%qv,vp%Dg
     &)

               else

                  io3d%bats_state = bats_state(i,j)   ! Moiz: Added
                  io3d%snow_state = snow_state(i,j)


                  call process_land(rhoair,rcp,
     &                  rm,soil(soiltyp(i,j)),
     &                  io3d,bal,pg3(i,j,iv),waterflx(i,j),   ! Moiz: Here iv =1 so, pg3(i,j,1) is used
     &                  Drw2d(i,j),Dgl2d(i,j,iv),GWstorage(i,j,iv),
     &                  diag3(i,j,iv),
     &                  frc(i,j)%psfc,frc(i,j)%qv,vp%Dg,
     &                  vegtyp(i,j,iv),
     &                  i2d,r3(i,j,iv),vp,
     &                  Ds2d(i,j),
     &                  ylat(i,j), pg3(i,j,0)%wetg,
     &                  store,
     &                  lai(i,j,iv),fpar(i,j,iv),
     &                  nroot(i,j,iv),ndeep(i,j,iv),vcover(i,j,iv)
     &)
                  bats_state(i,j) = io3d%bats_state    ! Moiz : Added
                  snow_state(i,j) = io3d%snow_state


               end if

                  diag2(i,j)%ptsflx = diag2(i,j)%ptsflx
     &                           + bal%hflux*vtypfrct(i,j,iv)
                  diag2(i,j)%qvsflx = diag2(i,j)%qvsflx
     &                           + bal%etmass*vtypfrct(i,j,iv)
                  diag2(i,j)%Rnflx = diag2(i,j)%Rnflx
     &             + (bal%radt(1)+bal%radt(2))*vtypfrct(i,j,iv)

                  diag2(i,j)%Ect_flx = diag2(i,j)%Ect_flx
     &                              + bal%ectmass*vtypfrct(i,j,iv)
                  diag2(i,j)%Eci_flx = diag2(i,j)%Eci_flx
     &                              + bal%ecimass*vtypfrct(i,j,iv)
                  diag2(i,j)%Egs_flx = diag2(i,j)%Egs_flx
     &                              + bal%egsmass*vtypfrct(i,j,iv)
                  diag2(i,j)%Egi_flx = diag2(i,j)%Egi_flx
     &                              + bal%egimass*vtypfrct(i,j,iv)

                  diag2(i,j)%tsfc = diag2(i,j)%tsfc
     &                              + io3d%tg*vtypfrct(i,j,iv)

                  runoff_inter(i,j)   = runoff_inter(i,j)
     &                 + bal%roffinter*vtypfrct (i,j,iv)
                  runoff_g(i,j)   = runoff_g(i,j)
     &                 + bal%roffg*vtypfrct (i,j,iv)
!RASMY
! WEB-S why we need to add extra_roff
                   Sst(i,j)  = Sst(i,j) + vtypfrct (i,j,iv)* 
     &                          (io3d%capac(2)+io3d%extra_sfc_roff)  

!                  Sst(i,j)  = Sst(i,j) + io3d%capac(2)*vtypfrct (i,j,iv)
                  Sst(i,j)  = amax1( 0., Sst(i,j) )

                  wetg_tmp(i,j) = wetg_tmp(i,j)
     &                        + pg3(i,j,iv)%wetg * vtypfrct(i,j,iv)
                  wetg_tmp(i,j) =  amax1( 0., wetg_tmp(i,j) ) !2007.10.30

               END DO
!
!#######################################################################
!
!     Flux unit convert from SiB2 to IGBHM
!     IGBHM                           SiB2
!     ptsflx= rhoair*bar( (w't') )      hflux=rhoair*cp*bar((w't'))
!     (K*kg/(m**2*s))                   (w/m2)
!
!     qvsflx=rhoair*bar( (w'qv') )      etmass=rhoair*bar((w'qv'))*dtlsm
!     (kg/(m**2*s))                      (kg/m**2)
!
!     runoff (m/s)                      roff (m)
!
!#######################################################################
!     delete   4/18/2007
!     ptsflx(i,j) = ptsflx(i,j) / cpair
!     qvsflx(i,j) = qvsflx(i,j)


               marki = i
               markj = j
                        
!            end if ! Moiz: Point simulation
 200           CONTINUE

            

! Parallelized loop type B
         end do   ! j
      end do   ! i

!!!! Parallelized loop type A
!!!            END DO              ! igrid
!!!         END DO                 ! iflow
!!!      END DO                    ! isub

!$OMP END DO

!      write(6,'(4i6)')  !, 3f12.2, 3f16.5)')
!     &      year,month,day,hour
!     &      io3d%tc,io3d%tg,rm%tm, ! hg/3600.0,hc/3600.0,
!     &      io3d%vlc(1), io3d%capac(2), io3d%snoww(2)

!$OMP END PARALLEL


      outflag = 0
!     Save soil and ground water status at end of the simulation
!     2007.10.29
      if ((inicon == 0) .and. (outflag == 1)) then
         if( curtim + dt_couple == tstop  ) then !03/10/2006
            do isub = inisub, finsub  !24, 40  for Agatsuma
               call save_soil_data(
     &               isub,vtypfrct,nroot,ndeep,pg3,Dgl2d,GWstorage)
            end do
         end if
      end if


!#######################################################################
!     Update prognostic variables on grid scale

      DO j=1,ny
         DO i=1,nx

            pg3(i,j,0)%tsfc    = 0.0
            pg3(i,j,0)%tsoil   = 0.0
            pg3(i,j,0)%tcanp   = 0.0
            pg3(i,j,0)%tland   = 0.0

            pg3(i,j,0)%vlcsfc  = 0.0
            pg3(i,j,0)%vlcrt   = 0.0
            pg3(i,j,0)%vlcdp   = 0.0

            pg3(i,j,0)%wetcanp = 0.0
            pg3(i,j,0)%wetg    = 0.0
            pg3(i,j,0)%snowc   = 0.0
            pg3(i,j,0)%snowg   = 0.0
            pg3(i,j,0)%rst2d   = 0.0

!     soil_hour (i,j,0) = 0.0    !8/14/2007
            diag2(i,j)%wsfc = 0.0
            diag2(i,j)%wrt = 0.0
            diag2(i,j)%wdp = 0.0
            Dgl2d (i,j,0) = 0.0
            diag2(i,j)%snowdens = 0.0 ! Moiz: Added
            diag2(i,j)%snowdep = 0.0 ! Moiz: Added


            DO iv = 1, nvtyps
               
               pg3(i,j,0)%tsfc = pg3(i,j,0)%tsfc +
     &                       pg3(i,j,iv)%tsfc*vtypfrct(i,j,iv)
               pg3(i,j,0)%tsoil = pg3(i,j,0)%tsoil +
     &                        pg3(i,j,iv)%tsoil*vtypfrct(i,j,iv)
               pg3(i,j,0)%tcanp = pg3(i,j,0)%tcanp +
     &                        pg3(i,j,iv)%tcanp*vtypfrct(i,j,iv)
               pg3(i,j,0)%tland = pg3(i,j,0)%tland +
     &                        pg3(i,j,iv)%tland*vtypfrct(i,j,iv)
               pg3(i,j,0)%vlcsfc = pg3(i,j,0)%vlcsfc +
     &                         pg3(i,j,iv)%vlcsfc*vtypfrct(i,j,iv)
               pg3(i,j,0)%vlcrt = pg3(i,j,0)%vlcrt +
     &                        pg3(i,j,iv)%vlcrt*vtypfrct(i,j,iv)
               pg3(i,j,0)%vlcdp = pg3(i,j,0)%vlcdp +
     &                        pg3(i,j,iv)%vlcdp*vtypfrct(i,j,iv)
               pg3(i,j,0)%wetcanp = pg3(i,j,0)%wetcanp +
     &                          pg3(i,j,iv)%wetcanp*vtypfrct(i,j,iv)
               pg3(i,j,0)%wetg = pg3(i,j,0)%wetg +
     &                       pg3(i,j,iv)%wetg*vtypfrct(i,j,iv)
               pg3(i,j,0)%snowc = pg3(i,j,0)%snowc +
     &                        pg3(i,j,iv)%snowc*vtypfrct(i,j,iv)
               pg3(i,j,0)%snowg = pg3(i,j,0)%snowg +
     &                        pg3(i,j,iv)%snowg*vtypfrct(i,j,iv)

! ASSERTION:
      if (     (pg3(i,j,0)%snowg < 0) .and.
     &         (pg3(i,j,0)%snowg .ne. -9999)) then ! (PRL)
         print *, '1228:snowg negative: ',pg3(i,j,0)%snowg,i,j,iv
         stop
      end if
      if (pg3(i,j,0)%snowg > 10) then ! (PRL)
         print *, '1232:snowg over: ',pg3(i,j,0)%snowg,':',i,':',j
         stop
      end if

               pg3(i,j,0)%rst2d = pg3(i,j,0)%rst2d +
     &                           pg3(i,j,iv)%rst2d *vtypfrct (i,j,iv)

               diag2(i,j)%wsfc = diag2(i,j)%wsfc +
     &                            diag3(i,j,iv)%wsfc * vtypfrct(i,j,iv)
               diag2(i,j)%wrt = diag2(i,j)%wrt +
     &                           diag3(i,j,iv)%wrt * vtypfrct(i,j,iv)
               diag2(i,j)%wdp = diag2(i,j)%wdp +
     &                           diag3(i,j,iv)%wdp * vtypfrct(i,j,iv)

               Dgl2d(i,j,0) = Dgl2d(i,j,0) +
     &                        Dgl2d(i,j,iv)*vtypfrct(i,j,iv)

            END DO

            diag2(i,j)%tsoil = pg3(i,j,0)%tsoil
            diag2(i,j)%snowc = pg3(i,j,0)%snowc
            diag2(i,j)%snowg = pg3(i,j,0)%snowg
            diag2(i,j)%wetcanp = pg3(i,j,0)%wetcanp
            diag2(i,j)%wetg = pg3(i,j,0)%wetg
            diag2(i,j)%roff_i = runoff_inter(i,j)
            diag2(i,j)%roff_g = runoff_g(i,j)
            diag2(i,j)%dgl = Dgl2d(i,j,0)
            diag2(i,j)%drw = Drw2d(i,j)
            diag2(i,j)%snowdens = snow_state(i,j)%snowdens  ! Moiz: Added
            diag2(i,j)%snowdep = snow_state(i,j)%snowdep    ! Moiz: Added
            



         END DO
      END DO

!      call write_tprlavg_grid_data(pg3,diag3)

      if (POINTSCALE_X == -1 .and. POINTSCALE_Y == -1) then
! This line handles all the 'fast' spatial variables
         call write_fast_sptl(1)

         call write_slow_sptl(wsfc_hdl,diag2%wsfc)
         call write_slow_sptl(wrt_hdl,diag2%wrt)
         call write_slow_sptl(wdp_hdl,diag2%wdp)

         call write_slow_sptl(dgl_hdl,diag2%dgl)
         call write_slow_sptl(swe_hdl,diag2%snowg)

         call write_slow_sptl(eci_hdl,diag2%Eci_flx)
         call write_slow_sptl(ect_hdl,diag2%Ect_flx)
         call write_slow_sptl(egi_hdl,diag2%Egi_flx)
         call write_slow_sptl(egs_hdl,diag2%Egs_flx)
         
         ! Moiz: Snow Vars
         call write_slow_sptl(snowdep_hdl,snow_state%snowdep)
         call write_slow_sptl(snowdens_hdl,snow_state%snowdens)

         ! Moiz: Precipitation
         call write_slow_sptl(rainh_hdl,rainh)

      end if

      RETURN
      END SUBROUTINE flux2d




      subroutine process_reservoir(rhoair,rcp,
     &                  rm,soil,io3d,bal,pg3,waterflx,Dgl2d,GWstorage,
     &                  diag3,psfc,qv,Dg)
         use point_scale, only: f_qvsat
         use mod_reservoir

         implicit none

! ===== INCLUDES
         include 'globcst.inc'

! ===== PARAMETERS
         real, intent(in) :: rhoair
         real, intent(in) :: rcp
         type(remainder), intent(in) :: rm
         type(soil_params), intent(in) :: soil

         type(ground_3d_vars), intent(inout) :: io3d
         type(balance_info), intent(out) :: bal
         type(prognostic3dgrid), intent(inout) :: pg3
         real, intent(out) :: waterflx
         real, intent(out) :: Dgl2d
         real, intent(out) :: GWstorage
         type(diagnostic3d), intent(out) :: diag3
         real, intent(in) :: psfc
         real, intent(in) :: qv
         real, intent(in) :: Dg


! ===== LOCAL VARIABLES
         real c_u, c_pt,c_uneu,c_ptneu
         real qvsfc


! ===== CODE
         io3d%tg   = pg3%tsfc
         qvsfc = f_qvsat( psfc, io3d%tg) !Saturation water vapor specific humidity (kg/kg)

         CALL cucptwtr1D(rm%zwind,rm%um,io3d%tg,rm%tmx,
     &               c_uneu,c_ptneu,c_u,c_pt)

         bal%hflux = rhoair * cpair * c_u * c_pt* rm%um*(io3d%tg-rm%tmx)

         bal%etmass = rhoair * c_u * c_pt *
     &               rm%um*(qvsfc-qv) * dtlsm !(mm): refer to notes
         bal%ectmass = 0.
         bal%egsmass = 0.
         bal%ecimass = 0.

!     qv:wanglei 11/25/2006
         bal%roffg      = 0.0
         bal%roffinter  = 0.0
         bal%roffsub    = 0.0

         waterflx = (rm%ppl + rm%ppc-bal%etmass)/1000.0
!     In the new manuscript, we calculate water surface temperature in current timestep
!     as the weighted mean of water surface temperature in last timestep (tsfc),
!     and air temperature in current timestep (tm )
         pg3%tsfc = (pg3%tsfc*23 + rm%tmx) /24.

         Dgl2d     = 0.0 !2007.11.05
         GWstorage = Dg*soil%GWcs !2007.11.05

         diag3%wsfc = 1.0
         diag3%wrt  = 1.0
         diag3%wdp  = 1.0

!> \todo The following code is completely arbitrary.  In the original WEBDHM
!> the value for overland flow from reservoirs would be the value from the
!> previous 'land' gridpoint.  Now just assume it is zero.
         io3d%capac(2) = 0

! RASMY
! WEB-S WHY we need this set up
         io3d%extra_sfc_roff = 0. ! Moiz: Added    

         return
      end subroutine process_reservoir



      
      subroutine process_land(rhoair,rcp,
     &                  rm,soil,io3d,bal,pg3,waterflx,
     &                  Drw2d,Dgl2d,GWstorage,
     &                  diag3,
     &                  psfc,qv,Dg,
     &                  vsib2,
     &               i2d,r3,vp,
     &               Ds2d,
     &               ylat, wetg,
     &               store,
     &               lai,fpar,
     &               nroot,ndeep,vcover_in
     &)
         use point_scale, only: sib2
         implicit none

! ===== INCLUDES
         include 'globcst.inc'
         include 'hydro.inc' ! (inicon)
         include 'incl_params.f' ! pointscale_x, pointscale_y


! ===== PARAMETERS
         real, intent(in) :: rhoair
         real, intent(in) :: rcp

         type(remainder), intent(inout) :: rm
         type(soil_params), intent(in) :: soil

         type(ground_3d_vars), intent(inout) :: io3d
         type(balance_info), intent(out) :: bal
         type(prognostic3dgrid), intent(inout) :: pg3
         real,    intent(out) :: waterflx

         real,    intent(inout) :: Drw2d
         real,    intent(inout) :: Dgl2d   !< depth to groundwater level in each GRID (m)
         real,    intent(out) :: GWstorage

         type(diagnostic3d), intent(out) :: diag3
         real,    intent(in) :: psfc
         real,    intent(in) :: qv
         real,    intent(in) :: Dg

         integer, intent(in) :: vsib2

         type(ground_2d_params), intent(inout) :: i2d
         type(remainder3dgrid), intent(inout) :: r3
         type(vertical_profile_in), intent(inout) :: vp

         real,    intent(in) :: Ds2d      !< depth of UZ soil of each GRID (m)
         real,    intent(in) :: ylat      !< latitude at current point
         real,    intent(in) :: wetg
!> \todo This is ridiculous (PRL)
         real    store             ! capac(2)  of previous step
         real,    intent(in) :: lai
         real,    intent(in) :: fpar
         integer, intent(in) :: nroot
         integer, intent(in) :: ndeep
         real,    intent(in) :: vcover_in

!     &               i2d,r3(i,j,iv),vp,
!     &               Ds2d(i,j),Dgl2d(i,j,iv),
!     &               ylat(i,j), pg3(i,j,0)%wetg,
!     &               store,
!     &               lai(i,j,iv),fpar(i,j,iv),
!     &               nroot(i,j,iv),ndeep(i,j,iv)



! ===== LOCAL VARIABLES
         real c_u, c_pt,c_uneu,c_ptneu
         real qvsfc
         type(remainder3dpoint) :: r3point
         type(vegetation) :: veg
         type(vegetation_link) :: vl
         type(vertical_profile_inout) :: vio
         integer k
         real wsat



! ===== ASSERTIONS
      if (isnan(pg3%tcanp)) stop "Tcanp == NaN"
      if (phy2(2)%anik <= 0) stop "ANIK2 <= 0"
      if (mrf2(vsib2)%laimax <= 0) stop "FLUX2D:LAIMAX=0"


! ===== CODE

         vio%Drw   = Drw2d

         call basin_to_pointscale_land(rhoair,
     &               i2d, io3d,veg,vl,rm,
     &               r3,
     &               soil,r3point,
     &               mrf2(vsib2),phy2(vsib2),opc(vsib2),vp,vio,
!     &               curtim,
!     &year,month,day,hour,minute,
     &               (curtim == tstart .and. inicon == 0),
     &               Ds2d,Dgl2d,
     &               rm%zwind,ylat, wetg,
     &               store, vp%Dg, r3%green2d,
     &               lai,fpar,GWstorage,
     &               pg3,
     &               nroot,ndeep,vcover_in
     &)

!#######################################################################
!
!     Call SiB2 to solve prognostic variables
!
!#######################################################################
!

! ASSERTION
         if (veg%z2 < r3point%ha) then
            print *, "ZXO:",veg%z2,r3point%ha
            stop "HEIGHT mismatcha"
         end if

! ASSERTION
      if (isnan(io3d%tc)) stop "Tc is NaN"

! ASSERTION QQQQQ
!      if ((vio%Dgl < vp%ds) .and. (vio%Dgl < 3.99)) then
!         print *, "RIDIC:process_land"
!         stop "RIDIC:process_land"
!      end if

         CALL SiB2(0, vsib2, soil,veg,
     &      i2d, r3point,io3d,vl,bal,rm,vp,vio,rhoair,rcp) ! SiB2 main module  !02/01/2006



         if ((i2d%prli == POINTSCALE_X) .and.
     &               (i2d%prlj == POINTSCALE_Y)) then
      write(123,'(7i6,52f15.6)')year,month,day,hour,vio%diag_isat,
     & vp%nrt,vp%ndp,
     & io3d%capac(2)*1000.,io3d%snoww(2)*1000.,io3d%tc,io3d%tg,io3d%td,
     & io3d%vlc(1),io3d%vlc(2),io3d%vlc(3),
     & vio%dgl, bal%hflux,bal%etmass,
     & bal%ectmass,bal%ecimass,bal%egsmass,bal%egimass,
     & bal%roffinter*1000.,bal%roffg*1000.,vio%diag_rech,
     & vio%diag_qsub_a*1000,
     & vio%gwst,
     & vio%lsm(1:vp%nrt+vp%ndp+1),
     & rm%ppl,rm%um,rm%tmx,rm%psur,
     & rm%radn(1,1),rm%radn(1,2),
     & rm%radn(2,1),rm%radn(2,2),
     & rm%radn(3,1),rm%radn(3,2),vio%Drw  ! swdown missing
         end if




!     if ipbl = 0, call dtcdtg.
!     ipbl = 0 , when we calculate the tc, tg;
!     ipbl = 1 , when we calculate the Atmospheric Boundary Layer conditions.
!
!#######################################################################
!
!     save SiB2 variables to arrays
!
!#######################################################################
!
!     01/24/2006
         Dgl2d =     vio%Dgl
         GWstorage = vio%GWst
         Drw2d =     vio%Drw

! ASSERTION:
         if ((vp%nrt+vp%ndp+1) > 30)
     &         print *, "NTOT too high: ", vp%nrt, vp%ndp

         DO k = 1, vp%nrt+vp%ndp+1
            pg3%vlcprf(k) = vio%lsm(k)
         END DO

         pg3%tcanp = io3d%tc
         pg3%tsfc =  io3d%tg
         pg3%tsoil = io3d%td
         pg3%tland = (r3point%vcover*io3d%tc**4.+
     &               (1-r3point%vcover)*io3d%tg**4.)**0.25     ! 2008.7.29


         pg3%vlcsfc = io3d%vlc(1)
         pg3%vlcrt  = io3d%vlc(2)
         pg3%vlcdp  = io3d%vlc(3)

!     Lei Wang   7/30/2007
!     soil_hour (i,j,iv)     =
!     & ( vlc(1)*thick(1)+vlc(2)*thick(2)+vlc(3)*thick(3) )/sodep/wsat  !8/7/2007

         wsat = soil%wsat

! Volumetric soil moisture
!         diag3%wsfc = pg3%vlcsfc
!         diag3%wrt = pg3%vlcrt
!         diag3%wdp = pg3%vlcdp

         diag3%wsfc = pg3%vlcsfc /wsat
         diag3%wrt =  pg3%vlcrt  /wsat
         diag3%wdp =  pg3%vlcdp  /wsat


         pg3%wetcanp = io3d%capac(1)
         pg3%wetg =    io3d%capac(2)
         pg3%snowc =   io3d%snoww(1)
         pg3%snowg =   io3d%snoww(2)
         pg3%rst2d =   io3d%rstcanp

         waterflx = 0.
      end subroutine process_land



      subroutine save_soil_data(
     &      isub,vtypfrct,nroot,ndeep,pg3,Dgl2d,GWstorage)
         implicit none


! ===== INCLUDES
         include 'globcst.inc'   ! simulation_dir
         include 'dims.inc'      ! nx,ny, etc.
         include 'hydro.inc'     ! sub_catch


! ===== PARAMETERS

         integer, intent(in) :: isub
         real, intent(in) :: vtypfrct(nx,ny,nvtyps)         !< Fraction of vegetation types
         integer, intent(in) :: nroot  (nx,ny,nvtyps)       !< number of layers in root zone
         integer, intent(in) :: ndeep  (nx,ny,nvtyps)       !< number of layers in deep soil
         type(prognostic3dgrid), intent(in) :: pg3(nx,ny,0:nvtyps)
         real, intent(in) :: Dgl2d(nx,ny,0:nvtyps)       !< depth to groundwater level in each GRID (m)
         real, intent(in) :: GWstorage(nx,ny,0:nvtyps)   !< depth to groundwater level in each GRID (m)


! ===== LOCAL VARIABLES

         integer smlunit
         integer iflow,ig,i,j,k,iv,ntotal_t


! ===== CODE

         call getunit(smlunit)

         open(smlunit,file=trim(simulation_dir)//
     &         sub_catch(isub)//'I_soil', status='unknown')
         do iflow = 1,num_flow(isub)
            do ig=1,ngrid(isub,iflow)
               i=grid_col(isub,iflow,ig)
               j=grid_row(isub,iflow,ig)
! TODO: ntotal -> nothing.  Check how this is read. (PRL)
! This variable (ntotal) is not valid at this point in the code.
!                     write(smlunit,*) iflow, i, j, ntotal, nvtyp
               write(smlunit,*) iflow, i, j, 0, nvtyps
               do iv = 1, nvtyps
                  if( vtypfrct(i,j,iv) > 0.0) then
                     ntotal_t = 1 + nroot(i,j,iv) + ndeep(i,j,iv)

                     write(smlunit,*)
     &                     (pg3(i,j,iv)%vlcprf(k), k=1, ntotal_t),
     &                     Dgl2d(i,j,iv),GWstorage(i,j,iv)
#ifdef SOILDATA_REV2
     &                     ,Sst(i,j),
     &                     wetcanp(i,j,iv),wetg(i,j,iv),
     &                     snowc(i,j,iv),snowg(i,j,iv) !2010.8.27
#endif
                  end if
               end do
            end do
         end do
         close(smlunit)
         call retunit(smlunit)

      end subroutine save_soil_data










       SUBROUTINE basin_to_pointscale_land(rhoair,
     &   i2d, io3d,veg,vl,rm,
!     &pg3,
     &   r3i,soil,r3o,
     &   mrf2,phy2,opc,vp,vio,
!     &   curtim,
     &is_nocond_start,
     &   Ds2d,Dgl2d,zwind_ini,ylat,wetg_currsum,store,Dg,green2d,
     &   lai,fpar,GWstorage,
     &   pg3,
!     &               tcanp,tsfc,tsoil,vlcsfc,vlcrt,vlcdp,vlcprf,
!     &               wetcanp,wetg,snowc,snowg,rst2d,
     &            nroot,ndeep,vcover_in
     &)
         use point_scale

         implicit none

! ===== INCLUDES
         include 'globcst.inc'


! ===== PARAMETERS

         real, intent(in) :: rhoair    !< Air pressure

         type(ground_2d_params), intent(inout) :: i2d

         type(ground_3d_vars), intent(inout) :: io3d
!         type(o3d_vars), intent(out) :: o3d
         type(vegetation), intent(out) :: veg
         type(vegetation_link), intent(out) :: vl

         type(remainder), intent(out) :: rm
!         type(prognostic3dgrid), intent(in) :: pg3
! TODO: Flip this back to IN only (PRL)
         type(remainder3dgrid), intent(inout) :: r3i
         type(soil_params), intent(in) :: soil
         type(remainder3dpoint), intent(out) :: r3o

         type(morphology2), intent(in) :: mrf2
         type(physiology2), intent(in) :: phy2
         type(optical), intent(in) :: opc

         type(vertical_profile_in), intent(inout) :: vp
         type(vertical_profile_inout), intent(inout) :: vio

!         integer, intent(in) :: curtim
!         integer, intent(in) :: year,month,day,hour,minute
         logical, intent(in) :: is_nocond_start

         real, intent(in) :: Ds2d
         real, intent(in) :: Dgl2d
         real, intent(in) :: zwind_ini
         real, intent(in) :: ylat
         real, intent(in) :: wetg_currsum

         real, intent(in) :: store
         real, intent(in) :: Dg

         real, intent(inout) :: green2d

         real, intent(in) :: lai    !< leaf area index
         real, intent(in) :: fpar
         real, intent(in) :: GWstorage


         type(prognostic3dgrid), intent(in) :: pg3

         integer, intent(in) :: nroot
         integer, intent(in) :: ndeep
         real, intent(in) :: vcover_in


! ===== LOCAL VARIABLES

         integer k
         logical IGCM !< Is the LSM connected to a GCM (online vs offline)

! ====== ASSERTIONS
!        assert(mrf2%zlw > 0)
      if (mrf2%leafw <= 0) green2d = 1/mrf2%leafw


! ===== CODE

!
!#######################################################################
!     For SiB2 vegetation type ( 1 ~ 9 )
!     Vegetation parameters initialized from SiB2 database
!#######################################################################
!         vsib2 = vegtyp(i,j,iv)

         veg%zc        = mrf2%zc !(vsib2) ! Canopy inflection height  (m)
         veg%z2        = mrf2%z2 !(vsib2) ! Canopy-top height  (m)
         veg%z1        = mrf2%z1 !(vsib2) ! Canopy-base height  (m)
         veg%chil      = mrf2%chil !(vsib2) ! Leaf angle distribution factor
!         rm%sodep         = Ds2d !(i,j)
         vl%rootd      = mrf2%rootd !(vsib2)
         vl%rootd      = amin1( vl%rootd, Ds2d*0.75 )
         r3o%thick(1) = mrf1%dsfc
         if (r3o%thick(1) == 0) r3o%thick(1) = 1 / r3o%thick(1)
         r3o%thick(2) = vl%rootd - mrf1%dsfc
         r3o%thick(3) = Ds2d - r3o%thick(1) - r3o%thick(2)

         rm%zwind    =  veg%z2 + zwind_ini !2008.04.04
! note:  zref is the reference height from the top height of the water body, the trees or the grassland.
         rm%zmet    =         veg%z2
! zwind and zmet should be modified due to input dataset
! In the dataset JRA25, zmet0 = 2.0m, Zwind0 =10.0m

         veg%tran (1,1) = opc%tranlv !(vsib2)
         veg%tran (1,2) = opc%trandv !(vsib2)
         veg%tran (2,1) = opc%tranln !(vsib2)
         veg%tran (2,2) = opc%trandn !(vsib2)
         veg%ref  (1,1) = opc%reflv  !(vsib2)
         veg%ref  (1,2) = opc%refdv  !(vsib2)
         veg%ref  (2,1) = opc%refln  !(vsib2)
         veg%ref  (2,2) = opc%refdn  !(vsib2)
         veg%vmax0      = phy2%vmax0  !(vsib2)
         veg%effcon     = phy2%effcon !(vsib2)
         veg%gradm      = phy2%gradm  !(vsib2)
         veg%binter     = phy2%binter !(vsib2)
         veg%atheta     = phy2%atheta !(vsib2)
         veg%hlti       = phy2%hlti   !(vsib2)
         veg%hhti       = phy2%hhti   !(vsib2)
         veg%phc        = phy2%phc    !(vsib2)
         veg%respcp     = phy2%respcp !(vsib2)

!         rm%shti           = shti_cst
!         rm%slti           = slti_cst
!         trda           = trda_cst
!         trdm           = trdm_cst
!         trop           = trop_cst
!         rm%btheta         = btheta_cst

         veg%soref (1)  = opc%sorefv !(vsib2)
         veg%soref (2)  = opc%sorefn !(vsib2)
         rm%g1             = mrf1%g1
         rm%g4             = mrf1%g4
         r3o%zlt            = lai      !(i,j,iv)
         
         ! Moiz: To match WEB-DHM-S 2.0.0
         if (snowmodel==1) then
         r3o%vcover = r3o%zlt / mrf2%laimax !(vsib2) ! Miller et al., GRL, 2006
         else
         r3o%vcover = vcover_in
         endif

         r3o%vcover = max(r3o%vcover,1.0e-7)
         r3o%vcover = min(r3o%vcover,1.0)
!#######################################################################
!
!     Derive aeodynamic and physiologic paramters once each day
!
!#######################################################################
         IF (hour == 0 .and. minute == 0) THEN !03/10/2007, 22/02/2013 added minute check (PRL)

!> \TODO : Move this flag (PRL)
            IGCM = .false.   !flag of GCM,  IGCM = 0 : offline

            CALL sibx(rhoair,rm%zwind,rm%zmet,
     1                    veg%z2,veg%zc,veg%z1,mrf1%z0g,
     1                    mrf2%leafl,mrf2%leafw, !leafl_v(vsib2),leafw_v(vsib2),
     1                    veg%chil,r3o%zlt,r3o%vcover,
     2                    rm%g1,rm%g4,IGCM,r3i%roufns,!(i,j,iv),
     2                    r3i%dd2d,
     3                    r3i%cc12d,r3i%cc2v2d,
     4                    r3i%cc2g2d,r3i%fug2d,
     5                    r3i%corb12d,r3i%corb22d,
     6                    r3i%ha2d,r3i%g22d,r3i%g32d)

! ASSERTION:
      if (r3i%dd2d <= 0) stop "SIBX:r3i%dd2d <= 0"

            CALL siby(year,month,day,ylat,!(i,j),
     &                     veg%chil,r3i%gmudmu2d) !gmudmu2d(i,j,iv))
            green2d =
     1                    amax1(r3o%zlt - mrf2%lais, 0.0) / !(vsib2)
     2                    amax1(mrf2%lais, r3o%zlt) !(vsib2)

         END IF
!#######################################################################
         r3o%ha         = r3i%ha2d     !(i,j,iv)
         r3o%z0d        = r3i%roufns   !(i,j,iv)
         r3o%dd         = r3i%dd2d     !(i,j,iv)
! ASSERTION:
      if (r3i%dd2d <= 0) stop "r3i%dd2d <= 0"
         r3o%g2         = r3i%g22d     !(i,j,iv)
         r3o%g3         = r3i%g32d     !(i,j,iv)
         r3o%cc1        = r3i%cc12d    !(i,j,iv)
         r3o%cc2v       = r3i%cc2v2d   !(i,j,iv)
         r3o%cc2g       = r3i%cc2g2d   !(i,j,iv)
         r3o%fug        = r3i%fug2d    !(i,j,iv)
         r3o%corb1      = r3i%corb12d / rhoair   !(i,j,iv) / rhoair !ask 02/22/2006
         r3o%corb2      = r3i%corb22d  !(i,j,iv)
         r3o%green      = green2d  !(i,j,iv)
         r3o%gmudmu     = r3i%gmudmu2d !(i,j,iv)
         r3o%fparc      = fpar !(i,j,iv)


!#######################################################################
!
!     Initialize prognostic variables           2010.8.27 modified by Lei Wang
!
!#######################################################################
!

! 456    ! Moiz: Turn off in case of snapshot
         !if (is_nocond_start) then   ! Moiz: This conditional was turned off for some reason
         !   vio%Dgl    =  Ds2d!(i,j) ! Moiz: testing to match inti condition of WEB-DHM-S 2.0.0
         !   vio%GWst   =  Dg * soil%GWcs
         !else
            vio%Dgl    =  Dgl2d!(i,j,iv) 
            vio%GWst   =  GWstorage!(i,j,iv)
        ! end if

         io3d%tc       = pg3%tcanp     !(i,j,iv)
         io3d%tg       = pg3%tsfc      !(i,j,iv)
         io3d%td       = pg3%tsoil     !(i,j,iv)
         io3d%vlc(1)   = pg3%vlcsfc    !(i,j,iv)
         io3d%vlc(2)   = pg3%vlcrt     !(i,j,iv)
         io3d%vlc(3)   = pg3%vlcdp     !(i,j,iv)

         vp%nrt = nroot!(i,j,iv)
         vp%ndp = ndeep!(i,j,iv)

         DO k = 1, (1 + vp%nrt + vp%ndp)
            vio%lsm(k) = pg3%vlcprf(k) !(i,j,k,iv)
         END DO

         io3d%capac(1) = pg3%wetcanp   !(i,j,iv)
         io3d%capac(2) = pg3%wetg * store / !wetg(i,j,iv)
     &                 max(1.0E-6,wetg_currsum) !note 2010.8.27
!     wetg(i,j,iv)/wetg(i,j,0) is the ratio that how much water in vegation(iv) in grid(i,j)
         io3d%snoww(1) = pg3%snowc     !(i,j,iv)
         io3d%snoww(2) = pg3%snowg     !(i,j,iv)


!******************************************************************
!
!     Updating snow cover
!     Using 8-daily MODIS/Terra Maximum_Snow_Extent
!
!     Lei Wang    4/18/2007
!******************************************************************
!     if(hour == 10) then
!     if(              snoww(1)+snoww(2)< 1.0e-6
!     &           .and.  snowe_8daily(i,j)==200
!     &           .and.  tm < tf                           )then
!     snoww(2) = 0.005
!     elseif(  snoww(1)+snoww(2)> 1.0e-3 .and.
!     &                          snowe_8daily(i,j)==25 ) then
!     snoww(1) = 0.0
!     snoww(2) = 0.0
!     endif
!     end if
!******************************************************************
         io3d%extra_sfc_roff = 0. ! Moiz: Added 
         io3d%rstcanp   = pg3%rst2d     !(i,j,iv)
      end SUBROUTINE basin_to_pointscale_land



      SUBROUTINE pointscale_to_basin_land()
      end SUBROUTINE pointscale_to_basin_land







!     ########################################################################
!     ########################################################################
!     #####                                                              #####
!     #####                                                              #####
!     #####          SUBROUTINE tatoesat                                 #####
!     #####                                                              #####
!     #####                                                              #####
!     ########################################################################
!     ########################################################################
!
!#######################################################################
!
!> \brief
!>    Using tetens formula to derive saturation vapor pressur with air
!>    temperature.
!
!#######################################################################

      pure SUBROUTINE tatoesat(ta,esat)
         implicit none

!#######################################################################

! ===== PARAMETERS
         real, intent(in) :: ta     !< Air temperature (K)
         real, intent(out) :: esat  !< saturation vapor pressur(Pa)


! ===== LOCAL VARIABLES
         real    t                  !< Air temperature (K)


! ===== CODE
         t=ta-273.15
         if(t.gt.0)then
            esat=6.1078*10**((7.5*t)/(237.7+t)) ! Tetens formula for water surface
         else
            esat=6.1078*10**((9.5*t)/(265.3+t)) ! Tetens formula for ice surface
         end if
         esat=esat*100             ! Convert hPa to Pa
         RETURN
      END SUBROUTINE tatoesat




!
!     ##################################################################
!     ##################################################################
!     ######                                                      ######
!     ######                SUBROUTINE radswpart                  ######
!     ######                                                      ######
!     ######                     Developed by                     ######
!     ######             NASA Goddard Space Flight Center         ######
!     ######                                                      ######
!     ##################################################################
!     ##################################################################
!
!#######################################################################
!
!> \brief
!>    Convert short wave radiation into four components: visible direct
!>    radn(1,1),visible diffuse radn(1,2),NIR direct radn(2,1), NIR
!>    diffuse radn(2,2)
!
!#######################################################################

      pure SUBROUTINE radswpart (swdown,sunang,cloud,radn)
         implicit none
!
!
!
!     AUTHOR: P. Sellers, J. Collatz, L. Bounoua
!     Last revision: dec. 15, 1993
!
!     MODIFICATION HISTORY:
!
!
!#######################################################################


! ===== PARAMETERS

         real, intent(in) :: swdown       !< incoming downward solar radiation (w/m2)
         real, intent(in) :: sunang       !< cos(zenith angle)
         real, intent(out) :: cloud       !< cloud fraction
         real, intent(out) :: radn(3,2)   !< radiation component


! ===== LOCAL VARIABLES

         real difrat                      !< diffuse fraction of total rad.
         real vnrat                       !< visisble fraction of total rad.


! ===== CODE

         IF (sunang <= 0.) THEN  ! it's nighttime or we're in shadow (PRL)
           cloud = 0.2  ! Moiz: changed from 0.58        ! TODO: review the physics for this code (PRL)
         ELSE
           cloud = (1160.*sunang - swdown) / (963. * sunang)
           cloud = amax1(cloud,0.)
           cloud = amin1(cloud,1.)
           cloud = amax1(0.2,cloud) ! Moiz: changed from 0.58
         ENDIF

         difrat = 0.0604 / max(0.01,( sunang-0.0223 )) + 0.0683
         IF ( difrat .lt. 0. ) difrat = 0.
         IF ( difrat .gt. 1. ) difrat = 1.

         difrat = difrat + ( 1. - difrat ) * cloud
         vnrat = ( 580. - cloud*464. ) / ( ( 580. - cloud*499. )
     &         + ( 580. - cloud*464. ) )

         radn(1,1) = (1.-difrat)*vnrat*swdown
         radn(1,2) = difrat*vnrat*swdown
         radn(2,1) = (1.-difrat)*(1.-vnrat)*swdown
         radn(2,2) = difrat*(1.-vnrat)*swdown

         RETURN

      END SUBROUTINE radswpart

      end module mod_flux2d

