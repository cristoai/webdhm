
      module hydro
         use time_and_fileio


         contains

!     ##################################################################
!     ##################################################################
!     ######                                                      ######
!     ######                 SUBROUTINE define_subcatchments      ######
!     ######                                                      ######
!     ######                     Developed by                     ######
!     ######     River and Environmental Engineering Laboratory   ######
!     ######                University of Tokyo                   ######
!     ######                                                      ######
!     ##################################################################
!     ##################################################################
!

      subroutine define_subcatchments
         implicit none


! ===== INCLUDES

         include 'globcst.inc' ! para_dir
         include 'hydro.inc'


! ===== LOCAL VARIABLES

         character*100 atmp
         data      kfs1/9*0/
         data      kfs2/81*0/
         data      kfs3/729*0/

         integer nlevel, level, l1, l2, l3, j, parunit

!-------------------------------------------------------------------
!     generate subcatchment name from kfs.dat
!-------------------------------------------------------------------
!

         call getunit(parunit)

         open(parunit,file=trim(para_dir)//'kfs.dat',status='old')
         read(parunit,*) atmp, nlevel
         read(parunit,*)
         do level = 1, nlevel-1
            if(level == 1) then   ! Level 1
               read(parunit,*)
               read(parunit,*) atmp, (kfs1(l1), l1=1,9)
            elseif(level == 2) then ! Level 2
               read(parunit,*)
               do l1=1,9
                  if(kfs1(l1) == 1) read(parunit,*)
     &                     atmp,(kfs2(l1,l2),l2=1,9)
               end do
            elseif(level == 3) then ! Level 3
               read(parunit,*)
               do l1=1,9
                  if(kfs1(l1) == 1) then
                     do l2=1,9
                        if(kfs2(l1,l2) == 1)
     &                     read(parunit,*) atmp,(kfs3(l1,l2,l3), l3=1,9)
                     end do
                  endif
               end do
            endif
         end do
         close(parunit)
         call retunit(parunit)

         num_sub=0
         do l1 = 9, 1, -1          ! Level 1
            if(kfs1(l1) == 0) then
               num_sub=num_sub+1
               sub_catch(num_sub)=
     &                     'ws'//char(48+l1)//char(48+0)//char(48+0)

            else
               do l2 = 9, 1, -1    ! Level 2
                  if(kfs2(l1,l2) .eq. 0) then
                     num_sub=num_sub+1
                     sub_catch(num_sub)=
     &                     'ws'//char(48+l1)//char(48+l2)//char(48+0)

                  else
                     do l3 = 9, 1, -1 ! Level 3
                        num_sub=num_sub+1
                        sub_catch(num_sub)=
     &                     'ws'//char(48+l1)//char(48+l2)//char(48+l3)

                     end do
                  endif
               end do
            endif
         end do

         return
      end subroutine define_subcatchments



      subroutine init_hydro
         implicit none
!
!************************************************************************
!     *
!     read sub-basin parameters                                  *
!     *
!************************************************************************
!     num_flow:   number of flow interval for each sub-catchment
!     ngrid:     number of grids in this flow interval
!     rdx:       length of flow intervals (m)
!     s0:        river-bed slope of the flow interval
!     b:         river width of the flow interval (m)
!     roughness: river roughness (Manning's coefficient) of the flow interval
!     Dr:        river depth of the flow interval (m)
!     grid_row:  row number of grid in this flow interval
!     grid_col:  column number of grid in this flow interval
!-----------------------------------------------------------
!     Modification history: 03/29/2006    WANG Lei
!-----------------------------------------------------------


! ===== INCLUDES

      include 'globcst.inc'
      include 'hydro.inc'
      include 'incl_params.f'


! ===== LOCAL VARIABLES

      character*200 river_file
      integer  i, ii, flunit, parunit, j

#if 0
      character*100 atmp
      data      kfs1/9*0/
      data      kfs2/81*0/
      data      kfs3/729*0/

      integer nlevel, level, l1, l2, l3, j

!-------------------------------------------------------------------
!     generate subcatchment name from kfs.dat
!-------------------------------------------------------------------
!
      call getunit(parunit)
      print *, 'parunit', parunit
      open(parunit,file=trim(para_dir)//'kfs.dat',status='old')
      read(parunit,*) atmp, nlevel
      read(parunit,*)
      do level = 1, nlevel-1
         if(level == 1) then   ! Level 1
            read(parunit,*)
            read(parunit,*) atmp, (kfs1(l1), l1=1,9)
         elseif(level == 2) then ! Level 2
            read(parunit,*)
            do l1=1,9
               if(kfs1(l1) == 1) read(parunit,*)
     &              atmp,(kfs2(l1,l2),l2=1,9)
            end do
         elseif(level == 3) then ! Level 3
            read(parunit,*)
            do l1=1,9
               if(kfs1(l1) == 1) then
                  do l2=1,9
                     if(kfs2(l1,l2) == 1)  read(parunit,*) atmp,
     &                    (kfs3(l1,l2,l3), l3=1,9)
                  end do
               endif
            end do
         endif
      end do
      close(parunit)
      call retunit(parunit)
!
      num_sub=0
      do l1 = 9, 1, -1          ! Level 1
         if(kfs1(l1) == 0) then
            num_sub=num_sub+1
            sub_catch(num_sub)='ws'//char(48+l1)//
     &           char(48+0)//char(48+0)
!     print *,num_sub,'   ',sub_catch(num_sub)
         else
            do l2 = 9, 1, -1    ! Level 2
               if(kfs2(l1,l2) == 0) then
                  num_sub=num_sub+1
                  sub_catch(num_sub)='ws'//char(48+l1)//char(48+l2)//
     &                 char(48+0)
!     print *,num_sub,'   ',sub_catch(num_sub)
               else
                  do l3 = 9, 1, -1 ! Level 3
                     num_sub=num_sub+1
                     sub_catch(num_sub)='ws'//char(48+l1)//char(48+l2)//
     &                    char(48+l3)
!     print *,num_sub,'   ',sub_catch(num_sub)
                  end do
               endif
            end do
         endif
      end do
      print *,'num_sub:',num_sub
#endif

!     read sub_catchment parameters
!-------------------------------------------------------------------
      print *,'Reading catchment parameters ...'
      do ii=1,num_sub
         river_file = sub_catch(ii)//'_river'
!         print *,'CTCH:',ii,sub_catch(ii)

!     change:02/20/2006
         CALL getunit(flunit)   !add
         open(flunit,file=trim(para_dir)//trim(river_file),
     &        status='old')
         read(flunit,*) num_flow(ii)

!     print *,ii,num_flow(ii)
         do i=1,num_flow(ii)
            read(flunit,*) ngrid(ii,i),rdx(ii,i),s0(ii,i),b(ii,i),
     &           roughness(ii,i),Drr(ii,i)
!***********************************************************
!     Calibrate for roughness, and river depth
!     0.03~0.05  (original input)
!***********************************************************

         if (RIVER_ROUGHNESS_OVERRIDE > 0) then
            roughness(ii,i) = RIVER_ROUGHNESS_OVERRIDE
         end if

!***********************************************************
            if(ngrid(ii,i) > 1000)
     &           print *,'more than 1000 grids:',ngrid(ii,i),ii,i
            read(flunit,*)
     &           (grid_row(ii,i,j),grid_col(ii,i,j), j=1,ngrid(ii,i))
            if(s0(ii,i).le.0.001) s0(ii,i)=0.001
         end do
         close(flunit)
         call retunit(flunit)   !add
      end do

      return
      END subroutine init_hydro


      subroutine init_river(Ds2d,Dr2d,Drw2d)
         include 'hydro.inc'
!         include '../dims.inc' ! nx,ny,nvtyps

! ===== PARAMETERS
         real, intent(in)  :: Ds2d(:,:)   !< depth of UZ soil of each GRID (m)
         real, intent(out) :: Dr2d(:,:)   !< river depth in each GRID (m)
         real, intent(out) :: Drw2d(:,:)  !< water depth in the river of each GRID (m)

! ===== LOCAL VARIABLES
         integer i,ii,ig,ix,iy

! ===== CODE
         do ii=1,num_sub
            do i=1,num_flow(ii)
               do ig=1,ngrid(ii,i)
                  iy=grid_row(ii,i,ig)
                  ix=grid_col(ii,i,ig)
                  Dr2d(ix,iy)=Drr(ii,i)

                  Drw2d(ix,iy) = Ds2d(ix,iy)  ! This probably has no meaning
               end do
            end do
         end do

      end subroutine init_river



!     ##################################################################
!     ##################################################################
!     ######                                                      ######
!     ######              SUBROUTINE lateral_routing              ######
!     ######                                                      ######
!     ######                     Developed by                     ######
!     ######     River and Environmental Engineering Laboratory   ######
!     ######                University of Tokyo                   ######
!     ######                                                      ######
!     ##################################################################
!     ##################################################################

!#######################################################################
!
!     Lateral redistribution of water moisture
!
!     River routing and Output results
!
!     AUTHOR: Wang Lei
!
!     Modification history
!     12/11/2006 (WangLei)
!#######################################################################
!
!>     (Module of Lateral Water Flow in sib-gbhm)
!
!>     (1) Using basin unit (Pfafstetter Basin Scheme)
!>     (2) Top-morphology parameterization (River-Hillslope Concept)
!>     (3) Kinematic water (surface and river) routing
!
!#######################################################################

! wetcanp,wetg,snowc,snowg is not used in this routine or its callees (PRL)
! URGENT: remove reliance on Dr2d
      subroutine lateral_routing(
     &     inbasin,area,slp2d,len2d,Dr2d,                ! static topographic input
     &     runoff_inter, runoff_g,Sst,Sstmax0,waterflx,  ! INPUT
     &     q2, Drw2d                                     ! OUTPUT
     &     )

!#######################################################################

         use river_rtg
         implicit none


! ===== INCLUDES

      include 'globcst.inc' ! curtim, tstart, prti.prtj
      include 'hydro.inc'
      include 'dims.inc' ! nx,ny,nvtyps


! ===== PARAMETERS

      integer, intent(in) ::  inbasin(nx,ny) !< flag to show the simulation boundary
      real,    intent(in) ::  area(nx,ny)    !< area of each GRID (m2)
      real,    intent(in) ::  slp2d(nx,ny)   !< slope of ground surface (ND)
      real,    intent(in) ::  len2d(nx,ny)   !< hillslope length of each GRID (m)
      real,    intent(in) ::  Dr2d(nx,ny)    !< river depth in each GRID (m)

!**************************************************************************
!     definition of variables for grids
!**************************************************************************

      real, intent(in)     :: runoff_inter(nx,ny)  !< inter flow <b>(m)</b>
      real, intent(in)     :: runoff_g(nx,ny)      !< exchange between groundwater and river <b>(m)</b>
      real, intent(inout)  :: Sst(nx,ny)
      real, intent(in)     :: Sstmax0(nx,ny)
      real, intent(in)     :: waterflx(nx,ny)

!********************************************
! transient status information for the GBHM run
!********************************************
      real, intent(out)    :: q2(nc,300)
      real, intent(out)    :: Drw2d(nx,ny)   !< water depth in the river of each GRID (m)
!      real, intent(out) :: hhr(100,300) !< only output for diagnostic purposes


! ===== LOCAL VARIABLES

      real qin(300)

      integer i, ig, ii, l1,l2,l3
      integer ix,iy

      character*100 atmp

      integer  flowunit
      integer  flowno(nx,ny)
!***************************************************************************
!     number of grids
!***************************************************************************
      integer m, num_grid(nc)   !the number of grids in each subcatchment
      integer ngrid_outlet      !the total number of grids for whole basin
      integer endflow           !for output Murakami-gauge's data 12/25/2006

!**************************************************************************
!     For output daily average values of each discharge gauges
!**************************************************************************
      integer ifg               !fg: flowgauge
!      integer code_start(nfg), code_end(nfg) ,gaugeloc(nfg)
      integer up_grid(nfg)      !total number of upstream grids

      real    up_rain(nfg, 366) !average daily rainfall of upstream
      real    interflow(nfg, 366)
      real    gwflow(nfg, 366)
      real    sfcflow(nfg, 366)
      integer ounit

!**************************************************************************
!     For output daily Land Surface Temperature
!**************************************************************************

      real    Tland_day(nx,ny,366)
      real    Tland_night(nx,ny,366)
      real    Tair_day(nx,ny,366)
      real    Tair_night(nx,ny,366)
      integer iLST
      real    rnum

!**************************************************************************
!     For output spatial distribution of whole basin
!**************************************************************************

      integer unit1
      character*4  ch4
      character*2  ch2a, ch2b
      character*3  ch3c
      integer      daynum, im, j

!**************************************************************************
      integer  outflag          ! control output or no output
!**************************************************************************

      integer dayinmonth(12)
      data    dayinmonth/31,28,31,30,31,30,31,31,30,31,30,31/


!     read sub_catchment parameters


      outflag = 1
      if(outflag == 1)        then !4/18/2007
         if(curtim == tstart) then
!            call getunit(flowunit)
!            open(flowunit,file='../output/other/flow.asc',
!     :           status='unknown')
!            do iy = 1, ny
!               do ix = 1, nx
!                  flowno(ix,iy) = -9999
!               end do
!            end do
!            write(flowunit,*) 'ncols        183'
!            write(flowunit,*) 'nrows        152'
!            write(flowunit,*) 'xllcorner     267228.28618242'
!            write(flowunit,*) 'yllcorner     4027010.4230879'
!            write(flowunit,*) 'cellsize      500'
!            write(flowunit,*) 'NODATA_value  -9999'
!            do iy = 1, ny
!               write(flowunit,'(183i8)')  (flowno(ix,iy),ix = 1, nx)
!            end do
!            close(flowunit)
!            call retunit(flowunit)
         end if
      end if                    !outflag

!*****************************************************************
!
!     sub-catchment loop
!
!*****************************************************************
!
!     Creat Sub-basin Model
!
!     ii>=24 .and. ii<=40                for Agatsuma Basin
!*****************************************************************

      ii=0
      do l1 = 9, 1, -1          ! Level 1
         if(kfs1(l1) == 0) then
            ii=ii+1

            if(ii>=inisub .and. ii<=finsub) then
               call surface_routing(ii,nx,ny,nvtyps,
     &              inbasin,area,slp2d,len2d,
     &              runoff_inter,runoff_g,
     &              Sst,Sstmax0,waterflx, qin)
               call river_routing(ii,1,l1,l2,l3,qin,q2)
            end if
         else
            do l2 = 9, 1, -1    ! Level 2
               if(kfs2(l1,l2) == 0) then
                  ii=ii+1
                  if(ii>=inisub .and. ii<=finsub) then
                     call surface_routing(ii,nx,ny,nvtyps,
     &                    inbasin,area,slp2d,len2d,
     &                    runoff_inter,runoff_g,
     &                    Sst,Sstmax0,waterflx, qin)
                     call river_routing(
     &                        ii,2,l1,l2,l3,qin, q2)
                  end if
               else
                  do l3 = 9, 1, -1 ! Level 3
                     ii=ii+1
                     if(ii>=inisub .and. ii<=finsub) then
                        call surface_routing(ii,nx,ny,nvtyps,
     &                       inbasin,area,slp2d,len2d,
     &                       runoff_inter,runoff_g,
     &                       Sst,Sstmax0,waterflx, qin)
                        call river_routing(
     &                     ii,3,l1,l2,l3,qin, q2)
                     end if
                     p3(l3)=ii
                  end do
               endif
               p2(l2)=ii
            end do
         endif
         p1(l1)=ii
      end do

      call distribute_river_water_depth(drw2d)

!#######################################################################################
!
!     Add by wanglei11/9 to calculate the average rainfall series of target subcatchment
!
!#######################################################################################
      if(curtim == tstart) then
         do ifg = 1,nfg
            do m = 1, 366
               up_rain  (ifg, m)  = 0.0
               interflow(ifg, m)  = 0.0
               gwflow   (ifg, m)  = 0.0
               sfcflow  (ifg, m)  = 0.0
            end do
         end do
      end if

      RETURN
      END subroutine lateral_routing
!


      subroutine distribute_river_water_depth(Drw2d)
         use river_rtg, only: hhr
         implicit none

         include 'hydro.inc'
         include 'dims.inc' ! nx,ny,nvtyps

 ! ===== PARAMETERS
         real, intent(out)    :: Drw2d(:,:)   !< water depth in the river of each GRID (m)

 ! ===== LOCAL VARIABLES
         integer ii,i,ig
         integer ix,iy

 ! ===== CODE
        do ii=1,num_sub
            do i=1,num_flow(ii)
               do ig=1,ngrid(ii,i)
                  iy=grid_row(ii,i,ig)
                  ix=grid_col(ii,i,ig)
                  Drw2d(ix,iy)=max(0.0, hhr(ii,i))

               end do
            end do
         end do

         return
      end subroutine distribute_river_water_depth




!     ##################################################################
!     ##################################################################
!     ######                                                      ######
!     ######              SUBROUTINE surface_routing              ######
!     ######                                                      ######
!     ######                     Developed by                     ######
!     ######     River and Environmental Engineering Laboratory   ######
!     ######                University of Tokyo                   ######
!     ######                                                      ######
!     ##################################################################
!     ##################################################################
!
! vtypfrct and vegtyp not used in this routine (PRL)
! qin() should be an output of this subroutine (PRL)
      subroutine surface_routing(ii,nx,ny,nvtyps,
     &     inbasin,area,slp2d,len2d,
     &     runoff_inter,runoff_g,
     &     Sst,Sstmax0,waterflx,qin)
         implicit none
!
!**********************************************************************
!     *
!     ***** Surface Routing on Hillslope  *****             *
!     *
!**********************************************************************
!     Variables Definition:
!     ii:     indicator of subcatchment loop
!     i:  indicator of flow interval loop
!     ig: indicator of grid loop within a flow interval
!
!
!     qg_hr:       subsurface runoff of one flow interval (m3/s/m)
!     q_hr:        surface runoff of one simulation unit (m/s, before flow into river)
!     q_hillslope: surface runoff of one simulation unit (m3/s/m, flow into river)
!
!     runoff_sub: subsurface runoff of a grid (m)
!     runoff_sur: surface runoff of a grid (m)
!     Sst:     surface storage of each grid (m)
!     Sst_max: maximum surface detension storage of each grid (m)
!
!     year:  calendar year
!     month: month in calendar
!     day:   day in calendar
!     hour:  hour of the reference day
!     minute:minute of the reference hour
!     second:second of the reference minute
!     jday:  Julian day starting from Jan. 1st of the year
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! ===== INCLUDES

      include 'globcst.inc' ! dthydro, JustStart
      include 'hydro.inc'
      include 'incl_params.f' ! ROUGHNESS_N_CALIB


! ===== PARAMETERS

      integer, intent(in) :: ii              !< indicator of subcatchment loop
      integer, intent(in) :: nx              !< number of gridpoints in x-dir.
      integer, intent(in) :: ny              !< number of gridpoints in y-dir.
      integer, intent(in) :: nvtyps          !< number of vegetation types

      integer, intent(in) :: inbasin(nx,ny)  !< flag to show the simulation boundary
      real, intent(in) :: area(nx,ny)        !< area of each GRID (m2)
      real, intent(in) :: slp2d(nx,ny)       !< slope of ground surface (ND)
      real, intent(in) :: len2d(nx,ny)       !< hillslope length of each GRID (m)

! These two are never used: (PRL)
!      real, intent(in) :: vtypfrct(nx,ny,nvtyps) !< fraction of vegetation types
!      integer, intent(in) :: vegtyp (nx,ny,nvtyps) !< Vegetation type at each point

      real, intent(in) :: runoff_inter(nx,ny)
      real, intent(in) :: runoff_g(nx,ny)
      real, intent(inout) :: Sst(nx,ny)      !< Liquid water on surface
      real, intent(in) :: Sstmax0(nx,ny)     !< Max. surface storage capacity
      real, intent(in) :: waterflx(nx,ny)
      real, intent(out) :: qin(:) !< lateral inflow of one flow interval (m3/sec/m)


! ===== LOCAL VARIABLES

      real dt           ! dt contains the grid title.
      integer i         ! indicator of flow interval loop
      integer ig        ! indicator of grid loop within a flow interval
      integer ic, ir
      integer ix, iy
      real q_surf, q_hillslope, water_depth
      real s0_slope
      real power


!      real depth_hr(nc,300)     ! water depth of river
      real basin_area, area_sub(nc) ! basin-area, subcatchment area
!
!      common /depth_rw1/     depth_hr ! depth of river water
      real  qin_inter           !inter flow from a flow interval
      real  qin_g               !exchanges between groundwater and river in a flow interval
      real  qin_sfc             !surface runoff from a flow interval
      real  qin_tmp             !totoal runoff from a flow interval

      real qin_sfc_tot, qin_inter_tot, qin_g_tot ! (PRL)

      real roughness_n


! ===== CODE

      dt = dthydro              ! second

!     ***********************************************************
!     model initiation
!     ***********************************************************

!     From testing, basin area is always 0 at this point (PRL)
!      print *,'basin area: ',basin_area
       basin_area = 0.0

      if(JustStart == 1) then

         if(ii == 1) then
            basin_area = 0.0
         endif
!
         area_sub(ii)=0.0
         do i=1,num_flow(ii)
            do ig=1,ngrid(ii,i)
               ir=grid_row(ii,i,ig)
               ic=grid_col(ii,i,ig)
               if(area(ic,ir) == -9999.0)
     $              print *,'wrong in grid-area:',ic,ir,area(ic,ir)
               basin_area=basin_area+area(ic,ir)
               area_sub(ii)=area_sub(ii)+area(ic,ir)
            end do
         end do
         if(ii == nc) print *,'total area:',basin_area
      endif

      if(JustStart == 1) return


!****************************************************************
!     start simulation
!****************************************************************
      do i=1,num_flow(ii)
         qin(i)  =  0.0         ! total lateral inflow (m3/s/m)
         qin_tmp   =  0.0       ! lateral inflow from present flow-interval

            qin_inter_tot =  0.0 ! (PRL)
            qin_g_tot     =  0.0
            qin_sfc_tot   =  0.0

         do ig=1,ngrid(ii,i)
            iy=grid_row(ii,i,ig)
            ix=grid_col(ii,i,ig)

            qin_inter =  0.0
            qin_g     =  0.0
            qin_sfc   =  0.0

!**************************************************
! for water body
!**************************************************
            if(waterflx(ix,iy) > 0.1E-10) then ! = prcp - evap on water surface (m)
               qin_sfc   = waterflx(ix,iy)*area(ix,iy)/dt ! m3/s  for waterbody
               qin_inter = 0.0
               qin_g     = 0.0
               qin_tmp   = qin_tmp + (qin_sfc + qin_inter + qin_g) ! m3/s

               qin_sfc_tot = qin_sfc_tot+qin_sfc ! (PRL)

!**************************************************
!     for other land use
!**************************************************
            else
!     1. subsurface runoff in (m2/s)
               qin_inter =  runoff_inter(ix,iy)
               qin_g     =  runoff_g(ix,iy)
!     2. surface runoff, in (m)
               q_surf=amax1(0.0,  Sst(ix,iy)-Sstmax0(ix,iy))
!     hillsope routing: steady constant sheet flow
               q_hillslope=0.0
               if(q_surf .le. 0.1E-8) goto 500 !note: important 2007.01.07
               water_depth = q_surf ! in meter
               Sst(ix,iy) =  Sst(ix,iy) - water_depth

               roughness_n = ROUGHNESS_N_CALIB
               if( len2d(ix,iy).le. 0.0)
     &              print *, 'wrong len2d:', ix,iy,len2d(ix,iy),
     &                 inbasin(ix,iy),area(ix,iy)
               s0_slope    = slp2d(ix,iy) +
     :              water_depth/len2d(ix,iy)
               if(s0_slope .lt. 0.1E-8) s0_slope = 0.1E-8
               power=1.6667
               q_hillslope = dt*sqrt(s0_slope)*
     :              water_depth**power/roughness_n
               q_hillslope = amin1(q_hillslope,water_depth*len2d(ix,iy))
               water_depth = water_depth - q_hillslope/len2d(ix,iy)
               Sst(ix,iy) =  Sst(ix,iy) + water_depth !update surface storage
               water_depth = 0.0
!     surface runoff
               qin_sfc   = q_hillslope/dt !m2/s

 500           qin_tmp = qin_tmp +
     &              (qin_sfc+qin_inter+qin_g)* area(ix,iy)/len2d(ix,iy) ! m3/s

               qin_sfc_tot =qin_sfc_tot+qin_sfc*area(ix,iy)/len2d(ix,iy) ! (PRL)
      qin_g_tot = qin_g_tot+
     &  q_surf*area(ix,iy)/len2d(ix,iy)
         qin_inter_tot =qin_inter_tot+qin_inter*area(ix,iy)/len2d(ix,iy)

            endif

         end do ! ig=1,ngrid(ii,i)
!
!     ----------------------------------------------------------------------
         qin(i)  = qin_tmp / rdx(ii,i) ! m^3/m/s, total lateral inflow

      end do ! i=1,num_flow(ii)


      return
      end subroutine surface_routing
!

      end module hydro
