
      module river_rtg
         use time_and_fileio
         use damops, only: damops_operate_dam => operate_dam, ! Moiz: Damops
     &                     damops_operate_hpp => operate_hpp,
     &                     use_damops,
     &                     idam, ihpp, dam, hpp, dam_s, hpp_s
         private


         include 'hydro.inc'


         public init_river_routing
         public river_routing

         public get_snapshot_size
         public save_snapshot
         public load_snapshot

         public hhr !< TODO: create an 'update river height' procedure to avoid this


!> State variables
         real, save :: q1(nc,300)      !< Outflow from each flow interval at the previous timestep
         real, save :: qlin1(nc,300)

!> For kinematic wave equation the Q-h curve at each flow interval is fixed (no hysteresis)
!> and so hhr is not a state variable.
         real, save :: hhr(nc,300)

      contains


!> If doing snapshots, be sure to load the snapshot after / instead of
!> river initialization
      subroutine init_river_routing(ii)
         implicit none

! ===== INCLUDES
         include 'globcst.inc'
         include 'timectrl.inc' ! rtinitdir
         include 'hydro.inc'    ! num_flow

! ===== PARAMETERS
         integer, intent(in) :: ii           !< subcatchment id
!         real, intent(out) :: hhr(100,300)   !< river height in each flow interval

! ===== LOCAL VARIABLES
         integer ldir_s
         integer smlunit
         integer itmp  ! only used to store unnecessary info

         integer i

! ===== CODE

         do i=1,num_flow(ii)
! This value will be corrected to be equal to qin when the
! processing routine is called during the first timestep.
! The reason for this two-step process is that qin is not yet known
! at this point in the initialization.
! However, the system should not overwrite the qlin1 value from
! the snapshot file if such is available.
            qlin1(ii,i) = -9999
         end do
!*********************************************************
!     specify initial condition
!     integer inicon   ! = 0, Give arbitrary value
!     ! = 1, Import from data file
!*********************************************************

!> \todo Why is initial flow being read from realtime directory?
         if(inicon == 1) then
            print *, 'reading initial I_flow from :',trim(rtinitdir)

            call getunit(smlunit)
            open(smlunit,file=trim(simulation_dir)//
     &           sub_catch(ii)//'I_flow',status='old')
            read (smlunit,*)(itmp,q1(ii,i), i=1,num_flow(ii))
            close(smlunit)
            call retunit(smlunit) !add:02/20/2006
         else

            q1(ii,1)=0.5        !initial waterflow in river channel of flow intervals: 5/10/2007
            do i=2,num_flow(ii)
               q1(ii,i)=q1(ii,i-1)+0.4 !initial flow
            end do
         endif

         call calc_kinematic_river_height(ii,hhr)
         return
      end subroutine




      subroutine calc_kinematic_river_height(ii,hhr)
         implicit none

! ===== INCLUDES
         include 'hydro.inc'    ! num_flow

! ===== PARAMETERS
         integer, intent(in) :: ii           !< subcatchment id
         real, intent(out) :: hhr(100,300)   !< river height in each flow interval

! ===== LOCAL VARIABLES
         real criterion, h1, h2, tmp, f, df
         integer k
         integer i

! ===== CODE
!     calculate initial river water depth
         do i=1,num_flow(ii)
            
            hhr(ii,i) = calc_kinematic_river_height_i(
     &                     q1(ii,i),b(ii,i),s0(ii,i),roughness(ii,i) )

         end do
         return
      end subroutine calc_kinematic_river_height



      function calc_kinematic_river_height_i(q1,b,s0,roughness)
         implicit none

! ===== PARAMETERS
         real :: calc_kinematic_river_height_i
         real, intent(in) :: q1              !< discharge
         real, intent(in) :: b               !< river width
         real, intent(in) :: s0              !< river slope
         real, intent(in) :: roughness       !< river roughness

! ===== LOCAL VARIABLES
         real criterion, h1, h2, tmp, f, df
         integer k

! ===== CODE
!     calculate river water depth
            criterion=0.01
         h1=q1/b
!     b: river width of the flow interval (m)

            k=1
         do while (k < 20)
            tmp=roughness*q1/sqrt(s0)
            f=b*h1-(tmp**0.6)*((b+2.0*h1)**0.4)
               if(k > 1 .and. abs(f) < criterion) exit
            df=b-0.8*((tmp/(b+2.0*h1))**0.6)
               h2=h1-f/df
               h1=h2
               k=k+1
            end do

         calc_kinematic_river_height_i = h2         ! depth of river water
         return
      end function calc_kinematic_river_height_i



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
! Begin snapshotting routines
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

      function get_snapshot_size()
         implicit none

         integer get_snapshot_size

         inquire(iolength=get_snapshot_size) q1,qlin1,hhr
         return
      end function get_snapshot_size


      subroutine save_snapshot(fileid,recid)
         implicit none

         integer, intent(in) :: fileid
         integer, intent(in) :: recid

         real :: hhrt(nc,300)
         integer :: ii,i

         write(fileid, rec=recid) q1,qlin1,hhr

         return
      end subroutine save_snapshot


      subroutine load_snapshot(fileid,recid)
         implicit none

! ===== PARAMETERS
         integer, intent(in) :: fileid
         integer, intent(in) :: recid

! ===== LOCAL VARIABLES
         integer ii

! ===== CODE
         read(fileid, rec=recid) q1,qlin1,hhr
         ! Moiz: These are the river snapshot vars
! Based on a request from Ohta-san, hhr is now stored as a snapshot field
! rather than recalculated from Q.
!         do ii=inisub,finsub
!            call calc_kinematic_river_height(ii,hhr)
!         end do
         return
      end subroutine load_snapshot

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! End snapshotting routines
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


!
!     ##################################################################
!     ##################################################################
!     ######                                                      ######
!     ######              SUBROUTINE river_routing                ######
!     ######                                                      ######
!     ######                     Developed by                     ######
!     ######     River and Environmental Engineering Laboratory   ######
!     ######                University of Tokyo                   ######
!     ######                                                      ######
!     ##################################################################
!     ##################################################################
!
!#######################################################################
!
!> \brief
!> Kinematic River Routing Model
!>
!**********************************************************************
!     Parameters:
!     q1:        discharge of last time step (m^3/s)
!     q2:        discharge of current time step (m^3/s)
!
!     Variables:
!
!     dt_couple: time step for coupling (s), from 'globcst.inc'
!     qin:       lateral inflow of one flow interval (m3/sec/m), from 'hydro.inc'
!
!     beta:      Manning's equation parameter
!     dt:        time step used in routing (s)
!     qlin1:     lateral inflow of last time step (m^3/m/s)
!     qlin2:     lateral inflow of current time step (m^3/m/s)
!     level:     Pfafstetter level
!     l1-3:      Pfafstetter number of local level (1-9)
!     p1:        basin name (number of whole system) of nine basin of level 1
!     p2:        basin name (number of whole system) of nine basin of level 2
!     p3:        basin name (number of whole system) of nine basin of level 3
!
!#######################################################################

      subroutine river_routing(ii,level,l1,l2,l3,qin,q2)
         use discharge_io
!        use damops, only:    ! Moiz: Commented for new damops
!     &                     damops_initialize => initialize,
!     &                     damops_pre_wave_solver => pre_wave_solver,
!     &                     damops_post_wave_solver => post_wave_solver,
!     &                     damops_post_timestep => post_timestep
         implicit none

!#######################################################################

! ===== INCLUDES

      include 'globcst.inc'   ! dthydro, curtim, tstart, month, day, hour, dt_couple, year
      include 'hydro.inc'
      include 'incl_params.f'


! ===== PARAMETERS

      integer, intent(in) :: ii           !< subcatchment id
      integer, intent(in) :: level
      integer, intent(in) :: l1
      integer, intent(in) :: l2
      integer, intent(in) :: l3
      real,    intent(in) :: qin(:)       !< lateral flow at current time step

!      real, intent(inout) :: q1(nc,300)  !< discharge of last time step (m^3/s)
      real, intent(inout) :: q2(nc,300)   !< discharge of current time step (m^3/s)
!      real, intent(inout) :: qlin1(nc,300)!< lateral inflow of last time step (m^3/m/s)
!      real, intent(out) :: hhr(100,300)   !< depth of water in the river (m)


! ===== LOCAL VARIABLES

      integer i

      real dt
      real qlin2(nc,300)

      real Qh(nc,8800)          ! hourly mean discharge (m3/s)
      real Qd(nc,366)           ! daily average discharge (m^3/s)
!     real Qm(nc,12)       ! monthly mean discharge (m3/s)
!     real Qday(nc,366)    ! daily average discharge for gauges or dams (m^3/s)

      save Qh                   ! hourly mean discharge (m3/s)
      save Qd                   ! daily average discharge (m^3/s)
!     save Qm       ! monthly mean discharge (m3/s)

      real q_outlet(366)       !wanglei 2005/10/1
      real qtest

      integer  smlunit, wsunit !add 02/20/2006

      integer iyy,imm,idd,ihh

      integer j, k , mtime,ii1,ii2
      real    h1, tmp, f,df,h2,h,q1_upper,q2_upper,rs,value
      integer outflag


! ===== CODE

! ### INJECT initialize ###
!         call damops_initialize ! Moiz: Commented for new damops


!***********************************************************************
!
!     define time interval
!
!***********************************************************************

      dt =  dthydro             ! in second

!***********************************************************************

      if( (curtim == tstart) .or.
     &     (month==1 .and. day==1 .and. hour==0) ) then
         do j=1,366
            Qd(ii,j)=0.0
         end do
      end if

!
!#######################################################################
!
!     lateral inflow
!
!#######################################################################

      if (curtim == tstart) then
         do i=1,num_flow(ii)
            if (qlin1(ii,i) == -9999) then
               qlin1(ii,i) = qin(i)
            end if
         end do
      end if
      
      do i=1,num_flow(ii)
         qlin2(ii,i)=qin(i)
      end do
!      if (ii == 1) print *,'IJ:',qlin2(ii,:) ! (PRL)
!
!#######################################################################
!
!     River routing
!
!#######################################################################
      mtime = int(dt_couple/dt)
      do 555 j = 1, mtime       !time loop
         do 333 i=1,num_flow(ii) !river segment loop
!
!     junction boundary condition: define the river network
            q1_upper=0.
            q2_upper=0.
            if(i == 1) then
               if(level == 1) then
                  if(l1 == 9 .or. mod(l1,2) == 0) then
                     q1_upper=0.
                     q2_upper=0.
                  else
                     ii1=p1(l1+1)
                     ii2=p1(l1+2)
                     q1_upper=q1(ii1,num_flow(ii1))+
     &                    q1(ii2,num_flow(ii2))
                     q2_upper=q2(ii1,num_flow(ii1))+
     &                    q2(ii2,num_flow(ii2))
                  endif
               elseif(level == 2) then
                  if(l2 == 9 .or. mod(l2,2) == 0) then
                     q1_upper=0.
                     q2_upper=0.
                     if(l2 == 9 .and. (l1.ne.9.and.mod(l1,2).ne.0)) then !taken from upper level
                        ii1=p1(l1+1)
                        ii2=p1(l1+2)
                        q1_upper=q1(ii1,num_flow(ii1))+
     &                       q1(ii2,num_flow(ii2))
                        q2_upper=q2(ii1,num_flow(ii1))+
     &                       q2(ii2,num_flow(ii2))
                     endif
                  else
                     ii1=p2(l2+1)
                     ii2=p2(l2+2)
                     q1_upper=q1(ii1,num_flow(ii1))+
     &                    q1(ii2,num_flow(ii2))
                     q2_upper=q2(ii1,num_flow(ii1))+
     &                    q2(ii2,num_flow(ii2))
                  endif
               elseif(level == 3) then
                  if(l3 == 9 .or. mod(l3,2) == 0) then
                     q1_upper=0.
                     q2_upper=0.
                     if(l3 == 9 .and. (l2 /= 9.and.mod(l2,2) /= 0)) then !taken from upper level
                        ii1=p2(l2+1)
                        ii2=p2(l2+2)
                        q1_upper=q1(ii1,num_flow(ii1))+
     &                       q1(ii2,num_flow(ii2))
                        q2_upper=q2(ii1,num_flow(ii1))+
     &                       q2(ii2,num_flow(ii2))
                     endif
                  else
                     ii1=p3(l3+1)
                     ii2=p3(l3+2)
                     q1_upper=q1(ii1,num_flow(ii1))+
     &                    q1(ii2,num_flow(ii2))
                     q2_upper=q2(ii1,num_flow(ii1))+
     &                    q2(ii2,num_flow(ii2))
                  endif
               endif
            else
               q1_upper=q1(ii,i-1)
               q2_upper=q2(ii,i-1)
            endif               !end of river network definition
!

! ### INJECT pre_wave_solver ###
!         call damops_pre_wave_solver(sub_catch(ii),i,q1_upper,q2_upper)


        rs=roughness(ii,i)

 111     call nkws(dt,rdx(ii,i),b(ii,i),s0(ii,i),rs,
     &           qlin1(ii,i),qlin2(ii,i),q1_upper,q1(ii,i),
     &           q2_upper,q2(ii,i),hhr(ii,i))
         
         ! Moiz: Damops
         if (use_damops > 0) then
            do idam=1,size(dam)
               if ((ii == dam(idam)%subcatch) .and.
     &             (i == dam(idam)%flowint)) then
                     !write(*,*) dam(idam)%damname
                     !write(*,*) dam_s(idam)
                     !write(*,*) 'Before:',q2(ii,i)                    
                     call damops_operate_dam(idam,q2(ii,i),use_damops)
                     !write(*,*) 'After:', q2(ii,i)
                     !write(*,*)
               endif
            enddo
            
            !do ihpp=1,size(hpp)
            !   if ((ii == hpp(ihpp)%subcatch) .and.
!     &             (i == hpp(ihpp)%flowint)) then
            !         call damops_operate_hpp(ihpp,q2(ii,i))
!               elseif (ihpp == 8) then  ! Moiz Temporary: Otozawa is out of bounds
!                                          ! This is called 9 times. Need to find a better
!                                          ! Conditional; the result however is correct
!                     call damops_operate_hpp(ihpp,q2(ii,i))
            !   endif
            !enddo
         endif
         ! Moiz: Damops


! ### INJECT post_wave_solver ###
!         q2(ii,i) = damops_post_wave_solver(
!     &                  sub_catch(ii),i,q1(ii,i),q2(ii,i)) ! Moiz: uncomment later

 333     continue   !river segment loop


! ***********************************
! Move data from current time step to previous time step
         do i=1,num_flow(ii)
            qlin1(ii,i)=qlin2(ii,i)
            q1(ii,i)=q2(ii,i)
         end do
! ***********************************

 555     continue  ! time loop



! ### INJECT post_timestep ###
!         call damops_post_timestep()


!> \todo: get rid of this debug info, rewrite as generalized (PRL)
      if(sub_catch(ii) == GAUGEPOINT_SUBCATCHMENT ) then          !murakami
         if (GAUGEPOINT_FLOWINT <= 0) then
            write(*,'(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,1f14.5)')  
     &      year,'-',month,'-',day,'-',hour,':',minute,':',second,
     &      q2(ii, num_flow(ii))
         else
            write(*,'(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,1f14.5)')  
     &      year,'-',month,'-',day,'-',hour,':',minute,':',second,
     &      q2(ii, max(1,GAUGEPOINT_FLOWINT))
         end if
      end if

!####################################################################72
!     24.0 if timestep = 1.0h;  48.0 if timestep = 0.5h:  02/06/2006
!      Qd(ii,jday)=Qd(ii,jday)+q2(ii,num_flow(ii))/24.0 !note by wanglei10/14

      ! if( (sub_catch(ii) == GAUGEPOINT_SUBCATCHMENT)
      !&                     .and. (hour == 23) ) then  !output daily discharge  Asif
      !   print *, year,month,day
      !end if

! TODO: Is it possible to pass a subarray with the second dimension varying
! without doing an array copy?
!      call save_discharge_info_bin(sub_catch(ii),q2(ii,:),num_flow(ii))


!####################################################################72
!
!     Save the status at end of the simulation

      outflag = 0               !SCE calibrate

      if ((inicon == 0) .and. (outflag == 1)) then
         if( hour == 23) then
            call save_discharge_info_text(
     &               sub_catch(ii), q2(ii,:),num_flow(ii))
         end if
      end if                    !outflag

!      call save_discharge_data()

      return
      end subroutine river_routing




!#######################################################################
!
!> \brief
!>     Sub-program of Nonlinear Kinematic Wave Scheme solving of river routing
!>     Calculate water depth using Newton's method
!
!     Definition of Variales:
!     q0 -- discharge of last time step;
!     q -- discharge of current time step;
!     qlin0 -- lateral inflow of last time step;
!     qlin -- lateral inflow of current time step;
!     h -- water depth;
!     p -- wetted perimeter;
!     b -- river width;
!     s0 -- river slope
!
!#######################################################################

      subroutine nkws(dt,dx,b,s0,roughness,qlin0,qlin,
     &     tmpq01,tmpq02,tmpq1,qtr,h)
        implicit none

!#######################################################################

! ===== PARAMETERS

         real, intent(in) :: dt     !< time delta
         real, intent(in) :: dx     !< streamlength delta
         real, intent(in) :: b      !< river width
         real, intent(in) :: s0     !< river slope
         real, intent(in) :: roughness
         real, intent(in) :: qlin0  !< lateral inflow of last time step
         real, intent(in) :: qlin   !< lateral inflow of current time step
         real, intent(in) :: tmpq01 !< discharge of last time step (upstream)
         real, intent(in) :: tmpq02 !< discharge of last time step (current point)
         real, intent(in) :: tmpq1  !< discharge of current time step (upstream)
         real, intent(out) :: qtr   !< discharge of current time step (current point)
         real, intent(out) :: h     !< water depth


! ===== LOCAL VARIABLES

         integer k
         real q0(2)     ! discharge of last time step
         real q(2)      ! discharge of current time step
         real beta
         real criterion
         real h1,h2
         real tmp,ctmp
         real f,df
         real p         ! wetted perimeter
         real alfa
         real aa,bb,cc
         real qq1,qq2


! ===== CODE

         beta=0.6
         q0(1)=tmpq01
         q0(2)=tmpq02
         q(1)=tmpq1
         criterion=0.005
         k=1
         h1=q0(2)/b
 15      tmp=roughness*q0(2)/sqrt(s0)
         f=b*h1-(tmp**0.6)*((b+2.0*h1)**0.4)
         if((k > 1) .and. (abs(f) < criterion)) goto 18
         df=b-0.8*((tmp/(b+2.0*h1))**0.6)
         h2=h1-f/df
         h1=h2
         if(k >= 20) goto 18
         k=k+1
         goto 15
 18      h=h2
         p=b+2.0*h


!#######################################################################
!     the initial discharge estimated using linear scheme

         alfa=(roughness*p**(2.0/3.0)/sqrt(s0))**0.6

         if((q0(2)+q(1)) .le. 0.0) then
            cc=0.0
            goto 22
         endif

         cc=(0.5*(q0(2)+q(1)))**(beta-1.0)
 22      aa=dt*q(1)/dx+alfa*beta*q0(2)*cc+0.5*(qlin0+qlin)*dt
         bb=dt/dx+alfa*beta*cc
         qq1=aa/bb

         if(qq1 .le. 0.1e-3) then
            qq2=0.0
            goto 30
         endif

!#######################################################################
!     Using Newton's method to calculate discharge
         criterion=0.005
         ctmp=dt*q(1)/dx+alfa*q0(2)**beta+0.5*dt*(qlin+qlin0)
         k=1
 20      f=dt*qq1/dx+alfa*qq1**beta-ctmp
         if((k > 1) .and. (abs(f) <= criterion)) goto 30
         df=dt/dx+alfa*beta*qq1**(beta-1.0)
         qq2=qq1-f/df

         if(qq2 <= 0.1e-3) then
            qq2=0.0
            goto 30
         endif

         qq1=qq2
         if(k >= 20) then
            print *,'NKWS: 20step limitA'
            goto 30
         end if
         k=k+1
         goto 20
 30      q(2)=qq2
         if(q(2) .le. 0.1e-3) q(2)=0.0
         qtr=q(2)

         criterion=0.01
         k=1
         h1=q(2)/b
 45      tmp=roughness*q(2)/sqrt(s0)
         f=b*h1-(tmp**0.6)*((b+2.0*h1)**0.4)
         if((k > 1) .and. (abs(f) < criterion)) goto 50
         df=b-0.8*((tmp/(b+2.0*h1))**0.6)
         h2=h1-f/df
         h1=h2
         if(k >= 20) then
            print *,'NKWS: 20step limitB'
            goto 50
         end if
         k=k+1
         goto 45
 50      h=h2
      return
      end subroutine nkws



      end module river_rtg
