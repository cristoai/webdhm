      module damops
        use time_and_fileio
      !====================================================================
      ! Purpose:
      ! --------
      !   This module contains all the variables, subroutines and functions
      !   needed for dam operations
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 17/11/2018          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      implicit none
	    save
      ! Data Dictionary:
      ! ----------------
      
      ! Real Variable Precision KINDs
      integer, parameter :: r4 = 4
      integer, parameter :: r8 = 8
      
      ! Flag for using damops
      integer :: use_damops = 0

      ! Temp iterator for use where different subroutines are called in WEB-DHM
      integer :: idam, ihpp
      

      ! Time varying data of each dam
      type dam_data_time
           real(r4), allocatable, dimension(:,:,:,:,:,:) :: inflow     ! Inflow (m3/s)
           real(r4), allocatable, dimension(:,:,:,:,:,:) :: outflow    ! Outflow (m3/s)
           real(r4), allocatable, dimension(:,:,:,:,:,:) :: spills     ! Spills (m3/s)
           real(r4), allocatable, dimension(:,:,:,:,:,:) :: t_outflow  ! Total Outflow (m3/s)
           real(r8), allocatable, dimension(:,:,:,:,:,:) :: wl         ! Water Level (m)
      end type dam_data_time

      ! Time invariant properties of each dam
      type dam_data
          ! Properties
          character(len=20)                     :: damname         ! Name
          integer                               :: subcatch        ! Subcatchment (code in subcatchment.dat)
          integer                               :: flowint         ! Flow Interval
          real(r4)                              :: h               ! Height (m)
          real(r8)                              :: gross_store     ! Gross Storage (m3)
          real(r8)                              :: eff_store       ! Effective Storage (m3)
          real(r8)                              :: lwl             ! Low Water Level (m)
          real(r8)                              :: pwl             ! Preliminary Release Water Level (m)
          real(r8)                              :: nwl             ! Normal (Constant) Water Level (m)
          real(r8)                              :: hwl             ! High Water Level (m)
          real(r4)                              :: qmax_spill      ! Maximum Spillway Discharge (m3/s)
          real(r4)                              :: qmin_spill      ! Minimum Spillway Discharge (m3/s)
          real(r4)                              :: qmax_out        ! Maximum Gates Discharge (m3/s)
          real(r4)                              :: qmin_out        ! Minimum Gates Discharge (m3/s)
          real(r4)                              :: flood_def       ! Flood Defintion in terms of discharge (m3/s)
          real(r4)                              :: des_flood       ! Design Flood Discharge (m3/s)        
          logical                               :: override_wl     ! Overrride Initial WaterLevel (m)
          real(r8)                              :: init_wl_override! Intial Water Level for Override (m el.)
          logical                               :: use_fixed_wl    ! Switch to Use Fixed Water Level
          real(r8)                              :: fixed_wl        ! Value of Fixed Water Level (m el.)          
          ! VH Curve
          real(r8),    allocatable, dimension (:,:) :: vh              ! VH Curve

          ! Time varying data
          type(dam_data_time) :: tvar ! Inflow, Outflow, Spills, Total Outflow Indexed by [YYYY:MM:DD:HH:MM:SS]

      end type dam_data
      
      ! Time varying data of each HPP
      type hpp_data_time
         real(r4), allocatable, dimension(:,:,:,:,:,:)    :: inflow         ! Inflow (m3/s)
         real(r4), allocatable, dimension(:,:,:,:,:,:)    :: outflow        ! Outflow (m3/s)
      end type hpp_data_time

      ! Time invariant properties of each HPP
      type hpp_data
        ! Properties
        character(len=20)                     :: hppname            ! Name
        integer                               :: subcatch           ! Subcatchment (code in subcatchment.dat)
        integer                               :: flowint            ! Flow Interval
        real(r4)                              :: nethead            ! Net Head (m)
        real(r4)                              :: qmax               ! Maximum Discharge (m3/s)
        real(r4)                              :: qmin               ! Minimum Discharge (m3/s)
        real(r4)                              :: installed_capacity ! Installed Capacity (MW)
        real(r4)                              :: efficiency         ! Efficiency of Powerplant
        real(r4)                              :: rate               ! Generation Rate (MW/(m3/s))
        logical                               :: use_ph             ! Switch to use PH Curves
        
        ! PH Curve
        real(r8),   allocatable, dimension(:,:) :: ph               ! PH Curve [Only Avaiable for Kurobe]
        
        ! Time varying data
        type(hpp_data_time) :: tvar ! Inflow, Outflow Indexed by [YYYY:MM:DD:HH:MM:SS]
        
      end type hpp_data


      
      ! Dam Reservoir State Variables
      type dam_state
        real(r4) :: q1in        ! Inflow (m3/s) at previous timestep
        real(r4) :: q2in        ! Inflow (m3/s) at current timestep
        real(r4) :: q1tout      ! Total Outflow = Outflow + Spills + Hydropower (m3/s) at previous timestep
        real(r4) :: q2tout      ! Total Outflow = Outflow + Spills + Hydropower (m3/s) at current timestep       
        real(r4) :: q1out       ! Outflow (m3/s) at previous timestep
        real(r4) :: q2out       ! Outflow (m3/s) at current timestep
        real(r4) :: q1hout      ! Outflow (m3/s) to all hydropower plants at previous timestep
        real(r4) :: q2hout      ! Outflow (m3/s) to all hydropower plants at current timestep
        real(r4) :: q1sout      ! Outflow (m3/s) through spillway at previous timestep
        real(r4) :: q2sout      ! Outflow (m3/s) through spillway at current timestep
        real(r8) :: wl1         ! Water Level (m) at previous timestep
        real(r8) :: wl2         ! Water Level (m) at current timestep
        real(r8) :: v1          ! Reservoir Storage (m3) at previous timestep
        real(r8) :: v2          ! Reservoir Storage (m3) at current timestep
      end type dam_state
      
      
      ! Hydropower Plant State Variables
      type hpp_state
        real(r4) :: power       ! Power Generation (MW) at current timestep
        real(r4) :: energy      ! Energy Production (kWh) at current timestep
        real(r4) :: nethead1    ! Net head (m) at previous timestep
        real(r4) :: nethead2    ! Net head (m) at current timestep
        real(r4) :: q1in        ! Inflow (m3/s) at previous timestep
        real(r4) :: q2in        ! Inflow (m3/s) at current timestep
        real(r4) :: q1out       ! Outflow (m3/s) at previous timestep
        real(r4) :: q2out       ! Outflow (m3/s) at current timestep
        ! Note: In the case of HPP, Inflow = Outflow, however, it is not clear if that is always true at complex junctions
        !       in the of the Kurobe Power System. Need to calrify this from Kansai Electric Power Company.
      end type hpp_state
      
      
      type(dam_data), allocatable, dimension(:) :: dam
      type(hpp_data), allocatable, dimension(:) :: hpp
      
      type(dam_state), allocatable, dimension(:) :: dam_s
      type(hpp_state), allocatable, dimension(:) :: hpp_s
      
      real(r8) :: test_wl = 14.5
      
      contains
      
      function get_snapshot_size()
      !
      ! Purpose:
      ! --------
      !   This function get the record length of dam and hpp state variables. This was
      !   developed for integration of Kurobe Dam System Model with WEB-DHM-S 4.0.0 snapshotting
      !   feature
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 26/09/2019          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------   
        implicit none

        integer get_snapshot_size

        inquire(iolength=get_snapshot_size)
     &          dam_s, hpp_s


      end function get_snapshot_size

      subroutine save_snapshot(fileid,recid)
        !
        ! Purpose:
        ! --------
        !   This subroutine saves the dam and hpp snapshot. This was
        !   developed for integration of Kurobe Dam System Model with WEB-DHM-S 4.0.0 snapshotting
        !   feature
        !
        ! Record of Revisions:
        ! --------------------
        !
        !   Date              Programmer                 Description of Change
        !   ----              ----------                 ---------------------
        ! 26/09/2019          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
        !====================================================================
        !
        ! Data Dictionary:
        ! ----------------   
          implicit none
  
          integer, intent(in) :: fileid
          integer, intent(in) :: recid
  
          write(fileid,rec=recid)
     &         dam_s, hpp_s
  
      end subroutine save_snapshot

      subroutine load_snapshot(fileid,recid)
        !
        ! Purpose:
        ! --------
        !   This subroutine loads in the dam and hpp snapshot. This was
        !   developed for integration of Kurobe Dam System Model with WEB-DHM-S 4.0.0 snapshotting
        !   feature
        !
        ! Record of Revisions:
        ! --------------------
        !
        !   Date              Programmer                 Description of Change
        !   ----              ----------                 ---------------------
        ! 26/09/2019          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
        !====================================================================
        !
        ! Data Dictionary:
        ! ----------------   
          implicit none
  
          integer, intent(in) :: fileid
          integer, intent(in) :: recid
          integer             :: i
  
          read(fileid,rec=recid)
     &         dam_s, hpp_s

          ! Water Level Override
          do i = 1, size(dam_s)
            if (dam(i)%override_wl == .true.) then
            dam_s(i)%wl1 = 
     &         store_to_wl(el_to_store(dam(i)%init_wl_override,
     &         dam(i)%vh),dam(i)%vh)
            dam_s(i)%v1 = el_to_store(dam(i)%init_wl_override,
     &                    dam(i)%vh)
            endif
          enddo
  
      end subroutine load_snapshot

      subroutine damops_load_snapshot
      !
      ! Purpose:
      ! --------
      !   This subroutine load the dam states and hpp states from snapshots
      !   if inicon == 1
      !   Note: This subroutine has been superseeded in WEB-DHM-S 4.0.0 so it is not used
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 28/11/2018          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------   
        implicit none
        include 'globcst.inc'
        
        integer :: idam, ihpp       ! Iterators
        
        integer           :: lu                       ! I/O Unit
        integer           :: io_status                ! I/O Status
        character(len=200):: io_message               ! I/O Message
            
            ! Dam State
            open(newunit=lu,file=trim(simulation_dir)//'dam_s.snap',
     &           iostat=io_status,iomsg=io_message,status='OLD')
            if (io_status/=0) then     
                write(*,*) io_message
            else
                do idam=1,size(dam_s)
                    read(lu,*) dam_s(idam)
                enddo
                close(lu)
            endif
            
            
            ! HPP State
            open(newunit=lu,file=trim(simulation_dir)//'hpp_s.snap',
     &           iostat=io_status,iomsg=io_message,status='OLD')
            if (io_status/=0) then     
                write(*,*) io_message
            else
                do ihpp=1,size(hpp_s)
                    read(lu,*) hpp_s(ihpp)
                enddo
                close(lu)
            endif
      
      
      end subroutine damops_load_snapshot
   
      
      subroutine damops_save_snapshot
      !
      ! Purpose:
      ! --------
      !   This subroutine saves the dam states and hpp states as snapshots
      !   if inicon == 0
      !   Note: This subroutine has been superseeded in WEB-DHM-S 4.0.0 so it is not used
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 28/11/2018          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------
        implicit none
        include 'globcst.inc'
        
        integer :: idam, ihpp       ! Iterators
        
        integer           :: lu                       ! I/O Unit
        integer           :: io_status                ! I/O Status
        character(len=200):: io_message               ! I/O Message
        
            ! Dam State
            open(newunit=lu,file=trim(simulation_dir)//'dam_s.snap',
     &           iostat=io_status,iomsg=io_message,status='UNKNOWN')
            if (io_status/=0) then     
                write(*,*) io_message
            else
                do idam=1,size(dam_s)
                    write(lu,*) dam_s(idam)
                enddo
                close(lu)
            endif
            
            ! HPP State
            open(newunit=lu,file=trim(simulation_dir)//'hpp_s.snap',
     &           iostat=io_status,iomsg=io_message,status='UNKNOWN')
            if (io_status/=0) then     
                write(*,*) io_message
            else
                do ihpp=1,size(hpp_s)
                    write(lu,*) hpp_s(ihpp)
                enddo
                close(lu)
            endif

      end subroutine damops_save_snapshot
      
      
      
      !=====================================================================
      subroutine damops_output
      !
      ! Purpose:
      ! --------
      !   This subroutines writes all the state variables for each dam and hpp
      !   at the end of each timestep
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 28/11/2018          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------
        implicit none
        include 'globcst.inc'

        integer :: idam, ihpp ! Iterators        
      
        ! Writing Dam State Variables
        do idam=1,size(dam)
            open(unit=(14000+idam),
     &           file='output/dam/'
     &           //trim(dam(idam)%damname)//'.hourly',
     &           status='unknown')
            write((14000+idam),
     &            '(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,14f20.5)')
     &  year,'-',month,'-',day,'-',hour,':',minute,':',second,
     &  dam_s(idam)    
        enddo
        
        ! Writing Powerplant State Variables
        do ihpp=1,size(hpp)
            open(unit=(15000+ihpp),
     &           file='output/powerplant/'
     &           //trim(hpp(ihpp)%hppname)//'.hourly',
     &           status='unknown')
            write((15000+ihpp),
     &            '(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,8f20.5)')
     &  year,'-',month,'-',day,'-',hour,':',minute,':',second,
     &  hpp_s(ihpp)     
        enddo
      
      end subroutine
      
      
      !=====================================================================
      subroutine operate_dam(idam,q2,use_damops)
      !
      ! Purpose:
      ! --------
      !   This is the main subroutine with operates the dam and updates
      !   all the necessary variables
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 23/11/2018          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------
      implicit none
      include 'globcst.inc'
      
      integer, intent(in)   :: idam
      real(r4),intent(inout):: q2
      integer               :: use_damops 
      
      
      ! Local Variable
      real(r4) :: ratio
      

      select case(use_damops)

      case(1)!================Based on Data==========================
      !++++++++++++++++++++++++++++Inflow++++++++++++++++++++++++++++
      dam_s(idam)%q2in = q2
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      !+++++++++++++++++++++++++++Ratio++++++++++++++++++++++++++++++
!      if (time_lookup_r4(dam(idam)%tvar%inflow,
!     &                             year,month,day,
!     &                             hour,minute,second) /= 0) then
!      ratio = dam_s(idam)%q2in/
!     &              time_lookup_r4(dam(idam)%tvar%inflow,
!     &                             year,month,day,hour,minute,second)
!     
!     else
!      ratio = 1
!      endif
      ratio = 1 ! Moiz: Override
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      ! ++++++++++++++++++++++++++++Outflow++++++++++++++++++++++++++
      dam_s(idam)%q2out = ratio*
     &              time_lookup_r4(dam(idam)%tvar%outflow,
     &                             year,month,day,hour,minute,second)
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      !+++++++++++++++++++++++++++++Spill++++++++++++++++++++++++++++
!      dam_s(idam)%q2sout = ratio*
!     &              time_lookup_r4(dam(idam)%tvar%spills,
!     &                             year,month,day,hour,minute,second)
      dam_s(idam)%q2sout = 0
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      
      ! +++++++++++++++++++++++Hydropower Outflow+++++++++++++++++++
      ! Calculation of hydropower outflow is different for each dam
      ! At this moment connectivity to hydropower plant is hardcoded
      ! Need to find a way to generalize this
      
      select case(idam)
      
      case(1) !-------------------Kurobe Dam------------------------
      !Kurobe 4
      hpp_s(1)%q2in = ratio*
     &              time_lookup_r4(hpp(1)%tvar%inflow,
     &                             year,month,day,hour,minute,second)
      hpp_s(1)%q2in = min(hpp_s(1)%q2in,hpp(1)%qmax)
      dam_s(idam)%q2hout = hpp_s(1)%q2in
      !-------------------------------------------------------------
      
      
      case(2) !-----------------Sennindani Dam----------------------
      !Shin-Kurobe 3
      hpp_s(2)%q2in = ratio*
     &              time_lookup_r4(hpp(2)%tvar%inflow,
     &                             year,month,day,hour,minute,second)
      hpp_s(2)%q2in = min(hpp_s(2)%q2in,hpp(2)%qmax)
 
      !Kurobe 3
      hpp_s(3)%q2in = ratio*
     &              time_lookup_r4(hpp(3)%tvar%inflow,
     &                             year,month,day,hour,minute,second)
      hpp_s(3)%q2in = min(hpp_s(3)%q2in,hpp(3)%qmax)
      
      ! Sennindani Dam = NK3 + K3
      dam_s(idam)%q2hout = hpp_s(2)%q2in + hpp_s(3)%q2in
     
      !Shin-Kurobe 2 ! This is not from any dam, directly connected to NK3
      hpp_s(4)%q2in = hpp_s(2)%q2in !NK3 --> NK2 --> River
      !--------------------------------------------------------------
     
      
      case(3) !----------------Koyadaira Dam------------------------
      ! Kurobe 2
      hpp_s(5)%q2in = ratio*
     &              time_lookup_r4(hpp(5)%tvar%inflow,
     &                             year,month,day,hour,minute,second)
      hpp_s(5)%q2in = min(hpp_s(5)%q2in,hpp(5)%qmax)
      
      dam_s(idam)%q2hout = hpp_s(5)%q2in
      !--------------------------------------------------------------
      
      
      case(4) !----------------Dashidaira Dam------------------------
      ! Dashidaira
      hpp_s(6)%q2in = ratio*
     &              time_lookup_r4(hpp(6)%tvar%inflow,
     &                             year,month,day,hour,minute,second)
      hpp_s(6)%q2in = min(hpp_s(6)%q2in,hpp(6)%qmax)
     
      ! Shin-Yana
      hpp_s(7)%q2in = ratio*
     &              time_lookup_r4(hpp(7)%tvar%inflow,
     &                             year,month,day,hour,minute,second)
      hpp_s(7)%q2in = min(hpp_s(7)%q2in,hpp(7)%qmax)     
     
      ! Otozawa
      hpp_s(8)%q2in = ratio*
     &              time_lookup_r4(hpp(8)%tvar%inflow,
     &                             year,month,day,hour,minute,second)
      hpp_s(8)%q2in = min(hpp_s(8)%q2in,hpp(8)%qmax)     
      
      dam_s(idam)%q2hout = hpp_s(6)%q2in +
     &                     hpp_s(7)%q2in +
     &                     hpp_s(8)%q2in
      !---------------------------------------------------------------
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      end select
      
      !++++++++++++++++++++++Total Outflow++++++++++++++++++++++++++++
      dam_s(idam)%q2tout = dam_s(idam)%q2out + 
     &                     dam_s(idam)%q2sout + 
     &                     dam_s(idam)%q2hout
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      case(2)!================Forecast Model==========================
      !++++++++++++++++++++++++++++Inflow++++++++++++++++++++++++++++
      dam_s(idam)%q2in = q2
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      ! ++++++++++++++++++++++++++++Outflow++++++++++++++++++++++++++
      dam_s(idam)%q2out = maintenance_release(dam(idam),
     &                    year,month,day,hour,0,0)
      
      !write(*,*) dam(idam)%damname, dam_s(idam)%q2out
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      !+++++++++++++++++++++++++++++Spill++++++++++++++++++++++++++++
      dam_s(idam)%q2sout = 0
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      
      ! +++++++++++++++++++++++Hydropower Outflow+++++++++++++++++++
      ! Calculation of hydropower outflow is different for each dam
      ! At this moment connectivity to hydropower plant is hardcoded
      ! Need to find a way to generalize this
      
      select case(idam)
      
      case(1) !-------------------Kurobe Dam------------------------
      !Kurobe 4
      hpp_s(1)%q2in =  time_lookup_r4(hpp(1)%tvar%inflow,
     &                 year,month,day,hour,minute,second)
      hpp_s(1)%q2in = min(hpp_s(1)%q2in,hpp(1)%qmax)
      dam_s(idam)%q2hout = hpp_s(1)%q2in
      !-------------------------------------------------------------
      
      
      case(2) !-----------------Sennindani Dam----------------------
      !Shin-Kurobe 3
      
      !call K4_to_NK3K3(hpp_s(1)%q2in,hpp_s(2)%q2in,dam_s(idam)%q2in)

      !hpp_s(3)%q2in = dam_s(idam)%q2in - dam_s(idam)%q2out
      
      !if (hpp_s(3)%q2in > hpp(3)%qmax) then
      !  hpp_s(3)%q2in = hpp(3)%qmax
      !endif
      

      ! Shin-Kurobe 3
      hpp_s(2)%q2in = dam_s(idam)%q2in - dam_s(idam)%q2out
      
      if (hpp_s(2)%q2in > hpp(2)%qmax) then
        hpp_s(3)%q2in = hpp_s(2)%q2in - hpp(2)%qmax
        hpp_s(2)%q2in = hpp(2)%qmax
      else
        hpp_s(3)%q2in = 0.0
      endif

      !Kurobe 3
      if (hpp_s(3)%q2in > hpp(3)%qmax) then
        !dam_s(idam)%q2sout = hpp_s(3)%q2in - hpp(3)%qmax
        hpp_s(3)%q2in = hpp(3)%qmax
      endif
      
      ! Sennindani Dam = NK3 + K3
      dam_s(idam)%q2hout = hpp_s(2)%q2in + hpp_s(3)%q2in
     
      !Shin-Kurobe 2 ! This is not from any dam, directly connected to NK3
      hpp_s(4)%q2in = hpp_s(2)%q2in !NK3 --> NK2 --> River
      !--------------------------------------------------------------
     
      
      case(3) !----------------Koyadaira Dam------------------------
      ! Kurobe 2
      hpp_s(5)%q2in = dam_s(idam)%q2in - dam_s(idam)%q2out

      if (hpp_s(5)%q2in > hpp(5)%qmax) then
        !dam_s(idam)%q2sout = hpp_s(5)%q2in - hpp(5)%qmax
        hpp_s(5)%q2in = hpp(5)%qmax
      endif
    
      dam_s(idam)%q2hout = hpp_s(5)%q2in
      !--------------------------------------------------------------
      
      
      case(4) !----------------Dashidaira Dam------------------------
      ! Release in made in order of Electric Water Ratio
      !Otozawa
      hpp_s(8)%q2in = dam_s(idam)%q2in - dam_s(idam)%q2out
      
      if (hpp_s(8)%q2in > hpp(8)%qmax) then
        hpp_s(7)%q2in = hpp_s(8)%q2in - hpp(8)%qmax
        hpp_s(8)%q2in = hpp(8)%qmax
      else
        hpp_s(7)%q2in = 0
        hpp_s(6)%q2in = 0
      endif

      !Shin-Yana
      if (hpp_s(7)%q2in > hpp(7)%qmax) then
        hpp_s(6)%q2in = hpp_s(7)%q2in - hpp(7)%qmax
        hpp_s(7)%q2in = hpp(7)%qmax
      else
        hpp_s(6)%q2in = 0
      endif

      !Dashidaira
      if (hpp_s(6)%q2in > hpp(6)%qmax) then
        !dam_s(idam)%q2sout = hpp_s(6)%q2in - hpp(6)%qmax
        hpp_s(6)%q2in = hpp(6)%qmax
      endif  
      
      dam_s(idam)%q2hout = hpp_s(6)%q2in +
     &                     hpp_s(7)%q2in +
     &                     hpp_s(8)%q2in
      !---------------------------------------------------------------
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      end select

      !++++++++++++++++++++++Total Outflow++++++++++++++++++++++++++++
      dam_s(idam)%q2tout = dam_s(idam)%q2out + 
     &                     dam_s(idam)%q2sout + 
     &                     dam_s(idam)%q2hout
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++       
      case(3)!======Outflow Override (Exact Data)=====================
      !++++++++++++++++++++++++++++Inflow++++++++++++++++++++++++++++
      dam_s(idam)%q2in = q2
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ratio=0
      ! ++++++++++++++++++++++++++++Outflow++++++++++++++++++++++++++
      dam_s(idam)%q2out = ratio*
     &              time_lookup_r4(dam(idam)%tvar%outflow,
     &                             year,month,day,hour,minute,second)
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      !+++++++++++++++++++++++++++++Spill++++++++++++++++++++++++++++
      dam_s(idam)%q2sout = ratio*
     &              time_lookup_r4(dam(idam)%tvar%spills,
     &                             year,month,day,hour,minute,second)

      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      
      ! +++++++++++++++++++++++Hydropower Outflow+++++++++++++++++++
      ! Calculation of hydropower outflow is different for each dam
      ! At this moment connectivity to hydropower plant is hardcoded
      ! Need to find a way to generalize this
      
      select case(idam)
      
      case(1) !-------------------Kurobe Dam------------------------
      !Kurobe 4
      hpp_s(1)%q2in = ratio*
     &              time_lookup_r4(hpp(1)%tvar%inflow,
     &                             year,month,day,hour,minute,second)
      hpp_s(1)%q2in = min(hpp_s(1)%q2in,hpp(1)%qmax)
      dam_s(idam)%q2hout = hpp_s(1)%q2in
      !-------------------------------------------------------------
      
      
      case(2) !-----------------Sennindani Dam----------------------
      !Shin-Kurobe 3
      hpp_s(2)%q2in = ratio*
     &              time_lookup_r4(hpp(2)%tvar%inflow,
     &                             year,month,day,hour,minute,second)
      hpp_s(2)%q2in = min(hpp_s(2)%q2in,hpp(2)%qmax)
 
      !Kurobe 3
      hpp_s(3)%q2in = ratio*
     &              time_lookup_r4(hpp(3)%tvar%inflow,
     &                             year,month,day,hour,minute,second)
      hpp_s(3)%q2in = min(hpp_s(3)%q2in,hpp(3)%qmax)
      
      ! Sennindani Dam = NK3 + K3
      dam_s(idam)%q2hout = hpp_s(2)%q2in + hpp_s(3)%q2in
     
      !Shin-Kurobe 2 ! This is not from any dam, directly connected to NK3
      hpp_s(4)%q2in = hpp_s(2)%q2in !NK3 --> NK2 --> River
      !--------------------------------------------------------------
     
      
      case(3) !----------------Koyadaira Dam------------------------
      ! Kurobe 2
      hpp_s(5)%q2in = ratio*
     &              time_lookup_r4(hpp(5)%tvar%inflow,
     &                             year,month,day,hour,minute,second)
      hpp_s(5)%q2in = min(hpp_s(5)%q2in,hpp(5)%qmax)
      
      dam_s(idam)%q2hout = hpp_s(5)%q2in
      !--------------------------------------------------------------
      
      
      case(4) !----------------Dashidaira Dam------------------------
      ! Dashidaira
      hpp_s(6)%q2in = ratio*
     &              time_lookup_r4(hpp(6)%tvar%inflow,
     &                             year,month,day,hour,minute,second)
      hpp_s(6)%q2in = min(hpp_s(6)%q2in,hpp(6)%qmax)
     
      ! Shin-Yana
      hpp_s(7)%q2in = ratio*
     &              time_lookup_r4(hpp(7)%tvar%inflow,
     &                             year,month,day,hour,minute,second)
      hpp_s(7)%q2in = min(hpp_s(7)%q2in,hpp(7)%qmax)     
     
      ! Otozawa
      hpp_s(8)%q2in = ratio*
     &              time_lookup_r4(hpp(8)%tvar%inflow,
     &                             year,month,day,hour,minute,second)
      hpp_s(8)%q2in = min(hpp_s(8)%q2in,hpp(8)%qmax)     
      
      dam_s(idam)%q2hout = hpp_s(6)%q2in +
     &                     hpp_s(7)%q2in +
     &                     hpp_s(8)%q2in
      !---------------------------------------------------------------
      case(5) !----------------Unazuki Dam------------------------
      ! Unazuki
      hpp_s(9)%q2in = ratio*
     &              time_lookup_r4(hpp(9)%tvar%inflow,
     &                             year,month,day,hour,minute,second)
      hpp_s(9)%q2in = min(hpp_s(9)%q2in,hpp(9)%qmax)

      ! Aimoto
      hpp_s(10)%q2in = ratio*
     &              time_lookup_r4(hpp(10)%tvar%inflow,
     &                             year,month,day,hour,minute,second)
      hpp_s(10)%q2in = min(hpp_s(10)%q2in,hpp(10)%qmax)      

      dam_s(idam)%q2hout = hpp_s(9)%q2in
      !--------------------------------------------------------------
      case(6) !----------------Kitamata Dam------------------------
      dam_s(idam)%q2hout = 0
      case(7) !----------------Futami Dam------------------------
      dam_s(idam)%q2hout = 0
      !--------------------------------------------------------------

      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      end select
      
      !++++++++++++++++++++++Total Outflow++++++++++++++++++++++++++++
      dam_s(idam)%q2tout = dam_s(idam)%q2out + 
     &                     dam_s(idam)%q2sout + 
     &                     dam_s(idam)%q2hout
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++       
      end select
      
      if (use_damops /= 3) then
        call update_dam_state(dam(idam),dam_s(idam),use_damops)
      endif
      
      q2 = dam_s(idam)%q2sout + dam_s(idam)%q2out

    !   write(*,*) dam(idam)%damname, 
    !  &           dam_s(idam)%q2out, 
    !  &           dam_s(idam)%q2sout,
    !  &           dam_s(idam)%q2hout
     
      end subroutine operate_dam
      !=====================================================================
     
      !=====================================================================
      subroutine operate_hpp(ihpp,q2)
      !
      ! Purpose:
      ! --------
      !   This is the main subroutine with operates the hpp and updates
      !   all the necessary variables
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 23/11/2018          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------
      implicit none
      include 'globcst.inc'
      include 'phycst.inc'
      
      integer, intent(in)   :: ihpp
      real(r4),intent(inout):: q2  ! Only for modification to the river model, not used to calculate any state variables for hpp
      real(r8)              ::wl_tmp

      if (hpp(ihpp)%use_ph) then
        select case(hpp(ihpp)%hppname)

        case('kurobe4')
          wl_tmp = 0.5*(dam_s(1)%wl1+dam_s(1)%wl2)
          hpp(ihpp)%rate = wl_to_rate(wl_tmp,hpp(ihpp)%ph)
          !write(*,*) wl_tmp,hpp(ihpp)%rate
        end select
      endif



      hpp_s(ihpp)%q2out = hpp_s(ihpp)%q2in

      if (hpp(ihpp)%rate /= -9999) then
        hpp_s(ihpp)%power = hpp_s(ihpp)%q2out*hpp(ihpp)%rate! MW
        hpp_s(ihpp)%energy = hpp_s(ihpp)%power*1E03              ! kWh
      else
        hpp_s(ihpp)%power = g*rhow*hpp(ihpp)%efficiency*1E-06*
     &                      hpp_s(ihpp)%q2out*hpp_s(ihpp)%nethead2

         hpp_s(ihpp)%energy = g*hpp(ihpp)%efficiency*1E-06*
     &                      hpp_s(ihpp)%q2out*hpp_s(ihpp)%nethead2*
     &                      dthydro/3600.0
      endif


!      hpp_s(ihpp)%power = g*rhow*hpp(ihpp)%efficiency*1E-06*
!     &                0.5*(hpp_s(ihpp)%q2out+hpp_s(ihpp)%q1out)*
!     &                0.5*(hpp_s(ihpp)%nethead2+hpp_s(ihpp)%nethead1)
         
!      hpp_s(ihpp)%energy = g*hpp(ihpp)%efficiency*
!     &                0.5*(hpp_s(ihpp)%q2out+hpp_s(ihpp)%q1out)*
!     &                0.5*(hpp_s(ihpp)%nethead2+hpp_s(ihpp)%nethead1)*
!     &                dthydro/3600.0
      
        if ((hpp(ihpp)%subcatch /= -9999) .and. ! Only if the outlet is within the model domain
     &     (hpp(ihpp)%flowint /= -9999)) then
            if (ihpp /= 2) then ! NK3 does not go to the river Moiz: Temporary
              !write(*,*) hpp(ihpp)%hppname
              !write(*,*) 'Before: ', q2
              if (ihpp == 9) then
                q2 = q2 + max(hpp_s(9)%q2out - hpp_s(10)%q2out,0.0)
              else
                q2 = q2 + hpp_s(ihpp)%q2out
              endif
              !write(*,*) 'After: ', q2
            endif
      endif
      

      end subroutine operate_hpp
      !=====================================================================
      
      
      !=====================================================================
      subroutine update_dam_state(dam,dam_s,use_damops)
      !
      ! Purpose:
      ! --------
      !   Updates the state (v,wl) of the dam based on simulated inflow &
      !   decided outflow using the follwoing equation:
      !
      !   (q1tout + q2tout)*0.5 - (q1in + q2in)*0.5 = (v2 - v1)/dt
      !
      !   We need to find v2.
      !   Known from previous timestep: q1tout,q1in,v1
      !   Known from current timestep: q2in
      !   Decided in current timestep: q2tout
      !
      !   Initializes the state variables of dam and hpp
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 22/11/2018          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------
      implicit none
      include 'globcst.inc'
      
      type(dam_data),   intent(in) :: dam
      type(dam_state),  intent(inout) :: dam_s
      integer,          intent(in)    :: use_damops
      
      ! Local Variable
      real(r4) :: diff
      real(r4) :: spill
      real(r4) :: deficit
      real(r8) :: total_hydro
      real(r8) :: lwl        ! Local Low Water Level
      real(r8) :: hwl        ! Local High Water Level


      dam_s%v2 = dam_s%v1 + 
     &           dthydro*(0.5*(dam_s%q1in+dam_s%q2in) - 
     &                    0.5*(dam_s%q1tout+dam_s%q2tout)) 

      if (dam%use_fixed_wl) then
        lwl = dam%fixed_wl
        hwl = dam%fixed_wl
      else
        lwl = dam%lwl
        hwl = dam%hwl
      endif

      if (dam_s%v2 > el_to_store(hwl,dam%vh)) then
        spill = (dam_s%q1in+dam_s%q2in) 
     &       - 2*((el_to_store(hwl,dam%vh) - dam_s%v1)/dthydro) 
     &       - dam_s%q1tout
     &       - (dam_s%q2hout + dam_s%q2out)
        
        if (spill < 0) then
          dam_s%v2 = dam_s%v1 + 
     &           dthydro*(0.5*(dam_s%q1in+dam_s%q2in) - 
     &           0.5*(dam_s%q1tout+(dam_s%q2hout + dam_s%q2out)))
          spill = 0
        else
          dam_s%v2 = el_to_store(hwl,dam%vh)
        endif

        dam_s%q2sout = dam_s%q2sout + spill
        dam_s%q2tout = dam_s%q2tout + spill
      
      elseif (dam_s%v2 < el_to_store(lwl,dam%vh)) then

        deficit = dam_s%q2tout - ((dam_s%q1in+dam_s%q2in) 
     &       - 2*((el_to_store(lwl,dam%vh) - dam_s%v1)/dthydro) 
     &       - dam_s%q1tout)

        dam_s%q2tout = (dam_s%q1in+dam_s%q2in) 
     &       - 2*((el_to_store(lwl,dam%vh) - dam_s%v1)/dthydro) 
     &       - dam_s%q1tout
        

        ! Priority 1: Reduce Hydropower
        select case(dam%damname)
          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          case('kurobe') !-------------Kurobe Dam----------------------
            ! Kurobe 4
            if (deficit>0) then
              diff = hpp_s(1)%q2in - deficit
              if (diff>=0) then
                deficit = 0
                hpp_s(1)%q2in = diff
              else
                deficit = deficit - hpp_s(1)%q2in
                hpp_s(1)%q2in = 0
              endif
            elseif (deficit<0) then
              hpp_s(1)%q2in = 0
            endif

            dam_s%q2hout = hpp_s(1)%q2in
          !------------------------------------------------------------
          case('sennindani') !----------Sennindani Dam-----------------
            ! Kurobe 3
            if (deficit>0) then
              diff = hpp_s(3)%q2in - deficit
              if (diff>=0) then
                deficit = 0
                hpp_s(3)%q2in = diff
              else
                deficit = deficit - hpp_s(3)%q2in
                hpp_s(3)%q2in = 0
              endif
            elseif (deficit<0) then
              hpp_s(3)%q2in = 0
            endif

            ! Shin-Kurobe 3
            if (deficit>0) then
              diff = hpp_s(2)%q2in - deficit
              if (diff>=0) then
                deficit = 0
                hpp_s(2)%q2in = diff
              else
                deficit = deficit - hpp_s(2)%q2in
                hpp_s(2)%q2in = 0
              endif
            elseif (deficit<0) then
              hpp_s(2)%q2in = 0
            endif

            dam_s%q2hout = hpp_s(2)%q2in + hpp_s(3)%q2in
            hpp_s(4)%q2in = hpp_s(2)%q2in !NK3 --> NK2 --> river
            !------------------------------------------------------------
          case('koyadaira') !----------Koyadaira Dam-----------------
            ! Kurobe 2
            if (deficit>0) then
              diff = hpp_s(5)%q2in - deficit
              if (diff>=0) then
                deficit = 0
                hpp_s(5)%q2in = diff
              else
                deficit = deficit - hpp_s(5)%q2in
                hpp_s(5)%q2in = 0
              endif
            elseif (deficit<0) then
              hpp_s(5)%q2in = 0
            endif

            dam_s%q2hout = hpp_s(5)%q2in
            !------------------------------------------------------------

          case('dashidaira') !----------Dashidaira Dam-----------------
            ! Dashdaira
            if (deficit>0) then
              diff = hpp_s(6)%q2in - deficit
              if (diff>=0) then
                deficit = 0
                hpp_s(6)%q2in = diff
              else
                deficit = deficit - hpp_s(6)%q2in
                hpp_s(6)%q2in = 0
              endif
            elseif (deficit<0) then
              hpp_s(6)%q2in = 0
            endif

            ! Shin-Yana
            if (deficit>0) then
              diff = hpp_s(7)%q2in - deficit
              if (diff>=0) then
                deficit = 0
                hpp_s(7)%q2in = diff
              else
                deficit = deficit - hpp_s(7)%q2in
                hpp_s(7)%q2in = 0
              endif
            elseif (deficit<0) then
              hpp_s(7)%q2in = 0
            endif

            ! Otozawa
            if (deficit>0) then
              diff = hpp_s(8)%q2in - deficit
              if (diff>=0) then
                deficit = 0
                hpp_s(8)%q2in = diff
              else
                deficit = deficit - hpp_s(8)%q2in
                hpp_s(8)%q2in = 0
              endif
            elseif (deficit<0) then
              hpp_s(8)%q2in = 0
            endif

            dam_s%q2hout = hpp_s(6)%q2in + 
     &                     hpp_s(7)%q2in +
     &                     hpp_s(8)%q2in
            !------------------------------------------------------------

        end select

        ! Priority 2: Reduce Outflow
        if (deficit>0) then
          diff = dam_s%q2out - deficit
          if (diff>=0) then
            deficit = 0
            dam_s%q2out = diff
          else
            deficit = deficit - dam_s%q2out
            dam_s%q2out = 0
          endif
        elseif (deficit<0) then
            dam_s%q2out = 0
        endif

          
         
         !++++++++++++++++++++++Total Outflow+++++++++++++++++++++++++
         dam_s%q2tout = dam_s%q2out + 
     &                  dam_s%q2sout + 
     &                  dam_s%q2hout
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         dam_s%v2 = el_to_store(lwl,dam%vh)
      endif
      
        dam_s%wl2 = store_to_wl(dam_s%v2,dam%vh)
      end subroutine update_dam_state
      !=====================================================================
      
      
      !=====================================================================
      subroutine prep_dam(idam)
      !
      ! Purpose:
      ! --------
      !   Initializes the 'previous timestep' for calculations in the
      !   next timestep.
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 23/11/2018          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------
      implicit none
      include 'globcst.inc'
      
      integer, intent(in)   :: idam
      
      dam_s(idam)%q1in = dam_s(idam)%q2in
      dam_s(idam)%q1tout = dam_s(idam)%q2tout
      dam_s(idam)%q1out = dam_s(idam)%q2out
      dam_s(idam)%q1sout = dam_s(idam)%q2sout
      dam_s(idam)%q1hout = dam_s(idam)%q2hout
      dam_s(idam)%wl1 = dam_s(idam)%wl2
      dam_s(idam)%v1 = dam_s(idam)%v2
      
      end subroutine prep_dam
      !=====================================================================
      
      !=====================================================================
      subroutine prep_hpp(ihpp)
      !
      ! Purpose:
      ! --------
      !   Initializes the 'previous timestep' for calculations in the
      !   next timestep.
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 23/11/2018          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------
      implicit none
      include 'globcst.inc'
      
      integer, intent(in)   :: ihpp
      
      hpp_s(ihpp)%q1in = hpp_s(ihpp)%q2in
      hpp_s(ihpp)%q1out = hpp_s(ihpp)%q2out
      
      end subroutine prep_hpp
      !=====================================================================
      
      

      
      
      !=====================================================================
      subroutine init_damops(use_damops,dam_snapshot)
      !
      ! Purpose:
      ! --------
      !   Reads in the names, properties and data of dams using the following
      !   subroutines in this module:
      !   1. read_dam
      !   2. read_hpp
      !
      !   Initializes the state variables of dam and hpp
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 22/11/2018          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------    
      implicit none
      include 'globcst.inc'
      include 'hydro.inc'
      
      ! Previous time step
      integer, intent(in) :: use_damops
      integer, intent(in) :: dam_snapshot
      integer :: pyear
      integer :: pmonth
      integer :: pday
      integer :: phour
      integer :: pminute
      integer :: psecond
      
      
      integer           :: i                        ! Iterator
      
      integer           :: lu                       ! I/O Unit
      integer           :: io_status                ! I/O Status
      character(len=200):: io_message               ! I/O Message

      if (use_damops == 0) return

      call read_dam(dam_dir,use_damops)
      call read_hpp(dam_dir,use_damops)
      
      ! Previous timestep for hydropower dam variable initialization
      call abss2ctim(abstini-int(dthydro),pyear,pmonth,pday,
     &                       phour,pminute,psecond)
      
      if (dam_snapshot == 1) then
        allocate(hpp_s(size(hpp)),stat=io_status, errmsg=io_message)
        allocate(dam_s(size(dam)),stat=io_status, errmsg=io_message)
        return
      endif
      

      ! ----------Allocate & Initialize HPP State Variables-------------
      allocate(hpp_s(size(hpp)),stat=io_status, errmsg=io_message)
      if (io_status /= 0) then
          write(*,*) io_message
      else
          do i=1,size(hpp_s)
          
            ! Power & Energy at Current Timestep
            ! These variables are assumed to be discrete
            hpp_s(i)%power = 0.0
            hpp_s(i)%energy = 0.0
            
            ! Nethead
            hpp_s(i)%nethead1 = hpp(i)%nethead ! It is assumed that head doesn't change with time
            hpp_s(i)%nethead2 = hpp(i)%nethead ! head due to headrace is much higher than that in the channel
            
            ! Inflow
            hpp_s(i)%q1in = time_lookup_r4(hpp(i)%tvar%inflow,
     &                                    pyear,pmonth,pday,
     &                                    phour,pminute,psecond)
            hpp_s(i)%q2in = 0.0
            
            ! Outflow (Outflow == Inflow)
            hpp_s(i)%q1out = time_lookup_r4(hpp(i)%tvar%inflow,
     &                                    pyear,pmonth,pday,
     &                                    phour,pminute,psecond)
            hpp_s(i)%q2out = 0.0
          enddo
          

      endif
      !---------------------------------------------------------------
      
      
      ! ----------Allocate & Initialize Dam State Variables-----------
      allocate(dam_s(size(dam)),stat=io_status, errmsg=io_message)
      if (io_status /= 0) then
          write(*,*) io_message
      else
          do i=1,size(dam_s)
             
            ! Initializing
            
            ! Inflow
            dam_s(i)%q1in = time_lookup_r4(dam(i)%tvar%inflow,
     &                                    pyear,pmonth,pday,
     &                                    phour,pminute,psecond)
            dam_s(i)%q2in = 0.0
            
            ! Outflow
            dam_s(i)%q1out = time_lookup_r4(dam(i)%tvar%outflow,
     &                                      pyear,pmonth,pday,
     &                                      phour,pminute,psecond)
            dam_s(i)%q2out = 0.0
            
            ! Spillway
            dam_s(i)%q1sout = time_lookup_r4(dam(i)%tvar%spills,
     &                                      pyear,pmonth,pday,
     &                                      phour,pminute,psecond)
            dam_s(i)%q2sout = 0.0
     
            ! Hydropower
            ! Connection of Dam to Hydropower Plants has been hardcoded
            ! Need to find a way to generalize this
            select case(i)
            
            case(1) ! Kurobe Dam
             dam_s(i)%q1hout = hpp_s(1)%q1in   ! Kurobe 4
            
            case(2) ! Sennindani Dam
             dam_s(i)%q1hout = hpp_s(2)%q1in + ! Shin-Kurobe 3 ----> Shin_Kurobe 2
     &                         hpp_s(3)%q1in   ! Kurobe 3
            
            case(3) ! Koyadaira Dam
             dam_s(i)%q1hout = hpp_s(5)%q1in   ! Kurobe 2
             
            case(4) ! Dashidaira Dam
             dam_s(i)%q1hout = hpp_s(6)%q1in + ! Dashidaira
     &                         hpp_s(7)%q1in + ! Shin-Yana
     &                         hpp_s(8)%q1in   ! Otozawa
            case(5) ! Unazuki Dam
             dam_s(i)%q1hout = hpp_s(9)%q1in   ! Unazuki ----> Aimoto 
            
            case default ! Other Dams
             dam_s(i)%q1hout = 0.0
             write(*,*) 'Connectivity for '//trim(dam(i)%damname)//
     &                  ' has not been specified. Setting q1hout=0.0'
             
            end select
             dam_s(i)%q2hout = 0.0
            
            
            ! Total Outflow = Outflow + Spillway + Hydropower
            dam_s(i)%q1tout = dam_s(i)%q1out + 
     &                        dam_s(i)%q1sout + 
     &                        dam_s(i)%q1hout
            dam_s(i)%q2tout = 0.0
            
            
            ! Water Level
            if (dam(i)%override_wl == .false.) then
            dam_s(i)%wl1 = time_lookup_r8(dam(i)%tvar%wl,
     &                                    pyear,pmonth,pday,
     &                                    phour,pminute,psecond)
            else
            dam_s(i)%wl1 = 
     &         store_to_wl(el_to_store(dam(i)%init_wl_override,
     &         dam(i)%vh),dam(i)%vh)
            endif

            dam_s(i)%wl2 = 0.0
            
            ! Volume
            dam_s(i)%v1 = wl_to_store(dam_s(i)%wl1,dam(i)%vh)
            dam_s(i)%v2 = 0.0
       
          enddo
      endif
      
      ! If inicon == 1, then read initial conditions from snapshot
      !if (inicon == 1) then
      !    write(*,*) 'Reading dam and hpp states from snapshot'
      !    call damops_load_snapshot  ! Moiz: commented temp 2019/09/26
      !endif
      !--------------------------------------------------------------------
      
      end subroutine init_damops
      !=====================================================================
      
      !=====================================================================
      subroutine read_dam(data_folder,use_damops)
      !
      ! Purpose:
      ! --------
      !   Reads in the names, properties and data of dams
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 19/11/2018          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------        
        implicit none
        
        character(len=*), intent(in)  ::   data_folder    ! Data Folder
        integer, intent(in)           ::   use_damops     ! Damops Configuration

        integer           :: i,j                      ! Iterators
		character(len=200):: in_file                  ! Input Dam Data
        
        integer           :: lu                       ! I/O Unit
		integer           :: io_status                ! I/O Status
		character(len=200):: io_message               ! I/O Message
        
        
        in_file = trim(data_folder)//'dam_para'
        write(*,*) in_file
        
        ! Datetime Index
		allocate(dam(csv_len(in_file)-1), 
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
        
        open(file=in_file,newunit=lu,iostat=io_status,iomsg=io_message,
     &       status='OLD')
        if (io_status/=0) then
            write(*,*) io_message
		else
          if (allocated(dam)) then
            read(lu,'(X)')
            
            do i = 1, size(dam)
            read(lu,*) dam(i)%damname,
     &                 dam(i)%subcatch,
     &                 dam(i)%flowint,
     &                 dam(i)%h,
     &                 dam(i)%gross_store,
     &                 dam(i)%eff_store,
     &                 dam(i)%lwl,
     &                 dam(i)%pwl,
     &                 dam(i)%nwl,
     &                 dam(i)%hwl,
     &                 dam(i)%qmax_spill,
     &                 dam(i)%qmin_spill,
     &                 dam(i)%qmax_out,
     &                 dam(i)%qmin_out,
     &                 dam(i)%flood_def,
     &                 dam(i)%des_flood,
     &                 dam(i)%override_wl,
     &                 dam(i)%init_wl_override,
     &                 dam(i)%use_fixed_wl,
     &                 dam(i)%fixed_wl
           
            call read_dam_data(trim(data_folder),use_damops,dam(i))

            enddo
     
          else
            write(*,*) 'Warning: Array is not allocated!'
          endif
          close(lu)
          
          
          
          
        endif  
            
        
        
      end subroutine read_dam
      !=====================================================================
      

      
      !=====================================================================
	    subroutine read_dam_data(data_folder,use_damops,dam)
	    !
      ! Purpose:
      ! --------
      !   This subroutine reads in the inflow, outflow, spills, water levels
      !   Total outflow is calculated as : Total Outflow = Outflow + Spills
      ! 
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 17/11/2018          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------
      implicit none     
      character(len=*), intent(in)  :: data_folder   ! Folder with Dam Inputs
      integer,          intent(in)  :: use_damops    ! Damops Configuration
      type(dam_data), intent(inout) :: dam           ! Dam TYPE Variable
      
      integer           :: i,j                      ! Iterators
      character(len=200):: in_file                  ! Input Dam Data
      character(len=30) :: date_format              ! Date Format YYYY-MM-DD HH:MM:SS

      integer           :: lu                       ! I/O Unit
      integer           :: io_status                ! I/O Status
      character(len=200):: io_message               ! I/O Message

      integer, allocatable, dimension(:)  :: year               ! Year [YYYY]
      integer  :: iyear,imonth,iday,ihour,imin,isec
      real(r4) :: r4_var
      real(r8) :: r8_var

          
        select case(use_damops)
        
        case(1) !------------------Observed Data-----------------------
        ! Input File
        in_file = trim(data_folder)//trim(dam%damname)//'_dam_in.csv'
        
        ! Reading Datetime
        allocate(year(csv_len(in_file)),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
        date_format = '(i4)'
        open(file=in_file,newunit=lu,iostat=io_status,iomsg=io_message,
     &       status='OLD')
        if (io_status/=0) then
            write(*,*) io_message
        else
          if (allocated(year)) then
            read(lu,date_format)
     &      ((year(i)),i=1,csv_len(in_file))
            close(lu)
          else
            write(*,*) 'Warning: Array is not allocated!'
          endif
        endif
	      
        ! Allocating Arrays
        
        ! Dam Inflow
        allocate(dam%tvar%inflow(year(1):year(size(year)),
     &                           1:12,
     &                           1:31,
     &                           0:23,
     &                           0:0,
     &                           0:0                ),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
        
        ! Dam Outflow
        allocate(dam%tvar%outflow(year(1):year(size(year)),
     &                           1:12,
     &                           1:31,
     &                           0:23,
     &                           0:0,
     &                           0:0                ),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
  
        
        ! Dam Spills
        allocate(dam%tvar%spills(year(1):year(size(year)),
     &                           1:12,
     &                           1:31,
     &                           0:23,
     &                           0:0,
     &                           0:0                ),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
        
        ! Total Outflows
        allocate(dam%tvar%t_outflow(year(1):year(size(year)),
     &                           1:12,
     &                           1:31,
     &                           0:23,
     &                           0:0,
     &                           0:0                ),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
        
        ! Water Level
        allocate(dam%tvar%wl(year(1):year(size(year)),
     &                           1:12,
     &                           1:31,
     &                           0:23,
     &                           0:0,
     &                           0:0                ),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
        
        
        ! VH Curve
        in_file = trim(data_folder)//trim(dam%damname)//'_dam_vh.csv'
        allocate(dam%vh(csv_len(in_file)-1,3),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
        
        
        
        ! Initializing all variables
        dam%tvar%inflow = 0.0
        dam%tvar%outflow = 0.0
        dam%tvar%spills = 0.0
        dam%tvar%t_outflow = 0.0
        dam%tvar%wl = 0.0
        dam%vh = 0.0

        
        ! Reading Inflow
        in_file = trim(data_folder)//trim(dam%damname)//'_dam_in.csv'
        open(file=in_file,newunit=lu,iostat=io_status,iomsg=io_message,
     &       status='OLD')
        if (io_status/=0) then
            write(*,*) io_message
		    else
          if (allocated(dam%tvar%inflow)) then
            write(*,*) 'Reading '//trim(in_file)
            do i = 1, csv_len(in_file)
               read(lu,'(i4,x,i2,x,i2,x,i2,x,i2,x,i2,x,f10.3)')
     &         iyear,imonth,iday,ihour,imin,isec,r4_var
               dam%tvar%inflow(iyear,imonth,iday,ihour,0,0) = 
     &         r4_var
            enddo
          else
            write(*,*) 'Warning: Array is not allocated!'
          endif
          close(lu)
        endif

        
        ! Reading Outflow
        in_file = trim(data_folder)//trim(dam%damname)//'_dam_out.csv'
        open(file=in_file,newunit=lu,iostat=io_status,iomsg=io_message,
     &       status='OLD')
        if (io_status/=0) then
            write(*,*) io_message
		    else
          if (allocated(dam%tvar%outflow)) then
            write(*,*) 'Reading '//trim(in_file)
            do i = 1, csv_len(in_file)
               read(lu,'(i4,x,i2,x,i2,x,i2,x,i2,x,i2,x,f10.3)')
     &         iyear,imonth,iday,ihour,imin,isec,r4_var
               dam%tvar%outflow(iyear,imonth,iday,ihour,0,0) = 
     &         r4_var
            enddo
          else
            write(*,*) 'Warning: Array is not allocated!'
          endif
          close(lu)
        endif
        
        
        ! Reading Spills
        in_file = trim(data_folder)//trim(dam%damname)//'_dam_spill.csv'
        open(file=in_file,newunit=lu,iostat=io_status,iomsg=io_message,
     &       status='OLD')
        if (io_status/=0) then
            write(*,*) io_message
		    else
          if (allocated(dam%tvar%spills)) then
            write(*,*) 'Reading '//trim(in_file)
            do i = 1, csv_len(in_file)
               read(lu,'(i4,x,i2,x,i2,x,i2,x,i2,x,i2,x,f10.3)')
     &         iyear,imonth,iday,ihour,imin,isec,r4_var
               dam%tvar%spills(iyear,imonth,iday,ihour,0,0) = 
     &         r4_var
            enddo
          else
            write(*,*) 'Warning: Array is not allocated!'
          endif
          close(lu)
        endif
        
        
        ! Reading Water Levels
        in_file = trim(data_folder)//trim(dam%damname)//'_dam_wl.csv'
        open(file=in_file,newunit=lu,iostat=io_status,iomsg=io_message,
     &       status='OLD')
        if (io_status/=0) then
            write(*,*) io_message
		    else
          if (allocated(dam%tvar%wl)) then
            write(*,*) 'Reading '//trim(in_file)
            do i = 1, csv_len(in_file)
               read(lu,'(i4,x,i2,x,i2,x,i2,x,i2,x,i2,x,f10.3)')
     &         iyear,imonth,iday,ihour,imin,isec,r8_var
               dam%tvar%wl(iyear,imonth,iday,ihour,0,0) = 
     &         r8_var
            enddo
          else
            write(*,*) 'Warning: Array is not allocated!'
          endif
          close(lu)
        endif
        
        ! Reading VH Curves
        in_file = trim(data_folder)//trim(dam%damname)//'_dam_vh.csv'
        open(file=in_file,newunit=lu,iostat=io_status,iomsg=io_message,
     &       status='OLD')
        if (io_status/=0) then
            write(*,*) io_message
		    else
          if (allocated(dam%vh)) then
            write(*,*) 'Reading '//trim(in_file)
            read(lu,'(X)') 
            read(lu,*)
     &      ((dam%vh(i,j), j=1,3), i = 1,(csv_len(in_file)-1))
          else
            write(*,*) 'Warning: Array is not allocated!'
          endif
          close(lu)
        endif

      case(2) !------------------Forecast Simulation-----------------------
        !Remove extra inputs if not needed
        ! Input File
        in_file = trim(data_folder)//trim(dam%damname)//
     &   '_dam_wl.csv'
        
        ! Reading Datetime
        allocate(year(csv_len(in_file)),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
        date_format = '(i4)'
        open(file=in_file,newunit=lu,iostat=io_status,iomsg=io_message,
     &       status='OLD')
        if (io_status/=0) then
            write(*,*) io_message
        else
          if (allocated(year)) then
            read(lu,date_format)
     &      ((year(i)),i=1,csv_len(in_file))
            close(lu)
          else
            write(*,*) 'Warning: Array is not allocated!'
          endif
        endif
	      
        ! Allocating Arrays
        
        ! Dam Inflow
        allocate(dam%tvar%inflow(year(1):year(size(year)),
     &                           1:12,
     &                           1:31,
     &                           0:23,
     &                           0:0,
     &                           0:0                ),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
        
        ! Dam Outflow
        allocate(dam%tvar%outflow(year(1):year(size(year)),
     &                           1:12,
     &                           1:31,
     &                           0:23,
     &                           0:0,
     &                           0:0                ),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
  
        
        ! Dam Spills
        allocate(dam%tvar%spills(year(1):year(size(year)),
     &                           1:12,
     &                           1:31,
     &                           0:23,
     &                           0:0,
     &                           0:0                ),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
        
        ! Total Outflows
        allocate(dam%tvar%t_outflow(year(1):year(size(year)),
     &                           1:12,
     &                           1:31,
     &                           0:23,
     &                           0:0,
     &                           0:0                ),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
        
        ! Water Level
        allocate(dam%tvar%wl(year(1):year(size(year)),
     &                           1:12,
     &                           1:31,
     &                           0:23,
     &                           0:0,
     &                           0:0                ),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
        
        
        ! VH Curve
        in_file = trim(data_folder)//trim(dam%damname)//'_dam_vh.csv'
        allocate(dam%vh(csv_len(in_file)-1,3),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
        
        ! Initializing all variables
        dam%tvar%inflow = 0.0
        dam%tvar%outflow = 0.0
        dam%tvar%spills = 0.0
        dam%tvar%t_outflow = 0.0
        dam%tvar%wl = 0.0
        dam%vh = 0.0
        
        ! Reading Water Levels
        in_file = trim(data_folder)//trim(dam%damname)//'_dam_wl.csv'
        open(file=in_file,newunit=lu,iostat=io_status,iomsg=io_message,
     &       status='OLD')
        if (io_status/=0) then
            write(*,*) io_message
		    else
          if (allocated(dam%tvar%wl)) then
            write(*,*) 'Reading '//trim(in_file)
            do i = 1, csv_len(in_file)
               read(lu,'(i4,x,i2,x,i2,x,i2,x,i2,x,i2,x,f10.3)')
     &         iyear,imonth,iday,ihour,imin,isec,r8_var
               dam%tvar%wl(iyear,imonth,iday,ihour,0,0) = 
     &         r8_var
            enddo
          else
            write(*,*) 'Warning: Array is not allocated!'
          endif
          close(lu)
        endif
        
        ! Reading VH Curves
        in_file = trim(data_folder)//trim(dam%damname)//'_dam_vh.csv'
        open(file=in_file,newunit=lu,iostat=io_status,iomsg=io_message,
     &       status='OLD')
        if (io_status/=0) then
            write(*,*) io_message
		    else
          if (allocated(dam%vh)) then
            write(*,*) 'Reading '//trim(in_file)
            read(lu,'(X)') 
            read(lu,*)
     &      ((dam%vh(i,j), j=1,3), i = 1,(csv_len(in_file)-1))
          else
            write(*,*) 'Warning: Array is not allocated!'
          endif
          close(lu)
        endif
      case(3) !------------------Outflow Override (Exact Data)-----------------------
        ! Input File
        in_file = trim(data_folder)//trim(dam%damname)//'_dam_in.csv'
        
        ! Reading Datetime
        allocate(year(csv_len(in_file)),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
        date_format = '(i4)'
        open(file=in_file,newunit=lu,iostat=io_status,iomsg=io_message,
     &       status='OLD')
        if (io_status/=0) then
            write(*,*) io_message
        else
          if (allocated(year)) then
            read(lu,date_format)
     &      ((year(i)),i=1,csv_len(in_file))
            close(lu)
          else
            write(*,*) 'Warning: Array is not allocated!'
          endif
        endif
	      
        ! Allocating Arrays
        
        ! Dam Inflow
        allocate(dam%tvar%inflow(year(1):year(size(year)),
     &                           1:12,
     &                           1:31,
     &                           0:23,
     &                           0:0,
     &                           0:0                ),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
        
        ! Dam Outflow
        allocate(dam%tvar%outflow(year(1):year(size(year)),
     &                           1:12,
     &                           1:31,
     &                           0:23,
     &                           0:0,
     &                           0:0                ),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
  
        
        ! Dam Spills
        allocate(dam%tvar%spills(year(1):year(size(year)),
     &                           1:12,
     &                           1:31,
     &                           0:23,
     &                           0:0,
     &                           0:0                ),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
        
        ! Total Outflows
        allocate(dam%tvar%t_outflow(year(1):year(size(year)),
     &                           1:12,
     &                           1:31,
     &                           0:23,
     &                           0:0,
     &                           0:0                ),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
        
        ! Water Level
        allocate(dam%tvar%wl(year(1):year(size(year)),
     &                           1:12,
     &                           1:31,
     &                           0:23,
     &                           0:0,
     &                           0:0                ),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
        
        
        ! VH Curve
        in_file = trim(data_folder)//trim(dam%damname)//'_dam_vh.csv'
        allocate(dam%vh(csv_len(in_file)-1,3),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
        
        
        
        ! Initializing all variables
        dam%tvar%inflow = 0.0
        dam%tvar%outflow = 0.0
        dam%tvar%spills = 0.0
        dam%tvar%t_outflow = 0.0
        dam%tvar%wl = 0.0
        dam%vh = 0.0

        
      end select
	  
	    end subroutine read_dam_data
      !====================================================================

      !=====================================================================
      subroutine read_hpp(data_folder,use_damops)
      !
      ! Purpose:
      ! --------
      !   Reads in the names, properties and data of hydropower plants
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 19/11/2018          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------        
        implicit none
        include 'phycst.inc'
        
        character(len=*), intent(in)  ::   data_folder    ! Data Folder
        integer,          intent(in)  ::   use_damops     ! Damops Configuration

        integer           :: i,j                      ! Iterators
		character(len=200):: in_file                  ! Input HPP Data
        
        integer           :: lu                       ! I/O Unit
		integer           :: io_status                ! I/O Status
		character(len=200):: io_message               ! I/O Message
        
        
        in_file = trim(data_folder)//'hpp_para'
        write(*,*) in_file
        
        ! Datetime Index
		allocate(hpp(csv_len(in_file)-1), 
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
        
        open(file=in_file,newunit=lu,iostat=io_status,iomsg=io_message,
     &       status='OLD')
        if (io_status/=0) then
            write(*,*) io_message
		else
          if (allocated(hpp)) then
            read(lu,'(X)')
            
            do i = 1, size(hpp)
            read(lu,*) hpp(i)%hppname,
     &                 hpp(i)%subcatch,
     &                 hpp(i)%flowint,
     &                 hpp(i)%nethead,
     &                 hpp(i)%qmax,
     &                 hpp(i)%qmin,
     &                 hpp(i)%installed_capacity,
     &                 hpp(i)%rate,
     &                 hpp(i)%use_ph
     
            
            hpp(i)%efficiency = (hpp(i)%installed_capacity*1E06)/
     &                          (rhow*g*hpp(i)%qmax*hpp(i)%nethead)
            
            call read_hpp_data(trim(data_folder),use_damops,hpp(i))
            
            enddo
     
          else
            write(*,*) 'Warning: Array is not allocated!'
          endif
          close(lu)
        endif         

      end subroutine read_hpp
      !=====================================================================
	  
      !=====================================================================
	  subroutine read_hpp_data(data_folder,use_damops,hpp)
	  !
      ! Purpose:
      ! --------
      !   This subroutine reads in the inflow for each Hydropower Plant
      !   Total outflow is calculated as : Total Outflow = Outflow + Spills
      ! 
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 17/11/2018          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------
	     implicit none
		   character(len=*), intent(in)  :: data_folder   ! Folder with HPP Inputs
       integer,          intent(in)  :: use_damops    ! Damops Configuration
       type(hpp_data), intent(inout) :: hpp           ! hpp TYPE Variable
		
		   integer           :: i,j                      ! Iterators
		   character(len=200):: in_file                  ! Input HPP Data
       character(len=30) :: date_format              ! Date Format YYYY-MM-DD HH:MM:SS

		   integer           :: lu                       ! I/O Unit
		   integer           :: io_status                ! I/O Status
       character(len=200):: io_message               ! I/O Message
       
       integer, allocatable, dimension(:)  :: year               ! Year [YYYY]
       integer  :: iyear,imonth,iday,ihour,imin,isec

       real(r4) :: r4_var
       real(r8) :: r8_var
       
       select case(use_damops)


       case(1)!---------------Observed Data----------------------------
       ! Input File
       in_file = trim(data_folder)//trim(hpp%hppname)//'_hpp_in.csv'
       
        ! Reading Datetime
        allocate(year(csv_len(in_file)),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
        date_format = '(i4)'
        open(file=in_file,newunit=lu,iostat=io_status,iomsg=io_message,
     &       status='OLD')
        if (io_status/=0) then
            write(*,*) io_message
        else
          if (allocated(year)) then
            read(lu,date_format)
     &      ((year(i)),i=1,csv_len(in_file))
            close(lu)
          else
            write(*,*) 'Warning: Array is not allocated!'
          endif
        endif
        
		
	      
        ! Allocating Arrays
        
        ! HPP Inflow
        allocate(hpp%tvar%inflow(year(1):year(size(year)),
     &                           1:12,
     &                           1:31,
     &                           0:23,
     &                           0:0,
     &                           0:0                ),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message        
        
        ! PH Curve
        if (hpp%use_ph) then
        in_file = trim(data_folder)//trim(hpp%hppname)//'_hpp_ph_s.csv'
        allocate(hpp%ph(csv_len(in_file)-1,2),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
        endif

        ! Initializing all variables
        hpp%tvar%inflow = 0.0
        hpp%ph = 0.0
        
        ! Reading Inflow
        in_file = trim(data_folder)//trim(hpp%hppname)//'_hpp_in.csv'
        open(file=in_file,newunit=lu,iostat=io_status,iomsg=io_message,
     &       status='OLD')
        if (io_status/=0) then
            write(*,*) io_message
		    else
          if (allocated(hpp%tvar%inflow)) then
            write(*,*) 'Reading '//trim(in_file)
            do i = 1, csv_len(in_file)
               read(lu,'(i4,x,i2,x,i2,x,i2,x,i2,x,i2,x,f10.3)')
     &         iyear,imonth,iday,ihour,imin,isec,r4_var
               hpp%tvar%inflow(iyear,imonth,iday,ihour,0,0) = 
     &         r4_var
            enddo
          else
            write(*,*) 'Warning: Array is not allocated!'
          endif
          close(lu)
        endif

        ! Reading PH Curves
        if (hpp%use_ph) then
        in_file = trim(data_folder)//trim(hpp%hppname)//'_hpp_ph_s.csv'
        open(file=in_file,newunit=lu,iostat=io_status,iomsg=io_message,
     &       status='OLD')
        if (io_status/=0) then
            write(*,*) io_message
		    else
          if (allocated(hpp%ph)) then
            write(*,*) 'Reading '//trim(in_file)
            read(lu,'(X)') 
            read(lu,*)
     &      ((hpp%ph(i,j), j=1,2), i = 1,(csv_len(in_file)-1))
          else
            write(*,*) 'Warning: Array is not allocated!'
          endif
          close(lu)
        endif
        
        endif
        
        !---------------------------------------------------------------------

       case(2)!---------------Forecast Simulation----------------------------
       ! Input File
       in_file = trim(data_folder)//trim(hpp%hppname)//
     &  '_hpp_in_opt.csv'
       
        ! Reading Datetime
        allocate(year(csv_len(in_file)),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
        date_format = '(i4)'
        open(file=in_file,newunit=lu,iostat=io_status,iomsg=io_message,
     &       status='OLD')
        if (io_status/=0) then
            write(*,*) io_message
        else
          if (allocated(year)) then
            read(lu,date_format)
     &      ((year(i)),i=1,csv_len(in_file))
            close(lu)
          else
            write(*,*) 'Warning: Array is not allocated!'
          endif
        endif
        

        ! Allocating Arrays
        
        ! HPP Inflow
        allocate(hpp%tvar%inflow(year(1):year(size(year)),
     &                           1:12,
     &                           1:31,
     &                           0:23,
     &                           0:0,
     &                           0:0                ),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message        
        
        ! PH Curve
        if (hpp%use_ph) then
        in_file = trim(data_folder)//trim(hpp%hppname)//'_hpp_ph_s.csv'
        allocate(hpp%ph(csv_len(in_file)-1,2),
     &           stat=io_status, errmsg=io_message)
        if (io_status /= 0) write(*,*) io_message
        endif

        ! Initializing all variables
        hpp%tvar%inflow = 0.0
        hpp%ph = 0.0
        
        ! Reading Inflow
        in_file = trim(data_folder)//trim(hpp%hppname)//
     &   '_hpp_in_opt.csv'
        open(file=in_file,newunit=lu,iostat=io_status,iomsg=io_message,
     &       status='OLD')
        if (io_status/=0) then
            write(*,*) io_message
		    else
          if (allocated(hpp%tvar%inflow)) then
            write(*,*) 'Reading '//trim(in_file)
            do i = 1, csv_len(in_file)
               read(lu,'(i4,x,i2,x,i2,x,i2,x,i2,x,i2,x,f10.3)')
     &         iyear,imonth,iday,ihour,imin,isec,r4_var
               hpp%tvar%inflow(iyear,imonth,iday,ihour,0,0) = 
     &         r4_var
            enddo
          else
            write(*,*) 'Warning: Array is not allocated!'
          endif
          close(lu)
        endif

        ! Reading PH Curves
        if (hpp%use_ph) then
        in_file = trim(data_folder)//trim(hpp%hppname)//'_hpp_ph_s.csv'
        open(file=in_file,newunit=lu,iostat=io_status,iomsg=io_message,
     &       status='OLD')
        if (io_status/=0) then
            write(*,*) io_message
		    else
          if (allocated(hpp%ph)) then
            write(*,*) 'Reading '//trim(in_file)
            read(lu,'(X)') 
            read(lu,*)
     &      ((hpp%ph(i,j), j=1,2), i = 1,(csv_len(in_file)-1))
          else
            write(*,*) 'Warning: Array is not allocated!'
          endif
          close(lu)
        endif
        endif
        !---------------------------------------------------------------------

      case(3)!======Outflow Override (Exact Data)=====================
        ! Input File
        in_file = trim(data_folder)//trim(hpp%hppname)//'_hpp_in.csv'
        
         ! Reading Datetime
         allocate(year(csv_len(in_file)),
     &           stat=io_status, errmsg=io_message)
         if (io_status /= 0) write(*,*) io_message
         date_format = '(i4)'
         open(file=in_file,newunit=lu,iostat=io_status,iomsg=io_message,
     &       status='OLD')
         if (io_status/=0) then
             write(*,*) io_message
         else
           if (allocated(year)) then
             read(lu,date_format)
     &      ((year(i)),i=1,csv_len(in_file))
             close(lu)
           else
             write(*,*) 'Warning: Array is not allocated!'
           endif
         endif
         
     
         
         ! Allocating Arrays
         
         ! HPP Inflow
         allocate(hpp%tvar%inflow(year(1):year(size(year)),
     &                           1:12,
     &                           1:31,
     &                           0:23,
     &                           0:0,
     &                           0:0                ),
     &           stat=io_status, errmsg=io_message)
         if (io_status /= 0) write(*,*) io_message        
         
         ! PH Curve
         if (hpp%use_ph) then
         in_file = trim(data_folder)//trim(hpp%hppname)//'_hpp_ph_s.csv'
         allocate(hpp%ph(csv_len(in_file)-1,2),
     &           stat=io_status, errmsg=io_message)
         if (io_status /= 0) write(*,*) io_message
         endif
 
         ! Initializing all variables
         hpp%tvar%inflow = 0.0
         hpp%ph = 0.0

         
         !---------------------------------------------------------------------

      end select
        

        ! Outflow  = Inflow (Reasonable Assumption)
        hpp%tvar%outflow = hpp%tvar%inflow

	  
	  end subroutine read_hpp_data
      !====================================================================
      
      
	  !====================================================================
      integer function csv_len(in_file)
      !
      ! Purpose:
      ! --------
      !   This function return the number of records in a csv file
      ! 
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 17/11/2018          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------
        character(len=*), intent(in) :: in_file       ! Input CSV file
        
        integer           :: lu            ! I/O Unit
        integer           :: io_status     ! I/O Status
        character(len=200):: io_message    ! I/O Error Message
		character(len=200):: dummy         ! Dummy Variable
      !--------------------------------------------------------------------
		
        open(file=in_file,newunit=lu,iostat=io_status,iomsg=io_message,
     &       status='OLD')
        
        if (io_status/=0) then
            write(*,*) io_message
            csv_len = 1
        
		else
		   csv_len = 0
		   do
            read(lu,'(a)',end=100) dummy
			csv_len = csv_len + 1
		   end do
           
100     continue
        endif
		
        close(lu)
        
      end function csv_len
	  !=====================================================================
      
      !====================================================================
      real(r8) function store_to_el(store,vh)
      !
      ! Purpose:
      ! --------
      !   Retrieves the value of reservoir volume (m3) against a given water level (m)
      !   using the reservoirs VH Curve
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 17/11/2018          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------      
        implicit none
        
        real(r8),dimension(:,:), intent(in) :: vh
        real(r8),                intent(in) :: store
        
        ! Local Variables
        real(r8), dimension(2) :: x, y  ! For linear interpolation
        integer                :: i, j  ! Iterators

           
        if (store >= vh(size(vh(:,3)),3)) then
            store_to_el = vh(size(vh(:,3)),2)
!            write(*,*) 
!     &      'Warning: Storage is higher than the maximum capacity'
        elseif (store <= vh(1,3)) then
            store_to_el = vh(1,2)
!            write(*,*) 
!     &      'Warning: Storage is lower than the minimum capacity'
        else
            do i = 1, size(vh(:,3))
               if (store <= vh(i,3)) then
                   
                 x(1) = vh(i-1,3)
                 x(2) = vh(i,3)
                 y(1) = vh(i-1,2)
                 y(2) = vh(i,2)
                   
                   
                 store_to_el=((y(2)-y(1))/(x(2)-x(1)))*(store-x(1))+y(1)
                   
                 exit
                   
               endif
            enddo
        endif

      end function store_to_el
      !=====================================================================
      
      
      !====================================================================
      real(r8) function el_to_store(el,vh)
      !
      ! Purpose:
      ! --------
      !   Retrieves the value of reservoir volume (m3) against a given water level (EL. m)
      !   using the reservoirs VH Curve
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 17/11/2018          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------      
        implicit none
        
        real(r8),dimension(:,:), intent(in) :: vh
        real(r8),                intent(in) :: el
        
        ! Local Variables
        real(r8), dimension(2) :: x, y  ! For linear interpolation
        integer                :: i, j  ! Iterators
        
        
           if (el >= vh(size(vh(:,2)),2)) then
              el_to_store = vh(size(vh(:,2)),3)
!              write(*,*) 
!     &         'Warning: Water Level is higher than the maximum level'
           elseif (el <= vh(1,2)) then
              el_to_store = vh(1,3)
!              write(*,*) 
!     &         'Warning: Water Level is lower than the minimum level'
           else
             do i = 1, size(vh(:,2))
                if (el <= vh(i,2)) then
                   
                   x(1) = vh(i-1,2)
                   x(2) = vh(i,2)
                   y(1) = vh(i-1,3)
                   y(2) = vh(i,3)
                   
                   el_to_store=((y(2)-y(1))/(x(2)-x(1)))*(el-x(1))+y(1)
                   
                   exit
                   
                endif
             enddo
           endif

      end function el_to_store
      !=====================================================================
      
      !====================================================================
      real(r8) function store_to_wl(store,vh)
      !
      ! Purpose:
      ! --------
      !   Retrieves the value of reservoir volume (m3) against a given water level (EL. m)
      !   using the reservoirs VH Curve
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 17/11/2018          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------      
        implicit none
        
        real(r8),dimension(:,:), intent(in) :: vh
        real(r8),                intent(in) :: store
        
        ! Local Variables
        real(r8), dimension(2) :: x, y  ! For linear interpolation
        integer                :: i, j  ! Iterators

           
        if (store >= vh(size(vh(:,3)),3)) then
            store_to_wl = vh(size(vh(:,3)),1)
!            write(*,*) 
!     &      'Warning: Storage is higher than the maximum capacity'
        elseif (store <= vh(1,3)) then
            store_to_wl = vh(1,1)
!            write(*,*) 
!     &      'Warning:Storage is lower than the minimum capacity'
        else
            do i = 1, size(vh(:,3))
               if (store <= vh(i,3)) then
                   
                 x(1) = vh(i-1,3)
                 x(2) = vh(i,3)
                 y(1) = vh(i-1,1)
                 y(2) = vh(i,1)
                  
                   
                 store_to_wl=((y(2)-y(1))/(x(2)-x(1)))*(store-x(1))+y(1)
                   
                 exit
                   
               endif
            enddo
        endif

      end function store_to_wl
      !=====================================================================
      
      
      
      !====================================================================
      real(r8) function wl_to_store(wl,vh)
      !
      ! Purpose:
      ! --------
      !   Retrieves the value of reservoir volume (m3) against a given water level (m)
      !   using the reservoirs VH Curve
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 17/11/2018          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------      
        implicit none
        
        real(r8),dimension(:,:), intent(in) :: vh
        real(r8),                intent(in) :: wl
        
        ! Local Variables
        real(r8), dimension(2) :: x, y  ! For linear interpolation
        integer                :: i, j  ! Iterators
        
        
        
           
           if (wl >= vh(size(vh(:,1)),1)) then
              wl_to_store = vh(size(vh(:,1)),3)
!              write(*,*) 
!     &         'Warning: Water Level is higher than the maximum level'
           elseif (wl <= vh(1,1)) then
              wl_to_store = vh(1,3)
!              write(*,*) 
!     &         'Warning: Water Level is lower than the minimum level'
           else
             do i = 1, size(vh(:,1))
                if (wl <= vh(i,1)) then
                   
                   x(1) = vh(i-1,1)
                   x(2) = vh(i,1)
                   y(1) = vh(i-1,3)
                   y(2) = vh(i,3)
                   
                   wl_to_store=((y(2)-y(1))/(x(2)-x(1)))*(wl-x(1))+y(1)
                   
                   exit
                   
                endif
             enddo
           endif

      end function wl_to_store
      !=====================================================================
      

      !====================================================================
      real(r8) function wl_to_rate_seasonal(wl,ph,summer)
      ! This is a bit complicated so this is not used. Instead the simpified version is used
      ! Purpose:
      ! --------
      !   Retrieves the value of energy rate against a given water level (m)
      !   using the reservoirs PH Curve
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 17/11/2018          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------      
        implicit none
        
        real(r8),dimension(:,:), intent(in) :: ph
        real(r8),                intent(in) :: wl
        logical,                 intent(in) :: summer
        integer                             :: ratio_col
        
        ! Local Variables
        real(r8), dimension(2) :: x, y  ! For linear interpolation
        integer                :: i, j  ! Iterators
        
        if (summer) then
          ratio_col = 3
        else
          ratio_col = 2
        endif
        
           
           if (wl >= ph(size(ph(:,1)),1)) then
              wl_to_rate_seasonal = ph(size(ph(:,1)),ratio_col)
!              write(*,*) 
!     &         'Warning: Water Level is higher than the maximum level'
           elseif (wl <= ph(1,1)) then
              wl_to_rate_seasonal = ph(1,ratio_col)
!              write(*,*) 
!     &         'Warning: Water Level is lower than the minimum level'
           else
             do i = 1, size(ph(:,1))
                if (wl <= ph(i,1)) then
                   
                   x(1) = ph(i-1,1)
                   x(2) = ph(i,1)
                   y(1) = ph(i-1,ratio_col)
                   y(2) = ph(i,ratio_col)
                   
                   wl_to_rate_seasonal=((y(2)-y(1))/(x(2)-x(1)))
     &              *(wl-x(1))+y(1)
                   
                   exit
                   
                endif
             enddo
           endif

      end function wl_to_rate_seasonal
      !=====================================================================      
      
      !====================================================================
      real(r8) function wl_to_rate(wl,ph)
      ! This is a bit complicated so this is not used. Instead the simpified version is used
      ! Purpose:
      ! --------
      !   Retrieves the value of energy rate against a given water level (m)
      !   using the reservoirs PH Curve
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 17/11/2018          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------      
        implicit none
        
        real(r8),dimension(:,:), intent(in) :: ph
        real(r8),                intent(in) :: wl
        
        ! Local Variables
        real(r8), dimension(2) :: x, y  ! For linear interpolation
        integer                :: i, j  ! Iterators
        
        
           
           if (wl >= ph(size(ph(:,1)),1)) then
              wl_to_rate = ph(size(ph(:,1)),2)
!              write(*,*) 
!     &         'Warning: Water Level is higher than the maximum level'
           elseif (wl <= ph(1,1)) then
              wl_to_rate = ph(1,2)
!              write(*,*) 
!     &         'Warning: Water Level is lower than the minimum level'
           else
             do i = 1, size(ph(:,1))
                if (wl <= ph(i,1)) then
                   wl_to_rate = ph(i,2)
                   exit
                   
                endif
             enddo
           endif

      end function wl_to_rate
      !=====================================================================      
      
      
      !====================================================================
      real(r4) function time_lookup_r4(var,
     &                         iyear,imonth,iday,ihour,iminute,isecond)
      !
      ! Purpose:
      ! --------
      !   Retrieves the value corresponding to a datetime from a dam variable 
      !   Only for single precision variables
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 17/11/2018          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------
        implicit none
        
        integer, intent(in) :: iyear, imonth, iday
        integer, intent(in) :: ihour, iminute, isecond
        real(r4),allocatable,dimension(:,:,:,:,:,:),   intent(in) :: var
        
        integer :: i

        time_lookup_r4 = var(iyear,imonth,iday,ihour,0,0)
        
      end function time_lookup_r4
      !=====================================================================
      
      !====================================================================
      real(r8) function time_lookup_r8(var,
     &                         iyear,imonth,iday,ihour,iminute,isecond)
      !
      ! Purpose:
      ! --------
      !   Retrieves the value corresponding to a datetime from a dam variable 
      !   Only for double precision variables
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 17/11/2018          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------
        implicit none
        
        integer, intent(in) :: iyear, imonth, iday
        integer, intent(in) :: ihour, iminute, isecond
        real(r8),allocatable,dimension(:,:,:,:,:,:),   intent(in) :: var
        
        integer :: i

        time_lookup_r8 = var(iyear,imonth,iday,ihour,0,0)
        
      end function time_lookup_r8
      !=====================================================================

      !====================================================================
      real(r4) function maintenance_release(dam,iyear,imonth,iday,
     &                                     ihour,iminute,isecond)
      !
      ! Purpose:
      ! --------
      !   Gets the value of maintenance flows for a given time and dam
      !   based on dam operation manual and observed data
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 29/01/2020          Abdul Moiz (REEL, UTOKYO)  Original Code (Developed for Kurobe)
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------
        implicit none
        
        type(dam_data), intent(in) :: dam
        integer, intent(in) :: iyear, imonth, iday
        integer, intent(in) :: ihour, iminute, isecond
        
        !----------------------------Kurobe Maintenance Flow----------------------
        if (dam%damname == 'kurobe') then
          if (((imonth == 6) .and. (iday >=26)) .or.
     &       (imonth == 7)) then
             if (ihour == 6) then
              maintenance_release = 6.4
             elseif ((ihour >= 7) .and. (ihour <=17)) then
              maintenance_release = 15.2
             elseif (ihour == 18) then
              maintenance_release = 7.9
             else
              maintenance_release = 0
             endif

          elseif ((imonth == 8) .and. (iday<=15)) then
              if (ihour == 6) then
                maintenance_release = 2.3
              elseif (ihour == 7) then
                maintenance_release = 11.7
              elseif ((ihour >= 8) .and. (ihour <= 17)) then
                maintenance_release = 15.2
              elseif (ihour == 18) then
                maintenance_release = 0.4
              else
                maintenance_release = 0
              endif

          elseif (((imonth == 8) .and. (iday>=16)) .or.
     &           ((imonth == 9) .and. (iday<=10))) then
              if (ihour == 6) then
                maintenance_release = 1.7
              elseif (ihour == 7) then
                maintenance_release = 8.8
              elseif ((ihour >=8) .and. (ihour<=17)) then
                maintenance_release = 10.2
              elseif (ihour == 18) then
                maintenance_release = 0.3
              else
                maintenance_release = 0
              endif

          elseif (((imonth == 9) .and. (iday >= 11)) .or.
     &           ((imonth == 10) .and. (iday <= 15))) then
              if (ihour == 7) then
                maintenance_release = 5.4
              elseif ((ihour >= 8) .and. (ihour <= 16)) then
                maintenance_release = 10.2
              elseif (ihour == 17) then
                maintenance_release = 5.3
              else
                maintenance_release = 0
              endif
          
          else
            maintenance_release = 0

          endif
          !-------------------------Sennindani Maintenance Flow---------------------
          elseif (dam%damname == 'sennindani') then
            if (((imonth >= 6) .and. (imonth <=10)) .or.
     &          ((imonth == 11) .and. (iday <=10))) then

                if ((((imonth >= 9) .and. (imonth <=10)) .or.
     &              ((imonth == 11) .and. (iday <=10))) .and.
     &              (dow(iyear,imonth,iday) == 0)) then
                  if ((ihour >= 9) .and. (ihour <=16)) then
                  maintenance_release = 2.78
                  else
                  maintenance_release = 0.3
                  endif
                else
                  if ((ihour >= 9) .and. (ihour <=16)) then
                  maintenance_release = 1.39
                  else
                  maintenance_release = 0.3
                  endif
                endif
            else
              maintenance_release = 0.3
            endif

          !-------------------------Koyadaira Maintenance Flow---------------------
          elseif (dam%damname == 'koyadaira') then
            if (((imonth >= 6) .and. (imonth <=10)) .or.
     &          ((imonth == 11) .and. (iday <=10))) then

                if ((((imonth >= 9) .and. (imonth <=10)) .or.
     &              ((imonth == 11) .and. (iday <=10))) .and.
     &              (dow(iyear,imonth,iday) == 0)) then
                  if ((ihour >= 9) .and. (ihour <=16)) then
                  maintenance_release = 2.78
                  else
                  maintenance_release = 0.45
                  endif
                else
                  if ((ihour >= 9) .and. (ihour <=16)) then
                  maintenance_release = 1.39
                  else
                  maintenance_release = 0.45
                  endif
                endif
            else
              maintenance_release = 0.45
            endif
            
          !-------------------------Dashidaira Maintenance Flow---------------------
          elseif (dam%damname == 'dashidaira') then
            if (((imonth >= 5) .and. (imonth <=10)) .or.
     &          ((imonth == 11) .and. (iday <=10))) then

                if ((((imonth >= 9) .and. (imonth <=10)) .or.
     &              ((imonth == 11) .and. (iday <=10))) .and.
     &              (dow(iyear,imonth,iday) == 0)) then
                  if ((ihour >= 7) .and. (ihour <=17)) then
                    maintenance_release = 3.65
                  else
                    maintenance_release = 0.5
                  endif
                else
                  if ((ihour >= 7) .and. (ihour <=17)) then
                    maintenance_release = 1.675
                  else
                    maintenance_release = 0.5
                  endif
                endif
            else
              maintenance_release = 0.5
            endif
          else
              maintenance_release = 0
          endif
      end function maintenance_release
      !=====================================================================

      !====================================================================
      !
      ! Purpose:
      ! --------
      !   Gets the day of the week
      ! Day_Of_Week: (0=Sunday,1=Monday...6=Saturday)
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 05/02/2003          J.D.Robertson              https://software.intel.com/en-us/forums/intel-fortran-compiler/topic/272033
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------
      integer function dow(yyyy,mm,dd)
      implicit none

      integer, intent(in) :: yyyy,mm,dd

      dow = mod((13*(mm+10-(mm+10)/13*12)-1)/5+dd+77 
     & +5*(yyyy+(mm-14)/12-(yyyy+(mm-14)/12)/100*100)/4 
     & +(yyyy+(mm-14)/12)/400-(yyyy+(mm-14)/12)/100*2,7)
      end function dow
      !=====================================================================

      !====================================================================
      !
      ! Purpose:
      ! --------
      !   Divides K4 discharge to K3 and NK3 discharge
      !
      ! Record of Revisions:
      ! --------------------
      !
      !   Date              Programmer                 Description of Change
      !   ----              ----------                 ---------------------
      ! 2020/01/31          Abdul Moiz                 Original Code
      !====================================================================
      !
      ! Data Dictionary:
      ! ----------------
      subroutine K4_to_NK3K3(k4,nk3,k3)
      implicit none

      real(r4), intent(in) :: k4
      real(r4), intent(out):: nk3
      real(r4), intent(out):: k3

      ! Local Variables
      real(r4), dimension(10,3) :: lookup

      real(r4), dimension(2) :: x, y  ! For linear interpolation
      integer                              :: i

      
      lookup(:,1) = [30.0,35.0,40.0,45.0,50.0,55.0,60.0,65.0,70.0,72.0]
      lookup(:,2) = [0.0,2.5,5.2,8.1,11.2,14.4,17.4,20.6,23.7,24.9]
      lookup(:,3) = [30.0,32.5,34.8,36.9,38.8,40.6,42.6,44.4,46.3,47.1]

      if (k4 >= lookup(size(lookup(:,1)),1)) then
        nk3 = lookup(size(lookup(:,3)),3)
        k3 =  k3+lookup(size(lookup(:,2)),2)

      elseif (k4 <= lookup(1,1)) then
        nk3 = min(lookup(1,3),k4)
        k3 = k3+0
      
      else
        do i = 1, size(lookup(:,1))
          if (k4 <= lookup(i,1)) then
            x(1) = lookup(i-1,1)
            x(2) = lookup(i,1)
            y(1) = lookup(i-1,3)/lookup(i-1,1)
            y(2) = lookup(i,3)/lookup(i,1)

            nk3 = k4*(((y(2)-y(1))/(x(2)-x(1)))*(k4-x(1))+y(1))
            k3 = k3+(k4-nk3)
            
            exit
          endif
        enddo
      endif

      end subroutine K4_to_NK3K3
      !=====================================================================



      end module damops
