      module discharge_io
         use time_and_fileio
         use sib2_public_vars


         contains


       subroutine save_discharge_info_bin(
     :    sub_catch, q2,num_flow)
         implicit none

! ===== INCLUDES
         include 'globcst.inc'  ! for simulation_dir
         include 'timectrl.inc' ! for ilead

! ===== PARAMETERS
         character, intent(in) :: sub_catch*(*)
         real, intent(in) :: q2(:)
         integer, intent(in) :: num_flow
!         integer, intent(in) :: step_num

! ===== LOCAL VARIABLES
         integer ldir_s
         integer smlunit
         integer arr_size, word_size
         integer i
!         character*4 iyear
         character(len=128) filepath

         integer step_num

! ===== CODE
         if (result_outflow_mode == 0) return

! add one to record id, records start from 1
! TODO: stepnum is probably already defined elsewhere (PRL)
         if (ilead < 0) return

         if (ilead == 0) then
           step_num = (((jday-1)*24+hour)*60+minute)*60+second
           step_num = (step_num / dt_couple)+1
           ! Convert the outflow path format into a true string
           write(filepath,outflow_path_format) sub_catch,year
         end if
         if (ilead > 0) then
           step_num = ilead
           write(filepath,outflow_path_format) sub_catch
         end if

!         write(iyear,  '(i4.4)') year

         ldir_s = len(simulation_dir)
         call strlnth(simulation_dir, ldir_s)
         call getunit(smlunit)

! Determine the units (words or bytes) for recl of unformatted IO
         inquire (IOLENGTH=arr_size) q2
         word_size = arr_size/size(q2)


!     &       simulation_dir(1:ldir_s)//
!!     &      sub_catch //'bI'//iyear//'_flow', status='unknown',
!     &      sub_catch //'_flow.bin', status='unknown',
!     &      form='unformatted', access='direct',recl=word_size*num_flow)

         open(smlunit,file=filepath, status='unknown',
     &      form='unformatted', access='direct',recl=word_size*num_flow)
!         print *, 'NOTES:',jday,',',hour,',',minute,'=',step_num
!        write(smlunit, rec=step_num) q2
         write(smlunit, rec=step_num) (q2(i), i=1,num_flow)
         close(smlunit)
         call retunit(smlunit)

       end subroutine

       subroutine save_discharge_info_text(sub_catch, q2,num_flow)
         implicit none

! ===== INCLUDES
         include 'globcst.inc'  ! for simulation_dir

! ===== PARAMETERS
         character, intent(in) :: sub_catch*(*)
         real, intent(in) :: q2(:)
         integer, intent(in) :: num_flow

! ===== LOCAL VARIABLES
         integer smlunit
         integer i

! ===== CODE
         call getunit(smlunit)
         open(smlunit,file=trim(simulation_dir)//
     &           sub_catch//'I_flow', status='unknown')
         do i=1,num_flow
            write(smlunit,'(2x,i4,f20.3)') i,q2(i)
         end do
         close(smlunit)
         call retunit(smlunit)

         return
      end subroutine save_discharge_info_text

#if 0
! I did not write this code (PRL)
!> \todo Get rid of this code, do something else instead
       subroutine save_discharge_info_text(
     :    ii,Qh)
          implicit none

          include 'globcst.inc'
          include 'hydro.inc'

          integer, intent(in) :: ii
          real, intent(in) :: Qh(nc,8800)          ! hourly mean discharge (m3/s)

          character dis_gauge*100, cyear*4
          integer j,k
          integer ldir_d, ldir_r, lfile
          integer obsunit,relunit
          integer iyy,imm,idd,ihh,value

c*************************************************
c     save river discharge at selected gauges
          if( month==12.and.day==31.and.hour==23
     &          .or.  (curtim + dt_couple .eq. tstop)  ) then
             ldir_d = len(data_dir)
             call strlnth(data_dir,ldir_d)
             ldir_r = len(result1_dir)
             call strlnth(result1_dir,ldir_r)

             write(cyear, '(i4.4)') year


c---------------------------------------------------------------------------
c    1. output discharge of gauges
c---------------------------------------------------------------------------
             dis_gauge='no_gauge'
             if(sub_catch(ii) == 'ws210') dis_gauge='murakami'
             if(sub_catch(ii) == 'ws910') dis_gauge='yuhara'
             if(sub_catch(ii) == 'ws800') dis_gauge='kosodebashi'
             if(sub_catch(ii) .eq. 'ws500') dis_gauge='yakatabara'
             if(sub_catch(ii) .eq. 'ws410') dis_gauge='kamikuya'
             if(sub_catch(ii) .eq. 'ws300') dis_gauge='iwamoto'
             if(sub_catch(ii) .eq. 'ws270') dis_gauge='iwashima'
             if(sub_catch(ii) .eq. 'ws100') dis_gauge='maebashi'

c     Dam inflow
             if(sub_catch(ii) .eq. 'ws970') dis_gauge='yagisawa_in'
             if(sub_catch(ii) .eq. 'ws960') dis_gauge='naramata_in'
             if(sub_catch(ii) .eq. 'ws800') dis_gauge='aimata_in'
             if(sub_catch(ii) .eq. 'ws430') dis_gauge='sonohara_in'

             lfile = len(dis_gauge)
             call strlnth(dis_gauge,lfile)
             if(trim(dis_gauge) == 'no_gauge') return


c     hourly discharge
             call getunit(relunit)
             call getunit(obsunit)
             open(relunit,file=trim(result1_dir)//trim(dis_gauge)
     &          //'.'//cyear//'.hourly',status='unknown')
             open(obsunit,file=trim(data_dir)//
     &          trim(dis_gauge)//'.hourly', status='old')
             k=0
 251         read(obsunit,*,end=351) iyy,imm,idd,ihh, value
             if(k == 0 .and. iyy > year) goto 351
             if(iyy == year .and. imm == 1 .and.
     &          idd == 1 .and. ihh == 1) then
                k=1
                j=0
             endif
             if(k == 1) then
                j=j+1
                write(relunit,'(4i6,2f15.3)')
     &             iyy,imm,idd,ihh,Qh(ii,j),value !wanglei12/01/2006
                if(iyy. eq. year .and.
     &           imm == 12 .and. idd == 31.and. ihh == 24) goto 101 !this part determines when it writes
                endif
             goto 251
 351         if(k == 0) write(relunit,451)
 451            format(1x,'No discharge data recorded for this year')
 101         close(relunit)
             close(obsunit)
             call retunit(relunit)
             call retunit(obsunit)

          end if
          return

       end subroutine
#endif

      end module discharge_io
