#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Main
!-----------------------------------------------------------------------------
! DESCRIPTION:
!     Main procedures for data assimilation
! 
! AUTHOR:
!     Lu Li, 12/2024
!     Zhilong Fan, Lu Li, 03/2024
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_DA_GRACE
   USE MOD_DA_SMAP
   IMPLICIT NONE
   SAVE


!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE init_DA ()
      
!-----------------------------------------------------------------------------
   IMPLICIT NONE
      
      IF (DEF_DA_GRACE) THEN
         CALL init_DA_GRACE ()
      ENDIF

      IF (DEF_DA_SMAP) THEN
         CALL init_DA_SMAP  ()
      ENDIF

   END SUBROUTINE init_DA



!-----------------------------------------------------------------------------

   SUBROUTINE run_DA (idate, deltim)
      
!-----------------------------------------------------------------------------
   IMPLICIT NONE

!---------------------Dummy arguments-----------------------------------------
   integer,  intent(in) :: idate(3)
   real(r8), intent(in) :: deltim   

!-----------------------------------------------------------------------------

      IF (DEF_DA_GRACE) THEN
         CALL run_DA_GRACE (idate, deltim)
      ENDIF

      IF (DEF_DA_SMAP) THEN
         CALL run_DA_SMAP  (idate, deltim)
      ENDIF
      
   END SUBROUTINE run_DA



!-----------------------------------------------------------------------------

   SUBROUTINE deallocate_DA_hist (idate, deltim, itstamp, etstamp) 

!-----------------------------------------------------------------------------
   USE MOD_DA_SMAP
   USE MOD_TimeManager
   IMPLICIT NONE

!---------------------Dummy arguments-----------------------------------------
      integer,  intent(in) :: idate(3)
      real(r8), intent(in) :: deltim
      type(timestamp), intent(in) :: itstamp
      type(timestamp), intent(in) :: etstamp

!-----------------------------------------------------------------------------
      logical :: lwrite

!-----------------------------------------------------------------------------

      SELECT CASE (trim(adjustl(DEF_HIST_FREQ)))
      CASE ('TIMESTEP')
         lwrite = .true.
      CASE ('HOURLY')
         lwrite = isendofhour (idate, deltim) .or. (.not. (itstamp < etstamp))
      CASE ('DAILY')
         lwrite = isendofday  (idate, deltim) .or. (.not. (itstamp < etstamp))
      CASE ('MONTHLY')
         lwrite = isendofmonth(idate, deltim) .or. (.not. (itstamp < etstamp))
      CASE ('YEARLY')
         lwrite = isendofyear (idate, deltim) .or. (.not. (itstamp < etstamp))
      CASE default
         lwrite = .false.
         write(*,*) 'Warning : Please USE one of TIMESTEP/HOURLY/DAILY/MONTHLY/YEARLY for history frequency.'
         write(*,*) '          Set to FALSE by default.                                                     '
      END SELECT

      IF (lwrite) THEN
         IF (allocated(pred_tb_h_out))   deallocate(pred_tb_h_out)
         IF (allocated(smap_tb_h_out))   deallocate(smap_tb_h_out)
         IF (allocated(pred_tb_v_out))   deallocate(pred_tb_v_out)
         IF (allocated(smap_tb_v_out))   deallocate(smap_tb_v_out)
      ENDIF

   END SUBROUTINE deallocate_DA_hist



!-----------------------------------------------------------------------------
   
   SUBROUTINE end_DA ()

!-----------------------------------------------------------------------------
   IMPLICIT NONE

      IF (DEF_DA_GRACE) THEN
         CALL end_DA_GRACE ()
      ENDIF

      IF (DEF_DA_SMAP) THEN
         CALL end_DA_SMAP  ()
      ENDIF

   END SUBROUTINE end_DA



!-----------------------------------------------------------------------------
END MODULE MOD_DA_Main
#endif
