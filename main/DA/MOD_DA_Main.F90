#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DataAssimilation

   USE MOD_Precision
   !USE MOD_DA_GRACE
   USE MOD_DA_SMAP
   IMPLICIT NONE

CONTAINS

   ! ----------
   SUBROUTINE init_DataAssimilation (idate, deltim)
      
   IMPLICIT NONE
   integer,  intent(in) :: idate(3)
   real(r8), intent(in) :: deltim
      
      !CALL init_DA_GRACE ()
      CALL allocate_SMAP  (idate, deltim)

   END SUBROUTINE init_DataAssimilation

   ! ----------
   SUBROUTINE do_DataAssimilation ()
      
   IMPLICIT NONE
      

      !CALL do_DA_GRACE (idate, deltim)
      CALL run_DA_SMAP  ()
      

   END SUBROUTINE do_DataAssimilation

   ! ---------
   SUBROUTINE final_DataAssimilation ()

   IMPLICIT NONE

      !CALL final_DA_GRACE ()
      CALL deallocate_SMAP  ()

   END SUBROUTINE final_DataAssimilation

END MODULE MOD_DataAssimilation
#endif
