#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_EnKF
!-----------------------------------------------------------------------
! DESCRIPTION:
!    ensemble Kalman filter (EnKF) 
! 
! AUTHOR:
! Lu Li, 12/2024: Initial version
!-----------------------------------------------------------------------
    USE MOD_Precision
    IMPLICIT NONE
    SAVE

! public functions
    PUBLIC :: letkf


!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

    SUBROUTINE letkf ( &
        num_ens, num_obs, num_state, &
        A, HA, y, R, loc_d, loc_r, infl, &
        A_analysis_mean, trans, A_analysis_ens)

!-----------------------------------------------------------------------
! Description:
!   local transform ensemble Kalman filter 
!
! Original author : 
!   Lu Li, 12/2024
!-----------------------------------------------------------------------
    USE MOD_Precision
    IMPLICIT NONE

!------------------------------dummy arguments----------------------------
    integer, intent(in)     :: num_ens                              ! ensemble size
    integer, intent(in)     :: num_obs                              ! number of observations
    integer, intent(in)     :: num_state                            ! number of state variables
    real(r8), intent(in)    :: A(num_state, num_ens)                ! ensemble matrix
    real(r8), intent(in)    :: HA(num_obs, num_ens)                 ! ensemble predicted observation matrix 
    real(r8), intent(in)    :: y(num_obs)                           ! observation vector
    real(r8), intent(in)    :: R(num_obs)                           ! observation error variance
    real(r8), intent(in)    :: loc_d(num_obs)                       ! localization distance
    real(r8), intent(in)    :: loc_r                                ! localization radius
    real(r8), intent(in)    :: infl                                 ! inflation factor
    real(r8), intent(out)   :: A_analysis_mean(num_state)           ! analysis mean
    real(r8), intent(out)   :: trans(num_ens, num_ens)              ! transform matrix (k x k)
    real(r8), intent(out)   :: A_analysis_ens(num_state, num_ens)   ! analysis ensemble (m x k)

!------------------------------local variables----------------------------
    real(r8) :: HA_mean(num_obs)                    ! mean of ensemble predicted observation (l)
    real(r8) :: dHA(num_obs, num_ens)               ! HA - mean(HA) (l x k)
    real(r8) :: A_mean(num_state)                   ! mean of ensemble (m)
    real(r8) :: dA(num_state, num_ens)              ! A - mean(A) (m x k)
    real(r8) :: dHA_t(num_ens, num_obs)             ! transpose of dHA (k x l)
    real(r8) :: C(num_ens, num_obs)                 ! C = (dHA)^T * (R)^-1 (k x l)
    real(r8) :: M1(num_ens, num_ens)                ! C * dHA (k x k)
    real(r8) :: pa_inv(num_ens, num_ens)            ! inverse of background error covariance matrix (k x k)
    real(r8) :: eigval(num_ens)                     ! eigenvalues of pa_inv (k)
    real(r8) :: eigvec(num_ens, num_ens)            ! eigenvectors of pa_inv (k x k)
    integer  :: lwork
    real(r8), allocatable :: work(:)
    integer  :: err
    real(r8) :: M2(num_ens, num_ens)                ! M2 = eigvec * eigval^-1 (k x k)
    real(r8) :: pa(num_ens, num_ens)                ! background error covariance matrix (k x k)
    real(r8) :: M3(num_ens, num_obs)                ! M3 = pa * C (k x l)
    real(r8) :: delta(num_obs)                      ! increment of observation (l)
    real(r8) :: w_avg(num_ens)                      ! weight (k)
    real(r8) :: A_pert(num_state)                   ! analysis mean perturbation (m)
    real(r8) :: M4(num_ens, num_ens)                ! M4 = eigvec * sqrt((k-1)/eigval) (k x k)
    real(r8) :: trans_pert(num_ens, num_ens)        ! perturbation transform matrix (k x k)
    real(r8) :: M5(num_state, num_ens)              ! M5 = dA * trans (m x k)
    integer  :: i, j

!-----------------------------------------------------------------------

        ! calculate observation space perturbation
        HA_mean = sum(HA, dim=2) / size(HA, dim=2)   !(lx1)
        DO j = 1, num_ens
            dHA(:,j) = HA(:,j) - HA_mean !(lxk)
        ENDDO   

        ! calculate background state perturbation
        A_mean = sum(A, dim=2) / size(A, dim=2)   !(mx1)
        DO j = 1, num_ens
            dA(:,j) = A(:,j) - A_mean !(mxk)
        ENDDO

        ! calculate C, intermediate matrix in localized observation 
        dHA_t = transpose(dHA) !(kxl)
        DO j = 1, num_obs
            C(:,j) = dHA_t(:,j) / (R(j))!*(exp((-loc_d(j)**2)/(2*loc_r**2)))) !(kxl) !//TODO: Lu Li: add localization
        ENDDO

        ! calculate C*dHA, intermediate matrix in background error M1 
        CALL dgemm('N', 'N', num_ens, num_ens, num_obs, 1.0d0, C, num_ens, dHA, num_ens, 0.0d0, M1, num_ens) !(kxk)

        ! calculate inverse of background error 
        DO i = 1, num_ens
            pa_inv(i,i) = M1(i,i) + (num_ens-1)*1.0d0/infl !(kxk)
        ENDDO

        ! eigenvalues and eigenvectors of inverse of background error
        lwork = 4 * num_ens
        allocate( work(lwork) )
        CALL dsyev('V', 'U', num_ens, pa_inv, num_ens, eigval, work, lwork, err)
        eigvec = pa_inv !(kxk)

        ! calculate background error covariance matrix pa = eigvec (eigval)^-1 eigvec^T
        DO i = 1, num_ens
            M2(:,i) = eigvec(:,i) / eigval(i) !(kxk)
        ENDDO
        CALL dgemm('N', 'T', num_ens, num_ens, num_obs, 1.0d0, M2, num_ens, eigvec, num_ens, 0.0d0, pa, num_ens) !(kxk)

        ! caculate pa * C, intermediate matrix in Kalman gain M3
        CALL dgemm('N', 'N', num_ens, num_obs, num_ens, 1.0d0, pa, num_ens, C, num_ens, 0.0d0, M3, num_ens) !(kxl)

        ! calculate weight
        delta = y - HA_mean
        CALL dgemm('N', 'N', num_ens, 1, num_obs, 1.0d0, M3, num_ens, delta, num_obs, 0.0d0, w_avg, num_ens) !(kx1)
        
        ! calculate mean of analysis ensemble 
        CALL dgemm('N', 'N', num_state, 1, num_ens, 1.0d0, dA, num_state, w_avg, num_ens, 0.0d0, A_pert, num_state) !(mx1)
        A_analysis_mean = A_mean + A_pert !(mx1)

        ! calculate pertubation transform matrix 
        DO j = 1, num_ens
            M4(:,j) = eigvec(:,j) * sqrt((num_ens-1) / eigval(j))  !(kxk)
        ENDDO
        CALL dgemm('N', 'T', num_ens, num_ens, num_ens, 1.0d0, M4, num_ens, eigvec, num_ens, 0.0d0, trans_pert, num_ens) !(kxk)

        ! calculate transform matrix
        DO j = 1, num_ens
            trans(:,j) = trans_pert(:,j) + w_avg(j) !(kxk)
        ENDDO

        ! caclulate analysis ensemble
        CALL dgemm('N', 'N', num_state, num_ens, num_ens, 1.0d0, dA, num_state, trans, num_ens, 0.0d0, M5, num_state) !(mxk)
        DO j = 1, num_ens
            A_analysis_ens(:,j) = A_mean(j) + M5(:,j) !(mxk)
        ENDDO

    END SUBROUTINE letkf



!-----------------------------------------------------------------------
END MODULE MOD_DA_EnKF
#endif