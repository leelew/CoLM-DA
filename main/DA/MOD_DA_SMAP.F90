#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_SMAP
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Data assimilation of brightness temperature 
! 
! AUTHOR:
!   Lu Li, 12/2024: Initial version, based on SMAP L1C TB data
!-----------------------------------------------------------------------------
    USE MOD_DataType
    USE MOD_SpatialMapping
    USE MOD_DA_ObsOperator
    USE MOD_DA_EnKF
    USE MOD_DA_Vars_TimeVariables
    USE MOD_Vars_Global, only : num_ens
    USE MOD_LandPatch
    IMPLICIT NONE
    SAVE

! PUBLIC MEMBER FUNCTIONS:
    PUBLIC :: allocate_SMAP
    PUBLIC :: run_DA_SMAP
    PUBLIC :: deallocate_SMAP

    PRIVATE

! local variables
    ! file path
    character(len=256) :: file_smap                           ! SMAP file path
    character(len=256) :: file_grid                           ! grid of SMAP file path

    ! grid
    type(grid_type) :: grid_smap                              ! grid of SMAP
    type(spatial_mapping_type) :: mg2p_smap                   ! mapping from grid to patch

    ! time
    integer :: month, mday, hour                              ! month, day, hour of current step
    character(len=256) :: yearstr, monthstr, daystr, hourstr  ! string of year, month, day, hour
    real(r8), allocatable :: smap_time (:)                    ! UTC time of each observation in current file
    real(r8), allocatable :: dt_begin  (:)
    real(r8), allocatable :: dt_end    (:)
    integer :: idate_begin, idate_end                         ! begin and end seconds of current step
    
    ! obs
    logical :: has_obs, has_file
    integer :: num_obs
    integer :: num_loc_obs
    real(r8), allocatable :: smap_lat  (:)                    ! latitude 
    real(r8), allocatable :: smap_lon  (:)                    ! longitude
    real(r8), allocatable :: smap_tb_h (:)                    ! H- polarized brightness temp
    real(r8), allocatable :: smap_tb_v (:)                    ! V- polarized brightness temp
    integer, allocatable :: smap_ii   (:)                    ! i-th grid in world map 
    integer, allocatable :: smap_jj   (:)                    ! j-th grid in world map
    integer, allocatable :: obs_index (:)                    ! index of avaliable obs for each patch  
    integer, allocatable :: obs_index_tmp (:)                    ! index of avaliable obs for each patch  
    real(r8), allocatable :: obs_tb_h (:)
    real(r8), allocatable :: obs_tb_v (:)
    real(r8), allocatable :: obs_lat (:)
    real(r8), allocatable :: obs_lon (:)
    real(r8), allocatable :: obs_d   (:)
    real(r8), allocatable :: obs_err (:)

    ! predict tb
    real(r8), allocatable :: pred_tb_obs_ens_h (:,:)          ! predicted H- polarized temp on obs grid
    real(r8), allocatable :: pred_tb_obs_ens_v (:,:)          ! predicted V- polarized temp on obs grid
    real(r8), allocatable :: pred_tb_pset_ens_h (:,:)         ! predicted H- polarized temp on patch
    real(r8), allocatable :: pred_tb_pset_ens_v (:,:)         ! predicted V- polarized temp on patch
    real(r8), allocatable :: pred_tb_loc_ens_h (:,:)          ! avaliable predicted tb for each patch
    real(r8), allocatable :: pred_tb_loc_ens_v (:,:)

    ! data assimilation
    real(r8), parameter   :: dres = 1.0                       ! search localization radius (deg)
    integer, parameter    :: num_state = 10                   ! number of state
    real(r8), parameter   :: loc_r = 1.0                      ! localization radius 
    real(r8), parameter   :: infl = 1.2                       ! inflation factor
    real(r8) :: trans(num_ens, num_ens)                       ! transform matrix on each patch


!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

    SUBROUTINE allocate_SMAP (idate, deltim)

!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Initialize the data assimilation of (SMAPL1C) TB data
!-----------------------------------------------------------------------------
    USE MOD_Spmd_Task
    USE MOD_Namelist, only : DEF_DA_obsdir
    USE MOD_Grid
    USE MOD_NetCDFSerial
    USE MOD_LandPatch
    USE MOD_Pixelset
    USE MOD_Vars_TimeInvariants, only : patchtype
    USE MOD_Forcing, only : forcmask_pch
    USE MOD_RangeCheck
    IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
    real(r8), intent(in) :: deltim
    integer, intent(in)  :: idate(3)

!-----------------------------------------------------------------------

        ! grid file path of EASE v2.0, 36km world grid
        file_grid = trim(DEF_DA_obsdir) // '/NSIDC0772_LatLon_EASE2_M36km_v1.0_pre.nc'

        ! read grid from file
        CALL grid_smap%define_from_file (file_grid, 'latitude', 'longitude')

        ! map SMAP grid to patch
        CALL mg2p_smap%build_arealweighted (grid_smap, landpatch)

        ! calculate year/month/day/hour of current step
        CALL julian2monthday(idate(1), idate(2), month, mday)
        hour = int(idate(3)/3600)

        ! define file path of SMAP of each day
        write (yearstr,  '(I4.4)') idate(1)
        write (monthstr, '(I2.2)') month
        write (daystr,   '(I2.2)') mday
        write (hourstr,  '(I2.2)') hour
        file_smap = trim(DEF_DA_obsdir) // '/SMAP_L1C_TB_D_' // &
            trim(yearstr) // '_' // trim(monthstr) // '_' // trim(daystr) // '_' // trim(hourstr) // '.nc'
        inquire (file=trim(file_smap), exist=has_file)

        ! whether have any obs in this time interval
        has_obs = .false.
        IF (has_file) THEN
            CALL ncio_read_bcast_serial (file_smap, 'time', smap_time)
            allocate ( dt_begin(size(smap_time)) )
            allocate ( dt_end  (size(smap_time)) )
            idate_begin = idate(3)
            idate_end = idate(3) + deltim
            dt_begin = smap_time - idate_begin
            dt_end = smap_time - idate_end 
            IF (any(dt_begin >= 0 .and. dt_end <= 0)) THEN
                has_obs = .true.
            ENDIF
            num_obs = size(smap_time)  ! all obs contains in current file (not constrain obs time)
        ELSE
            has_obs = .false.
        ENDIF

        ! allocate memory if have observations
        IF (has_obs) THEN
            IF (p_is_worker) THEN
                IF (numpatch > 0) THEN
                    ! observations
                    allocate (smap_lat   (num_obs))
                    allocate (smap_lon   (num_obs))
                    allocate (smap_tb_h  (num_obs))
                    allocate (smap_tb_v  (num_obs))
                    allocate (smap_ii    (num_obs))
                    allocate (smap_jj    (num_obs))
                    allocate (pred_tb_obs_ens_h  (num_obs, num_ens))
                    allocate (pred_tb_obs_ens_v  (num_obs, num_ens))
                    allocate (pred_tb_pset_ens_h (num_ens, numpatch))
                    allocate (pred_tb_pset_ens_v (num_ens, numpatch))
                ENDIF
            ENDIF
        ENDIF

    END SUBROUTINE allocate_SMAP



!-----------------------------------------------------------------------------

    SUBROUTINE run_DA_SMAP ()

!-----------------------------------------------------------------------------
    USE MOD_Spmd_task
    USE MOD_TimeManager
    USE MOD_NetCDFBlock
    USE MOD_Mesh
    USE MOD_LandElm
    USE MOD_LandPatch
    USE MOD_Vars_1DFluxes       
    USE MOD_Vars_TimeVariables
    USE MOD_Vars_TimeInvariants
    USE MOD_DA_Vars_TimeVariables
    USE MOD_RangeCheck 
    USE MOD_UserDefFun
    USE MOD_DA_EnKF
    IMPLICIT NONE

!------------------------ Local Variables ------------------------------
    real(r8) :: lat_p_n, lat_p_s, lon_p_w, lon_p_e
    integer :: ib, jb, il, jl, iens, iobs, np, i, n

    type(block_data_real8_3d) :: pred_tb_grid_ens_h
    type(block_data_real8_3d) :: pred_tb_grid_ens_v
    
!-----------------------------------------------------------------------------

        ! read obs 
        IF (has_obs) THEN
            
            CALL ncio_read_bcast_serial (file_smap, 'lat' , smap_lat)
            CALL ncio_read_bcast_serial (file_smap, 'lon' , smap_lon)
            CALL ncio_read_bcast_serial (file_smap, 'ii'  , smap_ii)
            CALL ncio_read_bcast_serial (file_smap, 'jj'  , smap_jj)
            CALL ncio_read_bcast_serial (file_smap, 'tb_h', smap_tb_h)
            CALL ncio_read_bcast_serial (file_smap, 'tb_v', smap_tb_v)

            ! forward model
            IF (p_is_worker) THEN
                DO iens = 1, num_ens
                    DO np = 1, numpatch
                        CALL forward ( &
                            patchtype(np), patchclass(np), lb_ens(iens,np), dz_soisno_ens(:,iens,np), &
                            forc_topo(np), &
                            tref_ens(iens,np), t_soisno_ens(:,iens,np), tleaf_ens(iens,np), &
                            wliq_soisno_ens(:,iens,np), wice_soisno_ens(:,iens,np), h2osoi_ens(:,iens,np), &
                            snowdp_ens(iens,np), lai_ens(iens,np), sai_ens(iens,np), &
                            vf_clay(:,np), vf_sand(:,np), BD_all(:,np), porsl(:,np), &  
                            pred_tb_pset_ens_h(iens,np), pred_tb_pset_ens_v(iens,np))
                    ENDDO
                ENDDO
            ENDIF
            CALL allocate_block_data (grid_smap, pred_tb_grid_ens_h, num_ens) 
            CALL mg2p_smap%pset2grid (pred_tb_pset_ens_h, pred_tb_grid_ens_h) ! map from patch to obs grid 

            ! crop the predicted obs to the same size as the obs
            IF (p_is_worker) THEN
                DO i = 1, num_obs
                    ib = grid_smap%xblk(smap_ii(i))
                    jb = grid_smap%yblk(smap_jj(i))
                    il = grid_smap%xloc(smap_ii(i))
                    jl = grid_smap%yloc(smap_jj(i)) 
                    pred_tb_obs_ens_h(i,:) = pred_tb_grid_ens_h%blk(ib, jb)%val(:, il, jl)
                    pred_tb_obs_ens_v(i,:) = pred_tb_grid_ens_v%blk(ib, jb)%val(:, il, jl)
                ENDDO
            ENDIF            
        
            ! data assimilation
            IF (p_is_worker) THEN
                DO np = 1, numpatch
                    ! find obs around 1 degree and time located in the interval
                    lat_p_n = patchlatr(np) + dres
                    lat_p_s = patchlatr(np) - dres
                    lon_p_w = patchlonr(np) - dres
                    lon_p_e = patchlonr(np) + dres

                    allocate (obs_index (0))
                    DO iobs = 1, num_obs
                        IF (smap_lat(iobs) < lat_p_n .and. smap_lat(iobs) > lat_p_s .and. &
                            smap_lon(iobs) > lon_p_w .and. smap_lon(iobs) < lon_p_e .and. &
                            smap_time(iobs) - idate_begin >= 0 .and. smap_time(iobs) - idate_end <= 0) THEN
                            
                            ! get the index of avaliable obs of target patch
                            n = size (obs_index)
                            IF (allocated(obs_index_tmp)) deallocate (obs_index_tmp)
                            allocate (obs_index_tmp(n+1))
                            obs_index_tmp(1:n) = obs_index(:)
                            obs_index_tmp(n+1) = iobs
                            IF (allocated(obs_index)) deallocate (obs_index)
                            allocate (obs_index(n+1))
                            obs_index(:) = obs_index_tmp(:)
                        ENDIF
                    ENDDO

                    ! crop avaliable obs and predicted obs 
                    num_loc_obs = size(obs_index)
                    allocate (obs_tb_h          (num_loc_obs))
                    allocate (obs_tb_v          (num_loc_obs))
                    allocate (obs_lat           (num_loc_obs))
                    allocate (obs_lon           (num_loc_obs))
                    allocate (obs_d             (num_loc_obs))
                    allocate (obs_err           (num_loc_obs))
                    allocate (pred_tb_loc_ens_h (num_loc_obs,num_ens))
                    allocate (pred_tb_loc_ens_v (num_loc_obs,num_ens))
                    obs_tb_h = smap_tb_h(obs_index)
                    obs_tb_v = smap_tb_v(obs_index)
                    obs_lat = smap_lat(obs_index)
                    obs_lon = smap_lon(obs_index)
                    obs_d = sqrt((obs_lat - patchlatr(np))**2 + (obs_lon - patchlonr(np))**2)
                    pred_tb_loc_ens_h = pred_tb_obs_ens_h(obs_index,:)
                    pred_tb_loc_ens_v = pred_tb_obs_ens_v(obs_index,:)
                    
                    ! data assimilation
                    obs_err(:) = 0.04   ! for SMAP
                    CALL letkf( &
                        num_ens, num_loc_obs, num_state, &
                        h2osoi_ens(:,:,np), pred_tb_loc_ens_h, obs_tb_h, obs_err, obs_d, loc_r, infl, &
                        h2osoi(:,np), trans, h2osoi_ens(:,:,np))
                ENDDO
            ENDIF
        ENDIF
        
    END SUBROUTINE run_DA_SMAP




!-----------------------------------------------------------------------------

    SUBROUTINE deallocate_SMAP ()

!-----------------------------------------------------------------------------
    IMPLICIT NONE

        ! time
        IF (allocated(smap_time))    deallocate(smap_time)
        IF (allocated(dt_begin))     deallocate(dt_begin)
        IF (allocated(dt_end))       deallocate(dt_end)

        ! obs
        IF (allocated(smap_lat))     deallocate(smap_lat)
        IF (allocated(smap_lon))     deallocate(smap_lon)
        IF (allocated(smap_tb_h))    deallocate(smap_tb_h)
        IF (allocated(smap_tb_v))    deallocate(smap_tb_v)
        IF (allocated(smap_ii))      deallocate(smap_ii)
        IF (allocated(smap_jj))      deallocate(smap_jj)
        IF (allocated(obs_index))    deallocate(obs_index)
        IF (allocated(obs_tb_h))     deallocate(obs_tb_h)
        IF (allocated(obs_tb_v))     deallocate(obs_tb_v)
        IF (allocated(obs_lat))      deallocate(obs_lat)
        IF (allocated(obs_lon))      deallocate(obs_lon)
        IF (allocated(obs_d))        deallocate(obs_d)
        IF (allocated(obs_err))      deallocate(obs_err)

        ! predicted
        IF (allocated(pred_tb_obs_ens_h))   deallocate(pred_tb_obs_ens_h)
        IF (allocated(pred_tb_obs_ens_v))   deallocate(pred_tb_obs_ens_v)
        IF (allocated(pred_tb_pset_ens_h))  deallocate(pred_tb_pset_ens_h)
        IF (allocated(pred_tb_pset_ens_v))  deallocate(pred_tb_pset_ens_v)
        IF (allocated(pred_tb_loc_ens_h))   deallocate(pred_tb_loc_ens_h)
        IF (allocated(pred_tb_loc_ens_v))   deallocate(pred_tb_loc_ens_v)
       
    END SUBROUTINE deallocate_SMAP


!-----------------------------------------------------------------------------
END MODULE MOD_DA_SMAP
#endif

    


