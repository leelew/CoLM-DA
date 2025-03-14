#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_SMAP
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Data assimilation of brightness temperature 
! 
! AUTHOR:
!   Lu Li, 12/2024: Initial version, based on SMAP L1C TB data
!   Zhilong Fan, Lu Li, 03/2024: Debug and clean codes
!-----------------------------------------------------------------------------
    USE MOD_DataType
    USE MOD_SpatialMapping
    USE MOD_DA_ObsOperator
    USE MOD_DA_EnKF
    USE MOD_DA_Vars_TimeVariables
    USE MOD_Vars_Global, only : num_ens, pi
    USE MOD_LandPatch
    IMPLICIT NONE
    SAVE

! public functions
    PUBLIC :: init_DA_SMAP
    PUBLIC :: run_DA_SMAP
    PUBLIC :: end_DA_SMAP

    ! use for save predicted obs and obs
    real(r8), allocatable, public :: pred_tb_h_out(:,:,:)    
    real(r8), allocatable, public :: smap_tb_h_out(:,:)       
    real(r8), allocatable, public :: pred_tb_v_out(:,:,:)    
    real(r8), allocatable, public :: smap_tb_v_out(:,:)       

    PRIVATE

! local variables
    ! file path
    character(len=256) :: file_smap                           ! SMAP file path
    character(len=256) :: file_grid                           ! SMAP world grid file path

    ! grid
    type(grid_type) :: grid_smap                              ! SMAP world grid
    type(spatial_mapping_type) :: mg2p_smap                   ! mapping between world grid to patch

    ! time (UTC) at current step
    integer :: month, mday, hour                              ! month, day, hour of current step
    character(len=256) :: yearstr, monthstr, daystr, hourstr  ! string of year, month, day, hour
    integer :: idate_b, idate_e                               ! begin & end seconds since begin of current day (UTC)

    ! time variables used to determine whether has obs
    real(r8), allocatable :: smap_time (:)                    ! seconds of all obs since begin of current day (UTC)
    real(r8), allocatable :: dt_b      (:)                    ! delta time between obs and begin seconds of current day
    real(r8), allocatable :: dt_e      (:)                    ! delta time between obs and end seconds of current day
    
    ! logical variables
    logical :: has_file                                       ! whether has file of SMAP at target day
    logical :: has_obs                                        ! whether has obs at current step

    ! observations (dimensions changes with time)
    integer :: num_obs                                        ! number of all obs in current file
    real(r8), allocatable :: smap_lat  (:)                    ! latitude of all obs 
    real(r8), allocatable :: smap_lon  (:)                    ! longitude of all obs
    real(r8), allocatable :: smap_tb_h (:)                    ! H- polarized brightness temperature of all obs ([K])
    real(r8), allocatable :: smap_tb_v (:)                    ! V- polarized brightness temperature of all obs ([K])
    integer,  allocatable :: smap_ii   (:)                    ! i-th lat grid in world map of all obs
    integer,  allocatable :: smap_jj   (:)                    ! j-th lon grid in world map of all obs

    ! ensemble predicted observations (at patch)
    real(r8), allocatable :: pred_tb_h_pset_ens (:,:)         ! predicted H- polarized temp on patch
    real(r8), allocatable :: pred_tb_v_pset_ens (:,:)         ! predicted V- polarized temp on patch
    real(r8), allocatable :: area_pset          (:)           ! area of each patch across world grid

    ! ensemble predicted observations (at world grid)
    type(block_data_real8_3d) :: pred_tb_h_wgrid_ens          ! predicted H- polarized temp on world grid
    type(block_data_real8_3d) :: pred_tb_v_wgrid_ens          ! predicted V- polarized temp on world grid
    type(block_data_real8_2d) :: area_wgrid                   ! area of each patch across world grid

    ! ensemble predicted observations (at obs grid)
    real(r8), allocatable :: pred_tb_h_ogrid_ens (:,:)        ! predicted H- polarized temp on obs grid
    real(r8), allocatable :: pred_tb_v_ogrid_ens (:,:)        ! predicted V- polarized temp on obs grid
    real(r8) :: area_wgrid_obs                                ! area of each patch across world grid for each obs

    ! observations around patch (dimensions changes with patch)  
    integer :: num_obs_loc                                    ! number of obs around each patch
    real(r8), parameter   :: dres = 2                         ! search localization radius (deg)
    integer,  allocatable :: obs_index     (:)                ! index of obs around each patch  
    integer,  allocatable :: obs_index_tmp (:)                ! temporary index of obs around each patch
    real(r8), allocatable :: smap_lat_loc  (:)                ! latitude of obs around each patch
    real(r8), allocatable :: smap_lon_loc  (:)                ! longitude of obs around each patch
    real(r8), allocatable :: smap_tb_h_loc (:)                ! H- polarized brightness temperature of obs around each patch ([K])
    real(r8), allocatable :: smap_tb_v_loc (:)                ! V- polarized brightness temperature of obs around each patch ([K])
    real(r8), allocatable :: d_loc         (:)                ! distance between obs and patch center

    ! ensemble predicted observations around patch
    real(r8), allocatable :: pred_tb_h_loc_ens (:,:)          ! predicted H- polarized temp around patch
    real(r8), allocatable :: pred_tb_v_loc_ens (:,:)          ! predicted V- polarized temp around patch

    ! observations around patch (after removing NaN)
    integer :: num_obs_loc_nonan                              ! number of obs around each patch after removing NaN
    real(r8), allocatable :: smap_tb_h_loc_nonan (:)          ! H- polarized brightness temperature of obs around each patch after removing NaN ([K])
    real(r8), allocatable :: smap_tb_v_loc_nonan (:)          ! V- polarized brightness temperature of obs around each patch after removing NaN ([K])
    real(r8), allocatable :: d_loc_nonan         (:)          ! distance between obs and patch center after removing NaN

    ! ensemble predicted observations around patch (after removing NaN)
    real(r8), allocatable :: pred_tb_h_loc_nonan_ens (:,:)    ! predicted H- polarized temp around patch after removing NaN
    real(r8), allocatable :: pred_tb_v_loc_nonan_ens (:,:)    ! predicted V- polarized temp around patch after removing NaN

    ! data assimilation
    real(r8), allocatable :: obs_err (:)                      ! observation error
    integer, parameter    :: num_state = 10                   ! number of state
    real(r8), parameter   :: loc_r = 1.0                      ! localization radius 
    real(r8), parameter   :: infl = 1.2                       ! inflation factor
    real(r8) :: trans(num_ens, num_ens)                       ! transform matrix on each patch

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

    SUBROUTINE init_DA_SMAP ()

!-----------------------------------------------------------------------------
    USE MOD_Spmd_Task
    USE MOD_Namelist, only : DEF_DA_obsdir
    USE MOD_Grid
    USE MOD_NetCDFSerial
    USE MOD_LandPatch
    USE MOD_Pixelset
    USE MOD_RangeCheck
    IMPLICIT NONE

!-----------------------------------------------------------------------

        ! grid file path of EASE v2.0, 36km world grid
        file_grid = trim(DEF_DA_obsdir) // '/NSIDC0772_LatLon_EASE2_M36km_v1.0_pre.nc'

        ! read grid from file
        CALL grid_smap%define_from_file (file_grid, 'latitude', 'longitude')

        ! map SMAP grid to patch
        CALL mg2p_smap%build_arealweighted (grid_smap, landpatch)

    END SUBROUTINE init_DA_SMAP



!-----------------------------------------------------------------------------

    SUBROUTINE run_DA_SMAP (idate, deltim)

!-----------------------------------------------------------------------------
    USE MOD_Spmd_task
    USE MOD_TimeManager
    USE MOD_NetCDFBlock
    USE MOD_Mesh
    USE MOD_LandElm
    USE MOD_LandPatch
    USE MOD_Vars_1DFluxes   
    USE MOD_Vars_1DForcing    
    USE MOD_Vars_TimeVariables
    USE MOD_Vars_TimeInvariants
    USE MOD_DA_Vars_TimeVariables
    USE MOD_RangeCheck 
    USE MOD_UserDefFun
    USE MOD_DA_EnKF
    IMPLICIT NONE

!------------------------ Dummy Arguments ---------------------------------
    integer,  intent(in) :: idate(3)
    real(r8), intent(in) :: deltim   

!------------------------ Local Variables ------------------------------
    real(r8) :: lat_p_n, lat_p_s, lon_p_w, lon_p_e
    integer :: ib, jb, il, jl, iens, iobs, np, i, n

!-----------------------------------------------------------------------------

        ! calculate year/month/day/hour of current step
        CALL julian2monthday(idate(1), idate(2), month, mday)
        hour = int(idate(3)/3600)

        ! whether has file of SMAP at target day 
        !//TODO: Lu Li: support ascending orbits
        write (yearstr,  '(I4.4)') idate(1)
        write (monthstr, '(I2.2)') month
        write (daystr,   '(I2.2)') mday
        write (hourstr,  '(I2.2)') hour
        file_smap = trim(DEF_DA_obsdir) // '/SMAP_L1C_TB_D_' // &
            trim(yearstr) // '_' // trim(monthstr) // '_' // trim(daystr) // '_' // trim(hourstr) // '.nc'
        inquire (file=trim(file_smap), exist=has_file)

        ! whether have obs at this time interval
        has_obs = .false.
        IF (has_file) THEN
            CALL ncio_read_bcast_serial (file_smap, 'time', smap_time)
            num_obs = size(smap_time)
            idate_b = idate(3)
            idate_e = idate(3) + deltim
            allocate (dt_b(num_obs))
            allocate (dt_e(num_obs))
            dt_b = smap_time - idate_b
            dt_e = smap_time - idate_e 
            IF (any(dt_b >= 0 .and. dt_e <= 0)) has_obs = .true.
            deallocate (dt_b)
            deallocate (dt_e)
        ELSE
            has_obs = .false.
        ENDIF

        ! allocate memory if have observations
        IF (has_obs) THEN
            IF (p_is_worker) THEN
                IF (numpatch > 0) THEN
                    IF (allocated (smap_lat ))           deallocate (smap_lat )
                    IF (allocated (smap_lon ))           deallocate (smap_lon )
                    IF (allocated (smap_tb_h))           deallocate (smap_tb_h)
                    IF (allocated (smap_tb_v))           deallocate (smap_tb_v)
                    IF (allocated (smap_ii  ))           deallocate (smap_ii  )
                    IF (allocated (smap_jj  ))           deallocate (smap_jj  )
                    IF (allocated (pred_tb_h_pset_ens))  deallocate (pred_tb_h_pset_ens)
                    IF (allocated (pred_tb_v_pset_ens))  deallocate (pred_tb_v_pset_ens)
                    IF (allocated (pred_tb_h_ogrid_ens)) deallocate (pred_tb_h_ogrid_ens)
                    IF (allocated (pred_tb_v_ogrid_ens)) deallocate (pred_tb_v_ogrid_ens)
                    IF (allocated (area_pset))           deallocate (area_pset)

                    allocate (smap_lat   (num_obs))
                    allocate (smap_lon   (num_obs))
                    allocate (smap_tb_h  (num_obs))
                    allocate (smap_tb_v  (num_obs))
                    allocate (smap_ii    (num_obs))
                    allocate (smap_jj    (num_obs))
                    allocate (pred_tb_h_pset_ens (num_ens, numpatch))
                    allocate (pred_tb_v_pset_ens (num_ens, numpatch))
                    allocate (pred_tb_h_ogrid_ens (num_obs, num_ens))
                    allocate (pred_tb_v_ogrid_ens (num_obs, num_ens))
                    allocate (area_pset (numpatch))
                ENDIF
            ENDIF
        ENDIF

        ! read obs 
        IF (has_obs) THEN
            ! read observations
            CALL ncio_read_bcast_serial (file_smap, 'lat' , smap_lat )
            CALL ncio_read_bcast_serial (file_smap, 'lon' , smap_lon )
            CALL ncio_read_bcast_serial (file_smap, 'ii'  , smap_ii  )
            CALL ncio_read_bcast_serial (file_smap, 'jj'  , smap_jj  )  
            CALL ncio_read_bcast_serial (file_smap, 'tb_h', smap_tb_h)  
            CALL ncio_read_bcast_serial (file_smap, 'tb_v', smap_tb_v)

            ! forward model
            IF (p_is_worker) THEN
                print *, month, mday, hour
                DO iens = 1, num_ens
                    DO np = 1, numpatch
                        CALL forward ( &
                            patchtype(np), patchclass(np), lb_ens(iens,np), dz_soisno_ens(:,iens,np), &
                            forc_topo(np), &
                            tref_ens(iens,np), t_soisno_ens(:,iens,np), tleaf_ens(iens,np), &
                            wliq_soisno_ens(:,iens,np), wice_soisno_ens(:,iens,np), h2osoi_ens(:,iens,np), &
                            snowdp_ens(iens,np), lai_ens(iens,np), sai_ens(iens,np), &
                            vf_clay(:,np), vf_sand(:,np), BD_all(:,np), porsl(:,np), &  
                            pred_tb_h_pset_ens(iens,np), pred_tb_v_pset_ens(iens,np))
                    ENDDO
                ENDDO
            ENDIF

            ! mapping predicted observations from patch to world grid
            CALL allocate_block_data (grid_smap, pred_tb_h_wgrid_ens, num_ens) 
            CALL allocate_block_data (grid_smap, pred_tb_v_wgrid_ens, num_ens)
            CALL mg2p_smap%pset2grid (pred_tb_h_pset_ens, pred_tb_h_wgrid_ens, spval)     
            CALL mg2p_smap%pset2grid (pred_tb_v_pset_ens, pred_tb_v_wgrid_ens, spval)     

            ! calculate area of each patch across world grid
            area_pset(:) = 1
            CALL allocate_block_data (grid_smap, area_wgrid)
            CALL mg2p_smap%pset2grid (area_pset, area_wgrid)

            ! crop the predicted observations from world grid to all obs grids
            IF (p_is_worker) THEN
                DO i = 1, num_obs
                    ib = grid_smap%xblk(smap_jj(i)+1) !//TODO: Lu Li: python index starts from 0
                    jb = grid_smap%yblk(smap_ii(i)+1)
                    il = grid_smap%xloc(smap_jj(i)+1)
                    jl = grid_smap%yloc(smap_ii(i)+1)
                    IF (ib == 0 .or. jb == 0) THEN
                        !//TODO: Lu Li: add explainations
                        pred_tb_h_ogrid_ens(i,:) = -9999.0
                        pred_tb_v_ogrid_ens(i,:) = -9999.0
                    ELSE
                        area_wgrid_obs = area_wgrid%blk(ib, jb)%val(il, jl)
                        IF (area_wgrid_obs == 0) THEN
                            pred_tb_h_ogrid_ens(i,:) = -9999.0
                            pred_tb_v_ogrid_ens(i,:) = -9999.0
                        ELSE
                            pred_tb_h_ogrid_ens(i,:) = (pred_tb_h_wgrid_ens%blk(ib, jb)%val(:, il, jl))/area_wgrid_obs
                            pred_tb_v_ogrid_ens(i,:) = (pred_tb_v_wgrid_ens%blk(ib, jb)%val(:, il, jl))/area_wgrid_obs
                        ENDIF
                    ENDIF
                ENDDO
            ENDIF
        
            ! data assimilation
            IF (p_is_worker) THEN
                !//TODO: the maximum number of obs around each patch depends on the size of selected regions
                allocate(pred_tb_h_out(5, num_ens, numpatch))   
                allocate(smap_tb_h_out(5, numpatch))             
                allocate(pred_tb_v_out(5, num_ens, numpatch))   
                allocate(smap_tb_v_out(5, numpatch))            
                pred_tb_h_out = -1.0d0
                smap_tb_h_out = -1.0d0
                pred_tb_v_out = -1.0d0
                smap_tb_v_out = -1.0d0

                DO np = 1, numpatch
                    ! regions info around target patch
                    lat_p_n = patchlatr(np)*180/pi + dres
                    lat_p_s = patchlatr(np)*180/pi - dres
                    lon_p_w = patchlonr(np)*180/pi - dres
                    lon_p_e = patchlonr(np)*180/pi + dres

                    ! get the index of obs around target patch and time at current step
                    allocate (obs_index (0))             
                    DO iobs = 1, num_obs
                        IF (smap_lat(iobs) < lat_p_n       .and. smap_lat(iobs) > lat_p_s       .and. &
                            smap_lon(iobs) > lon_p_w       .and. smap_lon(iobs) < lon_p_e       .and. &
                            smap_time(iobs) - idate_b >= 0 .and. smap_time(iobs) - idate_e <= 0) THEN
                            
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
                    num_obs_loc = size(obs_index)
                    allocate (smap_lat_loc      (num_obs_loc))
                    allocate (smap_lon_loc      (num_obs_loc))
                    allocate (smap_tb_h_loc     (num_obs_loc))
                    allocate (smap_tb_v_loc     (num_obs_loc))
                    allocate (pred_tb_h_loc_ens (num_obs_loc,num_ens))
                    allocate (pred_tb_v_loc_ens (num_obs_loc,num_ens))
                    allocate (d_loc             (num_obs_loc))
                    smap_tb_h_loc = smap_tb_h(obs_index)
                    smap_tb_v_loc = smap_tb_v(obs_index)
                    smap_lat_loc = smap_lat(obs_index)
                    smap_lon_loc = smap_lon(obs_index)
                    pred_tb_h_loc_ens = pred_tb_h_ogrid_ens(obs_index,:)
                    pred_tb_v_loc_ens = pred_tb_v_ogrid_ens(obs_index,:)

                    ! calculate distance between obs and patch center using haversine formula
                    d_loc = 2*6.3781e3*asin(sqrt(sin((smap_lat_loc*pi/180 - patchlatr(np))/2.0)**2 + &
                        cos(smap_lat_loc*pi/180)*cos(patchlatr(np))*sin((smap_lon_loc*pi/180 - patchlonr(np))/2.0)**2))
                    
                    ! remove NaN obs
                    num_obs_loc_nonan = count(pred_tb_h_loc_ens(:,1) > 0)
                    IF (num_obs_loc_nonan > 0) THEN
                        allocate (smap_tb_h_loc_nonan     (num_obs_loc_nonan))
                        allocate (smap_tb_v_loc_nonan     (num_obs_loc_nonan))
                        allocate (pred_tb_h_loc_nonan_ens (num_obs_loc_nonan, num_ens))
                        allocate (pred_tb_v_loc_nonan_ens (num_obs_loc_nonan, num_ens))
                        allocate (d_loc_nonan             (num_obs_loc_nonan))
                        allocate (obs_err                 (num_obs_loc_nonan))
                        smap_tb_h_loc_nonan = pack(smap_tb_h_loc, pred_tb_h_loc_ens(:,1) > 0)
                        smap_tb_v_loc_nonan = pack(smap_tb_v_loc, pred_tb_h_loc_ens(:,1) > 0)
                        DO i = 1, num_ens
                            pred_tb_h_loc_nonan_ens(:,i) = pack(pred_tb_h_loc_ens(:,i), pred_tb_h_loc_ens(:,1) > 0)
                            pred_tb_v_loc_nonan_ens(:,i) = pack(pred_tb_v_loc_ens(:,i), pred_tb_v_loc_ens(:,1) > 0)
                        ENDDO
                        d_loc_nonan = pack(d_loc, pred_tb_h_loc_ens(:,1) > 0)
                        obs_err(:) = 0.04   !//TODO: Lu Li: only for SMAP data

                        ! data assimilation
                        print *, 'Number of obs around patch ', np, ' is ', num_obs_loc_nonan
                        CALL letkf( &
                            num_ens, num_obs_loc_nonan, num_state, &
                            h2osoi_ens(:,:,np), pred_tb_h_loc_nonan_ens, smap_tb_h_loc_nonan, &
                            obs_err, d_loc_nonan, loc_r, infl, &
                            h2osoi(:,np), trans, h2osoi_ens(:,:,np))
        
                        IF (wliq_soisno_ens(1, 1, np) /= spval) then
                            DO i = 1, nl_soil
                                wliq_soisno_ens(i, :, np) = matmul(trans, wliq_soisno_ens(i, :, np))
                                wice_soisno_ens(i, :, np) = matmul(trans, wice_soisno_ens(i, :, np))
                            ENDDO
                        ENDIF
                    ELSE
                        print *, 'No avaliable obs for patch ', np
                    ENDIF

                    ! save predicted obs and obs
                    IF (num_obs_loc_nonan > 0) THEN    
                        pred_tb_h_out(1:num_obs_loc_nonan, :, np) = pred_tb_h_loc_nonan_ens(:, :)
                        smap_tb_h_out(1:num_obs_loc_nonan, np)    = smap_tb_h_loc_nonan(:)
                        pred_tb_v_out(1:num_obs_loc_nonan, :, np) = pred_tb_v_loc_nonan_ens(:, :)
                        smap_tb_v_out(1:num_obs_loc_nonan, np)    = smap_tb_v_loc_nonan(:)
                    ENDIF

                    !//TODO: Lu Li: adjust land surface variables based on the updated transform matrix

                    ! deallocate memory (cuz dimensions changes with patch)
                    IF (allocated(obs_index))               deallocate(obs_index)
                    IF (allocated(obs_index_tmp))           deallocate(obs_index_tmp)
                    IF (allocated(smap_lat_loc))            deallocate(smap_lat_loc)
                    IF (allocated(smap_lon_loc))            deallocate(smap_lon_loc)
                    IF (allocated(smap_tb_h_loc))           deallocate(smap_tb_h_loc)
                    IF (allocated(smap_tb_v_loc))           deallocate(smap_tb_v_loc)
                    IF (allocated(pred_tb_h_loc_ens))       deallocate(pred_tb_h_loc_ens)
                    IF (allocated(pred_tb_v_loc_ens))       deallocate(pred_tb_v_loc_ens)
                    IF (allocated(d_loc))                   deallocate(d_loc)
                    IF (allocated(smap_tb_h_loc_nonan))     deallocate(smap_tb_h_loc_nonan)
                    IF (allocated(smap_tb_v_loc_nonan))     deallocate(smap_tb_v_loc_nonan)
                    IF (allocated(pred_tb_h_loc_nonan_ens)) deallocate(pred_tb_h_loc_nonan_ens)
                    IF (allocated(pred_tb_v_loc_nonan_ens)) deallocate(pred_tb_v_loc_nonan_ens)
                    IF (allocated(d_loc_nonan))             deallocate(d_loc_nonan)
                    IF (allocated(obs_err))                 deallocate(obs_err)
                ENDDO
            ENDIF
        ENDIF
        
    END SUBROUTINE run_DA_SMAP

!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
    SUBROUTINE end_DA_SMAP ()
    IMPLICIT NONE

        IF (allocated(smap_time))               deallocate(smap_time)
        IF (allocated(dt_b))                    deallocate(dt_b)
        IF (allocated(dt_e))                    deallocate(dt_e)
        IF (allocated(smap_lat))                deallocate(smap_lat)
        IF (allocated(smap_lon))                deallocate(smap_lon)
        IF (allocated(smap_tb_h))               deallocate(smap_tb_h)
        IF (allocated(smap_tb_v))               deallocate(smap_tb_v)
        IF (allocated(smap_ii))                 deallocate(smap_ii)
        IF (allocated(smap_jj))                 deallocate(smap_jj)
        IF (allocated(pred_tb_h_pset_ens))      deallocate(pred_tb_h_pset_ens)
        IF (allocated(pred_tb_v_pset_ens))      deallocate(pred_tb_v_pset_ens)
        IF (allocated(area_pset))               deallocate(area_pset)
        IF (allocated(pred_tb_h_ogrid_ens))     deallocate(pred_tb_h_ogrid_ens)
        IF (allocated(pred_tb_v_ogrid_ens))     deallocate(pred_tb_v_ogrid_ens)
        IF (allocated(obs_index))               deallocate(obs_index)
        IF (allocated(obs_index_tmp))           deallocate(obs_index_tmp)
        IF (allocated(smap_lat_loc))            deallocate(smap_lat_loc)
        IF (allocated(smap_lon_loc))            deallocate(smap_lon_loc)
        IF (allocated(smap_tb_h_loc))           deallocate(smap_tb_h_loc)
        IF (allocated(smap_tb_v_loc))           deallocate(smap_tb_v_loc)
        IF (allocated(d_loc))                   deallocate(d_loc)
        IF (allocated(pred_tb_h_loc_ens))       deallocate(pred_tb_h_loc_ens)
        IF (allocated(pred_tb_v_loc_ens))       deallocate(pred_tb_v_loc_ens)
        IF (allocated(smap_tb_h_loc_nonan))     deallocate(smap_tb_h_loc_nonan)
        IF (allocated(smap_tb_v_loc_nonan))     deallocate(smap_tb_v_loc_nonan)
        IF (allocated(d_loc_nonan))             deallocate(d_loc_nonan)
        IF (allocated(pred_tb_h_loc_nonan_ens)) deallocate(pred_tb_h_loc_nonan_ens)
        IF (allocated(pred_tb_v_loc_nonan_ens)) deallocate(pred_tb_v_loc_nonan_ens)
        IF (allocated(obs_err))                 deallocate(obs_err)
       
    END SUBROUTINE end_DA_SMAP


!-----------------------------------------------------------------------------
END MODULE MOD_DA_SMAP
#endif

    


