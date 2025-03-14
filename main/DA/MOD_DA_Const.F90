#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Const
!-----------------------------------------------------------------------
! DESCRIPTION:
! Build look up table for constant used in radiative transfer model
! Convert CoLM land classification to ECOCLIMAP land classification (from L-MEB)
!
! REFERENCES:
!   [1] Drusch, M., Holmes, T., de Rosnay, P., Balsamo, G., 2009. 
!       Comparing ERA-40 based L-band brightness temperatures with Skylab observations:
!       a calibration/validationstudy using the Community Microwave Emission Model. J. Hydrometeorol.
!       doi: 10.1175/2008JHM964.1

!   [2] de Rosnay, P., Drusch, M., Boone, A., Balsamo, G., Decharme, B., Harris, P., Kerr, 
!       Y.,Pellarin, T., Polcher, J., Wigneron, J.-P., 2009a. Microwave land surface modellingevaluation 
!       against AMSR-E data over West Africa. The AMMA land surface modelintercomparison
!       experiment coupled to the community microwave emission model(ALMIP-MEM). 
!       J. Geophys. Res. 114. doi: 10.1029/2008JD010724

!   [3] Wang, S.; Wigneron, J.-P.; Jiang, L.M.; Parrens, M. Global-scale evaluation of roughness 
!       effects on C-Band AMSR-E observations. Remote Sens. 2015, 7, 5734–5757 doi: 10.3390/rs70505734
!
!   [4] Wigneron, J. P., Jackson, T. J., O'neill, P., De Lannoy, G., de Rosnay, P., Walker, 
!       J. P., ... & Kerr, Y. (2017). Modelling the passive microwave signature from land surfaces: 
!       A review of recent results and application to the L-band SMOS & SMAP soil moisture retrieval algorithms. 
!       Remote Sensing of Environment, 192, 238-262.
!
! AUTHOR:
!   Lu Li, 12/2024: Initial version
!-----------------------------------------------------------------------

    USE MOD_Precision
    USE MOD_Vars_Global, only : pi, N_land_classification

    IMPLICIT NONE
    SAVE

    ! Constant variables
    real(r8),    parameter :: theta = 40.0*pi/180   ! incidence angle (rad)
    real(r8),    parameter :: C = 2.998e8           ! speed of light (m/s)
    real(r8),    parameter :: fghz = 1.4            ! frequency (GHz)
    real(r8),    parameter :: f = 1.4e9             ! frequency (Hz)
    real(r8),    parameter :: omega = 2*pi*f        ! radian frequency (omega = 2. * pi * f), with f in Herz
    real(r8),    parameter :: mu0 = 4.*pi*1e-7       ! vacuum permeability (H/m)
    real(r8),    parameter :: eps0 = 8.854e-12      ! vacuum permittivity (Klein and Swift 1977) [Farads/meter]
    real(r8),    parameter :: z0 = sqrt(mu0/eps0)   ! impendace of free space (Ohm)
    real(r8),    parameter :: eps_w_inf = 4.9       ! dielectric constant at infinite frequency (Stogryn 1971),
    real(r8),    parameter :: eps_0 = 8.854e-12     ! dielectric constant of free space (Klein and Swift 1977) [Farads/meter]
    real(r8),    parameter :: rho_soil = 2.66       ! soil specific density (g/cm3)
    real(r8),    parameter :: f0w = 9.              ! relaxation frequency of liquid water (GHz)
    real(r8),    parameter :: lam = C/f             ! wavelength (m)
    real(r8),    parameter :: k = 2*pi/lam          ! wave number (rad/m)
    complex(r8), parameter :: jj = (0., 1.)         ! imaginary unit for complex number

    
!//TODO (Lu Li): support other land cover classification system
#ifdef LULC_IGBP 

    ! MODIS IGBP Land Use/Land Cover System Legend
    !---------------------------
    ! 0  Ocean
    ! 1  Evergreen Needleleaf Forests
    ! 2  Evergreen Broadleaf Forests
    ! 3  Deciduous Needleleaf Forests
    ! 4  Deciduous Broadleaf Forests
    ! 5  Mixed Forests
    ! 6  Closed Shrublands
    ! 7  Open Shrublands
    ! 8  Woody Savannas
    ! 9  Savannas
    !10  Grasslands
    !11  Permanent Wetlands
    !12  Croplands
    !13  Urban and Built-up Lands
    !14  Cropland/Natural Vegetation Mosaics
    !15  Permanent Snow and Ice
    !16  Barren
    !17  Water Bodies

    ! ECOCLIMAP Land Use/Land Cover System Legend
    !---------------------------
    ! 0 bare soil
    ! 1 decidious forests
    ! 2 coniferous forests
    ! 3 rain forests 
    ! 4 C3 grasslands
    ! 5 C4 grasslands
    ! 6 C3 crops
    ! 7 C4 crops

    ! index map IGBP to ECOCLIMAP
    integer, parameter, dimension(N_land_classification) :: igbp2eco &
        = (/2, 3, 2, 1, 1, 4, 4, 4, 4, 4, 0, 6, 0, 6, 0, 0, 0/)
    
    ! b parameters for Wigneron vegetation model
    real(r8), parameter, dimension(N_land_classification) :: b1 &
        = (/0.260, 0.226, 0.260, 0.226, 0.226, &
            0.0375, 0.0375, 0.0375, 0.0375, 0.0375, &
            0.0, 0.05, 0.0, 0.05, 0.0, 0.0, 0.0/)
    real(r8), parameter, dimension(N_land_classification) :: b2 &
        = (/0.006, 0.001, 0.006, 0.001, 0.001, &
            0.05, 0.05, 0.05, 0.05, 0.05, &
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    real(r8), parameter, dimension(N_land_classification) :: b3 &
        = (/0.69, 0.7, 0.69, 0.7, 0.7, &
            0.0, 0.0, 0.0, 0.0, 0.0, &
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)

    ! soil roughness parameters for H and V polarization
    ! (NOTE): we follow SMAP L2 algorithm (all 2.0)
    real(r8), parameter, dimension(N_land_classification) :: nrh &
        = (/2.00, 2.00, 2.00, 2.00, 2.00, &
            2.00, 2.00, 2.00, 2.00, 2.00, &
            2.00, 2.00, 2.00, 2.00, 2.00, &
            2.00, 2.00/)
    real(r8), parameter, dimension(N_land_classification) :: nrv &
        = (/2.00, 2.00, 2.00, 2.00, 2.00, &
            2.00, 2.00, 2.00, 2.00, 2.00, &
            2.00, 2.00, 2.00, 2.00, 2.00, &
            2.00, 2.00/)

    ! empirical parameters to account for incidence angle 
    ! in vegetation opacity calculation for H polarization
    ! (NOTE): from L-MEB model
    real(r8), parameter, dimension(N_land_classification) :: tth &
        = (/0.8, 1.0, 0.8, 0.49, 0.49, &
            1.0, 1.0, 1.0, 1.0, 1.0, &
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)
    real(r8), parameter, dimension(N_land_classification) :: ttv &
        = (/0.8, 1.0, 0.8, 0.46, 0.46, &
            1.0, 1.0, 1.0, 1.0, 1.0, &
            1.0, 2.0, 1.0, 2.0, 1.0, 1.0, 1.0/)
    
    ! empirical roughness parameters (Table 2 in Wigneron et al. 2017)
    real(r8), parameter, dimension(N_land_classification) :: hr &
        = (/0.160, 0.160, 0.160, 0.160, 0.160, &
            0.110, 0.110, 0.125, 0.156, 0.156, &
            0.100, 0.108, 0.000, 0.130, 0.000, &
            0.150, 0.000/)

    ! effective diffusion albedo (Table 3 in Wigneron et al. 2017)
    real(r8), parameter, dimension(N_land_classification) :: w &
        = (/0.050, 0.050, 0.050, 0.050, 0.050, &
            0.050, 0.050, 0.050, 0.080, 0.050, &
            0.050, 0.000, 0.065, 0.000, 0.000, &
            0.000, 0.000/)
#endif

END MODULE MOD_DA_Const
!-----------------------------------------------------------------------------
#endif
