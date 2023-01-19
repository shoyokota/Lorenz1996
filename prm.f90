module prm
  implicit none
  public
  integer, parameter :: r_sngl = 4     ! default precision for Single
  integer, parameter :: r_size = 8     ! default precision for REAL
  real(r_size), parameter :: pi = 3.1415926535897932384626433832795_r_size ! Pi
  !
  type obsdata
    integer :: id                      ! observation ID
    real(r_size) :: pos                ! observation position
    real(r_size) :: val                ! observation value
    real(r_size) :: err                ! observation error
  end type obsdata
  !
  integer, save :: i_detda             ! deterministic DA method (0:NoDA,1:KF,2:EnVar,3:ETLMVar,4:4DVar)
  integer, save :: i_ensda             ! ensemble DA method (0:NoDA,1-4:EDA,5:EnSRF,6:LETKF,7:LPF)
  !
  integer, save :: nx                  ! number of Lorenz1996 grid point
  real(r_size), save :: forcing        ! forcing parameter for Lorenz1996 (nature run)
  real(r_size), save :: dt             ! integration time step
  real(r_size), save :: ft_nature_pre  ! forecast time for spinup of nature run
  real(r_size), save :: ft_nature      ! forecast time for nature run
  !
  integer, save :: nslot               ! number of time slot in assimilation window
  integer, save :: aslot               ! analysis time slot in assimilation window
  integer, save :: niter               ! maximum number of iteration
  real(r_size), save :: dt_cycle       ! assimilation time interval
  real(r_size), save :: ft_cycle_pre   ! time for spinup of assimilation
  real(r_size), save :: ft_cycle       ! time for assimilation
  real(r_size), save :: ft_max         ! max forecast time for each analysis
  !
  integer, save :: nobs                ! number of observations in each timeslot
  integer, save :: obsthin             ! interval of observation thinning
  integer, save :: obsid               ! observation ID
  real(r_size), save :: obserr         ! observation error
  !
  integer, save :: nmem                ! number of ensemble members
  real(r_size), save :: hyb_betab      ! weight of climatological background error covariance
  real(r_size), save :: hyb_betae      ! weight of ensemble-based background error covariance
  real(r_size), save :: amp_bcli       ! amplitude of climatological background error covariance
  real(r_size), save :: sigma_bcli     ! exp(-1/2) scale of climatological background error covariance
  real(r_size), save :: sigma_sloc     ! exp(-1/2) scale for spatial localization
  real(r_size), save :: sigma_tloc     ! exp(-1/2) scale for temporal localization
  real(r_size), save :: sigma_sloce    ! exp(-1/2) scale for spatial localization in ETLM
  real(r_size), save :: sigma_tloce    ! exp(-1/2) scale for temporal localization in ETLM
  real(r_size), save :: infl           ! multiplicative inflation
  real(r_size), save :: rtpp           ! relaxation to prior purturbations
  !
  character(len=64), save :: truthfile ! file of truth
  character(len=64), save :: obsfile   ! file of observation
  character(len=64), save :: detdafile ! file of deterministic DA
  character(len=64), save :: ensdafile ! file of ensemble DA
  character(len=64), save :: biasfile  ! file of BIAS
  character(len=64), save :: rmsefile  ! file of RMSE
  character(len=64), save :: sprdfile  ! file of ensemble spread
  !
contains
  !
  subroutine read_namelist
    implicit none
    !
    namelist/namprm/ &
         & i_detda, i_ensda, &
         & nx, forcing, dt, ft_nature_pre, ft_nature, &
         & nslot, aslot, niter, dt_cycle, ft_cycle_pre, ft_cycle, ft_max, &
         & nobs, obsthin, obsid, obserr, &
         & nmem, hyb_betab, hyb_betae, amp_bcli, sigma_bcli, &
         & sigma_sloc, sigma_tloc, sigma_sloce, sigma_tloce, infl, rtpp
    namelist/namfile/ &
         & truthfile, obsfile, detdafile, ensdafile, biasfile, rmsefile, sprdfile
    !
    i_detda       = 0
    i_ensda       = 0
    !
    nx            = 40
    forcing       = 8.0_r_size
    dt            = 0.01_r_size
    ft_nature_pre = 0.2_r_size * 365.0_r_size
    ft_nature     = 0.2_r_size * 100.0_r_size
    !
    nslot         = 1
    aslot         = 1
    niter         = 20
    dt_cycle      = 0.2_r_size * 0.25_r_size
    ft_cycle_pre  = 0.2_r_size * 100.0_r_size
    ft_cycle      = 0.2_r_size * 100.0_r_size
    ft_max        = 0.2_r_size * 5.0_r_size
    !
    nobs          = 40
    obsthin       = 1
    obsid         = 1
    obserr        = 1.0_r_size
    !
    nmem          = 3
    hyb_betab     = 1.0_r_size
    hyb_betae     = 0.0_r_size
    amp_bcli      = 0.5_r_size
    sigma_bcli    = 2.0_r_size
    sigma_sloc    = 2.0_r_size
    sigma_tloc    = 0.2_r_size * 0.25_r_size
    sigma_sloce   = 2.0_r_size
    sigma_tloce   = 0.2_r_size * 0.25_r_size
    infl          = 0.1_r_size
    rtpp          = 0.0_r_size
    !
    truthfile = 'TRUTH.txt'
    obsfile   = 'OBS.txt'
    detdafile = 'DETDA.txt'
    ensdafile = 'ENSDA.txt'
    biasfile  = 'BIAS.txt'
    rmsefile  = 'RMSE.txt'
    sprdfile  = 'SPRD.txt'
    !
    read(5,namprm)
    read(5,namfile)
    !
    return
  end subroutine read_namelist
  !
end module prm
