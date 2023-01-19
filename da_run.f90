program da_run
  use prm, only : r_size, &
       & i_detda, i_ensda, &
       & nx, &
       & nslot, dt_cycle, ft_cycle_pre, ft_cycle, &
       & nobs, obsdata, &
       & nmem, amp_bcli, sigma_bcli, infl, rtpp, &
       & truthfile, obsfile, biasfile, rmsefile, sprdfile
  use prm, only : read_namelist
  use io, only : open_data, read_grid, read_obs, write_bias, write_rmse, write_sprd
  use random, only : init_random, get_gauss_random
  use lorenz1996, only : lorenz1996_run
  use matrix, only : matrix_gauss
  use kf, only : kf_run
  use envar, only : envar_run
  use etlmvar, only : etlmvar_run
  use fdvar, only : fdvar_run
  use ensrf, only : ensrf_run
!  use letkf, only : letkf_run
!  use lpf, only : lpf_run
  implicit none
  real(r_size),  allocatable :: x(:), xf(:), xt(:)
  real(r_size),  allocatable :: xe(:,:), xfe(:,:)
  real(r_size),  allocatable :: pb(:,:), pf(:,:)
  type(obsdata), allocatable :: obs(:,:), obsf(:,:)
  integer :: i, ii, j, k, kk
  integer :: iu_grid, iu_obs, iu_bias, iu_rmse, iu_sprd
  real(r_size) :: time
  !
  call read_namelist
  ! initial state
  allocate( x   (1:nx)           )
  allocate( xf  (1:nx)           )
  allocate( xt  (1:nx)           )
  allocate( xe  (1:nx,1:nmem)    )
  allocate( xfe (1:nx,1:nmem)    )
  allocate( pb  (1:nx,1:nx)      )
  allocate( pf  (1:nx,1:nx)      )
  allocate( obs (1:nslot,1:nobs) )
  allocate( obsf(1:nslot,1:nobs) )
  call init_random( 0 )
  call get_gauss_random( nx, x(1:nx) )
  call lorenz1996_run( 1, 1, ft_cycle_pre, x(1:nx) )
  do j = 1, nmem
     call get_gauss_random( nx, xe(1:nx,j) )
     xe(1:nx,j) = xe(1:nx,j) + x(1:nx)
  end do
  call matrix_gauss( nx, sigma_bcli, pb(1:nx,1:nx) )
  pb(1:nx,1:nx) = amp_bcli * pb(1:nx,1:nx)
  ! forecast-analysis cycle
  call open_data( iu_grid, truthfile )
  call open_data( iu_obs,  obsfile   )
  call open_data( iu_bias, biasfile  )
  call open_data( iu_rmse, rmsefile  )
  call open_data( iu_sprd, sprdfile  )
  cycle:do k = 1, nint(ft_cycle/dt_cycle)
     ! verification of deterministic forecast
     xf(1:nx) = x(1:nx)
     time = dt_cycle * ( k - 1 )
     do kk = 1, nslot
        call lorenz1996_run( 1, kk, dt_cycle/nslot, xf(1:nx) )
        time = time + dt_cycle / nslot * kk
        call read_grid ( iu_grid, time, xt(1:nx)           )
        call write_bias( iu_bias, time, xf(1:nx), xt(1:nx) )
        call write_rmse( iu_rmse, time, xf(1:nx), xt(1:nx) )
     end do
     ! read observations
     do kk = 1, nslot
        call read_obs( iu_obs, time, obs(kk,1:nobs) )
     end do
     if( i_detda >= 1 ) then
        ! deterministic forecast and analysis
        if( i_detda == 1 ) then      ! KF
           call kf_run( x(1:nx), pb(1:nx,1:nx), obs(1:nslot,1:nobs) )
        else if( i_detda == 2 ) then ! EnVar
           call envar_run( x(1:nx), xe(1:nx,1:nmem), pb(1:nx,1:nx), obs(1:nslot,1:nobs) )
        else if( i_detda == 3 ) then ! ETLMVar
           call etlmvar_run( x(1:nx), xe(1:nx,1:nmem), pb(1:nx,1:nx), obs(1:nslot,1:nobs) )
        else if( i_detda == 4 ) then ! 4DVar
           call fdvar_run( x(1:nx), xe(1:nx,1:nmem), pb(1:nx,1:nx), obs(1:nslot,1:nobs) )
        else
           write(*,*) 'i_detda error in da_run'
           stop 99
        end if
        call write_bias( iu_bias, time, x(1:nx), xt(1:nx) )
        call write_rmse( iu_rmse, time, x(1:nx), xt(1:nx) )
     end if
     if( i_ensda >= 1 ) then 
        ! verification of ensemble forecast
        xfe(1:nx,1:nmem) = xe(1:nx,1:nmem)
        time = dt_cycle * ( k - 1 )
        do kk = 1, nslot
           do j = 1, nmem
              call lorenz1996_run( 1, kk, dt_cycle/nslot, xfe(1:nx,j) )
           end do
           time = time + dt_cycle / nslot * kk
           call write_sprd( iu_sprd, time, xfe(1:nx,1:nmem) )
        end do
        ! ensemble forecast and analysis
        if( i_ensda == 1 ) then      ! KF-EDA
           pf(1:nx,1:nx) = pb(1:nx,1:nx)
           do j = 1, nmem
              obsf(1:nslot,1:nobs) = obs(1:nslot,1:nobs)
              call get_gauss_random( nobs*nslot, obsf(1:nslot,1:nobs)%val )
              obsf(1:nslot,1:nobs)%val = obs(1:nslot,1:nobs)%val &
                   & + obsf(1:nslot,1:nobs)%val * obsf(1:nslot,1:nobs)%err
              xe(1:nx,j) = x(1:nx) + ( xe(1:nx,j) - x(1:nx) ) * ( 1.0_r_size + infl )
              call kf_run( xe(1:nx,j), pf(1:nx,1:nx), obsf(1:nslot,1:nobs) )
           end do
        else if( i_ensda == 2 ) then ! EnVar-EDA
           do j = 1, nmem
              obsf(1:nslot,1:nobs) = obs(1:nslot,1:nobs)
              call get_gauss_random( nobs*nslot, obsf(1:nslot,1:nobs)%val )
              obsf(1:nslot,1:nobs)%val = obs(1:nslot,1:nobs)%val &
                   & + obsf(1:nslot,1:nobs)%val * obsf(1:nslot,1:nobs)%err
              xe(1:nx,j) = x(1:nx) + ( xe(1:nx,j) - x(1:nx) ) * ( 1.0_r_size + infl )
              call envar_run( xe(1:nx,j), xfe(1:nx,1:nmem), pb(1:nx,1:nx), obsf(1:nslot,1:nobs) )
           end do
        else if( i_ensda == 3 ) then ! ETLMVar-EDA
           do j = 1, nmem
              obsf(1:nslot,1:nobs) = obs(1:nslot,1:nobs)
              call get_gauss_random( nobs*nslot, obsf(1:nslot,1:nobs)%val )
              obsf(1:nslot,1:nobs)%val = obs(1:nslot,1:nobs)%val &
                   & + obsf(1:nslot,1:nobs)%val * obsf(1:nslot,1:nobs)%err
              xe(1:nx,j) = x(1:nx) + ( xe(1:nx,j) - x(1:nx) ) * ( 1.0_r_size + infl )
              call etlmvar_run( xe(1:nx,j), xfe(1:nx,1:nmem), pb(1:nx,1:nx), obsf(1:nslot,1:nobs) )
           end do
        else if( i_ensda == 4 ) then ! 4DVar-EDA
           do j = 1, nmem
              obsf(1:nslot,1:nobs) = obs(1:nslot,1:nobs)
              call get_gauss_random( nobs*nslot, obsf(1:nslot,1:nobs)%val )
              obsf(1:nslot,1:nobs)%val = obs(1:nslot,1:nobs)%val &
                   & + obsf(1:nslot,1:nobs)%val * obsf(1:nslot,1:nobs)%err
              xe(1:nx,j) = x(1:nx) + ( xe(1:nx,j) - x(1:nx) ) * ( 1.0_r_size + infl )
              call fdvar_run( xe(1:nx,j), xfe(1:nx,1:nmem), pb(1:nx,1:nx), obsf(1:nslot,1:nobs) )
           end do
        else if( i_ensda == 5 ) then ! EnSRF
           call ensrf_run( xe(1:nx,1:nmem), obs(1:nslot,1:nobs) )
        else if( i_ensda == 6 ) then ! LETKF
!           call letkf_run( xe(1:nx,1:nmem), obs(1:nslot,1:nobs) )
        else if( i_ensda == 7 ) then ! LPF
!           call lpf_run( xe(1:nx,1:nmem), obs(1:nslot,1:nobs) )
        else
           write(*,*) 'i_ensda error in da_run'
           stop 99
        end if
        ! recentering
        if( i_detda >= 1 ) then
           call recenter( xe(1:nx,1:nmem), x(1:nx) )
        else
           x(1:nx) = 0.0_r_size
           do j = 1, nmem
              x(1:nx) = x(1:nx) + xe(1:nx,j)
           end do
           x(1:nx) = x(1:nx) / nmem
           call write_bias( iu_bias, time, x(1:nx), xt(1:nx) )
           call write_rmse( iu_rmse, time, x(1:nx), xt(1:nx) )
        end if
        ! RTPP inflation
        do j = 1, nmem
           xe(1:nx,j) = x(1:nx) + rtpp * ( xfe(1:nx,j) - xf(1:nx) ) &
                & + ( 1.0_r_size - rtpp ) * ( xe(1:nx,j) - x(1:nx) )
        end do
        call write_sprd( iu_sprd, time, xe(1:nx,1:nmem) )
     end if
  end do cycle
  !
  stop
end program da_run
!
subroutine recenter( xe, x )
  use prm, only : r_size, nx, nmem
  implicit none
  real(r_size), intent(inout) :: xe(nx,nmem)
  real(r_size), intent(in) :: x(nx)
  real(r_size) :: xm(nx)
  integer :: j
  !
  xm(1:nx) = 0.0_r_size
  do j = 1, nmem
     xm(1:nx) = xm(1:nx) + xe(1:nx,j)
  end do
  xm(1:nx) = xm(1:nx) / nmem
  !
  do j = 1, nmem
     xe(1:nx,j) = xe(1:nx,j) - xm(1:nx)
     xe(1:nx,j) = xe(1:nx,j) + x(1:nx)
  end do
  !
  return
end subroutine recenter
