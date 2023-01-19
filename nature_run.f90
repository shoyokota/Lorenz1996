program nature_run
  use prm, only : r_size, &
       & nx, ft_nature_pre, ft_nature, &
       & nslot, dt_cycle, &
       & nobs, obsid, obserr, obsthin, obsdata, &
       & truthfile, obsfile
  use prm, only : read_namelist
  use io, only : open_data, write_grid, write_obs
  use random, only : init_random, get_gauss_random
  use lorenz1996, only : lorenz1996_run
  implicit none
  real(r_size), allocatable :: x(:), dx(:)
  type(obsdata), allocatable :: obs(:)
  integer :: i, j, k, nft
  integer :: iu_grid, iu_obs
  real(r_size) :: time
  !
  call read_namelist
  allocate( x(1:nx)     )
  allocate( dx(1:nx)    )
  allocate( obs(1:nobs) )
  ! 
  call init_random( 0 )
  call get_gauss_random( nx, x(1:nx) )
  call lorenz1996_run( 1, 1, ft_nature_pre, x(1:nx) )
  !
  call open_data( iu_grid, truthfile )
  call open_data( iu_obs,  obsfile   )
  time = 0.0_r_size
  !
  do k = 1, nint(ft_nature/dt_cycle/nslot)
     call lorenz1996_run( 1, 1, dt_cycle/nslot, x(1:nx) )     
     time = dt_cycle / nslot * k
     call write_grid( iu_grid, time, x(1:nx) )
     call get_gauss_random( nx, dx(1:nx) )
     j = 1
     do i = 1, nx, obsthin
        obs(j)%id  = obsid
        obs(j)%pos = real(i,r_size)
        obs(j)%val = x(i) + dx(i) * obserr
        obs(j)%err = obserr
        j = j + 1
        if( j > nobs ) exit
     end do
     call write_obs( iu_obs, time, obs(1:nobs) )
  end do
  !
  stop
end program nature_run
