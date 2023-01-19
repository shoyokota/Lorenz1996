module io
  use prm, only : r_size, nx, nobs, obsdata, nmem
  implicit none
  integer, save :: iu = 10
  character(64), save :: f_grid
  character(64), save :: f_obs
  character(64), save :: f_verif
  real(r_size), parameter :: eps = 1.0e-10
  private
  public :: open_data, close_data
  public :: read_grid, write_grid
  public :: read_obs, write_obs
  public :: write_bias, write_rmse, write_sprd
contains
  !
  subroutine open_data( iunit, filename )
    implicit none
    integer, intent(out) :: iunit
    character(*), intent(in) :: filename
    !
    iu = iu + 1
    iunit = iu
    open( iunit, file=filename, form='FORMATTED' )
    f_grid='(XXes20.10)'
    write( f_grid(2:3), '(i2)' ) nx + 1
    f_obs='(es20.10,i5,3es20.10)'
    f_verif='(2es20.10)'
    !
    return
  end subroutine open_data
  !
  subroutine close_data( iunit )
    implicit none
    integer, intent(in) :: iunit
    !
    close( iunit )
    !
    return
  end subroutine close_data
  !
  subroutine read_grid( iunit, time, x )
    implicit none
    integer, intent(in) :: iunit
    real(r_size), intent(in) :: time
    real(r_size), intent(out) :: x(nx)
    integer :: i
    real(r_size) :: t
    !
    read( iunit, f_grid ) t, ( x(i), i=1,nx )
    if( abs(t-time) > eps ) then
       write(*,*) 'time error in read_grid'
       stop 99
    end if
    !
    return
  end subroutine read_grid
  !
  subroutine write_grid( iunit, time, x )
    implicit none
    integer, intent(in) :: iunit
    real(r_size), intent(in) :: time
    real(r_size), intent(in) :: x(nx)
    integer :: i
    !
    write( iunit, f_grid ) time, ( x(i), i=1,nx )
    !
    return
  end subroutine write_grid
  !
  subroutine read_obs( iunit, time, obs )
    implicit none
    integer, intent(in) :: iunit
    real(r_size), intent(in) :: time
    type(obsdata), intent(out) :: obs(nobs)
    integer :: j
    real(r_size) :: t
    !
    do j = 1, nobs
       read( iunit, f_obs ) t, obs(j)%id, obs(j)%pos, obs(j)%val, obs(j)%err
       if( abs(t-time) > eps ) then
          write(*,*) 'time error in read_obs'
          stop 99
       end if
    end do
    !
    return
  end subroutine read_obs
  !
  subroutine write_obs( iunit, time, obs )
    implicit none
    integer, intent(in) :: iunit
    real(r_size), intent(in) :: time
    type(obsdata), intent(in) :: obs(nobs)
    integer :: j
    !
    do j = 1, nobs
       write( iunit, f_obs ) time, obs(j)%id, obs(j)%pos, obs(j)%val, obs(j)%err
    end do
    !
    return
  end subroutine write_obs
  !
  subroutine write_bias( iunit, time, x, xt )
    implicit none
    integer, intent(in) :: iunit
    real(r_size), intent(in) :: time
    real(r_size), intent(in) :: x(nx)
    real(r_size), intent(in) :: xt(nx)
    real(r_size) :: bias
    integer :: i
    !
    bias = 0.0_r_size
    do i = 1, nx
       bias = bias + ( x(i) - xt(i) )
    end do
    bias = bias / nx
    write( iunit, f_verif ) time, bias
    write( 6, * ) 't=', time, ', bias=', bias
    !
    return
  end subroutine write_bias
  !
  subroutine write_rmse( iunit, time, x, xt )
    implicit none
    integer, intent(in) :: iunit
    real(r_size), intent(in) :: time
    real(r_size), intent(in) :: x(nx)
    real(r_size), intent(in) :: xt(nx)
    real(r_size) :: rmse
    integer :: i
    !
    rmse = 0.0_r_size
    do i = 1, nx
       rmse = rmse + ( x(i) - xt(i) )**2
    end do
    rmse = rmse / nx
    rmse = sqrt( rmse )
    write( iunit, f_verif ) time, rmse
    write( 6, * ) 't=', time, ', rmse=', rmse
    !
    return
  end subroutine write_rmse
  !
  subroutine write_sprd( iunit, time, xe )
    implicit none
    integer, intent(in) :: iunit
    real(r_size), intent(in) :: time
    real(r_size), intent(in) :: xe(nx,nmem)
    real(r_size) :: mean, sprd
    integer :: i, j
    !
    sprd = 0.0_r_size
    do i = 1, nx
       mean = 0.0_r_size
       do j = 1, nmem
          mean = mean + xe(i,j)
       end do
       mean = mean / nmem
       do j = 1, nmem
          sprd = sprd + ( xe(i,j) - mean )**2
       end do
    end do
    sprd = sprd / ( nx * nmem )
    sprd = sqrt( sprd )
    write( iunit, f_verif ) time, sprd
    write( 6, * ) 't=', time, ', sprd=', sprd
    !
    return
  end subroutine write_sprd
  !
end module io
