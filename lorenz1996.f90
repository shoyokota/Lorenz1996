module lorenz1996
  use prm, only : r_size, nx, dt, forcing, nslot
  implicit none
  real(r_size), allocatable, save :: xb(:,:,:,:)
  private
  public :: lorenz1996_run
contains
  !
  subroutine lorenz1996_run( iflg, k, ft, x )
    implicit none
    integer, intent(in)         :: iflg  ! 1:non-linear,2:basic-state,3:tangent-linear,4:adjoint
    integer, intent(in)         :: k     ! time slot
    real(r_size), intent(in)    :: ft    ! forecast time
    real(r_size), intent(inout) :: x(nx) ! forecast variables
    integer :: kk, nt
    !
    nt = nint(ft/dt)
    if( iflg == 1 ) then      ! non-linear
       do kk = 1, nt
          call lorenz1996_1step( x(1:nx) )
       end do
    else if( iflg == 2 ) then ! basic-state
       if( k == 1 ) then
          if( allocated(xb) ) deallocate(xb)
          allocate( xb(1:nslot,1:nx,0:nt-1,4) )
       end if
       xb(k,1:nx,0:nt-1,4) = 0.0_r_size
       do kk = 1, nt
          call lorenz1996_1step( x(1:nx), xb(k,1:nx,k-1,1:4) )
       end do
    else if( iflg == 3 ) then ! tangent-linear
       do kk = 1, nt
          call lorenz1996_1step_tl( x(1:nx), xb(k,1:nx,k-1,1:4) )
       end do
    else if( iflg == 4 ) then ! adjoint
       do kk = nt, 1, -1
          call lorenz1996_1step_ad( x(1:nx), xb(k,1:nx,k-1,1:4) )
       end do
    end if
    !
    return
  end subroutine lorenz1996_run
  !
  subroutine lorenz1996_1step( x, xxb )
    implicit none
    real(r_size), intent(inout)           :: x(nx)     ! forecast variables
    real(r_size), intent(inout), optional :: xxb(nx,4) ! basic-state
    real(r_size) :: xx(nx,4)
    real(r_size) :: xu(nx,4)
    !
    xx(:,1) = x(:)
    call lorenz1996_core( xx(:,1), xu(:,1) )
    xx(:,2) = x(:) + xu(:,1) * dt / 2.0_r_size
    call lorenz1996_core( xx(:,2), xu(:,2) )
    xx(:,3) = x(:) + xu(:,2) * dt / 2.0_r_size
    call lorenz1996_core( xx(:,3), xu(:,3) )
    xx(:,4) = x(:) + xu(:,3) * dt
    call lorenz1996_core( xx(:,4), xu(:,4) )
    x(:) = x(:) + ( &
         & xu(:,1) + xu(:,2)*2.0_r_size + xu(:,3)*2.0_r_size + xu(:,4) &
         & ) * dt / 6.0_r_size
    if( present(xxb) ) xxb(:,:) = xx(:,:)
    !
    return
  end subroutine lorenz1996_1step
  !
  subroutine lorenz1996_1step_tl( x, xxb )
    implicit none
    real(r_size), intent(inout) :: x(nx)     ! forecast variables
    real(r_size), intent(in)    :: xxb(nx,4) ! basic-state
    real(r_size) :: xx(nx,4)
    real(r_size) :: xu(nx,4)
    !
    xx(:,1) = x(:)
    call lorenz1996_core_tl( xxb(:,1), xx(:,1), xu(:,1) )
    xx(:,2) = x(:) + xu(:,1) * dt / 2.0_r_size
    call lorenz1996_core_tl( xxb(:,2), xx(:,2), xu(:,2) )
    xx(:,3) = x(:) + xu(:,2) * dt / 2.0_r_size
    call lorenz1996_core_tl( xxb(:,3), xx(:,3), xu(:,3) )
    xx(:,4) = x(:) + xu(:,3) * dt
    call lorenz1996_core_tl( xxb(:,4), xx(:,4), xu(:,4) )
    x(:) = x(:) + ( &
         & xu(:,1) + xu(:,2)*2.0_r_size + xu(:,3)*2.0_r_size + xu(:,4) &
         & ) * dt / 6.0_r_size
    !
    return
  end subroutine lorenz1996_1step_tl
  !
  subroutine lorenz1996_1step_ad( x, xxb )
    implicit none
    real(r_size), intent(inout) :: x(nx)     ! forecast variables
    real(r_size), intent(in)    :: xxb(nx,4) ! basic-state
    real(r_size) :: xx(nx,4)
    real(r_size) :: xu(nx,4)
    !
    xu(:,4) = x(:) * dt / 6.0_r_size
    xu(:,3) = x(:) * dt / 6.0_r_size * 2.0_r_size
    xu(:,2) = x(:) * dt / 6.0_r_size * 2.0_r_size
    xu(:,1) = x(:) * dt / 6.0_r_size
    xx(:,4) = 0.0_r_size
    call lorenz1996_core_ad( xxb(:,4), xx(:,4), xu(:,4) )
    x(:) = x(:) + xx(:,4)
    xu(:,3) = xu(:,3) + xx(:,4) * dt
    xx(:,3) = 0.0_r_size
    call lorenz1996_core_ad( xxb(:,3), xx(:,3), xu(:,3) )
    x(:) = x(:) + xx(:,3)
    xu(:,2) = xu(:,2) + xx(:,3) * dt / 2.0_r_size
    xx(:,2) = 0.0_r_size
    call lorenz1996_core_ad( xxb(:,2), xx(:,2), xu(:,2) )
    x(:) = x(:) + xx(:,2)
    xu(:,1) = xu(:,1) + xx(:,2) * dt / 2.0_r_size
    xx(:,1) = 0.0_r_size
    call lorenz1996_core_ad( xxb(:,1), xx(:,1), xu(:,1) )
    x(:) = x(:) + xx(:,1)
    !
    return
  end subroutine lorenz1996_1step_ad
  !
  subroutine lorenz1996_core( xx, xu )
    implicit none
    real(r_size), intent(in)  :: xx(nx) ! forecast variables (current time)
    real(r_size), intent(out) :: xu(nx) ! forecast variables (1-step after)
    integer :: i, ip1, im1, im2
    !
    do i = 1, nx
       ip1 = i + 1
       im1 = i - 1
       im2 = i - 2
       if( ip1 > nx ) ip1 = ip1 - nx
       if( im1 <= 0 ) im1 = im1 + nx
       if( im2 <= 0 ) im2 = im2 + nx
       xu(i) = (xx(ip1) - xx(im2)) * xx(im1) - xx(i) + forcing
    end do
    !
    return
  end subroutine lorenz1996_core
  !
  subroutine lorenz1996_core_tl( xxb, xx, xu )
    implicit none
    real(r_size), intent(in)  :: xxb(nx) ! basic-state
    real(r_size), intent(in)  :: xx(nx)  ! forecast variables (current time)
    real(r_size), intent(out) :: xu(nx)  ! forecast variables (1-step after)
    integer :: i, ip1, im1, im2
    !
    do i = 1, nx
       ip1 = i + 1
       im1 = i - 1
       im2 = i - 2
       if( ip1 > nx ) ip1 = ip1 - nx
       if( im1 <= 0 ) im1 = im1 + nx
       if( im2 <= 0 ) im2 = im2 + nx
       xu(i) = (xx(ip1) - xx(im2)) * xxb(im1) &
            & + (xxb(ip1) - xxb(im2)) * xx(im1) - xx(i)
    end do
    !
    return
  end subroutine lorenz1996_core_tl
  !
  subroutine lorenz1996_core_ad( xxb, xx, xu )
    implicit none
    real(r_size), intent(in)    :: xxb(nx) ! basic-state
    real(r_size), intent(inout) :: xx(nx)  ! forecast variables (current time)
    real(r_size), intent(in)    :: xu(nx)  ! forecast variables (1-step after)
    integer :: i, ip1, im1, im2
    !
    do i = nx, 1, -1
       ip1 = i + 1
       im1 = i - 1
       im2 = i - 2
       if( ip1 > nx ) ip1 = ip1 - nx
       if( im1 <= 0 ) im1 = im1 + nx
       if( im2 <= 0 ) im2 = im2 + nx
       xx(i)   = xx(i)   - xu(i)
       xx(im1) = xx(im1) + xu(i) * (xxb(ip1) - xxb(im2))
       xx(im2) = xx(im2) - xu(i) * xxb(im1)
       xx(ip1) = xx(ip1) + xu(i) * xxb(im1)
    end do
    !
    return
  end subroutine lorenz1996_core_ad
  !
end module lorenz1996
