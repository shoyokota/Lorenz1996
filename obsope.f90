module obsope
  use prm, only : r_size, nx, nslot, obsdata
  implicit none
  real(r_size), allocatable, save :: hb(:,:)
  private
  public :: obsope_run
contains
  !
  subroutine obsope_run( iflg, k, n, obs, x, hx )
    implicit none
    integer, intent(in)         :: iflg   ! 1:non-linear,2:basic-state,3:tangent-linear,4:adjoint
    integer, intent(in)         :: k      ! time slot
    integer, intent(in)         :: n      ! number of observations
    type(obsdata), intent(in)   :: obs(n) ! observations
    real(r_size), intent(inout) :: x(nx)  ! grid-space variables
    real(r_size), intent(inout) :: hx(n)  ! observation-space variables
    integer :: l
    !
    if( iflg == 1 ) then      ! non-linear
       call obsope_nl( n, obs(1:n), x(1:nx), hx(1:n) )
    else if( iflg == 2 ) then ! basic-state
       if( k == 1 ) then
          if( allocated(hb) ) deallocate(hb)
          allocate( hb(1:nslot,1:n) )
       end if
       hb(k,1:n) = 0.0_r_size
       call obsope_nl( n, obs(1:n), x(1:nx), hx(1:n), hb(k,1:n) )
    else if( iflg == 3 ) then ! tangent-linear
       call obsope_tl( n, obs(1:n), x(1:nx), hx(1:n), hb(k,1:n) )
    else if( iflg == 4 ) then ! adjoint
       call obsope_ad( n, obs(1:n), x(1:nx), hx(1:n), hb(k,1:n) )
    end if
    !
    return
  end subroutine obsope_run
  !
  subroutine obsope_nl( n, obs, x, hx, hhb )
    implicit none
    integer, intent(in)                 :: n      ! number of observations
    type(obsdata), intent(in)           :: obs(n) ! observations
    real(r_size), intent(in)            :: x(nx)  ! grid-space variables
    real(r_size), intent(out)           :: hx(n)  ! observation-space variables
    real(r_size), intent(out), optional :: hhb(n) ! basic-state
    integer :: i, ip1, l
    integer :: ipos
    real(r_size) :: hhx
    !
    do l = 1, n
       do i = 1, nx
          ipos = int(obs(l)%pos)
          if( mod(ipos-1,nx)+1 == i ) then
             ip1 = i + 1
             if( ip1 > nx ) ip1 = ip1 - nx
             hhx = ( ipos + 1 - obs(l)%pos ) * x(i) &
                  & + ( obs(l)%pos - ipos ) * x(ip1)
             exit
          end if
       end do
       if( present(hhb) ) hhb(l) = hhx
       if( obs(l)%id == 1 ) then
          hx(l) = hhx
       else if( obs(l)%id == 2 ) then
          hx(l) = hhx**2
       else if( obs(l)%id == 3 ) then
          hx(l) = exp( hhx )
       else if( obs(l)%id == 4 ) then
          hx(l) = log( hhx )
       else
          write(*,*) 'obsid error in obsope_nl'
          stop 99
       end if
    end do
    !
    return
  end subroutine obsope_nl
  !
  subroutine obsope_tl( n, obs, x, hx, hhb )
    implicit none
    integer, intent(in)       :: n      ! number of observations
    type(obsdata), intent(in) :: obs(n) ! observations
    real(r_size), intent(in)  :: x(nx)  ! grid-space variables
    real(r_size), intent(out) :: hx(n)  ! observation-space variables
    real(r_size), intent(in)  :: hhb(n) ! basic-state
    integer :: i, ip1, l
    integer :: ipos
    real(r_size) :: hhx
    !
    do l = 1, n
       do i = 1, nx
          ipos = int(obs(l)%pos)
          if( mod(ipos-1,nx)+1 == i ) then
             ip1 = i + 1
             if( ip1 > nx ) ip1 = ip1 - nx
             hhx = ( ipos + 1 - obs(l)%pos ) * x(i) &
                  & + ( obs(l)%pos - ipos ) * x(ip1)
             exit
          end if
       end do
       if( obs(l)%id == 1 ) then
          hx(l) = hhx
       else if( obs(l)%id == 2 ) then
          hx(l) = 2.0_r_size * hhb(l) * hhx
       else if( obs(l)%id == 3 ) then
          hx(l) = exp( hhb(l) ) * hhx
       else if( obs(l)%id == 4 ) then
          hx(l) = hhx / hhb(l)
       else
          write(*,*) 'obsid error in obsope_tl'
          stop 99
       end if
    end do
    !
    return
  end subroutine obsope_tl
  !
  subroutine obsope_ad( n, obs, x, hx, hhb )
    implicit none
    integer, intent(in)         :: n      ! number of observations
    type(obsdata), intent(in)   :: obs(n) ! observations
    real(r_size), intent(inout) :: x(nx)  ! grid-space variables
    real(r_size), intent(in)    :: hx(n)  ! observation-space variables
    real(r_size), intent(in)    :: hhb(n) ! basic-state
    integer :: i, ip1, l
    integer :: ipos
    real(r_size) :: hhx
    !
    do l = n, 1, -1
       hhx = 0.0_r_size
       if( obs(l)%id == 1 ) then
          hhx = hhx + hx(l)
       else if( obs(l)%id == 2 ) then
          hhx = hhx + 2.0_r_size * hhb(l) * hx(l)
       else if( obs(l)%id == 3 ) then
          hhx = hhx + exp( hhb(l) ) * hx(l)
       else if( obs(l)%id == 4 ) then
          hhx = hhx + hx(l) / hhb(l)
       else
          write(*,*) 'obsid error in obsope_ad'
          stop 99
       end if
       do i = nx, 1, -1
          ipos = int(obs(l)%pos)
          if( mod(ipos-1,nx)+1 == i ) then
             ip1 = i + 1
             if( ip1 > nx ) ip1 = ip1 - nx
             x(i) = x(i) + ( ipos + 1 - obs(l)%pos ) * hhx
             x(ip1) = x(ip1) + ( obs(l)%pos - ipos ) * hhx
             exit
          end if
       end do
    end do
    !
    return
  end subroutine obsope_ad
  !
end module obsope
