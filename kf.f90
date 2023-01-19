module kf
  use prm, only : r_size, &
       & nx, &
       & nslot, dt_cycle, &
       & nobs, obsdata, &
       & infl
  use lorenz1996, only : lorenz1996_run
  use obsope, only : obsope_run
  use matrix, only : matrix_eigen, matrix_sym
  implicit none
  private
  public :: kf_run
contains
  !
  subroutine kf_run( x, pb, obs )
    implicit none
    real(r_size), intent(inout) :: x(nx)         ! analysis variables
    real(r_size), intent(inout) :: pb(nx,nx)     ! analysis error covariances
    type(obsdata), intent(in) :: obs(nslot,nobs) ! observations
    integer :: i, ii, k, l, ll
    real(r_size) :: amean
    real(r_size) :: sqinfl
    real(r_size) :: dep(nobs)
    real(r_size) :: hpb(nobs,nx)
    real(r_size) :: hpbh(nobs,nobs)
    real(r_size) :: kg(nx,nobs)
    !
    do k = 1, nslot
       ! forecast
       call lorenz1996_run( 2, k, dt_cycle/nslot, x(1:nx) )       ! x_f = M x_b
       do i = 1, nx
          call lorenz1996_run( 3, k, dt_cycle/nslot, pb(1:nx,i) ) ! M P_b
       end do
       do i = 1, nx
          call lorenz1996_run( 3, k, dt_cycle/nslot, pb(i,1:nx) ) ! P_f = M (M P_b)^T
       end do
       call matrix_sym( nx, pb )
       ! multiplicative inflation
       sqinfl = (1.0_r_size + infl)**2
       pb(1:nx,1:nx) = pb(1:nx,1:nx) * sqinfl
       ! observation operator
       call obsope_run( 2, k, nobs, obs(k,1:nobs), x(1:nx), dep(1:nobs) )           ! H x_f
       do i = 1, nx
          call obsope_run( 3, k, nobs, obs(k,1:nobs), pb(1:nx,i), hpb(1:nobs,i) )   ! H P_f
       end do
       do l = 1, nobs
          call obsope_run( 3, k, nobs, obs(k,1:nobs), hpb(l,1:nx), hpbh(1:nobs,l) ) ! H (H P_f)^T
       end do
       call matrix_sym( nobs, hpbh )
       dep(1:nobs) = obs(k,1:nobs)%val - dep(1:nobs)                          ! y - H x_f
       ! Kalman gain
       do l = 1, nobs
          hpbh(l,l) = hpbh(l,l) + obs(k,l)%err**2                 ! H P_f H^T + R
       end do
       call matrix_eigen( 3, nobs, hpbh(1:nobs,1:nobs) )          ! ( H P_f H^T + R )^(-1)
       kg(1:nx,1:nobs) = 0.0_r_size
       do l = 1, nobs
          do ll = 1, nobs
             kg(1:nx,l) = kg(1:nx,l) + hpb(ll,1:nx) * hpbh(ll,l)  ! K = (H P_f)^T ( H P_f H^T + R )^(-1)
          end do
       end do
       ! analysis
       do l = 1, nobs
          x(1:nx) = x(1:nx) + kg(1:nx,l) * dep(l)            ! x_a = x_f + K ( y - H x_f )
          do i = 1, nx
             pb(1:nx,i) = pb(1:nx,i) - kg(1:nx,l) * hpb(l,i) ! P_a = P_f - K (H P_f)
          end do
       end do
       call matrix_sym( nx, pb )
    end do

    return
  end subroutine kf_run
  !
end module kf
