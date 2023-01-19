module minimizer
  use prm, only : r_size, niter
  implicit none
  integer, save :: iter
  real(r_size), save :: line_coef
  real(r_size), save :: cost_pre
  real(r_size), allocatable, save :: hesinv(:,:)
  real(r_size), allocatable, save :: grad_pre(:,:)
  real(r_size), allocatable, save :: vctl_pre(:,:)
  real(r_size), parameter :: line_coef0 = 0.1
  real(r_size), parameter :: eps = 1.0e-10
  private
  public :: bfgs, lbfgs
contains
  !
  subroutine bfgs( n, cost, grad, vctl, iflg )
    implicit none
    integer, intent(in) :: n               ! dimension
    real(r_size), intent(in) :: cost       ! cost function
    real(r_size), intent(in) :: grad(n)    ! gradient of cost function
    real(r_size), intent(inout) :: vctl(n) ! control variables
    integer, intent(out) :: iflg           ! 0:finish,1:update,2:line-search
    integer :: i, j, k
    real(r_size) :: sty
    real(r_size) :: syt(n,n)
    real(r_size) :: yst(n,n)
    real(r_size) :: sst(n,n)
    real(r_size) :: work(n,n)
    !
    if( sum(grad(1:n)**2) < eps ) then
       iflg = 0
       return
    else if( sum(vctl(1:n)**2) < eps ) then
       iflg = 1
       line_coef = line_coef0
       cost_pre = cost
       ! Approximate inverse of hessian
       if( .not. allocated(hesinv) ) allocate( hesinv(1:n,1:n) )
       hesinv(1:n,1:n) = 0.0_r_size
       do i = 1, n
          hesinv(i,i) = 1.0_r_size
       end do
       ! update control variables
       if( .not. allocated(vctl_pre) ) allocate( vctl_pre(1:n,1) )
       if( .not. allocated(grad_pre) ) allocate( grad_pre(1:n,1) )
       grad_pre(1:n,1) = grad(1:n)
       vctl_pre(1:n,1) = vctl(1:n)
       do i = 1, n
          do j = 1, n
             vctl(i) = vctl(i) - hesinv(i,j) * grad(j) * line_coef
          end do
       end do
       return
    else if( cost_pre < cost ) then
       iflg = 2
       line_coef = line_coef * 0.5_r_size
       ! update control variables
       do i = 1, n
          do j = 1, n
             vctl(i) = vctl_pre(i,1) - hesinv(i,j) * grad_pre(j,1) * line_coef
          end do
       end do
       return
    else
       iflg = 1
       line_coef = line_coef0
       cost_pre = cost
    end if
    ! update inverse of hessian matrix
    sty = 0.0_r_size
    do i = 1, n
       sty = sty + ( vctl(i) - vctl_pre(i,1) ) * ( grad(i) - grad_pre(i,1) )
    end do
    sty = 1.0_r_size / sty
    syt(1:n,1:n) = 0.0_r_size
    yst(1:n,1:n) = 0.0_r_size
    do i = 1, n
       syt(i,i) = 1.0_r_size
       yst(i,i) = 1.0_r_size
       do j = 1, n
          syt(i,j) = syt(i,j) - sty * ( vctl(i) - vctl_pre(i,1) ) * ( grad(j) - grad_pre(j,1) )
          yst(i,j) = yst(i,j) - sty * ( vctl(j) - vctl_pre(j,1) ) * ( grad(i) - grad_pre(i,1) )
          sst(i,j) = sty * ( vctl(i) - vctl_pre(i,1) ) * ( vctl(j) - vctl_pre(j,1) )
       end do
    end do
    work(1:n,1:n) = 0.0_r_size
    do i = 1, n
       do j = 1, n
          do k = 1, n
             work(i,j) = work(i,j) + hesinv(i,k) * yst(k,j)
          end do
       end do
    end do
    hesinv(1:n,1:n) = sst(1:n,1:n)
    do i = 1, n
       do j = 1, n
          do k = 1, n
             hesinv(i,j) = hesinv(i,j) + syt(i,k) * work(k,j)
          end do
       end do
    end do
    ! update control variables
    grad_pre(1:n,1) = grad(1:n)
    vctl_pre(1:n,1) = vctl(1:n)
    do i = 1, n
       do j = 1, n
          vctl(i) = vctl(i) - hesinv(i,j) * grad(j) * line_coef
       end do
    end do
    !
    return
  end subroutine bfgs
  !
  subroutine lbfgs( n, cost, grad, vctl, iflg )
    implicit none
    integer, intent(in) :: n               ! dimension
    real(r_size), intent(in) :: cost       ! cost function
    real(r_size), intent(in) :: grad(n)    ! gradient of cost function
    real(r_size), intent(inout) :: vctl(n) ! control variables
    integer, intent(out) :: iflg           ! 0:finish,1:update,2:line-search
    integer :: i, j, k
    real(r_size) :: sty(niter)
    real(r_size) :: stw(niter)
    real(r_size) :: ytw(niter)
    real(r_size) :: work(n)
    !
    if( sum(grad(1:n)**2) < eps ) then
       iflg = 0
       return
    else if( sum(vctl(1:n)**2) < eps ) then
       iflg = 1
       iter = 1
       line_coef = line_coef0
       cost_pre = cost
       ! update control variables
       if( .not. allocated(grad_pre) ) allocate( grad_pre(1:n,1:niter) )
       if( .not. allocated(vctl_pre) ) allocate( vctl_pre(1:n,1:niter) )
       grad_pre(1:n,iter) = grad(1:n)
       vctl_pre(1:n,iter) = vctl(1:n)
       vctl(1:n) = vctl(1:n) - grad(1:n) * line_coef
       return
    else if( cost_pre < cost ) then
       iflg = 2
       line_coef = line_coef * 0.5_r_size
       vctl(1:n) = vctl_pre(1:n,iter)
       work(1:n) = grad_pre(1:n,iter)
    else
       iflg = 1
       iter = iter + 1
       line_coef = line_coef0
       cost_pre = cost
       grad_pre(1:n,iter) = grad(1:n)
       vctl_pre(1:n,iter) = vctl(1:n)
       work(1:n) = grad(1:n)
    end if
    ! update gradient of cost function
    do k = 1, iter-1
       sty(k) = 0.0_r_size
       do i = 1, n
          sty(k) = sty(k) + ( vctl_pre(i,k+1) - vctl_pre(i,k) ) * ( grad_pre(i,k+1) - grad_pre(i,k) )
       end do
       sty(k) = 1.0_r_size / sty(k)
    end do
    do k = iter-1, 1, -1
       stw(k) = 0.0_r_size
       do i = 1, n
          stw(k) = stw(k) + sty(k) * ( vctl_pre(i,k+1) - vctl_pre(i,k) ) * work(i)
       end do
       work(1:n) = work(1:n) - stw(k) * ( grad_pre(1:n,k+1) - grad_pre(1:n,k) )
    end do
    do k = 1, iter-1
       ytw(k) = 0.0_r_size
       do i = 1, n
          ytw(k) = ytw(k) + sty(k) * ( grad_pre(i,k+1) - grad_pre(i,k) ) * work(i)
       end do
       work(1:n) = work(1:n) + ( stw(k) - ytw(k) ) * ( vctl_pre(1:n,k+1) - vctl_pre(1:n,k) )
    end do
    ! update control variables
    vctl(1:n) = vctl(1:n) - work(1:n) * line_coef
    !
    return
  end subroutine lbfgs
  !
end module minimizer
