module etlmvar
  use prm, only : r_size, &
       & nx, &
       & nslot, aslot, niter, dt_cycle, &
       & nobs, obsdata, &
       & nmem, hyb_betab, hyb_betae, &
       & sigma_sloc, sigma_tloc, sigma_sloce, sigma_tloce, infl
  use lorenz1996, only : lorenz1996_run
  use obsope, only : obsope_run
  use matrix, only : matrix_eigen, matrix_gauss
  use minimizer, only : lbfgs
  implicit none
  private
  public :: etlmvar_run
contains
  !
  subroutine etlmvar_run( x, xe, pb, obs )
    implicit none
    real(r_size), intent(inout) :: x(nx)         ! analysis variables
    real(r_size), intent(in) :: xe(nx,nmem)      ! ensemble forecasts
    real(r_size), intent(in) :: pb(nx,nx)        ! static background error covariance
    type(obsdata), intent(in) :: obs(nslot,nobs) ! observations
    integer :: i, ii, j, k, l, iter
    integer :: iflg
    real(r_size) :: xf(nx)
    real(r_size) :: xfe(nx,nmem)
    real(r_size) :: xm(nx)
    real(r_size) :: dep(nslot,nobs)
    real(r_size) :: sqrt_pb(nx,nx)
    real(r_size) :: sqrt_pe(0:nslot,nx,nmem)
    real(r_size) :: pl(0:nslot,nx,nx)
    real(r_size) :: inv_pe(nx)
    real(r_size) :: etlm_loc(nslot,nx,nx)
    real(r_size) :: sqrt_loc(nslot,nx,nx)
    real(r_size) :: vctl(nx*(1+nmem))
    real(r_size) :: costb, costo
    real(r_size) :: grad(nx*(1+nmem))
    real(r_size) :: dx(nx)
    real(r_size) :: dxb(nx)
    real(r_size) :: mdxb(nx)
    real(r_size) :: dxe(nx)
    real(r_size) :: dxf(nx)
    real(r_size) :: dxp(nx)
    real(r_size) :: hdx(nobs)
    real(r_size) :: hdxmd
    real(r_size) :: rhdxmd(nslot,nobs)
    real(r_size) :: hrhdxmd(nx)
    real(r_size) :: hrhdxmde(nx)
    real(r_size) :: mhrhdxmdb(nx)
    real(r_size) :: hrhdxmdb(nx)
    !
    xf(1:nx) = x(1:nx)
    xfe(1:nx,1:nmem) = xe(1:nx,1:nmem)
    do k = 0, nslot
       if( k > 0 ) then
          ! forecast
          call lorenz1996_run( 1, k, dt_cycle/nslot, xf(1:nx) )       ! x_f = M x_b
          do j = 1, nmem
             call lorenz1996_run( 1, k, dt_cycle/nslot, xfe(1:nx,j) ) ! x_f_e = M x_e
          end do
          ! observation operator
          call obsope_run( 2, k, nobs, obs(k,1:nobs), xf(1:nx), dep(k,1:nobs) ) ! H x_f
          dep(k,1:nobs) = obs(k,1:nobs)%val - dep(k,1:nobs)                     ! d = y - H x_f
       end if
       ! P_e^(1/2)
       xm(1:nx) = 0.0_r_size
       do j = 1, nmem
          xm(1:nx) = xm(1:nx) + xfe(1:nx,j)
       end do
       xm(1:nx) = xm(1:nx) / nmem
       !xm(1:nx) = xm(1:nx) + ( xf(1:nx) - xm(1:nx) ) * 1.0e-3 ! for full-rank P_e
       do j = 1, nmem
          sqrt_pe(k,1:nx,j) = ( xfe(1:nx,j) - xm(1:nx) ) / ( nmem - 1 )**0.5
          sqrt_pe(k,1:nx,j) = sqrt_pe(k,1:nx,j) * (1.0_r_size + infl) ! multiplicative inflation
       end do
       ! first guess (at aslot)
       if( k == aslot ) then
          x(1:nx) = xf(1:nx)
          do i = 1, nx
             inv_pe(i) = 1.0_r_size / sum( sqrt_pe(0,i,:)**2 )
          end do
       end if
    enddo
    ! P_b^(1/2)
    sqrt_pb(1:nx,1:nx) = pb(1:nx,1:nx)
    call matrix_eigen( 4, nx, sqrt_pb(1:nx,1:nx) )
    ! L^(1/2)
    do k = 1, nslot
       call matrix_gauss( nx, sigma_sloc, sqrt_loc(k,1:nx,1:nx) )
       sqrt_loc(k,1:nx,1:nx) = sqrt_loc(k,1:nx,1:nx) * exp( -0.5_r_size * ( ( k - aslot ) * dt_cycle / nslot / sigma_tloc )**2 )
       call matrix_eigen( 4, nx, sqrt_loc(k,1:nx,1:nx) )
    end do
    do k = 1, nslot
       call matrix_gauss( nx, sigma_sloce, etlm_loc(k,1:nx,1:nx) )
       etlm_loc(k,1:nx,1:nx) = etlm_loc(k,1:nx,1:nx) * exp( -0.5_r_size * ( k * dt_cycle / nslot / sigma_tloce )**2 )
    end do
    ! output P_b and P_e
    if(hyb_betab==1.0) then
       pl=0.0_r_size
       do i=1,nx;do ii=1,nx;do j=1,nmem;do l=1,nx
          pl(aslot,i,ii)=pl(aslot,i,ii)+etlm_loc(aslot,i,l)*sqrt_pe(aslot,i,j)*sqrt_pe(0,l,j)*inv_pe(l)*sqrt_pb(l,ii)
          !pl(0,i,ii)=pl(0,i,ii)+etlm_loc(aslot,i,l)*sqrt_pe(0,i,j)*sqrt_pe(0,l,j)*inv_pe(l)*sqrt_pb(l,ii)
       enddo;enddo;enddo;enddo
       open(90,file='MAT_ETLMVar_MBcM.txt',position='append')
       write(90,'(i5,40es20.10)')0,(0.0_r_size,ii=1,nx)
       do i=1,nx
          write(90,'(i5,40es20.10)')i,(sum(pl(aslot,i,:)*pl(aslot,ii,:)),ii=1,nx)
          !write(90,'(i5,40es20.10)')i,(sum(pl(aslot,i,:)*pl(0,ii,:)),ii=1,nx)
       enddo
       write(*,*)
       close(90)
       pl=0.0_r_size
       do i=1,nx;do ii=1,nx;do j=1,nmem;do l=1,nx
          pl(aslot,i,ii)=pl(aslot,i,ii)+sqrt_loc(aslot,i,l)*sqrt_pe(aslot,i,j)*sqrt_pe(aslot,ii,j)*sqrt_loc(aslot,ii,l)
          !pl(aslot,i,ii)=pl(aslot,i,ii)+sqrt_loc(aslot,i,l)*sqrt_pe(aslot,i,j)*sqrt_pe(0,ii,j)*sqrt_loc(0,ii,l)
          !pl(aslot,i,ii)=pl(aslot,i,ii)+sqrt_loc(0,i,l)*sqrt_pe(0,i,j)*sqrt_pe(0,ii,j)*sqrt_loc(0,ii,l)
       enddo;enddo;enddo;enddo
       open(90,file='MAT_ETLMVar_MBeM.txt',position='append')
       write(90,'(i5,40es20.10)')0,(0.0_r_size,ii=1,nx)
       do i=1,nx
          write(90,'(i5,40es20.10)')i,(pl(aslot,i,ii),ii=1,nx)
       enddo
       write(*,*)
       close(90)
    endif
    ! minimization
    vctl(1:nx*(1+nmem)) = 0.0_r_size
    dx(1:nx) = 0.0_r_size
    iteration: do iter = 1, niter
       ! Jb
       costb = 0.0_r_size
       do i = 1, nx*(1+nmem)
          costb = costb + vctl(i)**2 ! 2 * Jb
          grad(i) = vctl(i)          ! gJb
       end do
       costb = costb * 0.5_r_size
       ! Jo
       costo = 0.0_r_size
       ! dx = beta_b M_e P_b^(1/2) v_b + beta_e P_e^(1/2) v_e
       dxb(1:nx) = 0.0_r_size
       do ii = 1, nx
          do i = 1, nx
             ! P_b^(1/2) v_b
             dxb(ii) = dxb(ii) + sqrt_pb(ii,i) * vctl(i)
          end do
       end do
       do k = 1, nslot
          mdxb(1:nx) = 0.0_r_size
          do ii = 1, nx
             do j = 1, nmem
                do i = 1, nx
                   ! M_e P_b^(1/2) v_b
                   mdxb(ii) = mdxb(ii) + etlm_loc(k,ii,i) * sqrt_pe(k,ii,j) * sqrt_pe(0,i,j) * inv_pe(i) * dxb(i)
                end do
             end do
          end do
          dxe(1:nx) = 0.0_r_size
          do ii = 1, nx
             do j = 1, nmem
                do i = 1, nx
                   ! P_e^(1/2) v_e
                   dxe(ii) = dxe(ii) + sqrt_loc(k,ii,i) * sqrt_pe(k,ii,j) * vctl(nx*j+i)
                   !dxe(ii) = dxe(ii) + sqrt_pb(ii,i) * vctl(nx*j+i)
                end do
             end do
          end do
          dxf(1:nx) = hyb_betab**0.5 * mdxb(1:nx) + hyb_betae**0.5 * dxe(1:nx)
          ! analysis increment (at aslot)
          if( k == aslot ) then
             dxp(1:nx) = dxf(1:nx)
          end if
          ! H dx
          call obsope_run( 3, k, nobs, obs(k,1:nobs), dxf(1:nx), hdx(1:nobs) )
          ! 2 * Jo = ( H dx - d )^T R^(-1) ( H dx - d )
          do l = 1, nobs
             hdxmd = hdx(l) - dep(k,l)              ! H dx - d
             rhdxmd(k,l) = hdxmd / obs(k,l)%err**2  ! R^(-1) ( H dx - d )
             costo = costo + hdxmd * rhdxmd(k,l)    ! 2 * Jo
          end do
       end do
       ! gJo = [ beta_b M_e P_b^(1/2), beta_e P_e^(1/2) ]^T H^T R^(-1) * ( H dx - d )
       mhrhdxmdb(1:nx) = 0.0_r_size
       do k = nslot, 1, -1
          ! H^T R^(-1) * ( H dx - d )
          hrhdxmd(1:nx) = 0.0_r_size
          call obsope_run( 4, k, nobs, obs(k,1:nobs), hrhdxmd(1:nx), rhdxmd(k,1:nobs) )
          hrhdxmdb(1:nx) = hyb_betab**0.5 * hrhdxmd(1:nx)
          hrhdxmde(1:nx) = hyb_betae**0.5 * hrhdxmd(1:nx)
          do ii = nx, 1, -1
             do j = nmem, 1, -1
                do i = nx, 1, -1
                   ! P_e^(T/2) H^T R^(-1) * ( H dx - d )
                   grad(nx*j+i) = grad(nx*j+i) + sqrt_loc(k,ii,i) * sqrt_pe(k,ii,j) * hrhdxmde(ii)
                   !grad(nx*j+i) = grad(nx*j+i) + sqrt_pb(ii,i) * hrhdxmde(ii)
                end do
             end do
          end do
          do ii = nx, 1, -1
             do j = nmem, 1, -1
                do i = nx, 1, -1
                   ! M_e^T H^T R^(-1) * ( H dx - d )
                   mhrhdxmdb(i) = mhrhdxmdb(i) + etlm_loc(k,ii,i) * sqrt_pe(k,ii,j) * sqrt_pe(0,i,j) * inv_pe(i) * hrhdxmdb(ii)
                end do
             end do
          end do
       end do
       do ii = nx, 1, -1
          do i = nx, 1, -1
             ! P_b^(T/2) M_e^T H^T R^(-1) * ( H dx - d )
             grad(i) = grad(i) + sqrt_pb(ii,i) * mhrhdxmdb(ii)
          end do
       end do
       costo = costo * 0.5_r_size
       ! L-BFGS
       call lbfgs( nx*(1+nmem), costb+costo, grad(1:nx*(1+nmem)), vctl(1:nx*(1+nmem)), iflg )
       write(6,'(2i5,4es20.10)') iter, iflg, costb+costo, costb, costo, sum(grad(:)**2)
       if( iflg <= 1 ) then
          dx(1:nx) = dxp(1:nx)
       end if
       if( iflg == 0 ) exit iteration
    end do iteration
    ! analysis (at aslot)
    x(1:nx) = x(1:nx) + dx(1:nx)
    ! forecast to the end of assimilation window
    do k = aslot+1, nslot
       call lorenz1996_run( 1, k, dt_cycle/nslot, x(1:nx) )
    end do
    !
    return
  end subroutine etlmvar_run
  !
end module etlmvar
