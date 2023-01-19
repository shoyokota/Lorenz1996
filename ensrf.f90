module ensrf
  use prm, only : r_size, &
       & nx, &
       & nslot, aslot, dt_cycle, &
       & nobs, obsdata, &
       & nmem, sigma_sloc, sigma_tloc, infl
  use lorenz1996, only : lorenz1996_run
  use obsope, only : obsope_run
  implicit none
  private
  public :: ensrf_run
contains
  !
  subroutine ensrf_run( xe, obs )
    implicit none
    real(r_size), intent(inout) :: xe(nx,nmem)   ! ensemble forecasts
    type(obsdata), intent(in) :: obs(nslot,nobs) ! observations
    integer :: i, j, k, kk, l
    integer :: iflg
    real(r_size) :: distance
    real(r_size) :: xfe(nx,nmem)
    real(r_size) :: xm(nslot,nx)
    real(r_size) :: sqrt_pe(nslot,nx,nmem)
    real(r_size) :: loc(nslot,nx,nobs)
    real(r_size) :: hxfe(nmem)
    real(r_size) :: hxm
    real(r_size) :: hxd(nmem)    
    real(r_size) :: hxhxr
    real(r_size) :: kg(nx)
    !
    xfe(1:nx,1:nmem) = xe(1:nx,1:nmem)
    do k = 1, nslot
       ! ensemble forecast
       do j = 1, nmem
          call lorenz1996_run( 1, k, dt_cycle/nslot, xfe(1:nx,j) ) 
       end do
       ! x_m
       xm(k,1:nx) = 0.0_r_size
       do j = 1, nmem
          xm(k,1:nx) = xm(k,1:nx) + xfe(1:nx,j)
       end do
       xm(k,1:nx) = xm(k,1:nx) / nmem
       ! P_e^(1/2)
       do j = 1, nmem
          sqrt_pe(k,1:nx,j) = ( xfe(1:nx,j) - xm(k,1:nx) ) / ( nmem - 1 )**0.5
          sqrt_pe(k,1:nx,j) = sqrt_pe(k,1:nx,j) * ( 1.0_r_size + infl ) ! multiplicative inflation
       end do
    enddo
    ! set localization
    do k = 1, nslot
       do l = 1, nobs
          do i = 1, nx
             distance = abs( obs(k,l)%pos - i )
             if( distance > nx/2 ) distance = nx - distance
             loc(k,i,l) = exp( -0.5_r_size * ( distance / sigma_sloc )**2 )
          end do
       end do
       loc(k,1:nx,1:nobs) = loc(k,1:nx,1:nobs) * exp( -0.5_r_size * ( ( k - aslot ) * dt_cycle / nslot / sigma_tloc )**2 )
    end do
    ! loop for observations
    do k = 1, nslot
       do l = 1, nobs
          ! H x_m, H P_e^(1/2)
          do j = 1, nmem
             xfe(1:nx,j) = xm(k,1:nx) + sqrt_pe(k,1:nx,j)
             call obsope_run( 1, k, 1, obs(k,l), xfe(1:nx,j), hxfe(j) )
          end do
          hxm = 0.0_r_size
          do j = 1, nmem
             hxm = hxm + hxfe(j)
          end do
          hxm = hxm / nmem
          hxd(1:nmem) = hxfe(1:nmem) - hxm 
          ! H P_e^(1/2) (H P_e^(1/2))^T + R
          hxhxr = obs(k,l)%err**2
          do j = 1, nmem
             hxhxr = hxhxr + hxd(j)**2
          end do
          do kk = 1, nslot
             ! K = P_e^(1/2) (H P_e^(1/2))^T [ H P_e^(1/2) (H P_e^(1/2))^T + R ]^(-1)
             kg(1:nx) = 0.0_r_size
             do j = 1, nmem
                kg(1:nx) = kg(1:nx) + sqrt_pe(kk,1:nx,j) * hxd(j)
             end do
             kg(1:nx) = kg(1:nx) / hxhxr
             ! x_m = x_m + K ( y - H x_m )
             xm(kk,1:nx) = xm(kk,1:nx) + kg(1:nx) * loc(kk,1:nx,l) * ( obs(k,l)%val - hxm )
             ! P_e^(1/2) = ( I - alpha K H ) P_e^(1/2)
             kg(1:nx) = kg(1:nx) / ( 1.0_r_size + obs(k,l)%err / sqrt(hxhxr) )
             do j = 1, nmem
                sqrt_pe(kk,1:nx,j) = sqrt_pe(kk,1:nx,j) - kg(1:nx) * loc(kk,1:nx,l) * hxd(j)
             end do
          end do
       end do
    end do
    do j = 1, nmem
       ! analysis (at aslot)
       xe(1:nx,j) = xm(aslot,1:nx) + sqrt_pe(aslot,1:nx,j) * ( nmem - 1 )**0.5
       ! forecast to the end of assimilation window
       do k = aslot+1, nslot
          call lorenz1996_run( 1, k, dt_cycle/nslot, xe(1:nx,j) )
       end do
    end do
    !
    return
  end subroutine ensrf_run
  !
end module ensrf
