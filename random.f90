module random
  use prm, only : r_size, pi
  implicit none
  integer, allocatable, save :: seed(:)
  private
  public :: init_random, get_random, get_gauss_random
contains
  !
  subroutine init_random( s )
    implicit none
    integer, intent(in) :: s
    integer :: nseed, i
    !
    call random_seed( size=nseed )
    allocate( seed(nseed) )
    do i = 1, nseed
       seed(i) = s + 37 * (i-1)
    end do
    call random_seed( put=seed )
    !
    return
  end subroutine init_random
  !
  subroutine get_random( n, x )
    implicit none
    integer, intent(in) :: n
    real(r_size), intent(out) :: x(n)
    !
    call random_number( x )
    !
    return
  end subroutine get_random
  !
  subroutine get_gauss_random( n, x )
    implicit none
    integer, intent(in) :: n
    real(r_size), intent(out) :: x(n)
    real(r_size) :: tmp(3)
    integer :: i
    !
    do i = 1, n
       call random_number( tmp )
       ! Gaussian random number obtained by Box-Muller
       if( tmp(3) >= 0.5_r_size ) then
          x(i) = sqrt( -2.0_r_size * log( tmp(1) ) ) * sin( 2.0_r_size * pi * tmp(2) )
       else
          x(i) = sqrt( -2.0_r_size * log( tmp(1) ) ) * cos( 2.0_r_size * pi * tmp(2) )
       end if
    end do
    !
    return
  end subroutine get_gauss_random
  !
end module random
