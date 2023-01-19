module matrix
  use prm, only : r_size
  implicit none
  real(r_size), parameter :: eps = 1.0e-10
  private
  public :: matrix_eigen, matrix_sym, matrix_gauss
contains
  !
  subroutine matrix_eigen( iflg, n, A )
    implicit none
    integer, intent(in) :: iflg           ! 1:eigenvalue,2:eigenvector,3:inverse,4:square-root,5:inverse-square-root
    integer, intent(in) :: n              ! dimension
    real(r_size), intent(inout) :: A(n,n) ! matrix
    integer :: i, j, k, info
    real(r_size) :: eigenvector(n,n)
    real(r_size) :: eigenvalue(n)
    real(r_size) :: work(n*66)
    !
    eigenvector(1:n,1:n) = A(1:n,1:n)
    call dsyev('V','U',n,eigenvector,n,eigenvalue,work,n*66,info)
    !
    A(1:n,1:n) = 0.0_r_size
    if( iflg == 1 ) then      ! eigenvalue
       do i = 1, n
          A(i,i) = eigenvalue(i)
       end do
    else if( iflg == 2 ) then ! eigenvector
       do i = 1, n
          do j = 1, n
             A(i,j) = eigenvector(i,j)
          end do
       end do
    else if( iflg == 3 ) then ! inverse
       do k = 1, n
          if( abs(eigenvalue(k)) > eps ) then
             do i = 1, n
                do j = 1, n
                   A(i,j) = A(i,j) + eigenvector(i,k) * eigenvector(j,k) / eigenvalue(k)
                end do
             end do
          end if
       end do
    else if( iflg == 4 ) then ! square-root
       do k = 1, n
          if( eigenvalue(k) > eps ) then
             do i = 1, n
                do j = 1, n
                   A(i,j) = A(i,j) + eigenvector(i,k) * eigenvector(j,k) * sqrt(eigenvalue(k))
                end do
             end do
          end if
       end do
    else if( iflg == 5 ) then ! inverse-square-root
       do k = 1, n
          if( eigenvalue(k) > eps ) then
             do i = 1, n
                do j = 1, n
                   A(i,j) = A(i,j) + eigenvector(i,k) * eigenvector(j,k) / sqrt(eigenvalue(k))
                end do
             end do
          end if
       end do
    end if
    !
    return
  end subroutine matrix_eigen
  !
  subroutine matrix_sym( n, A )
    implicit none
    integer, intent(in) :: n              ! dimension
    real(r_size), intent(inout) :: A(n,n) ! matrix
    integer :: i, j
    real(r_size) :: amean
    !
    do i = 1, n
       do j = i+1, n
          amean = ( A(i,j) + A(j,i) ) * 0.5_r_size
          A(i,j) = amean
          A(j,i) = amean
       end do
    end do
    !
    return
  end subroutine matrix_sym
  !
  subroutine matrix_gauss( n, sigma, A )
    implicit none
    integer, intent(in) :: n              ! dimension
    real(r_size), intent(in) :: sigma     ! length scale
    real(r_size), intent(inout) :: A(n,n) ! matrix
    integer :: i, j, d
    !
    do i = 1, n
       A(i,i) = 1.0_r_size
       do j = i+1, n
          d = min( j - i, i + n - j )
          A(i,j) = exp( -0.5_r_size * ( d / sigma )**2 )
          A(j,i) = A(i,j)
       end do
    end do
    !
    return
  end subroutine matrix_gauss
  !
end module matrix
