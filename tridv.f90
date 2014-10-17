subroutine tridv(bdar,ny,n)
  !**********************************************************************
  !
  !     purpose: tridiagonal matrix solver
  !              with multiple vectors
  !
  !     input:   n   size of a,b,d and r
  !              b   below diagonal elements (b(1) not used)
  !              d   diagonal elements
  !              a   above diagonal elements (a(n) not used)
  !              r   right hand side
  !              ny  number of input vectors
  !
  !     output:  r   solution vector
  !
  !**********************************************************************
  implicit none
  integer :: ny, n
  real, target  :: bdar(ny,n,4)
  real, pointer :: b(:,:), d(:,:), a(:,:), r(:,:)
  real :: fac
  integer :: i, j, ii, nm1

  b => bdar(:,:,1)
  d => bdar(:,:,2)
  a => bdar(:,:,3)
  r => bdar(:,:,4)
  
  if(n .le. 1 ) then
     do j=1,ny
        r(j,1) = r(j,1)/d(j,1)
     enddo
  else
     do j=1,ny
        d(j,1) = 1.0/d(j,1)
     enddo
     do i=2,n
        do j=1,ny
           fac = b(j,i)*d(j,i-1)
           d(j,i) = 1.0/(d(j,i) - fac*a(j,i-1))
           r(j,i) = r(j,i) - fac*r(j,i-1)
        enddo
     enddo
     do j=1,ny
        r(j,n) = r(j,n)*d(j,n)
     enddo
     nm1 = n - 1
     do ii=1,nm1
        i = n - ii
        do j=1,ny
           r(j,i) = d(j,i)*(r(j,i) - a(j,i)*r(j,i+1))
        enddo
     enddo
  endif
end subroutine tridv


