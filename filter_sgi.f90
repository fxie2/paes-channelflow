module fftstorage
  real scale
  real, pointer :: trigsx(:), trigsy(:), table(:)
  real, allocatable, target :: trigsxs(:), trigsys(:), tables(:)
  real, allocatable, target :: trigsxb(:), trigsyb(:), tableb(:)

  integer :: nx, ny, nxb, nyb
  real, pointer :: nxm1, nym1, nxym1
  real, target  :: nxm1s, nym1s, nxym1s
  real, target  :: nxm1b, nym1b, nxym1b
end module fftstorage






subroutine fft_init(nnx,nny)
  use fftstorage
  implicit none
  integer :: nnx,nny

  nx = nnx
  ny = nny
  nxm1s = 1./float(nnx)
  nym1s = 1./float(nny)
  nxym1s = 1./float(nnx)/float(nny)

  nxb = 3*nnx/2
  nyb = 3*nny/2
  nxm1b = 1./float(nxb)
  nym1b = 1./float(nyb)
  nxym1b = 1./float(nxb)/float(nyb)
  
  allocate(trigsxs(nx+15))
  allocate(trigsys(ny+15))
  allocate(tables((nx+15)+2*(ny+15)))

  allocate(trigsxb(nxb+15))
  allocate(trigsyb(nyb+15))
  allocate(tableb((nxb+15)+2*(nyb+15)))

  call scfftm1dui(nx,trigsxs)
  if (ny.ne.1) then
     call scfftm1dui(ny,trigsys)
     call scfft2dui(nx,ny,tables)
  endif
  
  call scfftm1dui(nxb,trigsxb)
  if (ny.ne.1) then
     call scfftm1dui(nyb,trigsyb)
     call scfft2dui(nxb,nyb,tableb)
  endif
end subroutine fft_init




      
      

subroutine fft_finalize()
  use fftstorage
  implicit none
 
  deallocate(trigsxs)
  deallocate(trigsys)
  deallocate(tables)
  deallocate(trigsxb)
  deallocate(trigsyb)
  deallocate(tableb)
end subroutine fft_finalize
 

      
      
      
 

subroutine ffthrz(fa,n1,n2,isign)
  use fftstorage
  implicit none
  integer :: isign, n1, n2
  real :: fa(n1+2,n2)

  if (ny.eq.1) then
     if (n1.eq.nx) then
        nxm1 => nxm1s
        trigsx => trigsxs
     else
        nxm1 => nxm1b
        trigsx => trigsxb
     endif
     
     if (isign.eq.-1) then
        call scfftm1du(-1,n1,n2,fa,1,n1+2,trigsx)
        scale = nxm1
     else
        call csfftm1du(+1,n1,n2,fa,1,n1+2,trigsx)
        fa(1:n1,:) = fa(1:n1,:)*scale
     endif
     
  else

     if (n1.eq.nx) then
        nxym1 => nxym1s
        table => tables
     else
        nxym1 => nxym1b
        table => tableb
     endif
     
     if (isign.eq.-1) then
        call scfft2du(-1,n1,n2,fa,n1+2,table)
        scale = nxym1
     else
        call csfft2du(+1,n1,n2,fa,n1+2,table)
        fa(1:n1,:) = fa(1:n1,:)*scale
     endif
  endif
end subroutine ffthrz






subroutine xderiv(ain,ax,xk,n1,n2)
  ! DETERMINES THE SPATIAL DERIVATIVE IN THE X-DIRECTION,
  ! BY TAKING THE FFT, MULTIPLYING BY WAVENUMBER, AND FINDING
  ! THE INVERSE FFT.  A REAL, TRIGONOMETRIC FFT IS USED.
  !
  ! INPUT:
  ! A  -  MATRIX TO DIFFERENTIATE.
  ! NX, NY  - DIMENSIONS OF A AND AX.
  ! DELX  - SPACING IN X-DIRECTION.
  !
  ! OUTPUT:
  ! AX -  DERIVATIVE OF A, IN THE X-DIRECTION.

  use fftstorage
  implicit none
  
  integer :: ix,iy,n1,n2
  real :: xk(n1)
  real :: coef
  real :: ax(n1,n2),ain(n1,n2)
  real :: temp(n1+2,n2)

  if (n1.eq.nx) then
     nxm1 => nxm1s
     trigsx => trigsxs
  else
     nxm1 => nxm1b
     trigsx => trigsxb
  endif

  temp(1:n1,:)=ain(:,:)
  call scfftm1du(-1, n1, n2, temp, 1, n1+2, trigsx)

  do iy=1,n2
     do ix=1,n1-1,2 
        coef = temp(ix,iy)
        temp(ix,iy)  = -xk((ix+1)/2)*temp(ix+1,iy)
        temp(ix+1,iy)=  xk((ix+1)/2)*coef
     enddo
  enddo

  temp(n1+1:n1+2,:)=0.0
  call csfftm1du(1, n1, n2, temp, 1, n1+2, trigsx)
  call sscalm1d(n1, n2, nxm1, temp, 1, n1+2)

  ax(:,:)=temp(1:n1,:)
end subroutine xderiv






subroutine yderiv(ain,ay,yk,n1,n2)
  ! Determines the spatial derivative in the y-direction,
  ! by taking the FFT, multiplying by wavenumber, and finding
  ! the inverse FFT.  A real, trigonometric FFT is used.
  !
  ! Input:
  ! A  -  Matrix to differentiate.
  ! NX, NY  - Dimensions of A and AY.
  ! DELY  - Spacing in y-direction.
  !
  ! Output:
  ! AY -  Derivative of A, in the y-direction.

  use fftstorage
  implicit none

  integer :: ix,iy,n1,n2
  real :: yk(n2)
  real :: coef
  real :: ay(n1,n2), ain(n1,n2)
  real :: temp(n1,n2+2)

  if (ny .eq. 1) then
     ay = 0.0
  else

  if (n1.eq.nx) then
     nym1 => nym1s
     trigsy => trigsys
  else
     nym1 => nym1b
     trigsy => trigsyb
  endif

  temp(:,1:n2)=ain(:,:)
  call scfftm1du(-1, n2, n1, temp, n1, 1, trigsy)

  do iy=1,n2-1,2
     do ix=1,n1
        coef = temp(ix,iy)
        temp(ix,iy)  =-yk((iy+1)/2)*temp(ix,iy+1)
        temp(ix,iy+1)= yk((iy+1)/2)*coef
     enddo
  enddo

  temp(:,n2+1:n2+2)=0.
  call csfftm1du(1, n2, n1, temp, n1, 1, trigsy)
  call sscalm1d(n2, n1, nym1, temp, n1, 1)
  
  ay(:,:)=temp(:,1:n2)

  endif
end subroutine yderiv








subroutine expand(xsmall, xlarge)
  use fftstorage
  implicit none

  real, intent(IN)  :: xsmall(nx,ny)
  real, intent(OUT) :: xlarge(nxb,nyb)
  real              :: zsmall(nx+2,ny), zlarge(nxb+2,nyb)
  
  zsmall(1:nx,:) = xsmall
  call ffthrz(zsmall,nx,ny,-1)

  if (ny .eq.1) then
     zlarge(1:nx,1) = zsmall(1:nx,1)
     zlarge(nx+1:nxb+2,1) = 0.0
  else
     zlarge(1:nx,1:ny/2) = zsmall(1:nx,1:ny/2)
     zlarge(1:nx,ny+2:nyb) = zsmall(1:nx,ny/2+2:ny)
     zlarge(nx+1:nxb+2,:) = 0.0
     zlarge(1:nx,ny/2+1:ny+1) = 0.0
  endif

  call ffthrz(zlarge,nxb,nyb,+1)
  xlarge = zlarge(1:nxb,:)

end subroutine expand
  




subroutine compress(xlarge, xsmall)
  use fftstorage
  implicit none

  real, intent(OUT) :: xsmall(nx,ny)
  real, intent(IN)  :: xlarge(nxb,nyb)
  real              :: zsmall(nx+2,ny), zlarge(nxb+2,nyb)

  zlarge(1:nxb,:) = xlarge
  call ffthrz(zlarge,nxb,nyb,-1)

  if (ny .eq. 1) then
     zsmall(1:nx,1) = zlarge(1:nx,1)
     zsmall(nx+1:nx+2,1) = 0.0
  else
     zsmall(1:nx,1:ny/2) = zlarge(1:nx,1:ny/2)
     zsmall(1:nx,ny/2+2:ny) = zlarge(1:nx,ny+2:nyb)
     zsmall(nx+1:nx+2,:) = 0.0
     zsmall(1:nx,ny/2+1) = 0.0
  endif
 
  call ffthrz(zsmall,nx,ny,+1)
  xsmall = zsmall(1:nx,:)
end subroutine compress





subroutine deal(xin,nnx,nny)
  implicit none

  integer, intent(IN) :: nnx, nny
  real, intent(INOUT) :: xin(nnx,nny)
  real                :: buf(nnx+2,nny)

  buf(1:nnx,:) = xin
  call ffthrz(buf,nnx,nny,-1)
  buf(nnx+1:nnx+2,:) = 0.0
  if (nny .ne. 1) buf(1:nnx,nny/2+1) = 0.0
  call ffthrz(buf,nnx,nny,+1)
  xin = buf(1:nnx,:)
end subroutine deal


 
