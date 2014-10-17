module fftstorage
  integer                   :: nx, ny, nxb, nyb
  integer*8, target           :: psf, psb, pbf, pbb
  integer*8, target           :: pxsf, pxsb, pxbf, pxbb
  integer*8, target           :: pysf, pysb, pybf, pybb
  integer*8, pointer          :: pf, pb
  real, target              :: nxym1s, nxym1b, nxm1s, nxm1b, nym1s, nym1b 
  real, pointer             :: nxym1, nxm1, nym1
end module fftstorage

!Links to /usr/local/lib/librfftw.a and /usr/local/libsfftw.a on raga




subroutine fft_init(nnx,nny)
  use fftstorage
  implicit none
  integer, intent(IN) :: nnx,nny

  integer FFTW_FORWARD,FFTW_BACKWARD
  parameter (FFTW_FORWARD=-1,FFTW_BACKWARD=1)
  
  integer FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL
  parameter (FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1)
     
  integer FFTW_ESTIMATE,FFTW_MEASURE
  parameter (FFTW_ESTIMATE=0,FFTW_MEASURE=1)
     
  integer FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM
  parameter (FFTW_OUT_OF_PLACE=0)
  parameter (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)
     
  integer FFTW_THREADSAFE
  parameter (FFTW_THREADSAFE=128)

  nx = nnx
  ny = nny
  nxm1s = 1./real(nnx)
  nym1s = 1./real(nny)
  nxym1s = 1./real(nnx)/real(nny)

  nxb = 3*nnx/2
  nyb = 3*nny/2
  nxm1b = 1./real(nxb)
  nym1b = 1./real(nyb)
  nxym1b = 1./real(nxb)/real(nyb)
 
  call rfftw2d_f77_create_plan(psf, nnx, nny, FFTW_REAL_TO_COMPLEX, &
       FFTW_ESTIMATE+FFTW_IN_PLACE)
  call rfftw2d_f77_create_plan(psb, nnx, nny, FFTW_COMPLEX_TO_REAL, &
       FFTW_ESTIMATE+FFTW_IN_PLACE)
  call rfftw2d_f77_create_plan(pbf, nxb, nyb, FFTW_REAL_TO_COMPLEX, &
       FFTW_ESTIMATE+FFTW_IN_PLACE)
  call rfftw2d_f77_create_plan(pbb, nxb, nyb, FFTW_COMPLEX_TO_REAL, &
       FFTW_ESTIMATE+FFTW_IN_PLACE)

  call rfftw_f77_create_plan(pxsf, nnx, FFTW_REAL_TO_COMPLEX, &
       FFTW_ESTIMATE+FFTW_IN_PLACE)
  call rfftw_f77_create_plan(pxsb, nnx, FFTW_COMPLEX_TO_REAL, &
       FFTW_ESTIMATE+FFTW_IN_PLACE)
  call rfftw_f77_create_plan(pxbf, nxb, FFTW_REAL_TO_COMPLEX, &
       FFTW_ESTIMATE+FFTW_IN_PLACE)
  call rfftw_f77_create_plan(pxbb, nxb, FFTW_COMPLEX_TO_REAL, &
       FFTW_ESTIMATE+FFTW_IN_PLACE)

  call rfftw_f77_create_plan(pysf, nny, FFTW_REAL_TO_COMPLEX, &
       FFTW_ESTIMATE+FFTW_IN_PLACE)
  call rfftw_f77_create_plan(pysb, nny, FFTW_COMPLEX_TO_REAL, &
       FFTW_ESTIMATE+FFTW_IN_PLACE)
  call rfftw_f77_create_plan(pybf, nyb, FFTW_REAL_TO_COMPLEX, &
       FFTW_ESTIMATE+FFTW_IN_PLACE)
  call rfftw_f77_create_plan(pybb, nyb, FFTW_COMPLEX_TO_REAL, &
       FFTW_ESTIMATE+FFTW_IN_PLACE)
end subroutine fft_init


      




subroutine fft_finalize()
  use fftstorage
  implicit none

  call rfftwnd_f77_destroy_plan(psf)
  call rfftwnd_f77_destroy_plan(psb)
  call rfftwnd_f77_destroy_plan(pbf)
  call rfftwnd_f77_destroy_plan(pbb)

  call rfftw_f77_destroy_plan(pxsf)
  call rfftw_f77_destroy_plan(pxsb)
  call rfftw_f77_destroy_plan(pxbf)
  call rfftw_f77_destroy_plan(pxbb)

  call rfftw_f77_destroy_plan(pysf)
  call rfftw_f77_destroy_plan(pysb)
  call rfftw_f77_destroy_plan(pybf)
  call rfftw_f77_destroy_plan(pybb)
end subroutine fft_finalize

      



subroutine ffthrz(fa, n1, n2, isign)
  use fftstorage
  implicit none

  integer, intent(IN) :: n1, n2, isign
  real, intent(INOUT) :: fa(n1+2,n2)
  
  if (n1.eq.nx) then
     pf => psf
     pb => psb
     nxym1 => nxym1s
  else
     pf => pbf
     pb => pbb
     nxym1 => nxym1b
  endif
     
  if (isign.ge.0) then
     call rfftwnd_f77_one_complex_to_real(pb,fa)
  else 
     call rfftwnd_f77_one_real_to_complex(pf,fa)
     fa = fa*nxym1
  endif
end subroutine ffthrz







SUBROUTINE XDERIV(AIN,AX,XK,N1,N2)

  ! Determines the spatial derivative in the x-direction,
  ! by taking the FFT, multiplying by wavenumber, and finding
  ! the inverse FFT.  A real, trigonometric FFT is used.

  ! Input:
  ! A  -  Matrix to differentiate.
  ! N1, N2  - Dimensions of A and AX.
  ! DELX  - Spacing in x-direction.

  ! Output:
  ! AX -  Derivative of A, in the x-direction.

  USE FFTSTORAGE
  IMPLICIT NONE

  integer, intent(IN) :: n1, n2
  real, intent(IN)    :: ain(n1,n2), xk(n1)
  real, intent(OUT)   :: ax(n1,n2)
  real                :: coef, scratch(n1)
  integer             :: ix,iy

  if (n1.eq.nx) then
     pf => pxsf
     pb => pxsb
     nxm1 => nxm1s
  else
     pf => pxbf
     pb => pxbb
     nxm1 => nxm1b
  endif

  ax = ain
  call rfftw_f77(pf, n2, ax, 1, n1, scratch, 0, 0)

  do iy=1,n2
     do ix=2,n1/2
        coef = ax(ix,iy)
        ax(ix,iy) = xk(n1-ix+2)*ax(n1-ix+2,iy)
        ax(n1-ix+2,iy) = xk(ix)*coef
     enddo
  enddo

  ax(1,:) = 0.0
  ax(n1/2+1,:) = 0.0
  call rfftw_f77(pb, n2, ax, 1, n1, scratch, 0, 0)
  ax = ax*nxm1

end subroutine xderiv






subroutine yderiv(ain,ay,yk,n1,n2)

  ! Determines the spatial derivative in the y-direction,
  ! by taking the FFT, multiplying by wavenumber, and finding
  ! the inverse FFT.  A real, trigonometric FFT is used.

  ! Input:
  ! A  -  Matrix to differentiate.
  ! N1, N2  - Dimensions of A and AY.
  ! DELY  - Spacing in y-direction.

  ! Output:
  ! AY -  Derivative of A, in the y-direction.

  use fftstorage
  implicit none

  integer, intent(in) :: n1, n2
  real, intent(in)    :: ain(n1,n2), yk(n2)
  real, intent(out)   :: ay(n1,n2)
  real                :: coef, scratch(n2)
  integer             :: ix, iy

  if (ny .eq. 1) then
     ay = 0.0
  else

  if (n1.eq.nx) then
     pf => pysf
     pb => pysb
     nym1 => nym1s
  else
     pf => pybf
     pb => pybb
     nym1 => nym1b
  endif

  ay = ain
  call rfftw_f77(pf, n1, ay, n1, 1, scratch, 0, 0)
      
  do iy=2,n2/2
     do ix=1,n1
        coef = ay(ix,iy)
        ay(ix,iy) = yk(n2-iy+2)*ay(ix,n2-iy+2)
        ay(ix,n2-iy+2) = yk(iy)*coef
     enddo
  enddo

  ay(:,1) = 0.0
  ay(:,n2/2+1) = 0.0
  call rfftw_f77(pb, n1, ay, n1, 1, scratch, 0, 0)
  ay = ay*nym1

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


 
