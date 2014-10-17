
! Change the inversion height for the comparions with Andren et al study.

subroutine initwave(istart, iend)
  ! INITIALIZES THE WAVENUMBER SPACE FOR THE RUN.  THE RESULTS
  ! ARE USED BY SUBROUTINE PRESS DURING THE PRESSURE CALCULATION.
  ! OUTPUT: XK, YK, XK2, WAVEXY
  
  use storage
  use parameters
  implicit none
  integer :: ix, iy, ii, istart, iend, iz
  real :: xnn, temp, temp1, ktop

  do ix=1,nnxb
     xkb(ix)=float(ix-1)*2.0*pi/xl
     if (ix.gt.nnxb/2+1) xkb(ix) = -float(nnxb-ix+1)*2.0*pi/xl
  enddo
  do iy=1,nnyb
     ykb(iy)=float(iy-1)*2.0*pi/yl
     if (iy.gt.nnyb/2+1) ykb(iy) = -float(nnyb-iy+1)*2.0*pi/yl
  enddo

  do ix=1,nnx
     xk(ix)=float(ix-1)*2.0*pi/xl
     if (ix.gt.ncx) xk(ix) = -float(nnx-ix+1)*2.0*pi/xl
  enddo
  do iy=1,nny
     yk(iy)=float(iy-1)*2.0*pi/yl
     if (iy.gt.ncy) yk(iy) = -float(nny-iy+1)*2.0*pi/yl
  enddo

  ii=-1
  do ix=1,ncx
     ii=ii+2
     temp=xk(ix)*xk(ix)
     do iy=1,nny
        temp1=temp+yk(iy)*yk(iy)
        xk2(ii,iy)=temp1
        xk2(ii+1,iy)=temp1
     enddo
  enddo

end subroutine initwave





subroutine initvel(istart,iend)
  ! INITIALIZES THE VELOCITY, TEMPERATURE AND ENERGY FIELDS AT ALL 
  ! LEVELS; THE LOWEST THREE LEVELS OF THE VELOCITY FIELD ARE 
  ! GIVEN A RANDOMIZED COMPONENT.  
  ! OUTPUT: U, V, W, T, E, RU, RV, RW, RT AND RE.

  use storage
  use parameters
  implicit none

  real :: zwall !Distance from wall
  integer :: istart, iend
  real :: ufree, ampv, ampt, vmax, facv, pmax, wind
  real :: psi(nnx,nny), psix(nnx,nny), psiy(nnx,nny), x, y
  integer :: iz, nzi, seed_size, ix, iy, kp, lp
  integer, allocatable :: seed(:)

  utau  = 2.0*vk/log(1/z0)  !Design utau
  ! SET FIELD VALUES BASED ON POSITION RELATIVE TO INVERSION HEIGHT.
  ampv=0.3

  do iz=istart,iend
     
      if(iz.gt.nnz/2) then
         zwall = (nnz-iz+0.5)*dz
      else
         zwall = (iz-0.5)*dz
      end if
      u(:,:,iz) = utau/vk*log(zwall/z0)
      v(:,:,iz) = 0.0
      w(:,:,iz)=0.0

      e(:,:,iz)=emin


      call random_seed(size=seed_size)
      allocate(seed(seed_size))
      seed=10972
      call random_seed(put=seed)
      deallocate(seed)
      
      call random_number(psi)
      psi = psi - sum(dble(psi))*xym
      psi = psi/maxval(abs(psi))
      
      call xderiv(psi,psix,xk,nnx,nny)
      call yderiv(psi,psiy,yk,nnx,nny)
      
      vmax=sqrt(maxval(psix**2+psiy**2))
      facv=ampv/vmax

      if (nny .eq.1) then
         u(:,:,iz) = u(:,:,iz) + psi(:,:)*ampv
      else
         u(:,:,iz) = u(:,:,iz) - psiy(:,:)*facv
         v(:,:,iz) = v(:,:,iz) + psix(:,:)*facv
      endif


        
!       e(:,:,iz) = 1.0
      
  enddo
       
end subroutine initvel

