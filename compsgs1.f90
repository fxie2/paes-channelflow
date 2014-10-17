subroutine compsgs1(kb,iz,j,jp1,k,km1,l,lm1,lp1,nm1,n,np1,np2,istart)
  
  ! MODELS THE SUBGRID-SCALE STRESS TERMS

  use storage
  use compbuf
  use parameters
  implicit none

  ! INPUT PARAMETERS:
  integer iz,j,jp1,k,km1,l,lp1,lm1,istart,kb,nm1,n,np1,np2
   real :: average

  ! LOCAL VARIABLES:
  real :: temp(nnxb,nnyb)

  ! COMPUTE LENGTH SCALES FOR EACH SLAB.
  ! THE LENGTH SCALE DEPENDS ON ENERGY AND STABILITY
  ! FOR STABLE CONDITIONS; OTHERWISE IT IS JUST THE
  ! LENGTH SCALE FOR A VOLUME ELEMENT.
  
  if (iz.eq.istart) then

     call complen(lm1, kb-1, nm1, n, iz-1)
     call complen(l, kb, n, np1,iz)
     tauww(:,:,j) = -(xkk(:,:,lm1)+xkk(:,:,l))*wz(:,:,n)

     if (iz.eq.1) then
        call sufsim(k,km1)
     else
        tauuw(:,:,km1) = -xkk(:,:,lm1)*(uz(:,:,lm1)+wx(:,:,lm1))
        tauvw(:,:,km1) = -xkk(:,:,lm1)*(vz(:,:,lm1)+wy(:,:,lm1))
     endif
  endif

  if (iz.ne.nnz) then
     call complen(lp1, kb+1, np1, np2, iz+1)
  endif

  ! FIND TURBULENT STRESSES ON THE CURRENT LEVEL.
  tauuu = -(xkk(:,:,lm1)+xkk(:,:,l))*ux(:,:,n)
  tauvv = -(xkk(:,:,lm1)+xkk(:,:,l))*vy(:,:,n)
  tauuv = -(xkk(:,:,lm1)+xkk(:,:,l))*(uy(:,:,n)+vx(:,:,n))*0.5
  if (iz.eq.nnz) then
     tauww(:,:,jp1) = 0.0
     call sufsim_upper(k,km1)
  else
     tauww(:,:,jp1) = -(xkk(:,:,l)+xkk(:,:,lp1))*wz(:,:,np1)
     tauuw(:,:,k) = -xkk(:,:,l)*(uz(:,:,l)+wx(:,:,l))
     tauvw(:,:,k) = -xkk(:,:,l)*(vz(:,:,l)+wy(:,:,l))
  endif
  
end subroutine compsgs1









