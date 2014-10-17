subroutine complen(l,kb,n,np1,iz)
  use storage
  use compbuf
  use parameters
  implicit none
  integer l, kb, n, np1, iz

  ! THIS ROUTINE COMPUTES THE LOCAL SHEAR AND SGS DIFFUSIVITIES. THE
  ! SGS LENGTH SCALE IS ALSO COMPUTED, BASED ON EITHER MOENG (1984) 
  ! OR MASON AND DERBYSHIRE (1990).

  ! FIRST COMPUTE THE SHEAR S_ij*S_ji AND STABILITY
  shear(:,:,l) = ux(:,:,n)**2+ux(:,:,np1)**2 + &
       vy(:,:,n)**2+vy(:,:,np1)**2 + wz(:,:,n)**2+wz(:,:,np1)**2 + &
       0.5*(uy(:,:,n)+vx(:,:,n))**2 + 0.5*(uy(:,:,np1)+vx(:,:,np1))**2 + &
       (uz(:,:,l)+wx(:,:,l))**2 + (vz(:,:,l)+wy(:,:,l))**2

  ! NOW COMPUTE THE SGS DIFFUSIVITIES FOR MOMENTUM AND HEAT FOR THIS SLAB
  SELECT CASE(imodel)
  CASE(3)	!Moeng 84 model: vt=0.1*l*e**0.5
  xkk(:,:,l) = min(0.1*alk(:,:,l)*sqrt(ebuf(:,:,kb)), xkmax)
  CASE(4)	!Smag model Csmag=0.1: vt=(cs*l)^2*(2S_ijS_ij)**0.5
  xkk(:,:,l) = min((0.1*dsl)**2*sqrt(shear(:,:,l)), xkmax)
  END SELECT

!  write(*,*) 'Eddy Visc at level iz = ', iz, ' = ', sum(xkk(:,:,l))*xymb

end subroutine complen







