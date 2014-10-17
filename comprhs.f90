subroutine comprhs(istart,iend)

  ! COMPUTES THE RHS OF THE POISSON EQ OF PRESSURE

  use storage
  use compbuf
  use parameters
  implicit none

  ! INPUT PARAMETERS
  integer istart, iend

  ! LOCAL VARIABLES
  real :: fnt1(nnx,nny), fnt2(nnx,nny)
  integer iz, itop

  do iz = istart,iend
     call xderiv(ru(:,:,iz)+u(:,:,iz)*dtabm,fnt1,xk,nnx,nny)
     if (nny .eq. 1) then
        fnt2 = 0.0
     else
        call yderiv(rv(:,:,iz)+v(:,:,iz)*dtabm,fnt2,yk,nnx,nny)
     endif
     if (iz.eq.1) then
        p(1:nnx,:,iz) = fnt1(:,:) + fnt2(:,:) + dzm*(w(:,:,iz) &
             *dtabm+rw(:,:,iz))
     elseif(iz.eq.nnz) then
        p(1:nnx,:,iz) = fnt1(:,:) + fnt2(:,:) - dzm*(w(:,:,iz-1) &
             *dtabm+rw(:,:,iz-1))
     else
        p(1:nnx,:,iz) = fnt1(:,:) + fnt2(:,:) + ((w(:,:,iz)-w(:,:,iz-1)) &
             *dtabm + rw(:,:,iz)-rw(:,:,iz-1))*dzm
     endif

  enddo

end subroutine comprhs







