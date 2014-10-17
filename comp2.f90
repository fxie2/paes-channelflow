subroutine comp2(istart,iend)

  ! CALCULATES THE PRESSURE TERM ON THE RHS OF THE STATE EQUATIONS,
  ! AND ADDS IT IN.  THE FIELDS AT THE NEXT TIME STEP ARE THEN
  ! COMPUTED.
  
  use storage
  use compbuf
  use parameters
  implicit none

  integer istart,iend,iz
  real :: fnt1(nnx,nny)

  do iz=istart,iend

     ! X-DIRECTION PGF (U-WIND COMPONENT)
     call xderiv(p(1:nnx,:,iz),fnt1,xk,nnx,nny)
     ru(:,:,iz) = ru(:,:,iz) - fnt1(:,:)
     
     ! Y-DIRECTION PGF (V-WIND COMPONENT)
     if (nny .ne. 1) then
        call yderiv(p(1:nnx,:,iz),fnt1,yk,nnx,nny)
        rv(:,:,iz) = rv(:,:,iz) - fnt1(:,:)
     endif

     ! Z-DIRECTION PGF (W-WIND COMPONENT)
     if (iz.ne.nnz) then
        rw(:,:,iz) = rw(:,:,iz) - (p(1:nnx,:,iz+1)-p(1:nnx,:,iz))*dzm
        rw(:,:,iz) = rw(:,:,iz) - sum(dble(rw(:,:,iz)))*xym
     endif
        
     ! TIME DERIVATIVE USING THE ADAMS-BASHFORTH SCHEME
     u(:,:,iz) = u(:,:,iz) + dtab1*ru(:,:,iz)
     if (nny.ne.1) v(:,:,iz) = v(:,:,iz) + dtab1*rv(:,:,iz)
     if (iz.ne.nnz) then
        w(:,:,iz) = w(:,:,iz) + dtab1*rw(:,:,iz)
        e(:,:,iz) = max(e(:,:,iz)+dtab1*re(:,:,iz),emin)
        w(:,:,iz) = w(:,:,iz) - sum(dble(w(:,:,iz)))*xym
     endif


  enddo
  
end subroutine comp2











