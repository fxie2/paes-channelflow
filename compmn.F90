subroutine compmn(istart,iend)
  
  use storage
  use compbuf
  use parameters
  implicit none
  
  real average
  integer iz,istart,iend,i
  real, pointer, dimension(:,:) :: ulocal, vlocal, wlocal
  real, dimension(nnx,nny) :: plocal, p1, p2, fnt1, fnt2

  ulocal => rubuf(1:nnx,1:nny)
  vlocal => rvbuf(1:nnx,1:nny)
  wlocal => rwbuf(1:nnx,1:nny)

  ! FIND THE HORIZONTAL-MEAN PROFILES.
  do iz=istart,iend
     i = iz-istart+1
     uxym(i) = average(u(:,:,iz))
     vxym(i) = average(v(:,:,iz))
     wxym(i) = average(w(:,:,iz)) ! Added by Ganesh
!     fnt1(:,:) = sqrt( (u(:,:,iz) + ugal)**2 + (v(:,:,iz) + vgal)**2 ) !zhou
!     usm(i) = average(fnt1(:,:))     !zhou
     if (isbgr) engsbz(i) = average(e(:,:,iz))
  enddo
     
  ! FIND MEAN PROFILES OF VARIOUS MOMENTS
  do iz=istart,iend
     i = iz-istart+1

     if (iz.ne.iend) then
        ulocal = 0.5*((u(:,:,iz)-uxym(i))+(u(:,:,iz+1))-uxym(i+1))
        vlocal = 0.5*((v(:,:,iz)-vxym(i))+(v(:,:,iz+1))-vxym(i+1))
     else
        ulocal = 0.5*((u(:,:,iz)-uxym(i))+(u(:,:,iz+1))-average(u(:,:,iz+1)))
        vlocal = 0.5*((v(:,:,iz)-vxym(i))+(v(:,:,iz+1))-average(v(:,:,iz+1)))
     endif
     
     p1 = p(1:nnx,:,iz) - (e(:,:,iz-1)+e(:,:,iz))/3.  !for div-form :change

     if (iz.eq.nnz) p(1:nnx,:,iz+1) = p(1:nnx,:,iz)

     p2 = p(1:nnx,:,iz+1) - (e(:,:,iz)+e(:,:,iz+1))/3. !for div-from:change

     p1 = p1 - average(p1)
     p2 = p2 - average(p2)
     plocal = 0.5*(p1 + p2)
     wlocal = 0.5*(w(:,:,iz-1)+w(:,:,iz))
     
     uule(i) = average(ulocal*ulocal)
     vvle(i) = average(vlocal*vlocal)
     wwle(i) = average(w(:,:,iz)*w(:,:,iz))
     uwle(i) = average(ulocal*w(:,:,iz))
     vwle(i) = average(vlocal*w(:,:,iz))
     uvle(i) = average(ulocal*vlocal)
     englez(i) = (uule(i)+vvle(i)+wwle(i))*0.5

!      if((iz.lt.3).or.(iz.gt.62)) then
!         write(*,*)iz, uule(i), vvle(i), englez(i)
!      end if
  enddo

end subroutine compmn



function average(x) result(result)
  
  use storage
  use parameters
  
  implicit none
  real :: x(nnx,nny), result
  result = sum(dble(x))*xym
end function average
  
function averageb(x) result(resultb)
  
  use storage
  use parameters
  
  implicit none
  real :: x(nnxb,nnyb), resultb
  resultb = sum(dble(x))*xymb
end function averageb
