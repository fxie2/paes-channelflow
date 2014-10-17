subroutine derivs(kb,iz,istart,l,lp1,np2)
  
  use storage
  use compbuf
  use parameters
  implicit none
  integer :: iz, istart, kb, l, lp1, np2

  if(iz .lt. (nnz-1))  then
     call xderiv(ubuf(1,1,kb+2),ux(1,1,np2),xkb,nnxb,nnyb)
     call yderiv(ubuf(1,1,kb+2),uy(1,1,np2),ykb,nnxb,nnyb)
     call xderiv(vbuf(1,1,kb+2),vx(1,1,np2),xkb,nnxb,nnyb)
     call yderiv(vbuf(1,1,kb+2),vy(1,1,np2),ykb,nnxb,nnyb)
     wz(:,:,np2) = (wbuf(:,:,kb+2)-wbuf(:,:,kb+1))*dzm
     uz(:,:,lp1) = (ubuf(:,:,kb+2)-ubuf(:,:,kb+1))*dzm
     vz(:,:,lp1) = (vbuf(:,:,kb+2)-vbuf(:,:,kb+1))*dzm
     call xderiv(wbuf(1,1,kb+1),wx(:,:,lp1),xkb,nnxb,nnyb)
     call yderiv(wbuf(1,1,kb+1),wy(:,:,lp1),ykb,nnxb,nnyb)

     call xderiv(ebuf(1,1,kb),ex,xkb,nnxb,nnyb)
     call yderiv(ebuf(1,1,kb),ey,ykb,nnxb,nnyb)
  else
     ux(:,:,np2) = 0.0
     uy(:,:,np2) = 0.0
     vx(:,:,np2) = 0.0
     vy(:,:,np2) = 0.0
     uz(:,:,lp1) = uz(:,:,l)
     vz(:,:,lp1) = vz(:,:,l)
     wz(:,:,np2) = 0.0
     wx(:,:,lp1) = 0.0
     wy(:,:,lp1) = 0.0
  end if

!   if(iz .gt. 60) then
!      write(*,*) 'iz      ux        uy       uz        vx      vy      vz     wx    wy     wz'
!      write(*,*) iz, sum(ux(:,:,np2))*xymb, sum(uy(:,:,np2))*xymb, &
!           sum(uz(:,:,lp1))*xymb, sum(vx(:,:,np2))*xymb, sum(vy(:,:,np2))*xymb, &
!           ' ', sum(vz(:,:,lp1))*xymb, sum(wx(:,:,lp1))*xymb, &
!           sum(wy(:,:,lp1))*xymb, sum(wz(:,:,np2))*xymb
!   end if

end subroutine derivs







