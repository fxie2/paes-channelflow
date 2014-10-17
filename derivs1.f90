subroutine derivsbot(kb,iz,lm1,l,nm1,n,np1)

  use storage
  use compbuf
  use parameters
  implicit none
  integer :: kb,iz,lm1,l,nm1,n,np1
  
  call xderiv(ubuf(1,1,kb),ux(1,1,n),xkb,nnxb,nnyb)
  call xderiv(ubuf(1,1,kb+1),ux(1,1,np1),xkb,nnxb,nnyb)
  
  call yderiv(ubuf(1,1,kb),uy(1,1,n),ykb,nnxb,nnyb)
  call yderiv(ubuf(1,1,kb+1),uy(1,1,np1),ykb,nnxb,nnyb)
  
  call xderiv(vbuf(1,1,kb),vx(1,1,n),xkb,nnxb,nnyb)
  call xderiv(vbuf(1,1,kb+1),vx(1,1,np1),xkb,nnxb,nnyb)
  
  call yderiv(vbuf(1,1,kb),vy(1,1,n),ykb,nnxb,nnyb)
  call yderiv(vbuf(1,1,kb+1),vy(1,1,np1),ykb,nnxb,nnyb)
  
  wz(:,:,n) = (wbuf(:,:,kb)-wbuf(:,:,kb-1))*dzm
  wz(:,:,np1) = (wbuf(:,:,kb+1)-wbuf(:,:,kb))*dzm

  uz(:,:,l)=(ubuf(:,:,kb+1)-ubuf(:,:,kb))*dzm
  vz(:,:,l)=(vbuf(:,:,kb+1)-vbuf(:,:,kb))*dzm
  call xderiv(wbuf(1,1,kb),wx(1,1,l),xkb,nnxb,nnyb)
  call yderiv(wbuf(1,1,kb),wy(1,1,l),ykb,nnxb,nnyb)
  
  if (iz.eq.1) then
     ux(:,:,nm1) = 0.0
     uy(:,:,nm1) = 0.0
     vx(:,:,nm1) = 0.0
     vy(:,:,nm1) = 0.0
     wz(:,:,nm1) = 0.0
     wx(:,:,lm1) = 0.0
     wy(:,:,lm1) = 0.0
     uz(:,:,lm1) = uz(:,:,l)
     vz(:,:,lm1) = vz(:,:,l)
  else
     call xderiv(ubuf(1,1,kb-1),ux(1,1,nm1),xkb,nnxb,nnyb)
     call yderiv(ubuf(1,1,kb-1),uy(1,1,nm1),ykb,nnxb,nnyb)
     call xderiv(vbuf(1,1,kb-1),vx(1,1,nm1),xkb,nnxb,nnyb)
     call yderiv(vbuf(1,1,kb-1),vy(1,1,nm1),ykb,nnxb,nnyb)
     wz(:,:,nm1) = (wbuf(:,:,kb-1)-wbuf(:,:,kb-2))*dzm
     uz(:,:,lm1) = (ubuf(:,:,kb)-ubuf(:,:,kb-1))*dzm
     vz(:,:,lm1) = (vbuf(:,:,kb)-vbuf(:,:,kb-1))*dzm
     call xderiv(wbuf(1,1,kb-1),wx(1,1,lm1),xkb,nnxb,nnyb)
     call yderiv(wbuf(1,1,kb-1),wy(1,1,lm1),ykb,nnxb,nnyb)
  endif


!   if(iz .lt. 4) then
!      write(*,*) 'iz      ux        uy       uz        vx      vy      vz     wx    wy     wz'
!      write(*,*) iz, sum(ux(:,:,n))*xymb, sum(uy(:,:,n))*xymb, &
!           sum(uz(:,:,l))*xymb, sum(vx(:,:,n))*xymb, sum(vy(:,:,n))*xymb, &
!           ' ', sum(vz(:,:,l))*xymb, sum(wx(:,:,l))*xymb, &
!           sum(wy(:,:,l))*xymb, sum(wz(:,:,n))*xymb

!      write(*,*) 'iz      ux        uy       uz        vx      vy      vz     wx    wy     wz'
!      write(*,*) iz-1, sum(ux(:,:,nm1))*xymb, sum(uy(:,:,nm1))*xymb, &
!           sum(uz(:,:,lm1))*xymb, sum(vx(:,:,nm1))*xymb, sum(vy(:,:,nm1))*xymb, &
!           ' ', sum(vz(:,:,lm1))*xymb, sum(wx(:,:,lm1))*xymb, &
!           sum(wy(:,:,lm1))*xymb, sum(wz(:,:,nm1))*xymb
!   end if
  
end subroutine derivsbot




  
