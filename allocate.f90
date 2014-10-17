subroutine allocate_arrays(istart,iend,np,nm)

  use storage
  use compbuf
  use parameters
  implicit none
  integer :: istart,iend,imin,imax,ivar
  integer :: np(4), nptot(4)
  integer :: nm(4), nmtot(4)

  imin = istart-1
  imax = iend+1
  
  allocate(U(nnx,nny,istart-nm(1):iend+np(1)))
  allocate(V(nnx,nny,istart-nm(2):iend+np(2)))
  allocate(W(nnx,nny,istart-nm(3):iend+np(3)))
     
  allocate(RU(nnx,nny,istart:iend))
  allocate(RV(nnx,nny,istart:iend))
  allocate(RW(nnx,nny,imin:iend))

  allocate(ubuf(nnxb,nnyb,-nm(1):np(1)))
  allocate(vbuf(nnxb,nnyb,-nm(2):np(2)))
  allocate(wbuf(nnxb,nnyb,-nm(3):np(3)))
  allocate(rubuf(nnxb,nnyb))
  allocate(rvbuf(nnxb,nnyb))
  allocate(rwbuf(nnxb,nnyb))

  allocate(P(nnx+2,nny,istart:imax))

  if (isbgr) then
     allocate(E(nnx,nny,istart-nm(4):iend+np(4)))
     allocate(RE(nnx,nny,istart:iend))
     allocate(ebuf(nnxb,nnyb,-nm(4):np(4)))
     allocate(rebuf(nnxb,nnyb))
  endif
 
  allocate(ug(istart:iend))
  allocate(vg(istart:iend))
  
  allocate(UX(nnxb,nnyb,4))
  allocate(UY(nnxb,nnyb,4))
  allocate(UZ(nnxb,nnyb,3))
  allocate(VX(nnxb,nnyb,4))
  allocate(VY(nnxb,nnyb,4))
  allocate(VZ(nnxb,nnyb,3))
  allocate(WX(nnxb,nnyb,3))
  allocate(WY(nnxb,nnyb,3))

  if (isbgr) then
     allocate(EX(nnxb,nnyb))
     allocate(EY(nnxb,nnyb))
     allocate(wz(nnxb,nnyb,4))
     allocate(tauuu(nnxb,nnyb))
     allocate(tauuv(nnxb,nnyb))
     allocate(tauvv(nnxb,nnyb))
     allocate(tauww(nnxb,nnyb,2))
     allocate(tauuw(nnxb,nnyb,2))
     allocate(tauvw(nnxb,nnyb,2))
     allocate(xkk(nnxb,nnyb,3))
     allocate(alk(NNXB,NNYB,3))
     allocate(shear(nnxb,nnyb,3))
  endif
  if (isbgr) then
     allocate(wxjp2(nnxb,nnyb))
     allocate(wyjp2(nnxb,nnyb))
     allocate(uzjp2(nnxb,nnyb))
     allocate(vzjp2(nnxb,nnyb))
  endif

  allocate(wavexy(nnx+2,nny))
  allocate(xk2(nnx+2,nny))
  allocate(xk(NNX))
  allocate(yk(NNY))
  allocate(xkb(NNXB))
  allocate(ykb(NNYB))

  ivar = 1

  mcount = 19 !Heavily modified from Martin's code
  if (istart.eq.1) then
     allocate(meanbuf(mcount,1:nnz))
     allocate(zt(nnz))
     allocate(zm(nnz))
  else
     allocate(meanbuf(mcount,istart:iend))
  endif
  meanbuf = 0.0
  englez => meanbuf(1,:)
  engsbz => meanbuf(2,:)
  divz => meanbuf(3,:)

  uxym => meanbuf(4,:)
  vxym => meanbuf(5,:)
  wxym => meanbuf(6,:)

  uule => meanbuf(7,:)
  vvle => meanbuf(8,:)
  wwle => meanbuf(9,:)
  uwle => meanbuf(10,:)
  vwle => meanbuf(11,:)
  uvle => meanbuf(12,:)
  
  pple => meanbuf(13,:)

  uusb => meanbuf(14,:)
  vvsb => meanbuf(15,:)
  uvsb => meanbuf(16,:)
  wwsb => meanbuf(17,:)
  uwsb => meanbuf(18,:)
  vwsb => meanbuf(19,:)


  spectcount = 3
  if (spectav) then
     allocate(spectbuf(ncx,spectcount,15)) !15 z levels
  else
     allocate(spectbuf(ncx,spectcount,15)) !15 z levels
  endif
  spectbuf = 0.0
  
  uspec2d => spectbuf(:,1,:)
  vspec2d => spectbuf(:,2,:)
  wspec2d => spectbuf(:,3,:)
  

end subroutine allocate_arrays




subroutine deallocate_arrays()

  use storage
  use compbuf
  use parameters
  implicit none

  deallocate(U)
  deallocate(V)
  deallocate(W)
     
  deallocate(RU)
  deallocate(RV)
  deallocate(RW)
  deallocate(P)

  deallocate(ubuf)
  deallocate(vbuf)
  deallocate(wbuf)
  deallocate(rubuf)
  deallocate(rvbuf)
  deallocate(rwbuf)

  if (isbgr) then
     deallocate(E)
     deallocate(RE)
     deallocate(ebuf)
     deallocate(rebuf)
  endif
 
  deallocate(ug)
  deallocate(vg)
  deallocate(meanbuf)
      
  deallocate(UX)
  deallocate(UY)
  deallocate(UZ)
  deallocate(VX)
  deallocate(VY)
  deallocate(VZ)
  deallocate(WX)
  deallocate(WY)
  deallocate(EX)
  deallocate(EY)
  if (isbgr) deallocate(wz)
      
  deallocate(tauuw)
  deallocate(tauvw)

  deallocate(xkk)
  deallocate(alk)

  deallocate(wavexy)
  deallocate(xk2)
  deallocate(xk)
  deallocate(yk)

  
end subroutine deallocate_arrays








