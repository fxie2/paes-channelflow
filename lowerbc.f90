subroutine lower(istart,np,nm)
  use parameters
  use storage
  use compbuf
  implicit none
  integer :: np(7), nm(7), iz, i, istart
  real :: fnt1(nnx,nny)
  
  if (istart.eq.1) then
     u(:,:,1-nm(1):0) = 0.0
     v(:,:,1-nm(2):0) = 0.0
     w(:,:,1-nm(3):0) = 0.0
     if (isbgr) then
        do iz=1-nm(4),0
           e(:,:,iz) = e(:,:,1)
        enddo
     endif
  endif

end subroutine lower




subroutine sufsim(k,km1)

  ! DETERMINES THE SURFACE STRESS BOUNDARY CONDITIONS USING 
  ! SIMILARITY SCALING.

  use storage
  use compbuf
  use parameters
  implicit none
  
  integer :: ix,iy,iter,k,km1,ierror
  logical :: iloop
  real :: wind(nnxb,nnyb),windm,tep,thta,thstar
  real :: ufree,vsfc,utau2,zeta,dzeta,zeta_a,z_old
  real :: x,y,u1xy,v1xy,t1xy,t10xy,dnom,zeta_min
  real :: qhstar,chstar,facq,facc,fact,facw
  real :: q1xy,q10xy,c1xy,c10xy,auuwm,auvwm
  real :: psim, psih, bus_psim, bus_psih

     
  ! FIND MEAN FIELD VALUES AT THIS LEVEL.
  u1xy=sum(dble(ubuf(:,:,0)))*xymb
  v1xy=sum(dble(vbuf(:,:,0)))*xymb

  wind = sqrt(ubuf(:,:,0)**2+vbuf(:,:,0)**2)
  windm = sum(dble(wind))*xymb
  vsfc = sqrt(u1xy*u1xy+v1xy*v1xy)
  windm = max(windm, ufree)
  windm = max(windm, 0.01)
  vsfc = max(vsfc, 0.01)

  utau = vsfc*vk/(zody)
 
  
  if (utau.gt.10.) then
     write(*,'(" stop because utau=",e12.3," windm=",e10.3)') utau,windm
     call lesabort(" ")
  endif

  ! NOTE ROUNDOFF PROBLEM IF ANGLES ARE CLOSE TO MULTIPLES OF PI.
  tep = u1xy/vsfc
  tep = min(tep,1.)
  tep = max(tep,-1.)
  thta=acos(tep)
  utau2=utau*utau
  auuwm=-utau2*cos(thta)
  auvwm=-utau2*sin(thta)*sign(1.,v1xy)
! facw = -utau2/max(windm*sqrt(u1xy*u1xy+v1xy*v1xy),1.e-12)
  facw = -utau2/(windm*vsfc)

  
  if (abs(facw) .lt. 1.e4) then
     tauuw(:,:,km1) = facw*(windm*(ubuf(:,:,0)+ugal-u1xy) + wind(:,:)*u1xy)
     tauvw(:,:,km1) = facw*(windm*(vbuf(:,:,0)+vgal-v1xy) + wind(:,:)*v1xy)
  else
     tauuw(:,:,km1) = auuwm
     tauvw(:,:,km1) = auvwm
  endif

end subroutine sufsim












