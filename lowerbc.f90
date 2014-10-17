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

  ! new variables added by Paulo 10/17/14
  real :: sfluc(nnxb,nnyb),sfluc2(nnxb,nnyb),sfluc2m
  real :: ufluc2(nnxb,nnyb),vfluc2(nnxb,nnyb),ufluc2m
  real :: vfluc2m,ufluc2mroot,vfluc2mroot,sfluc2mroot
  real :: ufluc(nnxb,nnyb),vfluc(nnxb,nnyb)
  real :: sufluc(nnxb,nnyb),svfluc(nnxb,nnyb),suflucm,svflucm
  real :: denom,checku,checkv,checkfluc
  real :: moengu(nnxb,nnyb),moengv(nnxb,nnyb),moengu2(nnxb,nnyb),moengv2(nnxb,nnyb),moengu2m,moengv2m,moengu2mroot,moengv2mroot
  integer :: sx,sy


     
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

  !new modelling approach by Paulo 10/17/14
  dudt = sum(dble(ru))/(nnx*nny*nnz)
  dvdt = sum(dble(rv))/(nnx*nny*nnz)
  utau_s = sqrt(zl/2*sqrt((-dvdt)**2+(0.016-dudt)**2))
  ! global quantities only on lower because is called first than upper 

 
  
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


! additional models by Paulo 10/17/14
if (isurface>2) then
  sfluc=wind(:,:)-windm
  ufluc=ubuf(:,:,0)+ugal-u1xy
  vfluc=vbuf(:,:,0)+vgal-v1xy
  do sx = 1,nnxb
    do sy = 1,nnyb
      sfluc2(sx,sy)=sfluc(sx,sy)**2
      ufluc2(sx,sy)=ufluc(sx,sy)**2
      vfluc2(sx,sy)=vfluc(sx,sy)**2
      sufluc(sx,sy)=ufluc(sx,sy)*sfluc(sx,sy)
      svfluc(sx,sy)=vfluc(sx,sy)*sfluc(sx,sy)
      moengu(sx,sy)=(sfluc(sx,sy)*u1xy+windm*ufluc(sx,sy))
      moengv(sx,sy)=(sfluc(sx,sy)*v1xy+windm*vfluc(sx,sy))
      moengu2(sx,sy)=moengu(sx,sy)**2
      moengv2(sx,sy)=moengv(sx,sy)**2
    end do
  end do
  sfluc2m=sum(dble(sfluc2))*xymb
  ufluc2m=sum(dble(ufluc2))*xymb
  vfluc2m=sum(dble(vfluc2))*xymb
  sfluc2mroot=sqrt(sfluc2m)
  ufluc2mroot=sqrt(ufluc2m)
  vfluc2mroot=sqrt(vfluc2m)
  suflucm=sum(dble(sufluc))*xymb
  svflucm=sum(dble(svfluc))*xymb
  moengu2m=sum(dble(moengu2))*xymb
  moengv2m=sum(dble(moengu2))*xymb
  moengu2mroot=sqrt(moengu2m)
  moengv2mroot=sqrt(moengv2m)
end if

SELECT CASE(isurface)

  CASE(1)	!Schumann's model
  tauuw(:,:,km1) = -utau2/(vsfc)*(ubuf(:,:,0)+ugal)
  tauvw(:,:,km1) = -utau2/(vsfc)*(vbuf(:,:,0)+vgal)

  CASE(2)	!Moeng's model
  tauuw(:,:,km1) = -utau2/(windm*vsfc)*(windm*ufluc + wind(:,:)*u1xy)
  tauvw(:,:,km1) = -utau2/(windm*vsfc)*(windm*vfluc + wind(:,:)*v1xy)

  CASE(3)	!Wei's model with Schumman's Fluctuation
  tauuw(:,:,km1) = auuwm+betaa*ufluc*utau2/ufluc2mroot
  tauvw(:,:,km1) = auvwm+betaa*vfluc*utau2/vfluc2mroot

  CASE(4)	!Moeng's model without linearization
  
  denom=sqrt((windm*u1xy+suflucm)**2+(windm*v1xy+svflucm)**2)

  tauuw(:,:,km1) = -utau2*(wind(:,:)*(ubuf(:,:,0)+ugal))/denom
  tauvw(:,:,km1) = -utau2*(wind(:,:)*(vbuf(:,:,0)+vgal))/denom

  CASE(5)	!Wei's model with Moeng's Fluctuation

  tauuw(:,:,km1) = auuwm+betaa*utau2*moengu/moengu2mroot
  tauvw(:,:,km1) = auvwm+betaa*utau2*moengv/moengv2mroot

  CASE(6)	! NEW APPROACH + Schumann's model

  auuwm=-zl/2*(0.016-dudt)
  auvwm=-zl/2*(-dvdt)

  tauuw(:,:,km1) = auuwm-utau*utau_s*ufluc/vsfc
  tauvw(:,:,km1) = auvwm-utau*utau_s*vfluc/vsfc


END SELECT



end subroutine sufsim












