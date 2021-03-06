subroutine upper(np)
  
  ! DETERMINES THE UPPER RADIATION BOUNDARY CONDITION.
  ! TO BE CALLED WHEN IZ=NNZ-1.

  use storage
  use compbuf
  use parameters
  implicit none
  integer :: np(7), iz
     
  do iz=nnz+1,nnz+np(1)
     u(:,:,iz) = 0.0
  enddo
  
  do iz=nnz+1,nnz+np(2)
     v(:,:,iz) = 0.0
  enddo

  if (isbgr) then
     do iz=nnz,nnz+np(4)
        e(:,:,iz) = e(:,:,nnz-1)
     enddo
  endif
  
  w(:,:,nnz:nnz+np(3)) = 0.0
  
  rw(:,:,nnz) = 0.0

end subroutine upper
  

subroutine sufsim_upper(k,km1)

  ! DETERMINES THE SURFACE STRESS BOUNDARY CONDITIONS FOR UPPER WALL USING 
  ! SIMILARITY SCALING.

  use storage
  use compbuf
  use parameters
  implicit none
  
  integer :: ix,iy,iter,k,km1,ierror
  logical :: iloop
  real :: wind(nnxb,nnyb),windm, tep, thta
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
  u1xy=sum(dble(ubuf(:,:,0)))*xymb+ugal
  v1xy=sum(dble(vbuf(:,:,0)))*xymb+vgal

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
  auuwm=utau2*cos(thta)
  auvwm=utau2*sin(thta)*sign(1.,v1xy)
! facw = -utau2/max(windm*sqrt(u1xy*u1xy+v1xy*v1xy),1.e-12)
  facw = utau2/(windm*vsfc)

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
  tauuw(:,:,k) = utau2/(vsfc)*(ubuf(:,:,0)+ugal)
  tauvw(:,:,k) = utau2/(vsfc)*(vbuf(:,:,0)+vgal)

  CASE(2)	!Moeng's model
  tauuw(:,:,k) = utau2/(windm*vsfc)*(windm*ufluc + wind(:,:)*u1xy)
  tauvw(:,:,k) = utau2/(windm*vsfc)*(windm*vfluc + wind(:,:)*v1xy)

  CASE(3)	!Wei's model
  tauuw(:,:,k) = auuwm-betaa*ufluc*utau2/ufluc2mroot
  tauvw(:,:,k) = auvwm-betaa*vfluc*utau2/vfluc2mroot

  CASE(4)	!Moeng's model without linearization
  
  denom=sqrt((windm*u1xy+suflucm)**2+(windm*v1xy+svflucm)**2)

  tauuw(:,:,k) = utau2*(wind(:,:)*(ubuf(:,:,0)+ugal))/denom
  tauvw(:,:,k) = utau2*(wind(:,:)*(vbuf(:,:,0)+vgal))/denom

  CASE(5)	!Wei's model with Moeng's Fluctuation

  tauuw(:,:,k) = auuwm-betaa*utau2*moengu/moengu2mroot
  tauvw(:,:,k) = auvwm-betaa*utau2*moengv/moengv2mroot

  CASE(6)	! NEW APPROACH + Schumann's model

  auuwm=-zl/2*(0.016-dudt)
  auvwm=-zl/2*(-dvdt)

  tauuw(:,:,k) = -auuwm+utau*utau_s*ufluc/vsfc
  tauvw(:,:,k) = -auvwm+utau*utau_s*vfluc/vsfc

END SELECT

end subroutine sufsim_upper







