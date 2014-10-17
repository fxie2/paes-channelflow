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

  if (abs(facw) .lt. 1.e4) then
     tauuw(:,:,k) = facw*(windm*(ubuf(:,:,0)+ugal-u1xy) + wind(:,:)*u1xy)
     tauvw(:,:,k) = facw*(windm*(vbuf(:,:,0)+vgal-v1xy) + wind(:,:)*v1xy)
  else
     tauuw(:,:,k) = auuwm
     tauvw(:,:,k) = auvwm
  endif

end subroutine sufsim_upper







