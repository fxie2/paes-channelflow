subroutine compsgs3(kb,iz,j,jp1,k,km1,l,lm1,lp1,n,np1,mnout,istage,istart)
  
  ! ADDS THE SUBGRID-SCALE STRESS CONTRIBUTION TO THE TENDENCIES

  use storage
  use compbuf
  use parameters
  implicit none

  ! INPUT PARAMETERS:
  integer iz,j,jp1,k,km1,l,lp1,lm1,istart,kb,n,np1,istage
  logical mnout

  ! LOCAL VARIABLES:
  real :: ump, um, vmp, vm
  real :: fnt1(nnxb,nnyb),fnt2(nnxb,nnyb),temp(nnxb,nnyb)
  integer i
  
  i = iz-istart+1
  if (iz.ne.nnz) then

     ! DETERMINE THE DIFFUSION SOURCE CONTRIBUTION TO SGS TKE.
     fnt1 = (2.*xkk(:,:,l))*ex
     fnt2 = (2.*xkk(:,:,l))*ey
     call xderiv((fnt1),fnt1,xkb,nnxb,nnyb)
     call yderiv((fnt2),fnt2,ykb,nnxb,nnyb)
     temp = fnt1 + fnt2 + &
          ((xkk(:,:,lp1)+xkk(:,:,l))*(ebuf(:,:,kb+1)-ebuf(:,:,kb))- &
          (xkk(:,:,l)+xkk(:,:,lm1))*(ebuf(:,:,kb)-ebuf(:,:,kb-1)))*dzm2
     rebuf = rebuf + temp
     
     ! DISSIPATION TERM IN SUB-GRID KE EQUATION.
     temp = -ceps*ebuf(:,:,kb)*sqrt(ebuf(:,:,kb))
     rebuf = rebuf + temp
  endif

!  write(*,*)'Finished evaluating sgs terms in e equation level ', iz

  ! SUBGRID STRESS DIVERGENCE IN U-EQUATION
  call xderiv(tauuu,fnt1,xkb,nnxb,nnyb)
  call yderiv(tauuv,fnt2,ykb,nnxb,nnyb)
  rubuf = rubuf - fnt1 - fnt2 - (tauuw(:,:,k)-tauuw(:,:,km1))*dzm
!  if((iz.le.3).or.(iz.gt.60)) then
!     write(*,*) 'rubuf (',iz,') = ', sum(fnt1(:,:))*xymb, sum(fnt2(:,:))*xymb, sum((tauuw(:,:,k)-tauuw(:,:,km1))*dzm)*xymb,  sum(rubuf(:,:))*xymb
!  end if
  
!  write(*,*)'Finished evaluating sgs terms in u equation level ', iz
  ! SUBGRID STRESS DIVERGENCE IN V-EQUATION
  if (nny .ne. 1) then
     call xderiv(tauuv,fnt1,xkb,nnxb,nnyb)
     call yderiv(tauvv,fnt2,ykb,nnxb,nnyb)
     rvbuf = rvbuf - fnt1 - fnt2 - (tauvw(:,:,k)-tauvw(:,:,km1))*dzm
  endif
!  write(*,*)'Finished evaluating sgs terms in v equation level ', iz
   ! SUBGRID STRESS DIVERGENCE IN W-EQUATION
  if (iz.ne.nnz) then
     call xderiv(tauuw(:,:,k),fnt1,xkb,nnxb,nnyb)
     call yderiv(tauvw(:,:,k),fnt2,ykb,nnxb,nnyb)
     rwbuf = rwbuf - fnt1 - fnt2 - (tauww(:,:,jp1)-tauww(:,:,j))*dzm
  else
     fnt1 = 0.0
     fnt2 = 0.0
     rwbuf = 0.0
  endif
!  if((iz.le.3).or.(iz.gt.60)) then
!     write(*,*) 'rwbuf (',iz,') = ', sum(fnt1(:,:))*xymb, sum(fnt2(:,:))*xymb, sum((tauww(:,:,jp1)-tauww(:,:,j))*dzm)*xymb,  sum(rwbuf(:,:))*xymb
!  end if

  
!  write(*,*)'Finished evaluating sgs terms in w equation level ', iz
  if (mnout .and. istage.eq.1) then
     uvsb(i) = sum(dble(tauuv))*xymb
     uusb(i) = sum(dble(tauuu))*xymb
     vvsb(i) = sum(dble(tauvv))*xymb
     wwsb(i) = sum(dble(tauww(:,:,j)))*xymb
     uwsb(i) = sum(dble(tauuw(:,:,k)))*xymb
     vwsb(i) = sum(dble(tauvw(:,:,k)))*xymb
  endif
!  write(*,*)'Finished evaluating means of sgs terms level ', iz
end subroutine compsgs3
