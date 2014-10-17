!12/13/2006: I changed the nonlinear terms to be divergence form
!
subroutine compadv1(kb,iz,lm1,l,n,np1,istart)

  ! COMPUTES THE ADVECTION TERMS FOR MOMENTUM in divergence form

  use storage
  use compbuf
  use parameters
  implicit none

  ! INPUT ARGUMENTS:
  integer :: iz, lm1, l, istart, kb, n, np1

  ! LOCAL VARIABLES:
  real :: fnt1(nnxb,nnyb), fnt2(nnxb,nnyb)
  
  ! ADVECTION OF MOMENTUM  
  CALL xderiv( ubuf(:,:,kb)*ubuf(:,:,kb), fnt1, xkb, nnxb, nnyb)
  CALL yderiv( ubuf(:,:,kb)*vbuf(:,:,kb), fnt2, ykb, nnxb, nnyb)
  
  rubuf(:,:) = -fnt1(:,:) - fnt2(:,:) - 		     &
      ( ( ubuf(:,:,kb+1)+ubuf(:,:,kb) )*wbuf(:,:,kb  ) -     &
        ( ubuf(:,:,kb-1)+ubuf(:,:,kb) )*wbuf(:,:,kb-1) )*todzm
  
!   if((iz.le.2).or.(iz.gt.62)) then
!      write(*,*) 'rubuf (',iz,') = ', sum(fnt1(:,:))*xymb, sum(fnt2(:,:))*xymb, sum(( ( ubuf(:,:,kb+1)+ubuf(:,:,kb) )*wbuf(:,:,kb  ) - ( ubuf(:,:,kb-1)+ubuf(:,:,kb) )*wbuf(:,:,kb-1) )*todzm )*xymb,  sum(rubuf(:,:))*xymb
!   end if
  
  
  if (nny .eq. 1) then
     rvbuf(:,:) = 0.0
  else
    CALL xderiv( ubuf(:,:,kb)*vbuf(:,:,kb), fnt1, xkb, nnxb, nnyb)
    CALL yderiv( vbuf(:,:,kb)*vbuf(:,:,kb), fnt2, ykb, nnxb, nnyb)      
     rvbuf(:,:) =-fnt1(:,:) - fnt2(:,:) - 		      &
      ( ( vbuf(:,:,kb+1)+vbuf(:,:,kb) )*wbuf(:,:,kb) -      &
        ( vbuf(:,:,kb-1)+vbuf(:,:,kb) )*wbuf(:,:,kb-1) )*todzm 
     
  endif

  if (iz .eq. nnz) then
     rwbuf(:,:) = 0.
  else
    CALL xderiv((ubuf(:,:,kb)+ubuf(:,:,kb+1))*wbuf(:,:,kb), &
    		fnt1,xkb,nnxb,nnyb)
    CALL yderiv((vbuf(:,:,kb)+vbuf(:,:,kb+1))*wbuf(:,:,kb), &
    		fnt2,ykb,nnxb,nnyb)     
     rwbuf(:,:) = -0.5*fnt1(:,:) - 0.5*fnt2(:,:)  - 		&
	  ( (0.5*(wbuf(:,:,kb+1)+wbuf(:,:,kb  )))**2 - 	        &
            (0.5*(wbuf(:,:,kb  )+wbuf(:,:,kb-1)))**2 )*dzm
  endif

!      if((iz.le.4).or.(iz.gt.60)) then
!         write(*,*) 'rw (',iz,') =  ', sum(rwbuf(:,:))*xymb
!      end if

  ! ADVECTION OF SGS TKE - FLUX FORM
  if (isbgr .and. iz.ne.nnz) then
     call xderiv(ebuf(:,:,kb)*(ubuf(:,:,kb)+ubuf(:,:,kb+1)),fnt1,xkb,nnxb,nnyb)
     call yderiv(ebuf(:,:,kb)*(vbuf(:,:,kb)+vbuf(:,:,kb+1)),fnt2,ykb,nnxb,nnyb)
     rebuf(:,:) = -0.5*(fnt1(:,:)+fnt2(:,:)) - &
          (wbuf(:,:,kb+1)*ebuf(:,:,kb+1)-wbuf(:,:,kb-1)*ebuf(:,:,kb-1))*todzm
  endif
  
end subroutine compadv1


