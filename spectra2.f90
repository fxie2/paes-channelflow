subroutine spectra(fnt1, zspec2d, zspecx)
 
  use storage
  implicit none
 
  real    :: zspec2d(nnx/2+1), zspecx(nnx/2+1) 
  real    :: fnt1(nnx,nny)
  complex :: wave1(nnx/2+1,1:nny)
  real    :: temp1(nnx+2,1:nny)
  real    :: var1, var2, var3, dist, denom
  integer :: m(nnx/2+1), iin, i, j
   
  temp1(1:nnx,:) = fnt1
  zspec2d = 0.0
  zspecx = 0.0
  m = 0
  
  call ffthrz(temp1, nnx, nny, -1)

  do j=1,nny
     wave1(:,j) = transfer(temp1(:,j),wave1(:,j))
  enddo
  
  do i=1,ncx
     do j=1,nny
        zspecx(i) = zspecx(i) + wave1(i,j)*conjg(wave1(i,j))
     enddo
  enddo
  
  do i=2,ncx-1
     zspecx(i) = zspecx(i)*2.
  enddo
  
  do j=1,nny
     do i=1,ncx
        dist=sqrt(xk(i)**2+yk(j)**2)
        if((i.eq.nnx/2+1).and.(j.eq.nny/2+1))then
           dist=sqrt((pi*nnx/xl)**2+(pi*nny/yl)**2)
        else
           if(i.eq.nnx/2+1)dist=sqrt((pi*nnx/xl)**2+yk(j)**2)
           if(j.eq.nny/2+1)dist=sqrt(xk(i)**2+(pi*nny/yl)**2)
        endif
        iin=nint(dist/(2.0*pi/xl))+1
        if((iin.gt.nnx/2+1).or.((i.eq.1).and.(j.gt.nny/2+1)))then
           continue
        else
           zspec2d(iin)=zspec2d(iin)+wave1(i,j)*conjg(wave1(i,j))
           m(iin)=m(iin)+1
        endif
     enddo
  enddo
  
  do i=1,ncx
     j=i-1
     if(i.eq.1)then
        denom =1.0 
     else 
        denom = 2.0*float(m(i))/(pi*(j+0.5)**2-pi*(j-0.5)**2)
     endif
     zspec2d(i)=2.0*zspec2d(i)/denom
  enddo
  
!  var1 = sum(fnt1*fnt1)*xym
!  var2 = sum(zspec2d)
!  var3 = sum(zspecx)
  
!  write(*,*) "physical space variance ", var1, "variance from 2d spec ", &
!       var2,  "variance from 1d spec ", var3
  
end subroutine spectra




subroutine cospec(fnt1, fnt2, cospec2d, quadsp2d, cospecx, quadspx)
 
  use storage
  implicit none
 
  real    :: cospec2d(nnx/2+1), cospecx(nnx/2+1) 
  real    :: quadsp2d(nnx/2+1), quadspx(nnx/2+1) 
  real    :: fnt1(nnx,nny), fnt2(nnx,nny)
  complex :: wave1(nnx/2+1,1:nny), wave2(nnx/2+1,1:nny)
  real    :: temp1(nnx+2,1:nny), temp2(nnx+2,1:nny)
  real    :: var1, var2, var3, dist, denom
  integer :: m(nnx/2+1), iin, i, j
   
  temp1(1:nnx,:) = fnt1
  temp2(1:nnx,:) = fnt2
  cospec2d = 0.0
  cospecx = 0.0
  quadsp2d = 0.0
  quadspx = 0.0
  m = 0
 
  call ffthrz(temp1, nnx, nny, -1)
  call ffthrz(temp2, nnx, nny, -1)
  
  do j=1,nny
     wave1(:,j) = transfer(temp1(:,j),wave1(:,j))
     wave2(:,j) = transfer(temp2(:,j),wave2(:,j))
  enddo  

  do i=1,ncx
     do j=1,nny
        cospecx(i) = cospecx(i) + real(wave1(i,j)*conjg(wave2(i,j)))
        quadspx(i) = quadspx(i) + aimag(wave1(i,j)*conjg(wave2(i,j)))
     enddo
  enddo
  
  do i=2,ncx-1
     cospecx(i) = cospecx(i)*2.
     quadspx(i) = quadspx(i)*2.
  enddo
  
  do j=1,nny
     do i=1,ncx
        dist=sqrt(xk(i)**2+yk(j)**2)
        if((i.eq.nnx/2+1).and.(j.eq.nny/2+1))then
           dist=sqrt((pi*nnx/xl)**2+(pi*nny/yl)**2)
        else
           if(i.eq.nnx/2+1)dist=sqrt((pi*nnx/xl)**2+yk(j)**2)
           if(j.eq.nny/2+1)dist=sqrt(xk(i)**2+(pi*nny/yl)**2)
        endif
        iin=nint(dist/(2.0*pi/xl))+1
        if((iin.gt.nnx/2+1).or.((i.eq.1).and.(j.gt.nny/2+1)))then
           continue
        else
           cospec2d(iin) = cospec2d(iin) + real(wave1(i,j)*conjg(wave2(i,j)))
           quadsp2d(iin) = quadsp2d(iin) + aimag(wave1(i,j)*conjg(wave2(i,j)))
           m(iin)=m(iin)+1
        endif
     enddo
  enddo
  
  do i=1,ncx
     j=i-1
     if(i.eq.1)then
        denom =1.0 
     else 
        denom = 2.0*float(m(i))/(pi*(j+0.5)**2-pi*(j-0.5)**2)
     endif
     cospec2d(i) = 2.0*cospec2d(i)/denom
     quadsp2d(i) = 2.0*quadsp2d(i)/denom
  enddo
  
!  var1 = sum(fnt1*fnt2)*xym
!  var2 = sum(cospec2d)
!  var3 = sum(cospecx)
!  
!  write(*,*) "physical space variance ", var1, "variance from 2d spec ", &
!       var2,  "variance from 1d spec ", var3
  
end subroutine cospec       

