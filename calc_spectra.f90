program read_write_3d
!This program reads the files that contain the 3d data (velocity, pressure, scalars) in the ".field" files and calculates the velocity spectra
implicit none
include "/usr/global/fftw/3.2.2/pgi/include/fftw3.f"


character *100 :: buffer !A character buffer to read input argument
integer :: i,j,k !Counter variables
integer :: iz
! GRID PARAMETERS
integer :: nnx, nny, nnz      ! x, y, and z grid dimensions 
integer :: nxy                ! number of grid points in a horizontal plane
real :: xl, yl, zl         ! x, y, and z domain sizes
real :: dx, dy, dz         ! x, y, and z grid lengths
real :: ugtop, ugbot, vgtop, vgbot
real :: divgls, fcor, utau
real :: time_start, dt, z0
real :: ugal, vgal
integer :: fcounter !Number of files used in the averaging process
real, allocatable, dimension(:,:,:) :: u,v,w,e,p
real, allocatable, dimension(:) :: umean, vmean
real, allocatable, dimension(:) :: uvar, vvar, wvar
real, allocatable, dimension(:) :: l111, l112, l221, l222, l331, l332

!FFT variables
double precision, allocatable, dimension(:,:) :: fft_in !Temporary array to store the input to the FFTW function
double complex, allocatable, dimension(:,:) :: fft_out !Temporary array to store the output from the FFTW function
real, allocatable, dimension(:,:) :: uhat, vhat, what !Array to store the spectrum at pre-defined z levels
real, allocatable, dimension(:) :: waveN ! Wavenumber
real, allocatable, dimension(:) :: nCount !Number of points in each wave number ring.
integer :: nLevels = 16 !Number of levels at which the spectrum is to be measured.
integer, dimension(16) :: zLevels = (/1,2,3,4,5,6,7,8,9,10,12,18,24,30,36,48/) !The levels at which the spectrum is to be measured
real :: deltaWaveN !The difference between two successive wave numbers
real :: pi = 3.14159265358979
integer :: r ! Distance of 2d wave number from orgin
integer*8 plan ; !Plan to create the FFT
integer :: stat



fcounter = 0
open(unit=75,file="fieldfiles")

read(75,*,iostat=stat)buffer

!call getarg(1,buffer)
open(unit=50,file=buffer,form="unformatted")
read(50) time_start, nnx, nny, nnz, xl, yl, zl
read(50) dt, z0, utau
read(50) divgls, fcor, ugtop, ugbot, vgtop, vgbot

allocate(u(nnx,nny,nnz))
allocate(v(nnx,nny,nnz))
allocate(w(nnx,nny,nnz))
allocate(e(nnx,nny,nnz))
allocate(p(nnx,nny,nnz))

allocate(umean(nnz))
allocate(vmean(nnz))

close(50)

write(*,*) 'nnx = ', nnx
write(*,*) 'nny = ', nny
write(*,*) 'nnz = ', nnz

dx = xl/nnx
dy = yl/nny
dz = zl/nnz

write(*,*) 'dx = ', dx
write(*,*) 'dy = ', dy
write(*,*) 'dz = ', dz

!Initialise FFT variables
!!Allocate complex arrays
allocate(waveN(nnx))
allocate(nCount(0:nnx/2))
allocate(fft_in(nnx,nny))
allocate(fft_out(nnx/2+1,nny))
allocate(uhat(0:nnx/2, nLevels))
allocate(vhat(0:nnx/2, nLevels))
allocate(what(0:nnx/2, nLevels))


allocate(uvar(nLevels))
allocate(vvar(nLevels))
allocate(wvar(nLevels))
allocate(l111(nLevels))
allocate(l112(nLevels))
allocate(l221(nLevels))
allocate(l222(nLevels))
allocate(l331(nLevels))
allocate(l332(nLevels))
l111 = 0
l112 = 0
l221 = 0
l222 = 0
l331 = 0
l332 = 0
uvar = 0
vvar = 0
wvar = 0

!!Wave number
deltaWaveN = 1/xl
do i =1,nnx/2+1
   waveN(i) = (i-1)
end do
do i=nnx/2+2,nnx
   waveN(i) = i-nnx-1
end do
nCount = 0
uhat = 0
vhat = 0
what = 0

!!Create FFTW plan
call dfftw_plan_dft_r2c_2d(plan, nnx, nny, fft_in, fft_out, FFTW_ESTIMATE)


do while (stat .eq. 0) 
fcounter = fcounter + 1
open(unit=50,file=buffer,form="unformatted")
read(50) time_start, nnx, nny, nnz, xl, yl, zl
read(50) dt, z0, utau
read(50) divgls, fcor, ugtop, ugbot, vgtop, vgbot

nxy = nnx*nny
ugal = (max(0.,max(ugtop,ugbot))+min(0.,min(ugtop,ugbot)))*0.5
vgal = (max(0.,max(vgtop,vgbot))+min(0.,min(vgtop,vgbot)))*0.5

do iz=1,nnz
     read(50) u(:,:,iz), v(:,:,iz), w(:,:,iz), e(:,:,iz), &
          p(1:nnx,1:nny,iz)
end do
close(50)

!Getting the means
do iz=1,nnz
   umean(iz) = sum(u(:,:,iz))/real(nxy)
   vmean(iz) = sum(v(:,:,iz))/real(nxy)
end do

do iz=1,nLevels
   !U spectra
   fft_in = 0.5*( (u(:,:,zLevels(iz)) - umean(zLevels(iz))) + ((u(:,:,zLevels(iz)+1) - umean(zLevels(iz)+1))) )
   call dfftw_execute_dft_r2c(plan, fft_in, fft_out)
   fft_out = fft_out/real(nxy)
   do i = 2,nnx/2
      do j = 1,nny
         r = nint(sqrt(waveN(i)*waveN(i) + waveN(j)*waveN(j)))
         if( r .le. nnx/2) then
            nCount(r) = nCount(r) + 2
            uhat(r,iz) = uhat(r,iz) + 2.0*fft_out(i,j)*conjg(fft_out(i,j))
         end if
      end do
   end do
   do i = 1,nnx/2+1,nnx/2
      do j = 1,nny
         r = nint(sqrt(waveN(i)*waveN(i) + waveN(j)*waveN(j)))
         if( r .le. nnx/2) then
            nCount(r) = nCount(r) + 1
            uhat(r,iz) = uhat(r,iz) + fft_out(i,j)*conjg(fft_out(i,j))
         end if
      end do
   end do
   
   l111(iz) = l111(iz) + sum(fft_out(1,:)*conjg(fft_out(1,:))) + sum(fft_out(1,:nnx/2)*conjg(fft_out(1,:nnx/2)))
   l112(iz) = l112(iz) + sum(fft_out(:,1)*conjg(fft_out(:,1)))
   uvar(iz) = uvar(iz) + sum(abs( fft_out(:,:) )**2 ) + sum(abs(fft_out(2:nnx/2,:))**2)


   !V spectra
   fft_in = 0.5*( (v(:,:,zLevels(iz)) - vmean(zLevels(iz))) + ((v(:,:,zLevels(iz)+1) - vmean(zLevels(iz)+1))) ) 
   call dfftw_execute_dft_r2c(plan, fft_in, fft_out)
   fft_out = fft_out/real(nxy)
   do i = 2,nnx/2
      do j = 1,nny
         r = nint(sqrt(waveN(i)*waveN(i) + waveN(j)*waveN(j)))
         if( r .le. nnx/2) then
            nCount(r) = nCount(r) + 2
            vhat(r,iz) = vhat(r,iz) + 2.0*fft_out(i,j)*conjg(fft_out(i,j))
         end if
      end do
   end do
   do i = 1,nnx/2+1,nnx/2
      do j = 1,nny
         r = nint(sqrt(waveN(i)*waveN(i) + waveN(j)*waveN(j)))
         if( r .le. nnx/2) then
            nCount(r) = nCount(r) + 1
            vhat(r,iz) = vhat(r,iz) + fft_out(i,j)*conjg(fft_out(i,j))
         end if
      end do
   end do

   l221(iz) = l221(iz) + sum(fft_out(1,:)*conjg(fft_out(1,:))) + sum(fft_out(1,:nnx/2)*conjg(fft_out(1,:nnx/2)))
   l222(iz) = l222(iz) + sum(fft_out(:,1)*conjg(fft_out(:,1)))
   vvar(iz) = vvar(iz) + sum(abs( fft_out(:,:) )**2 ) + sum(abs(fft_out(2:nnx/2,:))**2)

   !W spectra
   fft_in = w(:,:,zLevels(iz))
   call dfftw_execute_dft_r2c(plan, fft_in, fft_out)
   fft_out = fft_out/real(nxy)
   do i = 2,nnx/2
      do j = 1,nny
         r = nint(sqrt(waveN(i)*waveN(i) + waveN(j)*waveN(j)))
         if( r .le. nnx/2) then
            nCount(r) = nCount(r) + 2
            what(r,iz) = what(r,iz) + 2.0*fft_out(i,j)*conjg(fft_out(i,j))
         end if
      end do
   end do
   do i = 1,nnx/2+1,nnx/2
      do j = 1,nny
         r = nint(sqrt(waveN(i)*waveN(i) + waveN(j)*waveN(j)))
         if( r .le. nnx/2) then
            nCount(r) = nCount(r) + 1
            what(r,iz) = what(r,iz) + fft_out(i,j)*conjg(fft_out(i,j))
         end if
      end do
   end do
   
   l331(iz) = l331(iz) + sum(fft_out(1,:)*conjg(fft_out(1,:))) + sum(fft_out(1,:nnx/2)*conjg(fft_out(1,:nnx/2)))
   l332(iz) = l332(iz) + sum(fft_out(:,1)*conjg(fft_out(:,1)))
   wvar(iz) = wvar(iz) + sum(abs( fft_out(:,:) )**2 ) + sum(abs(fft_out(2:nnx/2,:))**2)
   
   !Verification of spectra by calculating variance in 3 different ways
!   write(*,*) 'Level ', zLevels(iz), sum((u(:,:,zLevels(iz)) - umean(zLevels(iz)))*(u(:,:,zLevels(iz)) - umean(zLevels(iz))))/real(nxy-1), sum(abs( fft_out(:,:) )**2 ) + sum(abs(fft_out(2:nnx/2,:))**2)&
!        , sum(uhat(:,iz))
end do


read(75,*,iostat=stat)buffer
end do
close(75)

!Averaging over all the field files read so far
do i=1,nnx/2
   uhat(i,:) = uhat(i,:) * pi*((i+0.5)**2 - (i-0.5)**2) *3*nLevels*fcounter/nCount(i)
   vhat(i,:) = vhat(i,:) * pi*((i+0.5)**2 - (i-0.5)**2) *3*nLevels*fcounter/nCount(i)
   what(i,:) = what(i,:) * pi*((i+0.5)**2 - (i-0.5)**2) *3*nLevels*fcounter/nCount(i)
end do

uvar = uvar/fcounter
vvar = vvar/fcounter
wvar = wvar/fcounter


write(*,*) uvar/0.016

l111 = l111/fcounter*pi/uvar
l112 = l112/fcounter*pi/uvar
l221 = l221/fcounter*pi/vvar
l222 = l222/fcounter*pi/vvar
l331 = l331/fcounter*pi/wvar
l332 = l332/fcounter*pi/wvar

!uhat = uhat/fcounter
!vhat = vhat/fcounter
!what = what/fcounter


!write(*,*) (nCount(j),j=0,nnx/2)
!write(*,*) (waveN(j),j=1,nnx)
call dfftw_destroy_plan(plan)

open(unit=15, file="u_spectra.dat")
write(15,"(A, 16I)", advance="no") '#k ', (zLevels(iz),iz=1,nLevels)
write(15,*)
do j=0,nnx/2
   write(15,"(I,16F)", advance="no") j, (uhat(j,iz),iz=1,nLevels)
   write(15,*)
end do
close(15)


open(unit=15, file="v_spectra.dat")
write(15,"(A, 16I)", advance="no") '#k ', (zLevels(iz),iz=1,nLevels)
write(15,*)
do j=0,nnx/2
   write(15,"(I,16F)", advance="no") j, (vhat(j,iz),iz=1,nLevels)
   write(15,*)
end do
close(15)


open(unit=15, file="w_spectra.dat")
write(15,"(A, 16I)", advance="no") '#k ', (zLevels(iz),iz=1,nLevels)
write(15,*)
do j=0,nnx/2
   write(15,"(I,16F)", advance="no") j, (what(j,iz),iz=1,nLevels)
   write(15,*)
end do
close(15)


open(unit=15, file="integral_ls.dat")
write(15,*) '#k    l111   l112   l221   l222   l331    l332 '
do i=1,nLevels
   write(15,"(I,6F)", advance="no") zLevels(i), l111(i), l112(i), l221(i), l222(i), l331(i), l332(i)
   write(15,*)
end do
close(15)


end program read_write_3d
