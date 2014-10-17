program read_write_3d
!This program reads the files that contain the 3d data (velocity, pressure, scalars) in the ".field" files and calculates the PDF
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

!PDF variables
real, allocatable, dimension(:,:) :: vtemp !Temporary variable to store interpolated velocity on a plane (if reqd)
real, allocatable, dimension(:,:) :: updf, vpdf, wpdf !Array to store the spectrum at pre-defined z levels
real, allocatable, dimension(:) :: waveN ! Wavenumber
integer :: nLevels = 16 !Number of levels at which the spectrum is to be measured.
integer, dimension(16) :: zLevels = (/1,2,3,4,5,6,7,8,9,10,12,18,24,30,36,48/) !The levels at which the PDF is to be calculated
real :: pi = 3.14159265358979
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

allocate(vtemp(nnx,nny))

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
allocate(updf( , nLevels))
allocate(vpdf( , nLevels))
allocate(wpdf( , nLevels))


allocate(uvar(nLevels))
allocate(vvar(nLevels))
allocate(wvar(nLevels))
uvar = 0
vvar = 0
wvar = 0

uhat = 0
vhat = 0
what = 0


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
  write(*,*) 'Finished reading level ', iz
end do
close(50)

!Getting the means
do iz=1,nnz
   umean(iz) = sum(u(:,:,iz))/real(nxy)
   vmean(iz) = sum(v(:,:,iz))/real(nxy)
end do

do iz=1,nLevels
   !U PDF
   vtemp = 0.5*( (u(:,:,zLevels(iz)) - umean(zLevels(iz))) + ((u(:,:,zLevels(iz)+1) - umean(zLevels(iz)+1))) )

   !V PDF
   vtemp = 0.5*( (v(:,:,zLevels(iz)) - vmean(zLevels(iz))) + ((v(:,:,zLevels(iz)+1) - vmean(zLevels(iz)+1))) ) 


   !W spectra
   vtemp = w(:,:,zLevels(iz))
   
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

l111 = l111/fcounter*pi*0.5/uvar
l112 = l112/fcounter*pi*0.5/uvar
l221 = l221/fcounter*pi*0.5/vvar
l222 = l222/fcounter*pi*0.5/vvar
l331 = l331/fcounter*pi*0.5/wvar
l332 = l332/fcounter*pi*0.5/wvar
!uhat = uhat/fcounter
!vhat = vhat/fcounter
!what = what/fcounter


write(*,*) (nCount(j),j=0,nnx/2)
write(*,*) (waveN(j),j=1,nnx)
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
