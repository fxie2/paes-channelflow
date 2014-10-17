program read_write_3d
!This program reads the files that contain the 3d data (velocity, pressure, scalars) in the ".field" files and calculates the integral length scales

implicit none

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
real, allocatable, dimension(:) :: umean, vmean, wmean

real, allocatable, dimension(:,:) :: cuux,cuuy,cvvx,cvvy,cwwx,cwwy !Cross correlation coefficient of u(x)u(x+lx) at various z levels
real, allocatable, dimension(:) :: luux !Integral length scale of cuux at various z levels


integer :: stat !To determine if EOF is reached while reading file

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
allocate(wmean(nnz))


allocate(cuux(nnx,nnz))
allocate(cvvx(nnx,nnz))
allocate(cwwx(nnx,nnz))
allocate(cuuy(nnx,nnz))
allocate(cvvy(nnx,nnz))
allocate(cwwy(nnx,nnz))

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

cuux = 0
cuuy = 0
cvvx = 0
cvvy = 0
cwwx = 0
cwwy = 0
luux = 0

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

do iz=1,nnz
   umean(iz) = sum(u(:,:,iz))/real(nxy)
   vmean(iz) = sum(v(:,:,iz))/real(nxy)
   wmean(iz) = sum(w(:,:,iz))/real(nxy)
end do

do iz=1,nnz
   do j=0,nnx-1
      do i=1,nnx-j
         cuux(j+1,iz) = cuux(j+1,iz) + sum( (u(i,:,iz)-umean(iz))* (u(i+j,:,iz)-umean(iz)) )
         cvvx(j+1,iz) = cvvx(j+1,iz) + sum( (v(i,:,iz)-vmean(iz))* (v(i+j,:,iz)-vmean(iz)) )
         cwwx(j+1,iz) = cwwx(j+1,iz) + sum( (w(i,:,iz)-wmean(iz))* (w(i+j,:,iz)-wmean(iz)) )
         cuuy(j+1,iz) = cuuy(j+1,iz) + sum( (u(:,i,iz)-umean(iz))* (u(:,i+j,iz)-umean(iz)) )
         cvvy(j+1,iz) = cvvy(j+1,iz) + sum( (v(:,i,iz)-vmean(iz))* (v(:,i+j,iz)-vmean(iz)) )
         cwwy(j+1,iz) = cwwy(j+1,iz) + sum( (w(:,i,iz)-wmean(iz))* (w(:,i+j,iz)-wmean(iz)) )
      end do
!      cuux(j+1,iz) = cuux(j+1,iz)/((nnx-j)*nny)
   end do
end do


read(75,*,iostat=stat)buffer
end do
close(75)

do j=0,nnx-1
   cuux(j+1,:) = cuux(j+1,:)/((nnx-j)*nny*fcounter)
   cvvx(j+1,:) = cvvx(j+1,:)/((nnx-j)*nny*fcounter)
   cwwx(j+1,:) = cwwx(j+1,:)/((nnx-j)*nny*fcounter)
   cuuy(j+1,:) = cuuy(j+1,:)/((nnx-j)*nny*fcounter)
   cvvy(j+1,:) = cvvy(j+1,:)/((nnx-j)*nny*fcounter)
   cwwy(j+1,:) = cwwy(j+1,:)/((nnx-j)*nny*fcounter)
end do



open(unit=15,file="integral_ls.dat")
write(15,"(A,96F)",advance="no") '#rx ', (iz*dz - dz*0.5,iz=1,nnz)
write(15,*)
do j=0,nnx-1
   write(15,"(97F)",advance="no") j*dx, (cuux(j+1,iz)/cuux(1,iz),iz=1,nnz)
   write(15,*)
end do
do j=0,nnx-1
   write(15,"(97F)",advance="no") j*dx, (cuuy(j+1,iz)/cuuy(1,iz),iz=1,nnz)
   write(15,*)
end do
do j=0,nnx-1
   write(15,"(97F)",advance="no") j*dx, (cvvx(j+1,iz)/cvvx(1,iz),iz=1,nnz)
   write(15,*)
end do
do j=0,nnx-1
   write(15,"(97F)",advance="no") j*dx, (cvvy(j+1,iz)/cvvy(1,iz),iz=1,nnz)
   write(15,*)
end do
do j=0,nnx-1
   write(15,"(97F)",advance="no") j*dx, (cwwx(j+1,iz)/cwwx(1,iz),iz=1,nnz)
   write(15,*)
end do
do j=0,nnx-1
   write(15,"(97F)",advance="no") j*dx, (cwwy(j+1,iz)/cwwy(1,iz),iz=1,nnz)
   write(15,*)
end do
close(15)

end program read_write_3d
