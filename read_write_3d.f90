program read_write_3d
!This program reads the files that contain the 3d data (velocity, pressure, scalars) in the ".field" files and puts them into a tecplot file for visualization.
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
real, allocatable, dimension(:,:,:) :: u,v,w,e,p
real, allocatable, dimension(:) :: umean, vmean


call getarg(1,buffer)
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

write(*,*) 'nnx = ', nnx
write(*,*) 'nny = ', nny
write(*,*) 'nnz = ', nnz


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
end do

dx = xl/nnx
dy = yl/nny
dz = zl/nnz

!Writing to tecplot file
!Grid file

open(unit=15,file="resolved_field_grid.dat")
write(15,*)"TITLE = ""Resolved field grid"" "
write(15,*)"FILETYPE = GRID"
write(15,*)"VARIABLES = ""X"" ""Y"" ""Z"" "
write(15,*)"ZONE"
write(15,*)"I = ", nnx , ", J = ", nny, " K = ", nnz
write(15,*)"ZONETYPE = Ordered, DATAPACKING = BLOCK"
write(*,*) 'Writing grid'
do k = 1,nnz
   do j= 1,nny
      write(15,*) (dx*real(i-1), i=1,nnx)
   end do
end do
write(*,*) 'Wrote x'
do k = 1,nnz
   do j= 1,nny
      write(15,*) (dy*real(j-1), i=1,nnx)
   end do
end do
write(*,*) 'Wrote y'
do k = 1,nnz
   do j= 1,nny
      write(15,*) (dz * (real(k)-0.5) , i=1,nnx)
   end do
end do
write(*,*) 'Wrote z'
close(15)

!Solution file

open(unit=15,file="resolved_field_solution.dat")
write(15,*)"TITLE = ""Resolved Field solution"" "
write(15,*)"FILETYPE = SOLUTION"
write(15,*)"VARIABLES = ""U"" ""V"" ""W"" ""E"" ""Uprime"" ""Vprime"" "
write(15,*)"ZONE"
write(15,*)"I = ", nnx , ", J = ", nny, " K = ", nnz
write(15,*)"ZONETYPE = Ordered, DATAPACKING = BLOCK"
do k = 1,nnz
   do j= 1,nny
      write(15,*) (U(i,j,k)+ugal, i=1,nnx)
   end do
end do
write(*,*) 'Wrote U'
do k = 1,nnz
   do j= 1,nny
      write(15,*) (V(i,j,k)+vgal, i=1,nnx)
   end do
end do
write(*,*) 'Wrote V'
do k = 1,nnz
   do j= 1,nny
      write(15,*) (W(i,j,k), i=1,nnx)
   end do
end do
write(*,*) 'Wrote W'
do k = 1,nnz
   do j= 1,nny
      write(15,*) (E(i,j,k), i=1,nnx)
   end do
end do
write(*,*) 'Wrote E'
do k = 1,nnz
   do j= 1,nny
      write(15,*) (U(i,j,k)-umean(k), i=1,nnx)
   end do
end do
write(*,*) 'Wrote Uprime'
do k = 1,nnz
   do j= 1,nny
      write(15,*) (V(i,j,k)-vmean(k), i=1,nnx)
   end do
end do
write(*,*) 'Wrote Vprime'
close(15)

end program read_write_3d
