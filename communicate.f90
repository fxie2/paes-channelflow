subroutine commun_arrays(istart, iend, nscal, abovenode, belownode, &
     mynode, comm1d, np, nptot, nm, nmtot)
  
  use storage
  use parameters
  implicit none
  include "mpif.h"
  
  integer :: np(4), nptot(4)
  integer :: nm(4), nmtot(4)
  integer istart, iend, mynode, comm1d
  integer nscal, abovenode, belownode
  integer ier, istat(MPI_STATUS_SIZE)

  real, allocatable :: buffer(:,:,:)
  
  allocate(buffer(nnx,nny,max(nptot(4),nmtot(4))))


  ! COMMUNICATE THE VARIABLES NEEDED FROM THE LEVEL BELOW
  if (iend .ne. nnz) then
     buffer(:,:,1:nmtot(1)) = u(:,:,iend-nm(1)+1:iend)
     buffer(:,:,nmtot(1)+1:nmtot(2)) = v(:,:,iend-nm(2)+1:iend)
     buffer(:,:,nmtot(2)+1:nmtot(3)) = w(:,:,iend-nm(3)+1:iend)
     if (isbgr) buffer(:,:,nmtot(3)+1:nmtot(4)) = e(:,:,iend-nm(4)+1:iend)

     call mpi_send(buffer, nmtot(4)*nxy, mpi_real, abovenode, mynode, &
          comm1d, ier)
     
!     write(*,*) 'MPI_SEND FROM LEVEL ', mynode , ' TO ', abovenode
  endif
  
  if (istart .ne. 1) then
     call mpi_recv(buffer, nmtot(4)*nxy, mpi_real, belownode, belownode, &
          comm1d, istat, ier)
!     write(*,*) 'MPI_RECV FROM LEVEL ', belownode , ' TO ', mynode
          
     u(:,:,istart-nm(1):istart-1) = buffer(:,:,1:nmtot(1))
     v(:,:,istart-nm(2):istart-1) = buffer(:,:,nmtot(1)+1:nmtot(2))
     w(:,:,istart-nm(3):istart-1) = buffer(:,:,nmtot(2)+1:nmtot(3))
     if (isbgr) e(:,:,istart-nm(4):istart-1) = buffer(:,:,nmtot(3)+1:nmtot(4))
     
     ! COMMUNICATE THE VARIABLES NEEDED FROM THE LEVEL ABOVE
     
     buffer(:,:,1:nptot(1)) = u(:,:,istart:istart+np(1)-1)
     buffer(:,:,nptot(1)+1:nptot(2)) = v(:,:,istart:istart+np(2)-1)
     buffer(:,:,nptot(2)+1:nptot(3)) = w(:,:,istart:istart+np(3)-1)
     if (isbgr) buffer(:,:,nptot(3)+1:nptot(4)) = e(:,:,istart:istart+np(4)-1)
    
     call mpi_send(buffer, nptot(4)*nxy, mpi_real, belownode, &
          mynode, comm1d, ier)
!     write(*,*) 'MPI_SEND FROM LEVEL ', mynode , ' TO ', belownode
  endif
  
  if (iend .ne. nnz) then
     call mpi_recv(buffer, nptot(4)*nxy, mpi_real, abovenode, &
          abovenode, comm1d, istat, ier)
!     write(*,*) 'MPI_RECV FROM LEVEL ',  abovenode , ' TO ', mynode
     u(:,:,iend+1:iend+np(1)) = buffer(:,:,1:nptot(1))
     v(:,:,iend+1:iend+np(2)) = buffer(:,:,nptot(1)+1:nptot(2))
     w(:,:,iend+1:iend+np(3)) = buffer(:,:,nptot(2)+1:nptot(3))
     if (isbgr) e(:,:,iend+1:iend+np(4)) = buffer(:,:,nptot(3)+1:nptot(4))

  endif

  deallocate(buffer)
  
end subroutine commun_arrays





subroutine commun_rw(istart, iend, nscal, abovenode, belownode, &
     mynode, comm1d)
  
  use storage
  use parameters
  implicit none
  include "mpif.h"
  
  integer istart, iend, mynode, comm1d
  integer nscal, abovenode, belownode
  integer ier, istat(MPI_STATUS_SIZE)
  real buffer(nnx,nny,2)
     
  if (iend .ne. nnz) then
     buffer(:,:,1) = rw(:,:,iend)
     buffer(:,:,2) = w(:,:,iend)
     call mpi_send(buffer, 2*nxy, mpi_real, abovenode, mynode, &
          comm1d, ier)
  endif
  if (istart .ne. 1) then
     call mpi_recv(buffer, 2*nxy, mpi_real, belownode, &
          belownode, comm1d, istat, ier)
     rw(:,:,istart-1)=buffer(:,:,1)
     w(:,:,istart-1)=buffer(:,:,2)
  endif
  
end subroutine commun_rw




subroutine commun_press(istart, iend, nscal, abovenode, belownode, &
     mynode, comm1d)
  
  use storage
  use parameters
  implicit none
  include "mpif.h"
  
  integer istart, iend, mynode, comm1d
  integer nscal, abovenode, belownode
  integer ier, istat(MPI_STATUS_SIZE)
  real buffer(nnx,nny)
  
  if (istart .ne. 1) then
     buffer(:,:) = p(1:nnx,:,istart)
     call mpi_send(buffer, nxy, mpi_real, belownode, mynode, comm1d, ier)
  endif
     
  if (iend .ne. nnz) then
     call mpi_recv(buffer, nxy, mpi_real, abovenode, abovenode, comm1d, &
          istat, ier)
     p(1:nnx,:,iend+1) = buffer(:,:)
  endif

end subroutine commun_press










