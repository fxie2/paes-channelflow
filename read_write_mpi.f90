subroutine write_field_mpi(procs,mynode,filename,istart,iend,comm1d)
  use parameters
  use storage
  implicit none
  integer procs, mynode, istart, iend, iz, j, ier, comm1d
  character*25 filename

  call mpi_barrier(comm1d, ier)
  do j=1,procs
     if (mynode.eq.j-1) then
        open (14, file=filename, form="unformatted", position="append")
        do iz=istart,iend
           call sav_field(iz)
        enddo
        close(14)
     endif
     call mpi_barrier(comm1d, ier)
  enddo
end subroutine write_field_mpi


subroutine read_restart_mpi(istart,iend,procs,mynode,comm1d,nvar, &
     qscal_sav,cscal_sav)

  use parameters
  use storage
  implicit none
  include "mpif.h"

  logical qscal_sav, cscal_sav
  integer procs, mynode, istart, iend, iz, j, ier, comm1d, nvar
  integer, allocatable :: ist(:), ien(:)
  integer istat(MPI_STATUS_SIZE)
  real, allocatable :: buffer(:,:,:)
  
  allocate(buffer(nnx,nny,nvar+1))
  call mpi_barrier(comm1d, ier)
  if (mynode.eq.procs-1) then
     allocate(ist(procs))
     allocate(ien(procs))
  endif
  
  call mpi_gather(istart,1,mpi_real,ist,1,mpi_real,procs-1,comm1d,ier)
  call mpi_gather(iend,1,mpi_real,ien,1,mpi_real,procs-1,comm1d,ier)
  
  if (mynode.eq.procs-1) then
     do j=1,procs
        do iz=ist(j),ien(j)
           read(50) buffer
           if (j.ne.procs) then
              call mpi_send(buffer, (nvar+1)*nxy, mpi_real, j-1, mynode, &
                   comm1d, ier)
           else
              u(:,:,iz) = buffer(:,:,1)
              v(:,:,iz) = buffer(:,:,2)
              w(:,:,iz) = buffer(:,:,3)
              e(:,:,iz) = buffer(:,:,4)
           endif
        enddo
     enddo
     close(50)
  else
     close(50)
     do iz=istart,iend
        call mpi_recv(buffer, (nvar+1)*nxy, mpi_real, procs-1, &
             procs-1, comm1d, istat, ier)
        u(:,:,iz) = buffer(:,:,1)
        v(:,:,iz) = buffer(:,:,2)
        w(:,:,iz) = buffer(:,:,3)
        e(:,:,iz) = buffer(:,:,4)
     enddo
  endif

  deallocate(buffer)
  if (mynode.eq.procs-1) then
     deallocate(ist)
     deallocate(ien)
  endif

end subroutine read_restart_mpi
  





