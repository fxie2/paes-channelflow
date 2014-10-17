subroutine lesabort(message)
  implicit none
  
#ifdef MPI  
  include "mpif.h"
#endif
  
  integer ier
  character*(*) :: message
  write(*,*) message
  
#ifdef MPI  
  call mpi_abort(MPI_COMM_WORLD,1,ier)
#endif
  stop

end subroutine lesabort

