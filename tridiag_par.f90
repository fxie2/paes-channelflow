subroutine trid_solver(proc_numb, nb, j_end, i_start, i_end, &
     communicator, abcd_matrices, e_matrix, bl_abcd_matrices)

  !=======================================================================
  !
  !     purpose : solution of tridiagonal system(s).
  !               Data assumed distributed on parallel processors.
  !               The basic equation could be written
  !
  !               a(i) * x(i-1) + b(i) * x(i) + c(i) * x(i+1) = d(i)
  !
  !               x = unknown. a, b, c, d = coefficients. i = .
  !
  !               We allow more than one variable to be found on one call.
  !               The number of variables is num_vars.
  !
  !     method  : parallel block reduction
  !               See UCRL-JC-114756 by Mattor, Williams, and Hewett
  !               for related methods.
  !               Actual method used here independently derived by
  !               Eltgroth at LLNL, who may someday provide notes.
  !               The basic idea is to eliminate interior variables from
  !               the block of unknowns owned by each processor.
  !               Then each sends the border information (coefficients)
  !               to some (other) processor.
  !               There the info forms a new tridiagonal system with 2 *
  !               (the number of blocks) unknowns.
  !               That system is solved anyway you want (c.f., 
  !               final_trid_solver), and border
  !               values are sent back to the original block processors.
  !               While or after the border solution takes place, the
  !               block processors prepare other coefficients to use the
  !               border values. When the border values are returned,
  !               each block processor computes its internal unknowns.
  !               A variation on the above is allowed (see variable
  !               full_exchange.) In this case, one complete exchange
  !               of border information is accomplished, and each
  !               processor solves the reduced system. No further
  !               communication is needed.
  !
  !               The basic assumption about data is that each block
  !               processor "owns" a contiguous set of unknowns and
  !               coefficients.
  !
  !      input  : proc_numb is the # of this processor, 1 <= proc_numb <= nb
  !
  !               nb (num_block_procs)  is the number of blocks (and
  !                                     processors) used to solve the
  !                                     block system.
  !
  !               reduce_proc_numb is the reduced system processor id.
  !
  !               NOTE: See declarations below (in code) for array sizes
  !                     and shapes.
  !
  !               abcd_matrices stores the a, b, c, and d coefficents
  !                   indexed (in order) by j, i, and
  !                   c_id coefficient id where 1->a, 2->b, 3->c, 4->d.
  !                   It is defined
  !                   with j dimensions 1 : j_end,
  !                   i dimensions i_start : i_end.
  !                   Inside a block, j runs from 1 to j_end;
  !                   i (the tridiagonal ) runs from i_start to i_end.
  !
  !               scratch arrays are spaces for real values allocated
  !                   by the user to be used by the solver.
  !
  !     author  : pete eltgroth
  !
  !     date    : 6 december 1993
  !
  !=======================================================================
  ! EXTERNAL VARIABLES.

  implicit none
  include 'mpif.h'

  integer proc_numb, nb
  integer j_end, i_start, i_end
  integer communicator
  real :: abcd_matrices(j_end, i_start:i_end, 4)
  real :: bl_abcd_matrices(j_end, 2*nb, 3)
  real :: e_matrix(j_end, i_start:i_end, 2)

  ! Processor(s) solving the reduced level must have additional space
  ! to solve for the border data on all the blocks.
  ! The space is taken up by the arrays starting with bl_.
  ! In order to arrange for this to be available, it is declared
  ! for each processor, but only has to be allocated and used
  ! by the reduced level processor(s).

  real :: coef_matrix(j_end,2,3)
  
  integer stat(MPI_STATUS_SIZE)
  integer num_block_procs
  integer reduce_proc_numb
  integer indx, indxp, m_len, msg_proc, tag
  integer i, j, ierr
  real factor

  num_block_procs = nb
  reduce_proc_numb = 0
  
  ! NORMALIZATION
  ! NOTE THIS CHANGES INPUT DATA!
  ! BECAUSE OF THE NORMALIZATION, WE HAVE THE B COEFFICIENTS
  ! EQUAL TO ONE AFTER THIS POINT.

  do i = i_start, i_end
     do j = 1, j_end
        factor = 1.0 / abcd_matrices(j, i, 2)
        abcd_matrices(j, i, 1) = abcd_matrices(j, i, 1) * factor
        abcd_matrices(j, i, 3) = abcd_matrices(j, i, 3) * factor
        abcd_matrices(j, i, 4) = abcd_matrices(j, i, 4) * factor
     enddo
  enddo

  ! OBTAIN COEFFICIENTS AT BORDERS FOR USE BY REDUCED SYSTEM SOLVER.
  
  coef_matrix = bl_abcd_matrices(:,1:2,:)
  
  call find_border_coeffs(j_end, i_start, i_end, abcd_matrices, &
       e_matrix, coef_matrix)
  
  ! THE BORDER COEFFICIENTS IN COEF_MATRIX MUST BE SENT TO THE
  ! SOLVER AT THE REDUCED LEVEL. 
  ! NOTE THE SPECIFIC ASSUMPTION THAT THE LOWEST NUMBER PROCESSOR
  ! HAS THE LEFT BOUNDARY AND THE HIGHEST NUMBERED PROCESSOR HAS
  ! THE RIGHT BOUNDARY.
  ! PROCESSOR NUMBERS INCREASE LEFT TO RIGHT.

  if (proc_numb .ne. reduce_proc_numb) then
     
     ! SEND BORDER INFORMATION (IF NOT THE REDUCED PROCESSOR.)
     
     m_len = j_end * 6
     tag = proc_numb
     call MPI_SEND(coef_matrix, m_len, MPI_REAL, reduce_proc_numb, &
          tag, communicator, ierr)
     if (ierr .ne. MPI_SUCCESS) call stopcode("Error sending border coefs")

     ! RECEIVE BORDER SOLUTION (IF NOT THE REDUCED PROCESSOR.)
     
     m_len = 2*j_end
     tag = nb+proc_numb

     call MPI_RECV(bl_abcd_matrices(1,1,3), m_len, MPI_REAL, &
          reduce_proc_numb, tag, communicator, stat, ierr)
     if (ierr .ne. MPI_SUCCESS) call stopcode("Error receiving border soln")

  else
     
     bl_abcd_matrices(:,1:2,:) = coef_matrix
     
     ! PREPARE THE REDUCED SOLUTION.
     ! FIRST RECEIVE BORDER COEFFICIENTS

     m_len = j_end * 6
     do i = 2, num_block_procs
        indx = 2 * i - 1
        indxp = indx + 1
        
        msg_proc = i-1
        tag = msg_proc
        call MPI_RECV (coef_matrix, m_len, MPI_REAL, msg_proc, &
             tag, communicator, stat, ierr)
        
        if (ierr .ne. MPI_SUCCESS) call stopcode("Error recvng border coefs")
        bl_abcd_matrices(:,indx:indxp,:) = coef_matrix
     enddo
     
     ! THE FOLLOWING TRIDIAGONAL SOLVER SHOULD BE OPTIMIZED
     ! SINCE THE BLOCK PROCESSORS ARE WAITING FOR RESULTS.

     call final_trid_solver(j_end, 2*num_block_procs, bl_abcd_matrices)

     ! NOW SEND THE REDUCED SOLUTION
     
     m_len = 2*j_end
     do i=2, nb
        indx = 2*i-1
        
        tag = nb+i-1
        call MPI_SEND (bl_abcd_matrices(1,indx,3), m_len, MPI_REAL, i-1, &
             tag, communicator, ierr)
        if (ierr .ne. MPI_SUCCESS) call stopcode("Error sending border soln")
     enddo
     
  endif

  ! PREPARE TO CALCULATE INTERIOR INFORMATION WHEN THE BORDER
  ! VALUES ARE RETURNED.
  ! THIS CAN OVERLAP THE BORDER SOLUTION CALCULATION AND DATA
  ! TRANSFER.

  e_matrix(:,i_start,1) = 0.0
  
  do i = i_start + 1, i_end - 1
     do j = 1, j_end
        e_matrix(j, i, 2) = 1.0 / (1.0 + e_matrix(j, i-1, 1) * &
             abcd_matrices(j, i, 1))
        e_matrix(j, i, 1) = -abcd_matrices(j, i, 3) * &
             e_matrix(j, i, 2)
     enddo
  enddo
  
  ! IN THIS IMPLEMENTATION THE BORDER SOLUTIONS ARE STORED
  ! DIRECTLY INTO THE SOLUTION.

  abcd_matrices(:,i_start,4) = bl_abcd_matrices(:,1,3)
  abcd_matrices(:,i_end,  4) = bl_abcd_matrices(:,2,3)

  ! CALCULATE INTERIOR SOLUTION

  e_matrix(:, i_start, 2) = abcd_matrices(:, i_start, 4)
     
  do i = i_start + 1, i_end - 1
     do j = 1, j_end
        e_matrix(j, i, 2) = (abcd_matrices(j, i, 4) - &
             abcd_matrices(j, i, 1) * e_matrix(j, i - 1, 2)) * &
             e_matrix(j, i, 2)
     enddo
  enddo

  do i = i_end - 1, i_start + 1, -1
     do j = 1, j_end
        abcd_matrices(j, i,4) = e_matrix(j, i, 1) * &
             abcd_matrices(j, i+1, 4) + e_matrix(j, i, 2)
     enddo
  enddo

end subroutine trid_solver








subroutine final_trid_solver(j_end, i_end, fabcd)
  !=======================================================================
  !
  !     purpose : solution of tridiagonal system(s).

  !               The basic equation could be written
  !               a(i) * x(i-1) + b(i) * x(i) + c(i) * x(i+1) = d(i)
  !               We allow more than one variable to be found on one call.
  !               The number of variables is num_vars.
  !
  !     data    : data structures are compatable with those used
  !               in the ice model of the ocean code at llnl.

  !               j runs from j_start to j_end.
  !               var_id runs from 1 to num_vars.
  !               i (the tridiagonal index) runs from i_start to i_end.

  !               fabcd_matrices stores the a, b, c, and d coefficents
  !               indexed (in order) by
  !               j latitude
  !               i longitude, and
  !               c_id coefficient id where 1->a, 2->b, 3->c, 4->d.

  !     output  : fsolution_vectr contains results for variable x in
  !               the given domain.
  !
  !     author  : pete eltgroth
  !
  !     date    : 6 december 1993
  !
  !=======================================================================
  implicit none
  integer j_end, i_start, i_end

  real :: fabcd(j_end, i_end, 3)
  integer :: i, j, im1, ii
  real :: fac

  do j=1,j_end
     fabcd(j,2,3) = fabcd(j,2,3) - fabcd(j,2,1)*fabcd(j,1,3)
     fabcd(j,2,1) = 1./(1.-fabcd(j,2,1)*fabcd(j,1,2))
  enddo

  do i=3,i_end
     do j=1,j_end
        fac = fabcd(j,i,1)*fabcd(j,i-1,1)
        fabcd(j,i,1) = 1./(1.-fac*fabcd(j,i-1,2))
        fabcd(j,i,3) = fabcd(j,i,3) - fac*fabcd(j,i-1,3)
     enddo
  enddo

  do j=1,j_end
     fabcd(j,i_end,3) = fabcd(j,i_end,3)*fabcd(j,i_end,1)
  enddo

  do i=i_end-1,2,-1
     do j=1,j_end
        fabcd(j,i,3) = fabcd(j,i,1)*(fabcd(j,i,3) - fabcd(j,i,2) &
             *fabcd(j,i+1,3))
     enddo
  enddo

   do j=1,j_end
      fabcd(j,1,3) = fabcd(j,1,3) - fabcd(j,1,2)*fabcd(j,2,3)
   enddo

end subroutine final_trid_solver
      
      
      
      
      
     
subroutine find_border_coeffs(j_end, i_start, i_end, &
     abcd_matrices, e_matrix, coef_matrix)

  ! Finds coefficients at borders of block.
  ! Normalize so that the new diagonal values will be one.
      
  implicit none
  integer j_end, i_start, i_end

  real abcd_matrices(j_end, i_start:i_end, 4)
  real e_matrix(j_end, i_start:i_end,2)
  real coef_matrix(j_end, 2, 3)

  ! In coef_matrix, the last index has the following meanings:
  !    1 => a; 2 => b; 3 => c; 4 => d.


  ! The variable which_way is +1 if finding values at the left
  ! and                       -1 if finding values to the right.

  integer which_way

  integer i, j, ipass, indx_border
  integer indx_other_border, c_id_border, c_id_other_border
  integer indx_next_to_border, indx_next_to_other_border
  integer coef_indx
  real norm_factor

  integer ci, cio
  
  ! We loop through the procedure twice, once for left border, &
  ! and once for the right border.

  do ipass = 1,2
     
     if (ipass.eq.1) then
        which_way = 1
        coef_indx = 1
        indx_border = i_start
        indx_other_border = i_end
        c_id_border = 1
        c_id_other_border = 3
        ci = 1
        cio = 2
     else
        which_way = -1
        coef_indx = 2
        indx_border = i_end
        indx_other_border = i_start
        c_id_border = 3
        c_id_other_border = 1
        ci = 2
        cio = 1
     endif

     indx_next_to_border = indx_border + which_way
     indx_next_to_other_border = indx_other_border - which_way

     do j = 1, j_end
        e_matrix(j, indx_other_border, 1) = 0.0
     enddo

     do i = indx_next_to_other_border, indx_next_to_border, -which_way
        do j = 1, j_end
           e_matrix(j, i, 1) = &
                -abcd_matrices(j, i-which_way, c_id_other_border) / &
                (1.0 + e_matrix(j, i+which_way, 1) * &
                abcd_matrices(j, i+which_way, c_id_border))
        enddo
     enddo

     ! calculation of the new coefficients for border of
     ! reduced system.

     do j = 1, j_end
        e_matrix(j, i_start,2) = 1.0
        coef_matrix(j, coef_indx, 3) = abcd_matrices(j, indx_border, 4)
     enddo

     do i = indx_next_to_border, indx_next_to_other_border, which_way
        do j = 1, j_end
           e_matrix(j, i_start, 2) = e_matrix(j, i_start, 2) * &
                e_matrix(j, i, 1)
           coef_matrix(j, coef_indx, 3) = coef_matrix(j, coef_indx, 3) &
                + e_matrix(j, i_start, 2) * abcd_matrices(j, i, 4)
        enddo
     enddo

     ! NORMALIZE

     do j = 1, j_end
        norm_factor = 1.0 / (1.0 + e_matrix(j, indx_next_to_border, 1) * &
             abcd_matrices(j, indx_next_to_border, c_id_border))

        coef_matrix(j, coef_indx, cio) = ( &
             e_matrix(j, i_start, 2) * &
             abcd_matrices(j, indx_next_to_other_border, &
             c_id_other_border)) * norm_factor

        coef_matrix(j, coef_indx, ci) = &
             abcd_matrices(j, indx_border, c_id_border) &
             * norm_factor

        coef_matrix(j, coef_indx, 3) = &
             coef_matrix(j, coef_indx, 3) * &
             norm_factor
     enddo

  enddo
end subroutine find_border_coeffs
      
      




      

subroutine stopcode(msg)

  character*(*) msg
  integer ier

  write (6,*) msg
  call mpi_abort(MPI_COMM_WORLD,1,ier)

end subroutine stopcode













