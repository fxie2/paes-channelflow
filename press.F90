subroutine press(istart,iend,procs,mynode,comm1d)

  use storage
  use compbuf
  use parameters
  implicit none

  integer :: istart,iend,procs,comm1d,mynode
  integer :: ibot, itop, iz, kp, lp
  real, target  :: bdar(nny,istart:iend,4)
  real :: dzdz
  
#ifdef MPI
! EXTRA STORAGE FOR PARALLEL TRIDIAGONAL SOLVER
  real :: sc1(nny,istart:iend,2)
  real :: sc2(nny,2*procs,3)
#endif

  dzdz=dz*dz

#ifdef MPI      
  if (istart.eq.1) then
     ibot=2
  else
     ibot=istart
  endif
  if (iend.eq.nnz) then
     itop=nnz-1
  else
     itop=iend
  endif
#else 
  ibot=2
  itop=nnz-1
#endif

  ! FOURIER ANALYZE THE RIGHT HAND SIDE AND UPPER B.C.
  
  do iz=istart,iend
     call ffthrz(p(1,1,iz),nnx,nny,-1)
  enddo
     
  do kp=1,nnx
     do iz=ibot,itop
        do lp=1,nny
           bdar(lp,iz,1)  = 1.0
           bdar(lp,iz,2) = -dzdz*xk2(kp,lp) - 2.0
           bdar(lp,iz,3)  = 1.0
           bdar(lp,iz,4) = dzdz*p(kp,lp,iz)
        enddo
     enddo

     ! LOWER BOUNDARY
#ifdef MPI
     if (istart.eq.1) then
#endif
        do lp=1,nny
           bdar(lp,1,1) = 1.0
           bdar(lp,1,2) = -dzdz*xk2(kp,lp) - 1.0
           bdar(lp,1,3) = 1.0
           bdar(lp,1,4) = dzdz*p(kp,lp,1)
        enddo
#ifdef MPI
     endif
#endif
     
     ! UPPER BOUNDARY
#ifdef MPI
     if (iend.eq.nnz) then
#endif
        do lp=1,nny
           bdar(lp,nnz,1) = 1.0
           bdar(lp,nnz,2) = -dzdz*xk2(kp,lp) - 1.0
           bdar(lp,nnz,3) = 1.0
           bdar(lp,nnz,4) = dzdz*p(kp,lp,nnz)
        enddo
#ifdef MPI
     endif
#endif
     ! SPECIAL SITUATION FOR ZEROTH MODE
     if (kp .eq. 1 .or. kp .eq. 2) then
        do iz=istart,iend
           bdar(1,iz,1) = 0.0
           bdar(1,iz,2) = 1.0
           bdar(1,iz,3) = 0.0
           bdar(1,iz,4) = 0.0
        enddo
     endif

     ! NOW SOLVE THE SYSTEM
#ifdef MPI
     if (procs.eq.1) then
#endif
        call tridv(bdar,nny,nnz)
#ifdef MPI            
     else
        call trid_solver(mynode, procs, nny, istart, &
             iend, comm1d, bdar, sc1, sc2)
     endif
#endif
        p(kp,:,istart:iend) = bdar(:,:,4) 
  enddo

  ! TRANSFER BACK TO THE PHYSICAL SPACE, P IS DEV FROM HORIZ MEAN
  ! THE VALUE OF P AT NNZ NO MATTER, BECAUSE W AT NNZ-1 NOT COMP.
  do iz=istart,iend
    p(1,1,iz)=0.0
    p(2,1,iz)=0.0
    p(nnx+1:nnx+2,:,iz)=0.0
    if (nny .ne. 1) p(:,nny/2+1,iz)=0.0
    call ffthrz(p(1,1,iz),nnx,nny,+1)
    p(nnx+1:nnx+2,:,iz)=0.0
  enddo

end subroutine press







