subroutine compspectra(istart,iend,comm1d,mynode,procs,spectfile)
!Modified by Tie Wei 02/28/06 to compute and output spectra
!at (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60 /)
!motivation: to study the near wall region matching that of Yong's code
!Changes:
!lev(5)--> lev(15), levels(15) not used
!change real :: eta(5) to real :: eta(15)
!skip the relation between lev(:) and levels(:) :(levels(:) not used)
!remove the spectav: the averaging process:
!change the output as: write(30) ncx, 15 !15 is ilev
!***Need to modify allocate.f90: allocate(spectbuf(ncx,spectcount,15))
!
  
  use storage
  use parameters
  use compbuf
  implicit none
  
# ifdef MPI
  include "mpif.h"
# endif
  
  integer :: istart, iend, comm1d, mynode, procs, i, ier, iz
  real, pointer, dimension(:,:) :: ulocal, vlocal, tlocal
  real, pointer, dimension(:,:) :: clocal, qlocal
  real, dimension(nnx,nny) :: plocal, p1, p2
  INTEGER, DIMENSION(15) :: lev=(/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, &
				  20, 30, 40, 50, 60 /) 
  INTEGER:: il, ilev=15
  real :: eta(15)
  integer :: nsarray(procs), nstart(procs), ns
  character(25) :: spectfile
  nsarray = 0
  nstart = 0
  
  ulocal => rubuf(1:nnx,1:nny)
  vlocal => rvbuf(1:nnx,1:nny)
  tlocal => rtbuf(1:nnx,1:nny)
  if (iqscal) qlocal => rqbuf(1:nnx,1:nny)
  if (icscal) clocal => rcbuf(1:nnx,1:nny)


  eta = real(lev)/real(izi(1))


!TWO DO-LOOPS
!FIRST LOOP AROUND ISTART:IEND, SECHOD LOOP AROUND IL=1:ILEV
!IF IZ=LEVELS(IL): PROCESS THE DATA ON THAT GRID LEVEL
!NOTE THE TRICK HERE IS: LEVELS(IL) HAVE GRID LEVELS IN DIFFERENT CPU
!YOU NEED TO FIND OUT HOW MANY LEVELS ARE IN THE CURRENT LEVEL (VIS NS)
!ulocal, vlocal and tlocal and plocal are interpolated to (w) level  
  ns = 0
  do iz=istart,iend
     do il=1,ilev
        if (iz.eq.lev(il)) then
           
           ns = ns+1
           i = iz-istart+1
           if (iz.ne.iend) then
              ulocal = 0.5*((u(:,:,iz)-uxym(i))+(u(:,:,iz+1))-uxym(i+1))
              vlocal = 0.5*((v(:,:,iz)-vxym(i))+(v(:,:,iz+1))-vxym(i+1))
              tlocal = 0.5*((t(:,:,iz)-txym(i))+(t(:,:,iz+1))-txym(i+1))
              if (iqscal) qlocal = 0.5*((q(:,:,iz)-qxym(i))+(q(:,:,iz+1)) &
                   -qxym(i+1))
              if (icscal) clocal = 0.5*((c(:,:,iz)-cxym(i))+(c(:,:,iz+1)) &
                   -cxym(i+1))
           else
              ulocal = 0.5*((u(:,:,iz)-uxym(i)) + &
                   (u(:,:,iz+1))-sum(dble(u(:,:,iz+1)))*xym)
              vlocal = 0.5*((v(:,:,iz)-vxym(i)) + &
                   (v(:,:,iz+1))-sum(dble(v(:,:,iz+1)))*xym)
              tlocal = 0.5*((t(:,:,iz)-txym(i)) + &
                   (t(:,:,iz+1))-sum(dble(t(:,:,iz+1)))*xym)
              if (iqscal) qlocal = 0.5*((q(:,:,iz)-qxym(i)) + &
                   (q(:,:,iz+1))-sum(dble(q(:,:,iz+1)))*xym)
              if (icscal) clocal = 0.5*((c(:,:,iz)-cxym(i)) + &
                   (c(:,:,iz+1))-sum(dble(c(:,:,iz+1)))*xym)
           endif
        
           p1 = p(1:nnx,:,iz) - (e(:,:,iz-1)+e(:,:,iz))/3. - &
                0.5*(u(:,:,iz)**2 + v(:,:,iz)**2 + &
                0.5*(w(:,:,iz-1)**2+w(:,:,iz)**2))
           
           if (iz.eq.nnz) p(1:nnx,:,iz+1) = p(1:nnx,:,iz)

           p2 = p(1:nnx,:,iz+1) - (e(:,:,iz)+e(:,:,iz+1))/3. - &
                0.5*(u(:,:,iz+1)**2 + v(:,:,iz+1)**2 + &
                0.5*(w(:,:,iz)**2+w(:,:,iz+1)**2))
     
           p1 = p1 - sum(dble(p1))*xym
           p2 = p2 - sum(dble(p2))*xym
           plocal = 0.5*(p1 + p2)

           ! COMPUTE SPECTRA
           call spectra(ulocal, uspec2d(:,ns), uspec1d(:,ns))
           call spectra(vlocal, vspec2d(:,ns), vspec1d(:,ns))
           call spectra(w(:,:,iz), wspec2d(:,ns), wspec1d(:,ns))
           call spectra(tlocal, tspec2d(:,ns), tspec1d(:,ns))
           call spectra(plocal, pspec2d(:,ns), pspec1d(:,ns))
           if (iqscal) call spectra(qlocal, qspec2d(:,ns), qspec1d(:,ns))
           if (icscal) call spectra(clocal, cspec2d(:,ns), cspec1d(:,ns))

           ! COMPUTE CO- AND QUADRATURE SPECTRA
           call cospec(ulocal, w(:,:,iz), uwcospec2d(:,ns), uwqspec2d(:,ns), &
                uwcospec1d(:,ns), uwqspec1d(:,ns))
           call cospec(vlocal, w(:,:,iz), vwcospec2d(:,ns), vwqspec2d(:,ns), &
                vwcospec1d(:,ns), vwqspec1d(:,ns))
           call cospec(tlocal, w(:,:,iz), twcospec2d(:,ns), twqspec2d(:,ns), &
                twcospec1d(:,ns), twqspec1d(:,ns))
           call cospec(plocal, w(:,:,iz), pwcospec2d(:,ns), pwqspec2d(:,ns), &
                pwcospec1d(:,ns), pwqspec1d(:,ns))
           if (iqscal) call cospec(qlocal, w(:,:,iz), qwcospec2d(:,ns), &
                qwqspec2d(:,ns), qwcospec1d(:,ns), qwqspec1d(:,ns))
           if (icscal) call cospec(clocal, w(:,:,iz), cwcospec2d(:,ns), &
                cwqspec2d(:,ns), cwcospec1d(:,ns), cwqspec1d(:,ns))
           call cospec(p1, (w(:,:,iz)-w(:,:,iz-1))*dzm, pintcospec2d(:,ns), &
                pintqspec2d(:,ns), pintcospec1d(:,ns), pintqspec1d(:,ns))
        endif
     enddo
  enddo

# ifdef MPI
  if (procs.gt.1) then
     call mpi_gather(ns, 1, mpi_integer, nsarray, 1, mpi_integer, 0, &
          comm1d, ier)
     
     nstart(1) = 0
     do i=2,procs
        nstart(i) = sum(nsarray(1:i-1))
     enddo
     
     call mpi_gatherv((spectbuf), spectcount*ncx*ns, mpi_real, &
          spectbuf, spectcount*ncx*nsarray, spectcount*ncx*nstart, &
          mpi_real, 0, comm1d, ier)
  endif
# endif
  
  
  if (mynode.eq.0) then
     open(unit=30, file=spectfile, form="unformatted")
     write(30) ncx, 15
     write(30) xk(1:ncx)
     do iz=1,15
        write(30) eta(iz), uspec2d(:,iz), vspec2d(:,iz), wspec2d(:,iz),  &
             tspec2d(:,iz), qspec2d(:,iz), cspec2d(:,iz), pspec2d(:,iz), &
             uwcospec2d(:,iz), vwcospec2d(:,iz), twcospec2d(:,iz), &
             qwcospec2d(:,iz), cwcospec2d(:,iz), pwcospec2d(:,iz), &
             pintcospec2d(:,iz), uwqspec2d(:,iz), vwqspec2d(:,iz), &
             twqspec2d(:,iz), qwqspec2d(:,iz), cwqspec2d(:,iz), &
             pwqspec2d(:,iz), pintqspec2d(:,iz), &
             uspec1d(:,iz), vspec1d(:,iz), wspec1d(:,iz),  tspec1d(:,iz), &
             qspec1d(:,iz), cspec1d(:,iz), pspec1d(:,iz), &
             uwcospec1d(:,iz), vwcospec1d(:,iz), twcospec1d(:,iz), &
             qwcospec1d(:,iz), cwcospec1d(:,iz), pwcospec1d(:,iz), &
             pintcospec1d(:,iz), uwqspec1d(:,iz), vwqspec1d(:,iz), &
             twqspec1d(:,iz), qwqspec1d(:,iz), cwqspec1d(:,iz), &
             pwqspec1d(:,iz), pintqspec1d(:,iz)
     enddo
     close(30)
  endif
  
end subroutine compspectra







