program les
  use storage
  use parameters
  implicit none
  
#ifdef MPI  
  include "mpif.h"
  integer :: belownode, abovenode, j
#endif
  
  integer :: it, time_end, iz, istage, ier
  integer :: time_array0(8), time_array1(8), nscal, nvar, i
  real :: wallclock0, wallclock1, seconds=0.0
  logical :: mnout, fdout, tnout
  integer :: meancount, fieldcount
  logical :: update_sparams, initq, initc
  logical :: exists, qscal_sav, cscal_sav
  integer :: nxt, nyt, nzt
  real :: xlt, ylt, zlt
  real :: time, time_start
  integer :: np(4), nptot(4) ! NUMBER LEVELS NEEDED ABOVE THE CURRENT LEVEL FOR EACH VARIABLE
  integer :: nm(4), nmtot(4) ! NUMBER LEVELS NEEDED BELOW THE CURRENT LEVEL FOR EACH VARIABLE
  integer :: istart, iend, procs, mynode, comm1d=0
  integer, pointer :: pstart(:), pend(:), psize(:)
  character(25) :: filename, fname, meanfile, spectfile, wtfile
  real :: zetas(3) = (/ 0.0, 17.0/60.0, 5.0/12.0 /)
  real :: gama(3) = (/ 8./15., 5.0/12.0, 3.0/4.0 /)
  real :: ucfl, vcfl, wcfl, cfl
      
  namelist /RESTARTING/ restart, restart_file
  namelist /TIMING/ case, timax, savfrq, meanfrq, isbgr, cfl, dt,&
       cb, skc, spectout, spectav, imodel, isurface ! added by Paulo 10/17/14
  namelist /GRIDSIZE/ nnx, nny, nnz
  namelist /PARAMS/ xl, yl, zl, z0
      
  mynode = 0
  procs = 1
!  cfl = 0.75

#ifdef MPI
  call mpi_init(ier) 
  call mpi_comm_size(MPI_COMM_WORLD, procs, ier)
  call mpi_cart_create(MPI_COMM_WORLD, 1, procs, .false., .true., comm1d, ier)
  call mpi_comm_rank(comm1d, mynode, ier)
  call mpi_cart_shift(comm1d, 0,  1, belownode, abovenode, ier)
#endif
  
  if (mynode.eq.0) then
     write(*,*)
     write(*,*) "NCAR/PENN STATE LARGE-EDDY SIMULATION"
     write(*,*)
     call date_and_time(values=time_array0)
     wallclock0 = time_array0(5)*3600. + time_array0(6)*60. + &
          time_array0(7) + 0.001*time_array0(8)
  endif
  
  ! SET DEFAULT NAMELIST VALUES
  call default()
  
  ! MAKE SURE NAMELIST FILE EXISTS
  inquire(file="nlist.input", exist=exists)
  if (.not. exists) call lesabort("Stopping; can't find nlist.input")
  
  ! READ NAMELIST INFORMATION
  open (10,file="nlist.input", action="read")
  read(10,nml=RESTARTING)
  read(10,nml=TIMING)

  if (restart) then

     ! MAKE SURE RESTART FILE EXISTS
     inquire(file=restart_file, exist=exists)
     if (.not. exists) then
        call lesabort("Stopping; can't find restart file " &
             //trim(restart_file))
     endif

     if (mynode.eq.0) write(*,*) "Reading parameters from restart file"
     open(unit=50, file=restart_file, form="unformatted")
     read(50) time_start, nnx, nny, nnz, xl, yl, zl
     read(50) dt, z0, utau
     read(50) divgls, fcor, ugtop, ugbot, vgtop, vgbot     
  else
     if (mynode.eq.0) write(*,*) "Reading parameters from namelist file"
     read(10,nml=GRIDSIZE)
     read(10,nml=PARAMS)
     time_start = 0.0
  endif
     
  close(10)
  
  pi=4.0*atan(1.0)
  c1=(sqrt(960.)*cb)/(7.*(1+cb)*skc)
  cs=sqrt(8.*(1+cb)/(27.*pi*pi))
  a=(27./8./pi)**(1./3.)
  ce=(8.*pi/27.)**(1./3.)*cs**(4./3.)
  cnlin=cs*cs*c1
  
  nnxb = 3*nnx/2
  nnyb = 3*nny/2
  xymb = 1.0/float(nnxb)/float(nnyb)
  emin = 1.e-6
  ncx=nnx/2+1
  ncy=nny/2+1
  nxy=nnx*nny
  xym=1./float(nnx*nny)
  dx=xl/float(nnx)
  dy=yl/float(nny)
  dz=zl/float(nnz)
  z1=dz/2.
  zody=log(z1/z0)
  dzm=1./dz
  dzm2=dzm**2
  todzm=0.5/dz
  if (nny .eq.1) then
     dsl = sqrt(dx*dz)
  else
     dsl=(dx*dy*dz)**(1./3.)
  endif
  ceps=0.93/dsl
  almin=0.0075*dsl
  
  istart = 1
  iend = nnz


#ifdef MPI  
  allocate(pstart(0:procs-1))
  allocate(pend(0:procs-1))
  allocate(psize(0:procs-1))
  if (mynode.eq.0) then
     do j=0,procs-1
        call decomp1d(nnz, procs, j, pstart(j), pend(j))
        write(*,*) "Node ",j, " Has Levels ",pstart(j)," Through ",pend(j)
     enddo
     istart = pstart(0)
     iend = pend(0)
     psize = pend-pstart+1
  else
     call decomp1d(nnz, procs, mynode, istart, iend)
     pstart = 0
     psize = 0
  endif
#endif

  call setcombuf(np, nptot, nm, nmtot)
  call allocate_arrays(istart,iend, np, nm)
  call initwave(istart,iend)
  call fft_init(nnx,nny)
 
  nscal = 0
  if (isbgr) nscal = nscal+1

  nvar = 4

  
  ! READ 3-D SAVED FIELDS FROM FILE IF RESTART
  if (restart) then
     if (procs.gt.1) then
#ifdef MPI
        call read_restart_mpi(istart,iend,procs,mynode,comm1d,nvar, &
             qscal_sav,cscal_sav)
#endif
     else
        do iz=istart,iend
           call read_field(iz, qscal_sav, cscal_sav)
        enddo
     endif
     close(50)
  endif
  
  ! INITIALIZE FIELDS THAT AREN'T FROM RESTART
  if (.not. restart) call initvel(istart,iend)

  ! HEIGHT ARRAYS DEFINED ONLY ON FIRST PROCESSOR
  if (mynode.eq.0) then
     zm(1) = dz/2
     zt(1) = dz
     do iz=2,nnz
        zm(iz) = zm(iz-1)+dz
        zt(iz) = zt(iz-1)+dz
     enddo
  endif
      
  do iz=istart,iend
     call deal(u(:,:,iz),nnx,nny)
     call deal(v(:,:,iz),nnx,nny)
     call deal(w(:,:,iz),nnx,nny)
     call deal(e(:,:,iz),nnx,nny)
  enddo
	

#ifdef MPI
     if (procs.gt.1) call commun_arrays(istart, iend, nscal, abovenode, &
          belownode, mynode, comm1d, np, nptot, nm, nmtot)
#endif

  it=0
  time = 0.0
  meancount=1
  fieldcount=1

  do while(time .lt. float(timax))
!  do while(time .lt. 0.05)
!    ucfl = maxval(abs(u(:,:,istart:iend)))*float(nnxb)/xl
!    vcfl = maxval(abs(v(:,:,istart:iend)))*float(nnyb)/yl
!    wcfl = maxval(abs(w(:,:,istart:iend)))/dz
!    dt = cfl/max(ucfl,vcfl,wcfl,0.5/dz)
!     dt = 1.0
    
! #    ifdef MPI     
!      if (procs.gt.1) then
!         call mpi_allreduce((dt),dt,1,mpi_integer,mpi_min,comm1d,ier)
!      endif
! #    endif
     xkmax=dz*dz/dt/5.
     time = time+dt
     if (mynode.eq.0) then
        call date_and_time(values=time_array1)
        wallclock1 = time_array1(5)*3600 + time_array1(6)*60 &
             + time_array1(7) + 0.001*time_array1(8)

        if (wallclock1 .ge. wallclock0) then
           seconds = seconds + (wallclock1 - wallclock0)
        else
           ! DAY MUST HAVE CHANGED
           seconds = seconds + (86400.0-wallclock0) + wallclock1
        endif
        
        wallclock0 = wallclock1

        write(*,*) "Iteration ", it+1, "; model time = ", time+time_start, &
             "; delta t = ", dt
        write(*,*) "Wall Clock Time = ", seconds, "seconds"
        write(*,*)
     endif
     
     if (time .ge. float(meancount)*meanfrq) then
        meancount=meancount+1
        mnout=.true.
     else if (it.eq.0) then
        mnout = .true.
     else
        mnout=.false.
     endif

     if (time .ge. float(fieldcount*savfrq)) then
        fieldcount=fieldcount+1
        fdout=.true.
     else
        fdout=.false.
     endif
     
     do istage=1,3
        dtab1 = dt*gama(istage)
        dtab2 = dt*zetas(istage)
        dtabm = 1.0/dtab1
        
        ! SET LOWER BOUNDARY VALUES
        if (istart.eq.1) call lower(istart,np,nm)
     
        ! SET UPPER BOUNDARY VALUES
        if (iend.eq.nnz) call upper(np)
        ! COMPUTE TENDENCIES
        call comp1(it, mnout, istart, iend, nm, np, istage)
        ! COMPUTE PRESSURE FIELD
#       ifdef MPI     
        if (procs.gt.1) call commun_rw(istart, iend, nscal, abovenode, &
             belownode, mynode, comm1d)
#       endif
        call comprhs(istart, iend)
        call press(istart,iend,procs,mynode,comm1d)
#       ifdef MPI     
        if (procs.gt.1) call commun_press(istart, iend, nscal, abovenode, &
             belownode, mynode, comm1d)
#       endif     

        ! ADVANCE FIELDS TO NEXT TIME STEP
        call comp2(istart,iend)

#       ifdef MPI
        if (procs.gt.1) call commun_arrays(istart, iend, nscal, abovenode, &
             belownode, mynode, comm1d, np, nptot, nm, nmtot)
#       endif
     enddo
        
     if (mnout .or. it.eq.0) then
        call compmn(istart,iend)
#ifdef MPI
        call mpi_gatherv((meanbuf), mcount*(iend-istart+1), mpi_real,& 
             meanbuf, psize*mcount, (pstart-1)*mcount, mpi_real,&
             0, comm1d,ier)
#endif

        if (mynode.eq.0) then
           write(filename,'(F9.2)') time + time_start !(meancount-1)*meanfrq+int(time_start)
           meanfile=trim(case)//trim(adjustl(filename))//".mean"

           write(*,'(a4,3a9,3a11)') "IZ", "Uxym ", "Vxym ", "pple ", &
			"engle", "divz  ", "engsbz "
           
           do iz=1,nnz
              write(*,'(I4,3f9.3,3e11.3)') iz, uxym(iz), vxym(iz), &
               	pple(iz), englez(iz), divz(iz), engsbz(iz)
           enddo
           
	   
	   CALL output_profiles(mynode, time, time_start, it)  !Tie
	   
	   
           if (mnout .and. (meancount.ge.1)) then
              open(unit=35, file=meanfile) !, form="unformatted")
              write(35,*) nnz, time+time_start
              write(35,*) utau, z0
              write(35,*) zm, zt
              write(35,*) englez, engsbz
              write(35,*) uxym+ugal, vxym+vgal, wxym
              write(35,*) uule, vvle, wwle, uwle, vwle, uvle
              write(35,*) pple
              write(35,*) uusb, vvsb, wwsb, uvsb, uwsb, vwsb
              close(35)
           endif
        endif
        
        ! spectfile = trim(case)//trim(adjustl(filename))//".spectra"
        ! if (mnout.and.spectout.and.(meancount.ge.1)) then
        !    call compspectra(istart,iend,comm1d,mynode,procs, spectfile)
        ! endif
     endif

     if (fdout) then
        write(filename,'(I10)') (fieldcount-1)*savfrq+int(time_start)
        filename=trim(case)//trim(adjustl(filename))//".field"
        if (procs.gt.1) then
#ifdef MPI
        if (mynode .eq. 0) then
           open (unit=14, file=filename,form="unformatted")
           call write_params(time+time_start)
           close(14)
        endif
        call write_field_mpi(procs,mynode,filename,istart,iend,comm1d)
#endif
        else
        open (unit=14, file=filename, form="unformatted")
        call write_params(time+time_start)
        do iz=istart,iend
           call sav_field(iz)
        enddo
        close(14)
        endif
     endif

     it=it+1
  enddo
      
  call fft_finalize()
  call deallocate_arrays()
      
  
#ifdef MPI
  call mpi_barrier(comm1d, ier)
#endif
  
  if (mynode.eq.0) then
     call date_and_time(values=time_array1)
     wallclock1 = time_array1(5)*3600. + time_array1(6)*60. &
          + time_array1(7) + 0.001*time_array1(8)

     if (wallclock1 .ge. wallclock0) then
        seconds = seconds + (wallclock1 - wallclock0)
     else
        ! DAY MUST HAVE CHANGED
        seconds = seconds + (86400.0-wallclock0) + wallclock1
     endif
        
     write(*,*) "Job took", seconds, "seconds to execute"
     write(*,*)
  endif

#ifdef MPI
  call mpi_finalize(comm1d, ier)
#endif
  
end program les




 SUBROUTINE output_profiles(mynode, time, time_start, it)
  USE storage
  USE parameters 
  IMPLICIT NONE
  REAL, INTENT(IN) :: time, time_start
  INTEGER, INTENT(IN) :: it	!iteration step
  INTEGER :: iz, mynode
  REAL :: wstar

	IF (mynode .EQ. 0) THEN
	OPEN(500, FILE='parameters.dat', POSITION='append')
	   IF (it .eq. 0) THEN
	   WRITE(500,'(2a12)') '#Time', 'utau'
	   END IF
	   WRITE(500,'(2f12.3)') time+time_start, utau
	CLOSE(500)  
        ENDIF
	
 END SUBROUTINE output_profiles
  







