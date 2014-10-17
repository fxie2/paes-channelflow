program expand_les
  
  implicit none
  character(25) :: fieldin, fieldout
  integer :: nxnew, nynew, nznew
  namelist /GRIDSIZE/ nxnew, nynew, nznew
  namelist /FILES/ fieldin, fieldout

  logical :: qscal, cscal, exists
  integer :: itsfc, iqsfc, icsfc
  real :: dt, z0, wtsfc, tsfcc, amonin, utau
  real :: wqsfc, qsfcc, wcsfc, csfcc
  real :: dtdzf, divgls, fcor, ugtop, ugbot, vgtop, vgbot
  real :: xl, yl, zl, time
  
  real, allocatable, dimension(:,:) :: ua, ub, va, vb, wa, wb, wm, pa, pb
  real, allocatable, dimension(:,:) :: ta, tb, ea, eb, em, qa, qb, ca, cb
  real, allocatable, dimension(:,:) :: uhold, vhold, whold, thold, ehold
  real, allocatable, dimension(:,:) :: qhold, chold, phold
  real, allocatable, dimension(:,:) :: unew, vnew, wnew, tnew, enew
  real, allocatable, dimension(:,:) :: qnew, cnew, pnew
  
  real zfold, dzold, zf, dznew, zh, zhold
  integer izo, izn, nxold, nyold, nzold
  
  inquire(file="redim.input", exist=exists)
  if (.not. exists) then
     write(*,*) "Stopping; can't find file redim.input"
     stop
  endif
  
  open (40, file="redim.input", action="read")
  write(*,*) "Reading namelist file"
  read(40,nml=GRIDSIZE)
  read(40,nml=FILES)
  close(40)
  
  inquire(file=fieldin, exist=exists)
  if (exists) then
     
     write(*,*) "Reading parameters from coarse field file"
     open(unit=50, file=fieldin, form="unformatted")
     read(50) time, nxold, nyold, nzold, xl, yl, zl, qscal, cscal
     read(50) dt, z0, itsfc, wtsfc, tsfcc, amonin, utau
     read(50) iqsfc, wqsfc, qsfcc, icsfc, wcsfc, csfcc
     read(50) dtdzf, divgls, fcor, ugtop, ugbot, vgtop, vgbot
     
     open(unit=60, file=fieldout, form="unformatted")
     write(60) time, nxnew, nynew, nznew, xl, yl, zl, qscal, cscal
     write(60) dt, z0, itsfc, wtsfc, tsfcc, amonin, utau
     write(60) iqsfc, wqsfc, qsfcc, icsfc, wcsfc, csfcc
     write(60) dtdzf, divgls, fcor, ugtop, ugbot, vgtop, vgbot

     write(*,*) "allocating space for field data"
     call allocate_arrays
     
     write(*,*) "interpolating field data"
     call interpol
     
     write(*,*) "deallocating field data"
     call deallocate_arrays
     
     close(50)
     close(60)

  endif
  
contains

  subroutine read_fieldin
    implicit none
    
    if (qscal .and. cscal) then
       read(50) ub, vb, wb, tb, eb, pb, qb, cb
    else if (qscal) then
       read(50) ub, vb, wb, tb, eb, pb, qb
    else if (cscal) then
       read(50) ub, vb, wb, tb, eb, pb, cb
    else
       read(50) ub, vb, wb, tb, eb, pb
    endif
  end subroutine read_fieldin



  subroutine write_fieldout
    implicit none
    
    if (qscal .and. cscal) then
       write(60) unew, vnew, wnew, tnew, enew, pnew, qnew, cnew
    else if (qscal) then
       write(60) unew, vnew, wnew, tnew, enew, pnew, qnew
    else if (cscal) then
       write(60) unew, vnew, wnew, tnew, enew, pnew, cnew
    else
       write(60) unew, vnew, wnew, tnew, enew, pnew
    endif

  end subroutine write_fieldout


  subroutine allocate_arrays
    implicit none
    
    allocate(ua(nxold,nyold))
    allocate(ub(nxold,nyold))
    allocate(va(nxold,nyold))
    allocate(vb(nxold,nyold))
    allocate(wm(nxold,nyold))
    allocate(wa(nxold,nyold))
    allocate(wb(nxold,nyold))
    allocate(ta(nxold,nyold))
    allocate(tb(nxold,nyold))
    allocate(em(nxold,nyold))
    allocate(ea(nxold,nyold))
    allocate(eb(nxold,nyold))

    allocate(uhold(nxold,nyold))
    allocate(vhold(nxold,nyold))
    allocate(whold(nxold,nyold))
    allocate(thold(nxold,nyold))
    allocate(ehold(nxold,nyold))

    allocate(unew(nxnew,nynew))
    allocate(vnew(nxnew,nynew))
    allocate(wnew(nxnew,nynew))
    allocate(tnew(nxnew,nynew))
    allocate(enew(nxnew,nynew))
    
    allocate(pa(nxold,nyold))
    allocate(pb(nxold,nyold))
    allocate(phold(nxold,nyold))
    allocate(pnew(nxnew,nynew))
    
    if (qscal) then
       allocate(qa(nxold,nyold))
       allocate(qb(nxold,nyold))
       allocate(qhold(nxold,nyold))
       allocate(qnew(nxnew,nynew))
    endif

    if (cscal) then
       allocate(ca(nxold,nyold))
       allocate(cb(nxold,nyold))
       allocate(chold(nxold,nyold))
       allocate(cnew(nxnew,nynew))
    endif
  end subroutine allocate_arrays



  subroutine deallocate_arrays
    implicit none
    
    deallocate(ua)
    deallocate(ub)
    deallocate(va)
    deallocate(vb)
    deallocate(wm)
    deallocate(wa)
    deallocate(wb)
    deallocate(ta)
    deallocate(tb)
    deallocate(em)
    deallocate(ea)
    deallocate(eb)

    deallocate(uhold)
    deallocate(vhold)
    deallocate(whold)
    deallocate(thold)
    deallocate(ehold)

    deallocate(unew)
    deallocate(vnew)
    deallocate(wnew)
    deallocate(tnew)
    deallocate(enew)
    
    deallocate(pa)
    deallocate(pb)
    deallocate(phold)
    deallocate(pnew)
    
    if (qscal) then
       deallocate(qa)
       deallocate(qb)
       deallocate(qhold)
       deallocate(qnew)
    endif

    if (cscal) then
       deallocate(ca)
       deallocate(cb)
       deallocate(chold)
       deallocate(cnew)
    endif
  end subroutine deallocate_arrays
  
  
  
  subroutine interpol
    
    dzold = zl/float(nzold)
    dznew = zl/float(nznew)

    call read_fieldin
    ua=ub ; va=vb ; wa=wb ; ta=tb ; ea=eb ; pa=pb
    if (qscal) qa=qb
    if (cscal) ca=cb
    call read_fieldin

    wm = 0.0
    em = ea
    izo = 2
    zhold = float(izo)*dzold - 0.5*dzold
    zfold = float(izo-1)*dzold

    do izn=1,nznew

       if (mod(izn,10).eq.0) write(*,*) "Level", izn
       
       zh = float(izn)*dznew-0.5*dznew
       zf = float(izn)*dznew

       do while (zhold .lt. zh .and. izo .lt. nzold)
          wm=wa ; em=ea
          ua=ub ; va=vb ; wa=wb ; ta=tb ; ea=eb ; pa=pb
          if (qscal) qa=qb
          if (cscal) ca=cb
          izo = izo+1
          call read_fieldin
          zhold = float(izo)*dzold-0.5*dzold
          zfold = float(izo-1)*dzold
       enddo

       uhold = ub + (ub-ua)*(zh-zhold)/dzold
       vhold = vb + (vb-va)*(zh-zhold)/dzold
       whold = wb + (wb-wa)*(zh-zhold)/dzold
       thold = tb + (tb-ta)*(zh-zhold)/dzold
       phold = pb + (pb-pa)*(zh-zhold)/dzold
       if (qscal) qhold = qb + (qb-qa)*(zh-zhold)/dzold
       if (cscal) chold = cb + (cb-ca)*(zh-zhold)/dzold

       if (zf.le.zfold) then
          whold = wa + (wa-wm)*(zf-zfold)/dzold
          ehold = ea + (ea-em)*(zf-zfold)/dzold
       else
          whold = wb + (wb-wa)*(zf-(zfold+dzold))/dzold
          ehold = eb + (eb-ea)*(zf-(zfold+dzold))/dzold
       endif

       call redim(uhold,unew,nxold,nyold,nxnew,nynew)
       call redim(vhold,vnew,nxold,nyold,nxnew,nynew)
       call redim(whold,wnew,nxold,nyold,nxnew,nynew)
       call redim(thold,tnew,nxold,nyold,nxnew,nynew)
       call redim(ehold,enew,nxold,nyold,nxnew,nynew)
       call redim(phold,pnew,nxold,nyold,nxnew,nynew)
       if (qscal) call redim(qhold,qnew,nxold,nyold,nxnew,nynew)
       if (cscal) call redim(chold,cnew,nxold,nyold,nxnew,nynew)

       enew = max(enew,1.e-10)     
       enew = min(maxval(ehold),enew)

       call write_fieldout

    enddo

  end subroutine interpol

end program expand_les





subroutine redim(xin,xout,nxin,nyin,nxout,nyout)

  implicit none
  
  integer :: nxin, nyin
  integer :: nxout, nyout
  integer :: nxm, nym
  
  real :: xin(nxin,nyin)
  real :: xout(nxout,nyout)

  real :: wavein(nxin+2,nyin)
  real :: waveout(nxout+2,nyout)

  if (nxin.eq.nxout .and. nyin.eq.nyout) then
     xout = xin
     return
  end if

  call fft_init(nxin, nyin)
  wavein(1:nxin,:) = xin
  call ffthrz(wavein, nxin, nyin, -1)
	
  wavein(nxin+1:nxin+2,:) = 0.0
  wavein(1:nxin,nyin/2+1) = 0.0
  
  nxm = min(nxin,nxout)
  nym = min(nyin,nyout)
  waveout = 0.0
  
  waveout(1:nxm+2,1:nym/2+1) = wavein(1:nxm+2,1:nym/2+1)
  waveout(1:nxm+2,nyout-nym/2+2:nyout) = wavein(1:nxm+2,nyin-nym/2+2:nyin)
  
  call fft_finalize()
  
  call fft_init(nxout, nyout)
  call ffthrz(waveout, nxout, nyout, 1)
  xout = waveout(1:nxout,:)
  
  call fft_finalize()
  
end subroutine redim








subroutine lesabort(message)
  implicit none
      
  integer ier
  character*(*) :: message
  write(*,*) message
      
  stop
end subroutine lesabort
      

