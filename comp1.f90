subroutine comp1(it, mnout, istart, iend, nm, np, istage)

  ! ANALYZES THE ENTIRE RIGHT-HAND-SIDE OF THE STATE EQUATIONS
  ! (THE TIME DERIVATIVES), EXCEPT FOR THE PRESSURE AND RADIATION 
  ! TERMS.  THOSE TERMS ARE ADDED IN BY COMP2, AFTER THE ROUTINE
  ! PRESS HAS COMPLETED CALCULATION OF THE PRESSURE TERM.
  !
  ! ONE VERTICAL INTERFACE IS HANDLED AT A TIME BY COMP1.  THE
  ! ROUTINE MUST BE CALLED ONCE FOR EACH IZ (THE INTERFACE INDEX), 
  ! IN INCREASING ORDER.  THE MAIN INPUT TO COMP1 ARE THE FIELDS
  ! U()..C() AND THEIR TIME DERIVATIVES (RIGHT-HAND-SIDES)
  ! RU()..RC() AT TIME STEP IT.  THESE FIELDS ARE PASSED AS 
  ! COMMON BLOCKS.  WHEN THE ROUTINE HAS FINISHED, RU()..RC()
  ! HAVE BEEN CALCULATED FOR TIME STEP IT+1, WITH THE EXCEPTION
  ! OF THE PRESSURE AND RADIATION TERMS.  U()..C() HAVE ALSO
  ! BEEN UPDATED BY THE "PAST" TERM IN THE ADAMS-BASHFORTH
  ! DIFFERENTIATION SCHEME.  THE "CURRENT" TERM WILL BE
  ! ADDED BY THE COMP2 ROUTINE, ONCE THE CALCULATION OF RU()..RC()
  ! IS COMPLETE.
  !
  ! COMP1 COPIES THE FIELD ARRAYS U()..C() AND THEIR DERIVATIVES
  ! INTO TEMPORARY "BUFFERS."  EACH BUFFER CONTAINS THE FIELD ARRAYS 
  ! AT A SMALL NUMBER OF VERTICAL INTERFACES, CENTERED ON THE CURRENT 
  ! LEVEL IZ.  TO SAVE COMPUTATION TIME, EACH BUFFER IS STORED AS A 
  ! COMMON BLOCK AND PASSED TO COMP1 WHEN CALLED FOR IZ+1.  THEN
  ! ONLY THE NEW TOP INTERFACE IN THE BUFFER IS COPIED, INTO THE
  ! LOCATION OCCUPIED BY THE PREVIOUS BOTTOM INTERFACE.
  !
  ! GENERALLY THE BUFFERS SHOULD HAVE HORIZONTAL DIMENSIONS 3/2
  ! TIMES THE INPUT FIELDS.  SINCE THE CALCULATIONS IN COMP1 ARE
  ! NONLINEAR, HIGHER-ORDER SPATIAL HARMONICS ARE GENERATED.
  ! THESE MUST BE FILTERED OUT AT THE END OF THE COMPUTATION
  ! IN ORDER TO AVOID ALIASING ON THE NEXT ITERATION. THE FILTERING
  ! IS PERFORMED BY CUTTING OUT THE TOP THIRD OF THE WAVENUMBER
  ! SPECTRUM.  BECAUSE THE ENSUING COMPUTATIONS IN PRESS AND COMP2
  ! ARE LINEAR, THE DEALIASING IS DONE AT THE END OF COMP1.
  !
  ! INPUT:
  ! U()..C()          FIELDS AT TIME STEP IT
  ! RU()..RC()        TIME DERIVATIVES AT TIME STEP IT
  !
  ! OUTPUT:
  ! U()..C()          FIELDS AT TIME STEP IT+1/2
  ! RU()..RC()        TIME DERIVATIVES AT TIME STEP IT+1/2

  use storage
  use compbuf
  use parameters
  implicit none

  integer :: ier
  integer :: iz                ! vertical plane (z) index
  integer :: j                 ! buffer index, current layer
  integer :: jp1               ! buffer index, layer above
  integer :: k                 ! buffer index, current layer
  integer :: km1               ! buffer index, layer below
  integer :: it                ! time step index
  logical :: mnout             ! true if means are calculated
  integer :: l, lp1, lm1
  integer :: istart, iend, kb=0, i, idif
  integer :: mm2, mm1, m, mp1, mp2
  integer :: nm1, n, np1, np2
  integer :: nm(4), np(4), istage
  real    :: fnct, tmean
  
  ! MAIN LOOP OVER EACH HORIZ. LEVEL
  do iz=istart,iend

     ! SET THE LAYER BUFFER POINTERS.
     j=mod(iz,2)+1
     jp1=mod(iz+1,2)+1

     k=mod(iz,2)+1
     km1=mod(iz+1,2)+1

     lm1=mod(iz+2,3)+1
     l=mod(iz,3)+1
     lp1=mod(iz+1,3)+1
     
     nm1=mod(iz+3,4)+1
     n=mod(iz,4)+1
     np1=mod(iz+1,4)+1
     np2=mod(iz+2,4)+1
     
     mm2=mod(iz+4,5)+1
     mm1=mod(iz,5)+1
     m=mod(iz+1,5)+1
     mp1=mod(iz+2,5)+1
     mp2=mod(iz+3,5)+1
      
     if (iz.eq.istart) then
        
        if (nny.eq.1) vbuf = 0.0

        do i=-nm(1),np(1)-1
           call expand(u(1,1,i+istart),ubuf(1,1,i))
        enddo
        if (nny .ne. 1) then
           do i=-nm(2),np(2)-1
              call expand(v(1,1,i+istart),vbuf(1,1,i))
           enddo
        endif
        do i=-nm(3),np(3)-1
           call expand(w(1,1,i+istart),wbuf(1,1,i))
        enddo
        if (isbgr) then
           do i=-nm(4),np(4)-1
              call expand(e(1,1,i+istart),ebuf(1,1,i))
              ebuf(:,:,i) = max(ebuf(:,:,i),emin)
           enddo
        endif
     endif
     
!     write(*,*) 'Proc ', istart, ' Finished expanding 1'

     call expand(u(1,1,iz+np(1)),ubuf(1,1,np(1)))
     if (nny .ne. 1) call expand(v(1,1,iz+np(2)),vbuf(1,1,np(2)))
     call expand(w(1,1,iz+np(3)),wbuf(1,1,np(3)))
     if (isbgr)  call expand(e(1,1,iz+np(4)),ebuf(1,1,np(4)))

     ebuf(:,:,np(4)) = max(ebuf(:,:,np(4)),emin)
!     write(*,*) 'Proc ', istart, ' Finished expanding 2'
     
     ! MODIFY PROGRAM VARIABLES BY -1/2*RHS(N-1).
     ! (FIRST STEP IN ADAMS-BASHFORTH TIME DIFFERENTIATION.)
     if (istage.ne.1) then
        u(:,:,iz) = u(:,:,iz) - dtab2*ru(:,:,iz)
        if (nny .ne. 1) v(:,:,iz) = v(:,:,iz) - dtab2*rv(:,:,iz)
        if (iz .ne. nnz) then
           w(:,:,iz) = w(:,:,iz) - dtab2*rw(:,:,iz)
           if (isbgr)  e(:,:,iz) = e(:,:,iz) - dtab2*re(:,:,iz)
        endif
     endif
!     write(*,*) 'Proc ', istart, ' Finished Adams Bashforth'
     ! CALCULATE NEEDED DERIVATIVES 
     if (iz.eq.istart) call derivsbot(kb,iz,lm1,l,nm1,n,np1)
     if(iz.ne.nnz) then
        call derivs(kb,iz,istart,l,lp1,np2)
     end if
     

 !    write(*,*) 'Proc ', istart, ' Finished derivs'
     ! CALCULATE DIVERGENCE
!     if (mnout .and. istage.eq.1) divz(iz-istart+1) = sum(dble( &
     divz(iz-istart+1) = sum(dble((ux(:,:,n)+vy(:,:,n)+wz(:,:,n))**2))*xymb
     
     ! DETERMINE THE ADVECTIVE TERMS
     call compadv1(kb,iz,lm1,l,n,np1,istart)
!     write(*,*) 'Proc ', istart, ' Finished CompAdv1'

     !Imposed pressure gradient
     rubuf(:,:) = rubuf(:,:) + 0.016
     
     ! SUB-GRID EDDY TERMS (INCLUDING SURFACE FORCING)
     call compsgs1(kb,iz,j,jp1,k,km1,l,lm1,lp1,nm1,n,np1,np2,istart)
!     write(*,*) 'Level ', iz, ' Finished compsgs1'
     call compsgs3(kb,iz,j,jp1,k,km1,l,lm1,lp1,n,np1,mnout,istage,istart)
!     write(*,*) 'Level ', iz, ' Finished compsgs3'

     ! DONE WITH NONLINEAR CALCULATIONS; DEALIAS
     if (nny .eq. 1) then
        rv(:,:,iz) = 0.0
     else
        call compress(rvbuf,rv(1,1,iz))
     endif
     call compress(rubuf,ru(1,1,iz))
     if (iz.ne.nnz) call compress(rwbuf,rw(1,1,iz))
     if (isbgr .and. iz.ne.nnz)  call compress(rebuf,re(1,1,iz))
 
!     write(*,*) 'Level ', iz, ' Finished Dealias'
     ! SUBTRACT OUT THE HORIZONTAL MEAN PART TO ASSURE WMEAN=0.
     if (iz.ne.nnz) rw(:,:,iz) = rw(:,:,iz) - sum(dble(rw(:,:,iz)))*xym

     if (iz.eq.nnz) then
        rw(:,:,iz) = 0.0
        re(:,:,iz) = 0.0
     endif

     if (iz.ne.iend) then
        ubuf(:,:,-nm(1):np(1)-1) = ubuf(:,:,-nm(1)+1:np(1))
        if (nny .ne. 1) vbuf(:,:,-nm(2):np(2)-1) = vbuf(:,:,-nm(2)+1:np(2))
        wbuf(:,:,-nm(3):np(3)-1) = wbuf(:,:,-nm(3)+1:np(3))
         if (isbgr)  ebuf(:,:,-nm(4):np(4)-1) = ebuf(:,:,-nm(4)+1:np(4))
      endif

  enddo

end subroutine comp1











