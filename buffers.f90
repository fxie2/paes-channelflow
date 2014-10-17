subroutine setcombuf(np, nptot, nm, nmtot)
  
  use parameters
  implicit none
  integer i
  
  ! NUMBER LEVELS NEEDED ABOVE THE CURRENT LEVEL FOR EACH VARIABLE
  integer :: np(4), nptot(4)
  
  ! NUMBER LEVELS NEEDED BELOW THE CURRENT LEVEL FOR EACH VARIABLE
  integer :: nm(4), nmtot(4)
  
  ! NUMBER OF LEVELS NEEDED ABOVE CURRENT LEVEL IN BASE MODEL
  np(1)=1; np(2)=1; np(3)=1;
  np(4)=merge(1,0,isbgr)
  
  ! NUMBER OF LEVELS NEEDED BELOW CURRENT LEVEL IN BASE MODEL
  nm(1)=1; nm(2)=1; nm(3)=1; 
  nm(4)=merge(1,0,isbgr)
  
  ! EXTRA LEVELS NEEDED IF USING SUBGRID MODEL
  if (isbgr) then
     np(1)=max(np(1),2)
     np(2)=max(np(2),2)
     np(3)=max(np(3),2)
     nm(3)=max(nm(3),2)
  endif
  
  do i=1,4
     nptot(i) = sum(np(1:i))
     nmtot(i) = sum(nm(1:i))
  enddo
  
end subroutine setcombuf





