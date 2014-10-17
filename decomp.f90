subroutine decomp1d( n, numprocs, myid, s, e )
  integer n, numprocs, myid, s, e
  integer nlocal
  integer deficit

  nlocal  = n / numprocs
  s	      = myid * nlocal + 1
  deficit = mod(n,numprocs)
  s	      = s + min(myid,deficit)
  if (myid .lt. deficit) then
     nlocal = nlocal + 1
  endif
  e = s + nlocal - 1
  if (e .gt. n .or. myid .eq. numprocs-1) e = n
end subroutine decomp1d
