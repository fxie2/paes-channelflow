subroutine default
  use storage
  use parameters
  implicit none

  restart = .false.
  restart_file = "restart.bin" 
  update_grids = .false.
  spectout = .false.
  spectav = .true.
  case    = "les"
  timax   = 100
  savfrq  = 100
  meanfrq = 100
  isbgr = .true.
  nnx = 32
  nny = 32
  nnz = 32
  xl = 5000.0 
  yl = 5000.0 
  zl = 2000.0 
  z0 = 0.16     
  fcor = 0.0001
  dudzf = 0.0
  dvdzf = 0.0
  ugtop = 15.0   
  ugbot = 15.0   
  vgtop = 0.0 
  vgbot = 0.0 
  divgls = 0.0 
      
end subroutine default
