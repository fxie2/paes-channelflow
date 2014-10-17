module parameters

  implicit none

! MODEL CONSTANTS
  real, parameter :: grav=9.81
  real, parameter :: vk=0.4
  real :: emin
  real :: cb
  real :: skc
      
  logical :: isbgr
  logical :: restart, update_grids, update_tsurf, update_scals
  logical :: spectout, spectav
  integer timax, savfrq
  real meanfrq
  character(100) restart_file, case
      
! MODEL OPTIONS
  logical, parameter :: ICOF=.TRUE.
  logical, parameter :: IPGF=.TRUE.
  logical, parameter :: ISUBS=.FALSE. ! 1 IF LARGE SCALE SUBSIDENCE
  logical, parameter :: ICPMN=.FALSE.
  logical, parameter :: IALIAS=.TRUE.  ! TRUE TO DEALIAS;  
  
!MODELs
! imodel=0: RSFS model
! imodel=1: RANS eddy-viscosity model
! imodel=2: SMM94 model
! imodel=3: M84
! imodel=4: Smag

  INTEGER :: imodel  
end module parameters









