module storage

  ! GRID PARAMETERS
  integer :: nnx, nny, nnz      ! x, y, and z grid dimensions 
  integer :: ncx, ncy           ! num. unaliased points in x and y axis
  integer :: nxy                ! number of grid points in a horizontal plane
  integer :: nnxb, nnyb         ! size of expanded buffers
  real :: xl, yl, zl         ! x, y, and z domain sizes
  real :: dx, dy, dz         ! x, y, and z grid lengths
  real :: z1                 ! offset of vert. velocity grid, z-dir.
  real :: zi                 ! inversion height
  real :: xym, dzm, todzm    ! 1/(nnx*nny), 1/dz, 0.5*dz
  real :: dsl, xkmax, almin
  real :: dzm2, ceps
  real c1, cs, a, ce, cnlin  ! params in Branko's sgs model
  real pi, xymb
  real dtab1, dtab2, dtabm, zody
  integer izi(1)

  ! FIELD STORAGE ARRAYS:
  real, allocatable ::  u(:,:,:), v(:,:,:), w(:,:,:)
  real, allocatable ::  e(:,:,:), p(:,:,:) 

  ! RIGHT-HAND SIDES OF GOVERNING EQUATIONS (TIME DERIVATIVES)
  real, allocatable ::  ru(:,:,:), rv(:,:,:), rw(:,:,:)  
  real, allocatable ::  re(:,:,:) 

  ! SURFACE LAYER AND LOWER HORIZONTAL PLANE VARIABLES:     
  real :: utau               ! friction velocity
  real :: uwsfc              ! surface flux of long. velocity
  real :: vwsfc              ! surface flux of lateral velocity
  real :: phim               ! stability function for momentum
  

  ! MEAN VERTICAL PROFILE VARIABLES:    
  integer mcount
  real, allocatable, target   :: meanbuf(:,:)
  real, pointer, dimension(:) :: engsbz, englez, divz
  real, pointer, dimension(:) :: uxym, vxym, wxym
  real, pointer, dimension(:) :: uule, vvle, wwle, uwle, vwle, uvle
  real, pointer, dimension(:) :: pple, wple
  real, pointer, dimension(:) :: uusb, vvsb, wwsb, uvsb, uwsb, vwsb

  ! SPECTRUM VARIABLES:    
  integer spectcount
  real, allocatable, target   :: spectbuf(:,:,:)
  real, pointer, dimension(:,:) :: uspec2d, vspec2d, wspec2d

  ! WAVENUMBER SPACE.
  real, allocatable :: wavexy(:,:), xk2(:,:), xk(:), yk(:), xkb(:), ykb(:)

  ! INPUT VARIABLES      
  real :: dt, z0, fcor, dudzf, dvdzf, dtdzf
  real :: ugtop, vgtop, ugbot, vgbot, ugal, vgal, divgls
  real, allocatable :: ug(:), vg(:)
  real, allocatable :: zt(:), zm(:)
      
end module storage






