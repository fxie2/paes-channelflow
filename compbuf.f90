module compbuf
  real, allocatable :: ubuf(:,:,:), vbuf(:,:,:), wbuf(:,:,:)
  real, allocatable :: ebuf(:,:,:)
  real, allocatable, target :: rubuf(:,:), rvbuf(:,:), rwbuf(:,:)  
  real, allocatable, target :: rebuf(:,:) 

  real, allocatable :: UX(:,:,:)  ! X-DERIV. OF U
  real, allocatable :: UY(:,:,:)  ! Y-DERIV. OF U
  real, allocatable :: UZ(:,:,:)  ! Z-DERIV. OF U
  real, allocatable :: VX(:,:,:)  ! X-DERIV. OF V
  real, allocatable :: VY(:,:,:)  ! Y-DERIV. OF V
  real, allocatable :: VZ(:,:,:)  ! Z-DERIV. OF V
  real, allocatable :: WX(:,:,:)  ! X-DERIV. OF W
  real, allocatable :: WY(:,:,:)  ! Y-DERIV. OF W
  real, allocatable :: WZ(:,:,:)  ! Z-DERIV. OF W
  real, allocatable :: EX(:,:)    ! X-DERIV. OF E
  real, allocatable :: EY(:,:)    ! Y-DERIV. OF E
  real, allocatable :: tauuu(:,:)
  real, allocatable :: tauuv(:,:)
  real, allocatable :: tauvv(:,:)
  real, allocatable :: tauuw(:,:,:)
  real, allocatable :: tauvw(:,:,:)
  real, allocatable :: tauww(:,:,:)
  REAL, allocatable :: XKK(:,:,:)    ! SUBGRID EDDY COEF.
  REAL, allocatable :: ALK(:,:,:)    ! SUBGRID LENGTH SCALE
  real, allocatable :: shear(:,:,:)
      
  real, pointer :: ap(:,:,:,:), ttd(:,:,:,:), cplus(:,:,:)
  real, pointer :: wxjp2(:,:), wyjp2(:,:)
  real, pointer :: uzjp2(:,:), vzjp2(:,:)
      
end module compbuf
