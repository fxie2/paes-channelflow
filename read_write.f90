subroutine read_field(iz, qscal_sav, cscal_sav)
  use parameters
  use storage
  implicit none
  logical qscal_sav, cscal_sav
  integer iz

  read(50) u(:,:,iz), v(:,:,iz), w(:,:,iz), e(:,:,iz), &
       p(1:nnx,1:nny,iz)

end subroutine read_field




subroutine write_params(time)
  use parameters
  use storage
  implicit none
  real :: time
  
  write(14) time, nnx, nny, nnz, xl, yl, zl
  write(14) dt, z0, utau
  write(14) divgls, fcor, ugtop, ugbot, vgtop, vgbot
end subroutine write_params





subroutine sav_field(iz)
  use parameters
  use storage
  implicit none
  integer iz
  
  write(14) u(:,:,iz), v(:,:,iz), w(:,:,iz) , e(:,:,iz), &
       p(1:nnx,1:nny,iz)

end subroutine sav_field




