MODULE force_m
  implicit none
  public :: set_force,force
  character(*),allocatable,save :: calc_force
! this module calculates force
contains
subroutine set_force()
      character(*),intent(IN) :: typ_force
      if (
end subroutine set_force


END MODULE force_m


