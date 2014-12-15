MODULE vector
  implicit none
contains
  real(8) function norm(q)
  real(8),intent(IN) :: q(:)
  integer	     :: i

  norm=0.0D0

  DO i=1,size(q)
   norm=norm+q(i)*q(i)
  END DO

  norm=dsqrt(norm)

  end function norm

  real function inpro(q)
  real(8),intent(IN) :: q(:)
  integer	     :: i

  inpro=0.0D0

  DO i=1,size(q)
   inpro=inpro+q(i)*q(i)
  END DO
  end function inpro

