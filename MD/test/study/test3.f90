MODULE calc
implicit none
private
public :: push,writeAAA,dot
real(8),save :: AAA 
contains
subroutine push(a)
real(8) :: a
 AAA=a
end subroutine push

subroutine dot(A,B,C)
real(8),intent(IN) :: A(:),B(:)
real(8),intent(OUT):: C
integer :: i

C=0.0D0
DO i=1,size(A)
 C=C+A(i)*B(i)
END DO

end subroutine dot

subroutine writeAAA()

write(*,*) AAA

end subroutine writeAAA


END MODULE calc

PROGRAM test3
use calc
implicit none

call push(10.0D0)
call writeAAA()

END PROGRAM test3
