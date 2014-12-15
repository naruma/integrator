PROGRAM test_input
  use input
  implicit none
 real(8) :: x(9)
 real(8) :: r(3,3)
 integer :: i,j
 
 DO i=1,9
   x(i)=dble(i)
 END DO
 
 DO i=1,9
 write(*,*) x(i),atom(i)
 END DO
 call q_r(x,r)
 DO j=1,3
 DO i=1,3
 write(*,*) r(i,j)
 END DO
 END DO

END PROGRAM test_input
