MODULE input
 implicit none 
contains
 subroutine q_r(q,r)
! this routine convert generalized coordinate q 
! to cartesian r(3,Natom).

  real(8),intent(IN)	 :: q(:)
  real(8),intent(OUT) 	 :: r(3,int(size(q)/3))
  integer    		 :: i,Natom

  Natom=int(size(q)/3)
  r=0.0D0

! assume q is the array of [x1 y1 z1,x2 y2 z2,x3 y3 z3,...]
  DO i=1,size(q)
  if (mod(i,3) == 1) then
  r(1,atom(i))=q(i)
  else if (mod(i,3) == 2) then
  r(2,atom(i))=q(i)
  else if (mod(i,3) == 0) then
  r(3,atom(i))=q(i)
  end if
  END DO
 end subroutine q_r

 integer function atom(i)
! this function shows to which atom one degree of freedom in 3N
! belongs.
 integer :: i
! for the array [x1 y1 z1,x2 y2 z2,x3 y3 z3,...]
 atom=int((i-1)/3)+1
 end function atom 

END MODULE input
