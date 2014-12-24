MODULE force
 implicit none
  interface
 subroutine force(x_i,F_i)
  real(8),intent(IN) :: x_i(:)
  real(8),intent(OUT):: F_i(size(x_i))
 end subroutine
 end interface
END MODULE force 

MODULE calc_EOM
 implicit none 
 interface
   subroutine calc_EOM(q,p,m,dt,force)
   use force
   real(8),intent(INOUT) :: q(:),p(:)
   real(8),intent(OUT)   :: m(:),dt
   end subroutine calc_EOM
 end interface
END MODULE calc_EOM


