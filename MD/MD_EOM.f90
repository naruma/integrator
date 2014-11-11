MODULE MD_EOM
  implicit none
contains
SUBROUTINE VV
! Velocity velet method
 real(8),intent(INout) :: q(:),p(:)
 real(8),intent(IN) :: m(:),dt
 real(8),dimension(size(q)) :: F,F_new,q_new,p_new

 call force(q,F)
 q_new=q+dt*p/m+0.5D0*(dt**2)*F_new/m
 call force(q_new,F_new)
 p_new=p+0.5D0*dt*(F_new-F)/m

 q=q_new
 p=p_new

END SUBROUTINE VV

!SUBROUTINE GEAR
! Gear method

!END SUBROUTINE GEAR      

!SUBROUTINE LAI
! LAI method

!END SUBROUTINE LAI      

END MODULE EOM


