MODULE MD_EOM
        ! this routine numerically calculate the q,p after dt 
        ! by solving Netwon equation of motion.
        ! The input: 
        ! q:coordinate, p:momentum, dt:time interval
        ! m:mass
        ! The output:
        ! q: new coordinate, p:new momentum
        !
        implicit none
contains
!------------------------------------------------------------
! Velet method
SUBROUTINE VV(q,p,m,dt)
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

SUBROUTINE PV(q,p,m,dt)
! position velet method
 ! Velocity velet method
 real(8),intent(INout) :: q(:),p(:)
 real(8),intent(IN) :: m(:),dt
 real(8),dimension(size(q)) :: F,F_new,q_new,p_new,temp_q

 temp_q=(q+0.5d0*dt*p/m)**(size(q)
 call force(temp_q,F)
 p_new=p+dt*F
 q_new=q+0.5D0*dt*(p_new+p)/m**2

 q=q_new
 p=p_new
END SUBROUTINE PV

!-------------------------------------------------------------
SUBROUTINE GEAR(q,p,dt,m)
! Gear method
! dt will be optimized during the multi-value method
! set initial condition and dt by  

END SUBROUTINE GEAR      

!SUBROUTINE LAI
! LAI method

!END SUBROUTINE LAI      

END MODULE EOM


