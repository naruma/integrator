PROGRAM test
   use ODE
   implicit none
   real(8) :: dt=0.01,y(1)=1.0D0
   integer :: i

 DO i=1,1000
   call RK4(dt,y)
   write(1,*) dt*dble(i),y
 END DO
END PROGRAM test

   
