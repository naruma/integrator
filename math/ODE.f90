MODULE ODE
  implicit none
contains
SUBROUTINE RK4(dt,y)
 real(8),intent(IN) :: dt
 real(8),intent(INOUT) :: y(:)
 real(8),dimension(size(y)) :: k1,k2,k3,k4,temp_f,k

! write(*,*) "foo"
 call force(y,temp_f)
 k1(:)=dt*temp_f(:)

 call force(y+0.50D0*k1,temp_f)
 k2(:)=dt*temp_f(:)

 call force(y+0.50D0*k2,temp_f)
 k3(:)=dt*temp_f(:)

 call force(y+k3,temp_f)
 k4(:)=dt*temp_f(:)

! write(*,*) "k1,k2,k3,k4"
! write(*,*) k1,k2,k3,k4
 k=(k1+2.0D0*k2+2.0D0*k3+k4)/6.0D0

 y(:)=y(:)+k(:)

END SUBROUTINE RK4

SUBROUTINE Gear(dt,y)
        use gear_parameter
        ! warning! This Gear method needs r derivatives of
        ! y as an input. Namely,
        ! y(:,1) = y_n(:)
        ! y(:,2) = d(y_n(:))/dt
        ! y(:,3) = d^2(y_n(:))/dt^2
        ! ...

END SUBROUTINE Gear

subroutine force(y,f)
   real(8),intent(IN) :: y(:)
   real(8),intent(OUT) :: f(size(y))
   f=-0.5D0*y
 END subroutine

END MODULE ODE
