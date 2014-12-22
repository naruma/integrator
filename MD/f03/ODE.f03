MODULE ODE
  implicit none
contains
! For first order
SUBROUTINE RK4(dt,y,func)
        ! solve dy/dt=f(y,t) by 4th-order Runge-Kutta
 real(8),intent(IN) :: dt
 real(8),intent(INOUT) :: y(:)
 real(8),dimension(size(y)) :: k1,k2,k3,k4,temp_f,k
 interface
  subroutine  func(y,F)
  real(8),intent(IN)  :: y(:)
  real(8),intent(OUT) :: F(size(y))
  end subroutine func
 end interface

 ! write(*,*) "foo"
 call func(y,temp_f)
 k1(:)=dt*temp_f(:)

 call func(y+0.50D0*k1,temp_f)
 k2(:)=dt*temp_f(:)

 call func(y+0.50D0*k2,temp_f)
 k3(:)=dt*temp_f(:)

 call func(y+k3,temp_f)
 k4(:)=dt*temp_f(:)

! write(*,*) "k1,k2,k3,k4"
! write(*,*) k1,k2,k3,k4
 k=(k1+2.0D0*k2+2.0D0*k3+k4)/6.0D0

 y(:)=y(:)+k(:)
END SUBROUTINE RK4

SUBROUTINE Gear(dt,n,y,func)
        use gear_parameter
        ! This algorithm calculates the dy/dt=f 
        ! The gear method needs r derivatives of y
        ! as an input.
        real(8),intent(IN) :: dt
        integer,intent(IN) :: n
        real(8),intent(INOUT) :: y(:,:)
        integer               :: r,y_dim,i,count
        real(8),allocatable   :: y_np1(:,:),new_y_np1(:,:),B(:,:),C(:),f_np1(:)
        real(8)               :: epsilon=1.0D-5
	interface
	  subroutine func(y,F)
	  real(8),intent(IN) :: y(:)
	  real(8)	     :: F(size(y))
	  end subroutine func
	end interface

        y_dim=size(y(:,1))
        r=size(y(1,:))

        ! make Y_N vector
        DO i=0,r-1
                y(:,i)=(dt**i)*y(:,i)/fac(i)
        END DO
        
        allocate(y_np1(N,0:r-1),new_y_np1(N,0:r-1),B(r,r),C(r),f_np1(N))

        ! Predictor
        call get_gear_parameter(r,n,B,C)
        call dgemm('N','T',N,r,r,1.0D0,y(:,0:r-1),N,B,N,0.0D0,y_np1(:,0:r-1))
       
        !corrector
        count=0
        call func(y_np1(:,0),f_np1)
        do while ( maxval(f_np1(:)-y_np1(:,1)) > epsilon)

        if (n==1) then
                DO i=1,y_dim
                new_y_np1(i,:)=y_np1(i,:)+(f_np1(i)-y_np1(i,1))*dt*C(:)
                END DO
        else if(n==2) then
                DO i=1,y_dim
                new_y_np1(i,:)=y_np1(i,:)+0.5d0*(f_np1(i)-y_np1(i,1))*(dt**2)*C(:)
                END DO
        else
                write(*,*) "the order of ODE has to be 1 or 2"
                stop
        end if

        call func(new_y_np1(:,0),f_np1)
        y_np1=new_y_np1
        count=count+1
        end do
END SUBROUTINE Gear

!   function f(y)
!   real(8),intent(IN) :: y(:)
!   real(8)            :: f(size(y))
!   f=-0.5D0*y
!   end function f

   integer function fac(r)
        integer,intent(IN) :: r
        integer            :: i

 if (r==0) then
         fac=0
 else if(r .lt. 0) then
         write(*,*) "n! has to be n>0"
         stop
 else
         fac=1
         DO i=1,r
         fac=fac*i
         END DO
 end if
   end function fac


END MODULE ODE
