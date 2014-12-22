MODULE MD_EOM
        ! this routine numerically calculate the q,p after dt 
        ! by solving Netwon equation of motion.
        ! The input: 
        ! q:coordinate, p:momentum, dt:time interval
        ! m:mass
        ! The output:
        ! q: new coordinate, p:new momentum
        !
        ! the keyword of calc_EOM:
        ! VV    : 'velocity verlet method' Dec 12
        ! PV    : 'Position verlet method' Dec 12
        ! GEAR  : 'Gear method' Dec 12
        ! RK4   : '4th Runge-kutta method' Dec 12

        implicit none
contains 
SUBROUTINE EOM(q_init,p_init,m,dt,   &  !initial q, initial p, mass, time interval  
	       calc_EOM,force,       & !calculation method, force(name of subroutine)
               NSTEP               )    ! N step 
 real(8),intent(IN) :: q_init(:),p_init(:),m(:)
 real(8),intent(IN) :: dt
 integer,intent(IN) :: NSTEP
 character,intent(IN) :: calc_EOM ! calc_EOM
! force 
 interface
      subroutine force(x_i,F_i)
         real(8),intent(IN) :: x_i(:)
         real(8),intent(OUT):: F_i(size(x_i))
      end subroutine
 end interface
! 
  integer          :: i
  real(8)          :: q(size(q_init)),p(size(p_init))
  

  if (calc_EOM == "VV" ) then
   DO i=1,NSTEP
    call VV(q,p,m,dt,force)
   END DO
  else if( calc_EOM == "PV" ) then
   DO i=1,NSTEP
    call PV(q,p,m,dt,force)
   END DO
  else if (calc_EOM == "GEAR" ) then
    call gear_EOM(q_init,p_init,m,dt,Nstep,force)
  else if (calc_EOM == "RK4" ) then
    call RK4_EOM(q_init,p_init,m,dt,Nstep,force)
  else
    write(*,*) "invalid calculation method of EOM"
    STOP
  end if

END SUBROUTINE EOM
!-----------------------------------------------------------
!            calculation of Equation of motion
!-----------------------------------------------------------
! Velet method
SUBROUTINE VV(q,p,m,dt,force)
! Velocity velet method
 real(8),intent(INOUT) :: q(:),p(:)
 real(8),intent(IN) :: m(:),dt
 real(8),dimension(size(q)) :: F,F_new,q_new,p_new
 interface
 subroutine force(x_i,F_i)
  real(8),intent(IN) :: x_i(:)
  real(8),intent(OUT):: F_i(size(x_i))
 end subroutine
 end interface

 call force(q,F)
 q_new=q+dt*p/m+0.5D0*(dt**2)*F_new/m
 call force(q_new,F_new)
 p_new=p+0.5D0*dt*(F_new-F)/m

 q=q_new
 p=p_new

END SUBROUTINE VV

SUBROUTINE PV(q,p,m,dt,force)
! position velet method
 real(8),intent(INOUT) :: q(:),p(:)
 real(8),intent(IN) :: m(:),dt
 real(8),dimension(size(q)) :: F,F_new,q_new,p_new,temp_q
 interface
 subroutine force(x_i,F_i)
  real(8),intent(IN) :: x_i(:)
  real(8),intent(OUT):: F_i(size(x_i))
 end subroutine
 end interface

 temp_q=(q+0.5d0*dt*p/m)**(size(q))
 call force(temp_q,F)
 p_new=p+dt*F
 q_new=q+0.5D0*dt*(p_new+p)/m**2

 q=q_new
 p=p_new

END SUBROUTINE PV
!---------------------------------------------------
SUBROUTINE gear_EOM(q_init,p_init,m,dt,Nstep,force)
 use ODE
 real(8),intent(IN)         :: q_init(:),p_init(:),dt 
 real(8),intent(IN)         :: m(:)
 integer,intent(IN)         :: Nstep 
 interface
 subroutine force(x_i,F_i)
  real(8),intent(IN) :: x_i(:)
  real(8),intent(OUT):: F_i(size(x_i))
 end subroutine
 end interface
 integer               :: Ndim ! degree of Gear method
 integer               :: Ndrv ! solve Ndrv-th order ODE
 integer               :: i,N
 real(8),allocatable   :: y(:,:)

 !initial condition(Default)
 Ndim=6                 
 Ndrv=2                 

 N=size(q_init)

if (Ndrv==2) then
! y(:,:)= [q^(0) q^(1) q'^(2),...] 

 allocate(y(N,0:Ndim-1))
 
 y=0.0D0

 y(:,0)=q_init(:)
 y(:,1)=p_init/m

 DO i=1,Nstep
   call Gear(dt,2,y,force_m)
 END DO

else if (Ndrv == 1) then
! y(:,:)=[q q^(1) q^(2) q^(3) ... ]
!        [p p^(1) p^(2) p^(3) ... ]
!
  allocate(y(2*N,0:Ndim-1))

  y=0.0D0

  DO i=1,N
  y(i,0)=q(i)
  y(i,1)=p(i)/m(i)
  END DO
  DO i=1,N
  y(N+i,0)=p(i)
  END
  call force(y(1:N,0),y(N+1:2*N,1))

  DO i=1,Nstep
   call Gear(dt,1,y,f_first)
  end do
  
end if
contains 

subroutine force_m(x_i,F_m_i)
! calculate F/m
   real(8),intent(IN)  :: x_i(:)
   real(8),intent(OUT) :: F_m_i(size(x_i))
   real(8)             :: F_i(size(x_i))

   call force(x_i,F_i)
   
   F_m_i=F_i/m

end subroutine force_m


subroutine f_first(x_i,f_i)
! function for first order ODE
    real(8),intent(IN) :: x_i(:)
    real(8),intent(OUT) :: f_i(:)
    integer             :: i
    integer             :: N

    if (mod(N,2) \= 0 ) then
        write(*,*) 'invalid dimension of array for the first order ODE'
        stop
    end if
    DO i=1,N/2
     f_i(i)=x_i(i+N/2)/m(i)
    END DO    
    DO i=1,N/2
     call force(x_i(i),f_i(N/2+i))
    END DO
    END if
end subroutine f_first

END SUBROUTINE gear_EOM
!---------------------------------------------
!SUBROUTINE LAI
! LAI method

!END SUBROUTINE LAI      

END MODULE MD_EOM


