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
        real(8),save   :: m(:)
        interface
         subroutine force(x_i,F_i)
          real(8),intent(IN) :: x_i(:)
          real(8),intent(OUT):: F_i(size(x_i))
         end subroutine
        end interface

contains 
SUBROUTINE EOM(q_init,p_init,m,dt,   &  !initial q, initial p, mass, time interval  
	       calc_EOM,force,       & !calculation method, force(name of subroutine)
               NSTEP               )    ! N step 
 real(8),intent(IN) :: q_init(:),p_init(:),m(:)
 real(8),intent(IN) :: dt
 integer,intent(IN) :: NSTEP
 character(*),intent(IN) :: calc_EOM ! calc_EOM
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

   q=q_init
   p=p_init

!  write(*,*) calc_EOM
  if (calc_EOM == "VV" ) then
   DO i=1,NSTEP
    call VV(q,p,m,dt,force)
    write(1,*) i*dt,q(1)               !debug
   END DO
  else if( calc_EOM == "PV" ) then
   DO i=1,NSTEP
    call PV(q,p,m,dt,force)
    write(2,*) i*dt,q(1)              !debug
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
! real(8),intent(IN) :: m(:)
 real(8),intent(IN) :: dt
 real(8),dimension(size(q)) :: F,F_new,q_new,p_new
 interface
 subroutine force(x_i,F_i)
  real(8),intent(IN) :: x_i(:)
  real(8),intent(OUT):: F_i(size(x_i))
 end subroutine
 end interface

 call force(q,F)
 q_new=q+dt*p/m+0.5D0*(dt**2)*F/m
 call force(q_new,F_new)
 p_new=p+0.5D0*dt*(F_new+F)*m

 q=q_new
 p=p_new

END SUBROUTINE VV

SUBROUTINE PV(q,p,m,dt,force)
! position velet method
 real(8),intent(INOUT) :: q(:),p(:)
 real(8),intent(IN) :: dt
! real(8),intent(IN) :: m(:)
 real(8),dimension(size(q)) :: F,q_new,p_new,temp_q
 interface
 subroutine force(x_i,F_i)
  real(8),intent(IN) :: x_i(:)
  real(8),intent(OUT):: F_i(size(x_i))
 end subroutine
 end interface

 temp_q=(q+0.5d0*dt*p/m)
 call force(temp_q,F)
 p_new=p+dt*F
 q_new=q+0.5D0*dt*(p_new+p)/m**2

! write(*,*) q(1),p(1)
 q=q_new
 p=p_new

END SUBROUTINE PV
!---------------------------------------------------
SUBROUTINE gear_EOM(q_init,p_init,m,dt,Nstep,force)
 use ODE
 real(8),intent(IN)         :: q_init(:),p_init(:),dt 
! real(8),intent(IN),save         :: m(:)
 integer,intent(IN)         :: Nstep 
 interface,save
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
  y(i,0)=q_init(i)
  y(i,1)=p_init(i)/m(i)
  END DO
  DO i=1,N
  y(N+i,0)=p_init(i)
  END DO
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
! interface
! subroutine force(x_i,F_i)
!  real(8),intent(IN) :: x_i(:)
!  real(8),intent(OUT):: F_i(size(x_i))
! end subroutine
! end interface

   call force(x_i,F_i)
   
   F_m_i=F_i/m

end subroutine force_m


subroutine f_first(x_i,f_i)
! function for first order ODE
    real(8),intent(IN) :: x_i(:)
    real(8),intent(OUT) :: f_i(:)
    integer             :: i
    integer             :: N
! interface
! subroutine force(x_i,F_i)
!  real(8),intent(IN) :: x_i(:)
!  real(8),intent(OUT):: F_i(size(x_i))
! end subroutine
! end interface

    if (mod(N,2) .ne. 0 ) then
        write(*,*) 'invalid dimension of array for the first order ODE'
        stop
    end if
    DO i=1,N/2
     f_i(i)=x_i(i+N/2)/m(i)
    END DO    
     call force(x_i(1:N/2),f_i(N/2+1:N))
end subroutine f_first

END SUBROUTINE gear_EOM
!---------------------------------------------
subroutine RK4_EOM(q_init,p_init,m,dt,Nstep,force)
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
 integer               :: i,N
 real(8),allocatable   :: y(:)


 N=size(q_init)
! y(:,:)=[q q^(1) q^(2) q^(3) ... ]
!        [p p^(1) p^(2) p^(3) ... ]
!
  allocate(y(2*N))

  y=0.0D0

  DO i=1,N
  y(i)=q_init(i)
  END DO
  DO i=1,N
  y(N+i)=p_init(i)
  END DO

  DO i=1,Nstep
   call RK4(dt,y,f_first)
   write(3,*) i*dt, y(1)  !debug
  end do

 
contains 
  subroutine f_first(x_i,f_i,force)
! function for first order ODE
    real(8),intent(IN) :: x_i(:)
    real(8),intent(OUT) :: f_i(:)
    integer             :: i
    integer             :: N
 interface
 subroutine force(x_i,F_i)
  real(8),intent(IN) :: x_i(:)
  real(8),intent(OUT):: F_i(size(x_i))
 end subroutine
 end interface

    if (mod(N,2) .ne. 0 ) then
        write(*,*) 'invalid dimension of array for the first order ODE'
        stop
    end if
    DO i=1,N/2
     f_i(i)=x_i(i+N/2)/m(i)
    END DO    
     call force(x_i(1:N/2),f_i(N/2+1:N))
  end subroutine f_first


end subroutine RK4_EOM
!SUBROUTINE LAI
! LAI method

!END SUBROUTINE LAI      

END MODULE MD_EOM


