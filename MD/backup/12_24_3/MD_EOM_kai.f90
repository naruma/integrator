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
        private VV, PV,gear_EOM,RK4_EOM
        public EOM_INI, EOM,EOM_END
        interface
                subroutine force(x_i,F_i)
                 real(8),intent(IN) :: x_i(:)
                 real(8),intent(OUT):: F_i(size(x_i))
                end subroutine
        end interface

        procedure(force),pointer,save   :: MD_force

        type MD_t
        integer                    :: N
        real(8),allocatable        :: m(:)
        real(8),allocatable        :: q_init(:),p_init(:)
        real(8)                    :: dt
        integer                    :: NSTEP
        character(64)               :: calc_EOM ! calc_EOM
	end type MD_t
        type(MD_t),private,save	    :: MD
contains
 
SUBROUTINE EOM_INI(q_init,p_init,m,dt,   &  !initial q, initial p, mass, time interval  
	           calc_EOM,force,       & !calculation method, MD_force(name of subroutine)
                   NSTEP              )    ! N step,N dim 
 real(8),intent(IN) :: q_init(:),p_init(:),m(:)
 real(8),intent(IN) :: dt
 integer,intent(IN) :: NSTEP
 character(*),intent(IN) :: calc_EOM ! calc_EOM
 interface
      subroutine force(x_i,F_i)
         real(8),intent(IN) :: x_i(:)
         real(8),intent(OUT):: F_i(size(x_i))
      end subroutine
 end interface

 MD%N=size(q_init)
 allocate(MD%q_init(MD%N),MD%p_init(MD%N),MD%m(MD%N))
 MD%q_init=q_init
 MD%p_init=p_init
 MD%m=m
 MD%dt=dt
 MD%NSTEP=NSTEP
 MD%calc_EOM=calc_EOM
 MD_force => force

END SUBROUTINE EOM_ini
! 
SUBROUTINE EOM()
  integer          :: i
  real(8)          :: q(MD%N),p(MD%N)

   q=MD%q_init
   p=MD%p_init

!  write(*,*) calc_EOM
  if (trim(MD%calc_EOM) == "VV" ) then
   DO i=1,MD%NSTEP
    call VV(q,p)
    write(1,*) i*MD%dt,q(1)               !debug
   END DO
  else if( trim(MD%calc_EOM) == "PV" ) then
   DO i=1,MD%NSTEP
    call PV(q,p)
    write(2,*) i*MD%dt,q(1)              !debug
   END DO
  else if (trim(MD%calc_EOM) == "GEAR" ) then
    call gear_EOM()
  else if (trim(MD%calc_EOM) == "RK4" ) then
    call RK4_EOM()
  else
    write(*,*) "invalid calculation method of EOM"
    STOP
  end if

END SUBROUTINE EOM

SUBROUTINE EOM_END()
 deallocate(MD%q_init,MD%p_init,MD%m)
 MD%dt=0.0D0
 MD%NSTEP=0
 MD%calc_EOM=' '
 nullify(MD_force)
END SUBROUTINE EOM_END
!-----------------------------------------------------------
!            calculation of Equation of motion
!-----------------------------------------------------------
! Velet method
SUBROUTINE VV(q,p)
! Velocity velet method
 real(8),intent(INOUT) :: q(:),p(:)
 real(8),dimension(size(q)) :: F,F_new,q_new,p_new

 call MD_force(q,F)
 q_new=q+MD%dt*p/MD%m+0.5D0*(MD%dt**2)*F/MD%m
 call MD_force(q_new,F_new)
 p_new=p+0.5D0*MD%dt*(F_new+F)*MD%m

 q=q_new
 p=p_new

END SUBROUTINE VV

SUBROUTINE PV(q,p)
! position velet method
 real(8),intent(INOUT) :: q(:),p(:)
 real(8),dimension(size(q)) :: F,q_new,p_new,temp_q

 temp_q=(q+0.5d0*MD%dt*p/MD%m)
 call MD_force(temp_q,F)
 p_new=p+MD%dt*F
 q_new=q+0.5D0*MD%dt*(p_new+p)/MD%m**2

! write(*,*) q(1),p(1)
 q=q_new
 p=p_new

END SUBROUTINE PV
!---------------------------------------------------
SUBROUTINE gear_EOM()
 use ODE

 integer               :: Ndim ! degree of Gear method
 integer               :: Ndrv ! solve Ndrv-th order ODE
 integer               :: i
 real(8),allocatable   :: y(:,:)

 !initial condition(Default)
 Ndim=6                 
 Ndrv=2                 

if (Ndrv==2) then
! y(:,:)= [q^(0) q^(1) q'^(2),...] 

 allocate(y(MD%N,0:Ndim-1))
 
 y=0.0D0

 y(:,0)=MD%q_init(:)
 y(:,1)=MD%p_init/MD%m
 call MD_force(y(:,0),y(:,2))

 write(*,*) "y for gear is", y(1,0)

 DO i=1,MD%NSTEP
   call Gear(MD%dt,2,y(:,0:),MD_force_m)
   write(4,*) i*MD%dt, y(1,0)  !debug
 END DO

else if (Ndrv == 1) then
! y(:,:)=[q q^(1) q^(2) q^(3) ... ]
!        [p p^(1) p^(2) p^(3) ... ]
!
  allocate(y(2*MD%N,0:Ndim-1))

  y=0.0D0

  DO i=1,MD%N
  y(i,0)=MD%q_init(i)
  y(i,1)=MD%p_init(i)/MD%m(i)
  END DO
  DO i=1,MD%N
  y(MD%N+i,0)=MD%p_init(i)
  END DO
  call MD_force(y(1:MD%N,0),y(MD%N+1:2*MD%N,1))

  write(*,*) "y for gear is", y(1,0)
  DO i=1,MD%NSTEP
   call Gear(MD%dt,1,y(:,0:),f_first)
   write(5,*) i*MD%dt, y(1,0)  !debug
  end do
  
end if

deallocate(y)

END SUBROUTINE gear_EOM
!---------------------------------------------
subroutine RK4_EOM()
 use ODE
 integer               :: i
 real(8),allocatable   :: y(:)


! y(:,:)=[q q^(1) q^(2) q^(3) ... ]
!        [p p^(1) p^(2) p^(3) ... ]
!
  allocate(y(2*MD%N))

  y=0.0D0

  DO i=1,MD%N
  y(i)=MD%q_init(i)
  END DO
  DO i=1,MD%N
  y(MD%N+i)=MD%p_init(i)
  END DO

  DO i=1,MD%NSTEP
   call RK4(MD%dt,y,f_first)
   write(3,*) i*MD%dt, y(1)  !debug
  end do

  deallocate(y) 

end subroutine RK4_EOM

subroutine MD_force_m(x_i,F_m_i)
! calculate F/m
   real(8),intent(IN)  :: x_i(:)
   real(8),intent(OUT) :: F_m_i(size(x_i))
   real(8)             :: F_i(size(x_i))
! interface
! subroutine MD_force(x_i,F_i)
!  real(8),intent(IN) :: x_i(:)
!  real(8),intent(OUT):: F_i(size(x_i))
! end subroutine
! end interface

   call MD_force(x_i,F_i)
   
   F_m_i=F_i/MD%m

end subroutine MD_force_m

subroutine f_first(x_i,f_i)
! function for first order ODE
    real(8),intent(IN) :: x_i(:)
    real(8),intent(OUT) :: f_i(size(x_i))
    integer             :: i,S
    
     S=size(x_i)

    if (mod(S,2) .ne. 0 ) then
        write(*,*) 'invalid dimension of array for the first order ODE'
        stop
    end if
    DO i=1,S/2
     f_i(i)=x_i(i+S/2)/MD%m(i)
    END DO    

     call MD_force(x_i(1:S/2),f_i(S/2+1:S))

end subroutine f_first


!SUBROUTINE LAI
! LAI method

!END SUBROUTINE LAI      

END MODULE MD_EOM


