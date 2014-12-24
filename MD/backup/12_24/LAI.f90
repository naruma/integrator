MODULE LAI

Implicit none
     private
     save
! N_dim: degree of freedom, N_int: the number of intervals, 

! N_exp: degree of expansion
     integer  :: N_dim,N_int,N_exp
! for convergence condition
      real(8) :: eps       
! for derivative option
      logical :: anal_flag
      integer :: H_flag
! for time interval
      real(8) :: dt
contains

SUBROUTINE LAI_EOM(q_ini,p_ini,              &   ! input
                   N_dim,N_int,N_exp,eps,    &   ! global variables
                   anal_flag,H_flag,dt,      &   
                   q,p)                          ! output
! This program calculates time evolution of q and p
! using locally analytic integrator.
! for global 
     integer                ::     N_dim,N_int,N_exp
     real(8)                ::     eps       
     logical                ::     anal_flag
     integer                ::     H_flag
     real(8)                ::     dt

! q_ini,p_ini : initial condition of q and p
      real(8),intent(IN)    ::     q_ini(N_dim),p_ini(N_dim)
! q,p : q(t),p(t)
      real(8),intent(OUT)   ::     q(N_dim,N_int),p(N_dim,N_int)
      integer               ::     i,j,k
      real(8)               ::     a(N_dim,N_exp),a_new(N_dim,N_exp)
! for potential calculation
      real(8)               ::     dvpot_ini(N_dim),dvpot(N_dim,N_int),&
                                   norma



! create initial guess of a
     call guess_a(a)
       
! calculate derivative of potential dvpot at t=t_0
     call dv_pot(q_ini,dvpot_ini)

! self-consistent procedure
     DO k=1,N_dim
       DO j=1,N_exp
          a_new(k,j)=0.1D0+a(k,j)
       END DO
     END DO

     call norm_a(a,a_new,norma)
     do while (norma > eps)

! {q_k(t_i)}_{k,i} -> dvpot at t=t_i -> a at t-t_0
        call calc_q(q_ini,p_ini,dvpot_ini,a,q)
       DO i=1,N_int
        call dv_pot(q(:,i),dvpot(:,i))
       END DO
        call calc_a(dvpot_ini,dvpot,a_new)
     end do

END SUBROUTINE LAI_EOM

!---------------------------------------------------------

SUBROUTINE guess_a(a)
     real(8),intent(out) :: a(N_dim,N_exp)
     integer             :: k,j
     DO k=1,N_dim
        DO j=1,N_exp
           a(k,j)=0.0D0
        END DO
     END DO
END SUBROUTINE guess_a

!----------------------------------------------------------

SUBROUTINE dv_pot(q,dvpot)
     use hamiltonian
     real(8),intent(IN)               ::  q(N_dim)
     real(8),intent(OUT)              ::  dvpot(N_dim)
     logical                          ::  anal_OK

     if (anal_flag) then
        call check_anal_cal(H_flag,anal_OK)
        if (anal_OK) then
                call dH_dq_anal(H_flag,q,dvpot)
        else 
                call dH_dq_num(H_flag,q,dvpot)
        end if
     else
        call dH_dq_num(H_flag,q,dvpot)
     end if
END SUBROUTINE dv_pot

!----------------------------------------------------------        

SUBROUTINE norm_a(a,a_new,norma)
      real(8),intent(IN)  :: a(N_dim,N_exp),a_new(N_dim,N_exp)  
      real(8),intent(OUT) :: norma
      integer             :: k,j
      real(8)             :: temp

      norma=0.0D0
      DO k=1,N_dim
        DO j=1,N_exp
          temp=abs(a_new(k,j)-a(k,j))
          if (norma > temp) then
                  norma=temp
          end if
        END DO
      END DO

END SUBROUTINE norm_a

!----------------------------------------------------------
SUBROUTINE calc_q(q_ini,p_ini,dvpot_ini,a,q)
        real(8),intent(IN) :: q_ini(N_dim),p_ini(N_dim),   &
                              dvpot_ini(N_dim),a(N_dim,N_exp)
        real(8),intent(OUT):: q(N_dim,N_int)
        integer            :: k,i,j

        DO k=1,N_dim
          DO i=1,N_int
            q(k,i)=q_ini(k)+p_ini(k)*(dt*dble(i))+0.5D0   &
                  *dvpot_ini(k)*(dt*dble(i))**2
            DO j=1,N_exp
               q(k,i)=q(k,i)+a(k,j)*(dt*dble(i))**(j+2)
            END DO
          END DO
        END DO
END SUBROUTINE calc_q

!----------------------------------------------------------

SUBROUTINE calc_a(dvpot_ini,dvpot,a_new)
! algorithm of polynomial fitting
        real(8),intent(IN) :: dvpot_ini(N_dim),dvpot(N_dim,N_int) 
        real(8),intent(OUT):: a_new(N_dim,N_exp)
        integer            :: i,j,k
        real(8)            :: deltaf(N_dim,N_int),polycoeff(N_exp)
! for polcoe routine, the dimension of input and output has to be
! equal.
        if (N_int .ne. N_exp) then
          write(*,*)"ERROR: The number of sampling time and degree of polynomial &
                     has to be equal."
          STOP
        end if 
        DO k=1,N_dim
          DO i=1,N_int
             deltaf(k,i)=dvpot_ini(k)-dvpot(k,i)
          END DO
            call polcoe(deltaf(k,:),polycoeff(:))
          DO j=1,N_exp
            a_new(k,j)=polycoeff(j)/((j+1)*(j+2))
          END DO
        END DO
END SUBROUTINE
 

END MODULE LAI
