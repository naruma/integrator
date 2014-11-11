MODULE time_evolution
 implicit none
 real(8)            :: hbar=1.0D0
 complex(kind(0d0))  :: i=(0.0D0,1.0D0)
 real(8),parameter :: pi=3.141592653589793238
contains

 SUBROUTINE CN(psi_init,H,dt,psi)
! Crank-Nicholson method
! A psi=A* psi_init
 complex(kind(0d0)),intent(IN) :: psi_init(:),H(:,:)
 real(8),intent(IN) :: dt
 complex(kind(0d0)),intent(OUT) ::psi(size(psi_init))
 complex(kind(0d0))  :: temp_H(size(psi_init),size(psi_init)), &
                        temp_psi(size(psi_init))
 integer            :: j,k,N
 N=size(psi_init)
 !get hamiltonian

 !obtain temp_psi=temp_H* \dot psi_init
 temp_psi=0.0D0
! temp_H(:,:)=1.0D0+i*dt*H(:,:)/hbar

 DO j=1,N
   DO k=1,N
     temp_psi(j)=temp_psi(j)+(1.0D0-0.5D0*i*dt*H(j,k)/hbar)*(psi_init(k))
   END DO
 END DO
 write(*,*) temp_psi
 
 call zgemv('N',N,N,(1.0D0,0.0D0),temp_H,N,psi_init,N,0.0D0,temp_psi,N)
 write(*,*) temp_psi

!{1-i*H+dt/hbar}*psi_init
temp_H=0.0D0
 DO j=1,N
   DO k=1,N
    temp_H(j,k)=1+i*dt*H(j,k)/hbar
   END DO
 END DO

 !obtain psi= temp_H^{-1} temp_psi
 psi=0.0D0

 END SUBROUTINE CN

 
 SUBROUTINE CB(psi_init,H,dt,M,psi)
 complex(kind(0d0)),intent(IN) :: psi_init(:),H(:,:)
 real(8), intent(IN)    :: dt
 integer :: N
 integer :: M
 complex(kind(0d0)),intent(OUT) :: psi(size(psi_init))
 real(8) :: E_MAX,E_MIN,a,b,alpha
 complex(kind(0d0)) :: temp_H(size(psi_init),size(psi_init)), &
                       v(size(psi_init),0:M)    
 integer ::j,k,l
 ! parameter
 alpha=0.01D0
 N=size(psi_init)
 !a,b 
 ! get E_MAX,E_MIN
 E_MAX=0.0D0; E_MIN=0.0D0
 
 a=(0.5D0+alpha)*(E_max-E_min)
 b=0.5D0*(E_max+E_min)
 ! shift hamiltonian for chebyshev scheme
 
 temp_H=(H-b)/a
 ! calculate vector v(:,k) k=1,...M

 v(:,0)=psi_init(:)
! v(:,1)=H(:,:)*v(:,0)
!DO k=1,N-1
! v(:,k+1)=2.0D0*H(:,:)*v(:,k)-v(:,k-1)

 ! calculate ck
 END SUBROUTINE CB

END MODULE time_evolution
