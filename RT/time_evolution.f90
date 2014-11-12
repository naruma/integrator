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
 complex(kind(0d0)),intent(OUT),allocatable ::psi(:)
 complex(kind(0d0)),allocatable  :: temp_H(:,:),temp_psi(:)
 integer            :: j,k,N
! for lapack
 integer,allocatable :: IPIV(:)
 integer :: INFO
 N=size(psi_init)
 allocate(psi(N),temp_H(N,N),temp_psi(N))
 !get hamiltonian

 !obtain temp_psi=temp_H* \dot psi_init
! temp_H(:,:)=1.0D0+i*dt*H(:,:)/hbar
temp_psi=0.0D0
! DO j=1,N
!   DO k=1,N
!     temp_psi(j)=temp_psi(j)+(1.0D0-0.5D0*i*dt*H(j,k)/hbar)*(psi_init(k))
!   END DO
! END DO
! write(*,*) temp_psi(1),temp_psi(2)

 DO j=1,N
  DO k=1,N
     temp_H(j,k)=1.0D0-0.5D0*i*dt*H(j,k)/hbar
  END DO
 END DO

 call zgemv('N',N,N,(1.0D0,0.0D0),temp_H,N,psi_init,1,(0.0D0,0.0D0),temp_psi,1)
! write(*,*) temp_psi(1),temp_psi(2)

 temp_H=conjg(temp_H)

!{1-i*H+dt/hbar}*psi_init
!obtain psi= temp_H^{-1} temp_psi

allocate (IPIV(N))
 call zgesv(N,1,temp_H,N,IPIV,temp_psi,N,INFO)
 psi(:)=temp_psi(:)

 deallocate(temp_psi,temp_H,IPIV)
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
 call zgemv('N',N,N,(1.0D0,0.0D0),2.0D0*temp_H,N,v(:,0),1,(0.0D0,0.0D0),v(:,1),1)
 DO k=1,N-1
 call zgemv('N',N,N,(1.0D0,0.0D0),2.0D0*temp_H,N,v(:,k),1,(0.0D0,0.0D0),v(:,k+1),1)
 v(:,k+1)=v(:,k+1)-v(:,k-1)
 END DO

! v(:,1)=H(:,:)*v(:,0)
!DO k=1,N-1
! v(:,k+1)=2.0D0*H(:,:)*v(:,k)-v(:,k-1)

 ! calculate ck
 END SUBROUTINE CB

END MODULE time_evolution
