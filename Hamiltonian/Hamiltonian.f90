MODULE Hamiltonian
! This module contains the sets of Hamiltonian
! and its analytical derivative 
      implicit none
!      interface Hamiltonian
contains 
      SUBROUTINE Harmonic(q,p,H)
!       include "Harmonic.h"
       real(8),intent(IN)  :: q(:),p(:)
       real(8),intent(OUT) :: H
       real(8),allocatable :: k(:,:),r0(:,:),m(:) 
       integer             :: N,i,j

       N=ubound(q,1)
       allocate(k(N,N),r0(N,N),m(N))
! get parameter k,r0,m
!  test
       k  = 1.0D0
       r0 = 1.0D0
       m  = 1.0D0
!

       H=0.0D0

       DO i=1,N
        H=H+p(i)**2/(2.0D0*m(i))
       END DO
       DO i=1,N-1
         DO j=i,N
           H=H+0.5D0*k(i,j)*(dabs(q(i)-q(j))-r0(i,j))**2
         END DO
       END DO     

      deallocate(k,r0,m)
      END SUBROUTINE Harmonic
      
      SUBROUTINE LJ(q,p,H)
       real(8),intent(IN)  :: q(:),p(:)
       real(8),intent(OUT) :: H
       real(8),allocatable :: eps(:,:),r0(:,:),sigma(:,:),m(:)
       integer             :: N,i,j

       N=ubound(q,1)
       allocate(eps(N,N),r0(N,N),sigma(N,N),m(N))
! get parameter eps,r0,m,sigma
!  test
       eps   = 1.0D0
       r0    = 1.0D0
       m     = 1.0D0
       sigma = 1.0D0
!
       H=0.0D0
       DO i=1,N
        H=H+p(i)**2/(2.0D0*m(i))
       END DO
       DO i=1,N-1
         DO j=i,N
           H=H+4.0D0*eps(i,j)*((sigma(i,j)/(dabs(q(i)-q(j))-r0(i,j)))**12 &
             -(sigma(i,j)/(dabs(q(i)-q(j))-r0(i,j)))**6)
         END DO
       END DO     

       deallocate(eps,r0,sigma,m)

       END SUBROUTINE LJ

       SUBROUTINE Morse(q,p,H)
       real(8),intent(IN)  :: q(:),p(:)
       real(8),intent(OUT) :: H
       real(8),allocatable :: D(:,:),r0(:,:),a(:,:),m(:)
       integer             :: N,i,j


       N=ubound(q,1)
       allocate(D(N,N),r0(N,N),a(N,N),m(N))
! get parameter eps,r0,m,sigma
!  test
       D     = 1.0D0
       r0    = 1.0D0
       m     = 1.0D0
       a     = 1.0D0
!
       H=0.0D0
       DO i=1,N
        H=H+p(i)**2/(2.0D0*m(i))
       END DO
       DO i=1,N-1
         DO j=i,N
           H=H+D(i,j)*(1-dexp(-a(i,j)*(dabs(q(i)-q(j))-r0(i,j))))**2
         END DO
       END DO     

       deallocate(D,r0,a,m)

       END SUBROUTINE Morse
END MODULE Hamiltonian
