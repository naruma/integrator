MODULE Hamiltonian
! This module contains the sets of Hamiltonian
! and its analytical derivative 
! VPOT=1 : Harmonic potential
!      2 : LJ potential
!      3 : Morse potential

      implicit none
!      interface Hamiltonian
contains 
      SUBROUTINE Hamiltonian(q,p,H)
      use input
      real(8),intent(IN)  :: q(:),p(:)
      real(8),intent(OUT) :: H,vpot
      integer		  :: N,i,j      
      real(8) 		  :: r(3,int(size(q)/3))
      N=ubound(q,1)

      H=0.0D0
! sum up kinetic energy
       DO i=1,N
        H=H+p(i)**2/(2.0D0*m(i))
       END DO
! add potential energy
     call q_r(q,r)
      if (VPOT==1) then
      call Harmonic(r,vpot)
      else if (VPOT==2) then
      call LJ(r,vpot)
      else if (VPOT==3) then
      call Morse(r,vpot)
      else
      write(*,*) "Invalid selection of potential."
      stop
      end if

      H=H+vpot

      END SUBROUTINE Hamiltonian

! potential term
      SUBROUTINE Harmonic(r,vpot)
!       include "Harmonic.h"
       real(8),intent(IN)  :: r(:,:)
       real(8),intent(OUT) :: vpot
       real(8),allocatable :: k(:,:),r0(:,:),m(:) 
       integer             :: N,i,j
!       real(8)             :: qa(:,:) ! coordinate with atom

       N=ubound(r,2)
       allocate(k(N,N),r0(N,N),m(N),qa(N,N))
! get parameter k,r0,m
!  test
       k  = 1.0D0
       r0 = 1.0D0
       m  = 1.0D0
!

       vpot=0.0D0

       DO i=1,N-1
         DO j=i,N
           vpot=vpot+0.5D0*k(i,j)*(norm(r(:,i)-r(:,j))-r0(i,j))**2
         END DO
       END DO     

      deallocate(k,r0,m)
      END SUBROUTINE Harmonic
      
      SUBROUTINE LJ(r,vpot)
       real(8),intent(IN)  :: r(:,:)
       real(8),intent(OUT) :: vpot
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
      vpot=0.0d0

      DO i=1,N-1
         DO j=i,N
           vpot=vpot+4.0D0*eps(i,j)*((sigma(i,j)/(norm(r(:,i)-r(:,j))-r0(i,j)))**12 &
             -(sigma(i,j)/(norm(r(:,i)-r(:,j))-r0(i,j)))**6)
         END DO
       END DO     

       deallocate(eps,r0,sigma,m)

       END SUBROUTINE LJ

       SUBROUTINE Morse(r,vpot)
       real(8),intent(IN)  :: r(:,:)
       real(8),intent(OUT) :: vpot
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
      vpot=0.0d0

      DO i=1,N-1
         DO j=i,N
           vpot=vpot+D(i,j)*(1-dexp(-a(i,j)*(norm(r(:,i)-r(:,j))-r0(i,j))))**2
         END DO
       END DO     

       deallocate(D,r0,a,m)

       END SUBROUTINE Morse
END MODULE Hamiltonian
