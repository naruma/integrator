MODULE derivative_H
! This module contains analytical derivative of
! Hamiltonian.
! VPOT=1 : Harmonic potential
!      2 : LJ potential
!      3 : Morse potential
      implicit none

contains 
      SUBROUTINE dH_dp(q,p,m,dH_dp)
      real(8),intent(IN)  :: q(:),p(:),m(:)
      real(8),intent(OUT) :: dH_dp(size(q))
  
      dH_dp=0.0D0

      dH_dp=p/m
      ! for the calculation of dv_dp
!     ! if V is dependent to p, then add conditional
      ! sentence here.
!     if (VPOT== ?? ) then
!     call d(POTENTIAL)_dq(q,p,m,dH_dp)
!     end if

      dv_dp=0.0D0      
      dH_dp=dH_dp+dv_dp

      END SUBROUTINE dH_dp


      SUBROUTINE dH_dq(q,p,m,dH_dq)
      real(8),intent(IN)  :: q(:),p(:),m(:)
      real(8),intent(OUT) :: dH_dq(size(q))
      integer		  :: N,i,j
      
      N=ubound(q,1)

      dH_dq=0.0D0

! add potential energy
      if (VPOT==1) then
      call dHarmonic(q,dH_dq)
      else if (VPOT==2) then
      call dLJ(q,dH_dq)
      else if (VPOT==3) then
      call dMorse(q,dH_dq)
      else
      write(*,*) "Invalid selection of potential."
      stop
      end if

      END SUBROUTINE dH_dq

! potential term
      SUBROUTINE dHarmonic(q,dv_dq)
!       include "Harmonic.h"
       real(8),intent(IN)  :: q(:)
       real(8),intent(OUT) :: dv_dq(size(q))
       real(8),allocatable :: k(:,:),r0(:,:),m(:) 
       integer             :: N,i,j
       real(8),allocatable :: qa(N,N)

       N=ubound(q,1)
       allocate(k(N,N),r0(N,N),m(N),qa(N,N))
! get parameter k:force constant,r0:,m:mass
!  test
       k  = 1.0D0
       r0 = 1.0D0
       m  = 1.0D0
!
       dv_dq=0.0D0

! get atom number of q -> r(3, # of atom)
       call 
       DO i=1,N-1
         DO j=i,N
           vpot=vpot+0.5D0*k(i,j)*(dabs(q(i)-q(j))-r0(i,j))**2
         END DO
       END DO     

      deallocate(k,r0,m)
      END SUBROUTINE dHarmonic
      
      SUBROUTINE dLJ(q,vpot)
       real(8),intent(IN)  :: q(:)
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
           vpot=vpot+4.0D0*eps(i,j)*((sigma(i,j)/(dabs(q(i)-q(j))-r0(i,j)))**12 &
             -(sigma(i,j)/(dabs(q(i)-q(j))-r0(i,j)))**6)
         END DO
       END DO     

       deallocate(eps,r0,sigma,m)

       END SUBROUTINE dLJ

       SUBROUTINE dMorse(q,vpot)
       real(8),intent(IN)  :: q(:)
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
           vpot=vpot+D(i,j)*(1-dexp(-a(i,j)*(dabs(q(i)-q(j))-r0(i,j))))**2
         END DO
       END DO     

       deallocate(D,r0,a,m)

       END SUBROUTINE dMorse
END MODULE derivative_H
