MODULE calc_k
 implicit none
!       real(8) :: A(3),B(3),C
!       integer :: i
        interface
          subroutine k(a,b,c)
	  real(8) :: a(3),b(3),c
	  end subroutine
	end interface 
END MODULE calc_k

MODULE calc_test
  implicit none
contains 
       subroutine calc_tot(A,B,C,k)
 real(8) :: A(3),B(3),C
          interface
                  subroutine k(a,b,c)
        	  real(8) :: a(3),b(3),c
                  end subroutine
          end interface
       call k(A,B,C)
       end subroutine calc_tot

       subroutine dot(A,B,C)
 real(8) :: A(3),B(3),C
               integer :: i
        C=0.0d0
       do i=1,3
        C=C+A(i)*B(i)
       end do
       end subroutine dot

       subroutine sum(A,B,C)
  real(8) :: A(3),B(3),C
               integer :: i
        C=0.0d0
       do i=1,3
        C=C+A(i)+B(i)
       end do
       end subroutine sum
END MODULE calc_test


PROGRAM test
      use calc_test
      implicit none
      real(8) :: A(3),B(3),C
      A=1.0d0
      B=2.0d0
      call calc_tot(A,B,C,dot)
      write(*,*) C
      call calc_tot(A,B,C,sum)
      write(*,*) C
END PROGRAM test

