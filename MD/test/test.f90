PROGRAM test
      implicit none
      call calc(dot)
contains 
       subroutine calc(k)
       real(8) :: A(3),B(3),C
        interface
	  subroutine k(a,b,c)
	  real(8) :: a(3),b(3),c
	  end subroutine
	end interface 

       A=1.0d0
       B=2.0d0
       call k(A,B,C)
       write(*,*) C
       end subroutine calc

       subroutine dot(A,B,C)
       real(8) :: A(3),B(3),C
       integer :: i
        C=0.0d0
       do i=1,3
        C=C+A(i)*B(i)
       end do
       end subroutine dot
END PROGRAM test
