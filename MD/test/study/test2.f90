PROGRAM test3
      implicit none
      real(8) :: a,b,c
      real(8),save :: m

a=5.0D0
b=1.0D0
m=3.0D0
      call calc(a,b,c)
      write(*,*) c
contains
      subroutine calc(A,B,C)
              real(8),intent(IN) :: A,B
              real(8),intent(OUT) :: C
        write(*,*) m
              C=A/B*m
      end subroutine calc
end program test3
