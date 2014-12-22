MODULE B_G
 implicit none
contains
 real(8) function B(A,k)
 real(8) :: A
 real(8),save :: m=5.0D0
 interface 
         real(8) function k(A)
                 real(8) A
         end function k
 end interface
 B=G(A)-1.0D0
 end function B

 real(8) function G(A)
  real(8) :: A
  G=2.0D0*A
 end function G
END MODULE B_G

PROGRAM test1
        use B_G
        implicit none
        real(8) :: A,C
        A=1.0D0

        C=B(A,G)
        write(*,*) C
END PROGRAM test1
