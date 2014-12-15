MODULE derivative
    IMPLICIT NONE
contains
  SUBROUTINE threep(func,x,h,df)
    real(8),intent(IN) :: h,x
    real(8),intent(OUT):: df
    external           :: func
    real(8)            :: f1,f_1 

    call func(x+h,f1)
    call func(x-h,f_1)

    df=(f1-f_1)/(2.0D0*h)

  END SUBROUTINE threep

  SUBROUTINE fivep(func,x,h,df)
    real(8),intent(IN) :: h,x
    real(8),intent(OUT):: df
    external           :: func
    real(8)            :: f2,f1,f_1,f_2

    call func(x+2.0D0*h,f2)
    call func(x+h,f1)
    call func(x-h,f_1)
    call func(x-2.0D0*h,f_2)

    df=(-f2+8.0D0*f1-8.0D0*f_1+f_2)/(12.0D0*h)

  END SUBROUTINE fivep

END MODULE derivative
