MODULE gear_parameter
      implicit none
contains
        subroutine get_gear_parameter(r,n,B,C)
 integer,intent(IN) :: r              ! r-1th derivative
 integer,intent(IN) :: n              ! dy^(n)/dt^(n)=f
 real(8),intent(OUT) :: C(r),B(r,r)
 real(8),parameter ::                                                   &
 C4_1(4)=(/3.0D0/8.0D0, 1.0D0, 0.750D0, 1.0D0/6.0D0 /),                 &
 C5_1(5)=(/251.0D0/720.0D0,1.0D0,11.0D0/12.D0,1.0D0/3.0D0,1.0D0/24.0D0/),&
 C6_1(6)=(/ 95.0D0/288.0D0,1.0D0, 25.0D0/24.0D0, 35.0D0/72.0D0,          &
            5.0D0/48.0D0, 1.0D0/120.0D0 /),                              &
 C4_2(4)=(/ 1.0D0/6.0D0, 5.0D0/6.0D0, 1.0D0, 1.0D0/3.0D0 /),             &
 C5_2(5)=(/ 19.0D0/120.0D0, 3.0D0/4.0D0, 1.0D0, 0.50D0, 1.0D0/12.0D0 /),  &
 C6_2(6)=(/ 3.0D0/20.0D0, 251.0D0/360.0D0, 1.0D0, 11.0D0/18.0D0,         &
          1.0D0/6.0D0, 1.0D0/60.0D0 /),                                &
! C7(7)=(/ 863.0D0/6048.0D0, 665.0D0/1008.0D0, 1.0D0, 25.0D0/36.0D0,    &
!          35.0D0/144.0D0, 1.0D0/24.0D0, 1.0D0/360.0D0 /),              &
! C8(8)=(/ 275.0D0/2016.0D0, 19087.0D0/30240.0D0, 1.0D0,137.0D0/180.0D0,&
!         5.0D0/16.0D0, 17.0D0/240.0D0, 1.0D0/120.0D0, 1.0D0/2520.0D0/),&
 
  B4(4,4)=reshape((/1.0D0, 0.0D0, 0.0D0, 0.0D0,                        &
                    1.0D0, 1.0D0, 0.0D0, 0.0D0,                        &
                    1.0D0, 2.0D0, 1.0D0, 0.0D0,                        &
                    1.0D0, 3.0D0, 3.0D0, 1.0D0/), shape(B4)), &
  B5(5,5)=reshape((/1.0D0,0.0D0,0.0D0,0.0D0,0.0D0,                     &
                    1.0D0,1.0D0,0.0D0,0.0D0,0.0D0,                     &
                    1.0D0,2.0D0,1.0D0,0.0D0,0.0D0,                     &
                    1.0D0,3.0D0,3.0D0,1.0D0,0.0D0,                     &
                    1.0D0,4.0D0,6.0D0,4.0D0,1.0D0/),shape(B5)),        &
  B6(6,6)=reshape((/1.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,               &
                    1.0D0,1.0D0,0.0D0,0.0D0,0.0D0,0.0D0,               &
                    1.0D0,2.0D0,1.0D0,0.0D0,0.0D0,0.0D0,               &
                    1.0D0,3.0D0,3.0D0,1.0D0,0.0D0,0.0D0,               &
                    1.0D0,4.0D0,6.0D0,4.0D0,1.0D0,0.0D0,               &
                    1.0D0,5.0D0,10.0D0,10.0D0,5.0D0,1.0D0/),shape(B6))
!  B7(7,7)
!  B8(8,8)
 if (n==1) then
     if (r==4) then
             B=B4
             C=C4_1
     else if (r==5) then
             B=B5
             C=C5_1
     else if (r==6) then
             B=B6
             C=C6_1
     else
           write(*,*) "The current Gear method only supports 4,5,6th order of &
                       derivatives. The current order is",r
           stop
     end if
 else if (n==2) then
     if (r==4) then
             B=B4
             C=C4_2
     else if (r==5) then
             B=B5
             C=C5_2
     else if (r==6) then
             B=B6
             C=C6_2
     else
           write(*,*) "The current Gear method only supports 4,5,6th order of &
                       derivatives. The current order is",r
           stop
     end if
  else
           write(*,*) "The current Gear method only supports 1st and 2nd order of &
                       ODE. The current order is",n
           stop
 end if
 end subroutine get_gear_parameter

END MODULE gear_parameter

