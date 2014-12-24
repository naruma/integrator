PROGRAM test_MD
     use MD_EOM
     implicit none
     real(8) :: q_init(2),p_init(2),m(2)
     
     q_init=1.0D0
     p_init=1.0D0
     m=1.0D0

     call EOM_INI(q_init,p_init,m,0.01d0,"VV",harmonic,1000)
     call EOM()
     call EOM_END()
     call EOM_ini(q_init,p_init,m,0.01d0,"PV",harmonic,1000)
     call EOM()
     call EOM_END()
     call EOM_ini(q_init,p_init,m,0.01d0,"RK4",harmonic,1000)
     call EOM()
     call EOM_END()
     call EOM_ini(q_init,p_init,m,0.01d0,"GEAR",harmonic,1000)
     call EOM()
     call EOM_END()

      contains
              subroutine harmonic(q,F)
              real(8),intent(IN) :: q(:)
              real(8),intent(OUT):: F(size(q))
              real(8) :: k
              
              k=1.0D0

              F=-k*q

              end subroutine harmonic

END PROGRAM test_MD
