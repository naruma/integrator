module m_gl  
  use OpenGL_GL
  use OpenGL_GLUT
  implicit none  
  real, save :: xrot = 0.0, yrot = 0.0
contains
  subroutine SetupRC() bind(c)
    ! Black background
    call glClearColor(0.0, 0.0, 0.6, 1.0 )
    ! Set drawing color to green
    call glColor3f(1.0, 0.0, 0.0)
  end subroutine SetupRC

  subroutine KeyPressFunc(key, x, y) bind(C)
    integer(GLbyte), intent(in), value :: key
    integer(GLint ), intent(in), value  :: x, y
    if ( key == 27 ) stop ! ESC
  end subroutine KeyPressFunc

  subroutine KeySpecialFunc(key, x, y) bind(C)
    integer(GLint), intent(in), value  :: key, x, y
    if ( key == GLUT_KEY_UP    ) xRot = xrot - 5.0
    if ( key == GLUT_KEY_DOWN  ) xRot = xrot + 5.0
    if ( key == GLUT_KEY_LEFT  ) yRot = yrot - 5.0
    if ( key == GLUT_KEY_RIGHT ) yRot = yrot + 5.0 
    if ( xrot > 356.0 )then
      xRot = 0.0
    else if ( xrot < -1.0 )then
      xRot = 355.0
    endif
    if ( yrot > 356.0 )then
      yRot = 0.0
    else if ( yrot < -1.0 )then
      yRot = 355.0
    endif
    ! Refresh the Window
    call glutPostRedisplay
  end subroutine KeySpecialFunc

  subroutine ChangeSize(win, hin) bind(C)
    integer(GLcint), intent(IN), value :: win, hin
    integer(GLcint) :: w, h
    real(GLdouble)  :: Zero, One, Range = 100.0, Aspect
    w = win
    h = hin
    ! Prevent a divide by zero, when window is too short
    ! (you cant make a window of zero width).
    if( h == 0 ) h = 1
    ! Set the viewport to be the entire window
    call glViewport(0, 0, w, h)
    ! Reset coordinate system
    call glMatrixMode(GL_PROJECTION)
    call glLoadIdentity
    ! Establish clipping volume (left, right, bottom, top, near, far)
    aspect = float(w)/float(h)
    ! Keep the square square
    if( w <= h )then
      call glOrtho( -Range, Range, &
                    -Range / Aspect, Range / Aspect,&
                    -Range, Range  )
    else
      call glOrtho( -Range * Aspect, Range * Aspect, &
                    -Range, Range, &
                    -Range, Range  )
    endif
    call glMatrixMode(GL_MODELVIEW)
    call glLoadIdentity
  end subroutine ChangeSize 

  subroutine RenderScene() bind(C)
    integer, parameter :: n = 10**4
    real(glfloat), save :: dx(n - 1), dy(n - 1), dz(n - 1)
    logical, save :: first = .true.
    integer :: i
    ! make Lorenz attractor data        

    if (first) then
      block  
        integer(8) :: m = 2_8**31
        integer(8) :: ix 
        real(glfloat) :: x(n), y(n), z(n)
        first = .false.
        ix = 1
        do i = 1, n
          ix = mod(65539 * ix, m)  ! ;  print *, ix
          x(i) = real(ix) / m 
          ix = mod(65539 * ix, m)  ! ;  print *, ix
          y(i) = real(ix) / m
          ix = mod(65539 * ix, m)  ! ;  print *, ix
          z(i) = real(ix) / m
        end do
        dx = (x(2:) - x(:n - 1)) * 80.0_glfloat ! scale data 
        dy = (y(2:) - y(:n - 1)) * 80.0_glfloat
        dz = (z(2:) - z(:n - 1)) * 80.0_glfloat
      end block  
    end if  
    ! Clear the window with current clearing color
    call glClear(GL_COLOR_BUFFER_BIT)
    ! Save matrix state and do the rotation
    call glPushMatrix
    call glRotatef(xRot, 1.0, 0.0, 0.0)
    call glRotatef(yRot, 0.0, 1.0, 0.0)
    ! draw Lorenz attractor
    do i = 1, n - 1
      call glTranslatef( dx(i), dy(i), dz(i) )
      call glutSolidCube(0.5_gldouble )
    end do
    !
    call glPopMatrix
    ! Flush drawing commands
    call glutSwapBuffers
  end subroutine RenderScene
end module m_gl

program RANDU
   use m_gl
   implicit none
   integer :: iwin
   call glutInit()
   call glutInitDisplayMode(ior(GLUT_DOUBLE,ior(GLUT_RGB,GLUT_DEPTH)))
   call glutInitWindowSize(800, 600)
   iwin = glutCreateWindow('RANDU'//char(0))
   call glutReshapeFunc ( ChangeSize     )
   call glutKeyboardFunc( KeyPressFunc   )
   call glutSpecialFunc ( KeySpecialFunc )
   call glutDisplayFunc ( RenderScene    )
   call SetupRC()
   call glutMainLoop()
end program RANDU

