      program Bisection

      implicit none
      real*8 a, b, c, nc, error, tolerance, test,f
      integer i
      print*,'introducir tolerancia'
      read*, tolerance
      print*,'introducir a'
      read*,a
      print*,'introducir b'
      read*,b
      c=0.
      nc =0.
      i=0
      error = 10.
c     Initialization of some variables. 
      do while(error.GT.tolerance)
        nc = (a+b)/2
        test = f(a)*f(nc)
        IF(test.GT.0) THEN !if they have different signs, update a. If
c       same sign, update b.
          a = nc
        ELSE
          b=nc
        END IF
        error = ABS(nc-c)
        c = nc
        i = i+1
      end do
      PRINT*, 'la raiz es ', c
c     we could give the explicit interval, but the midpoint is good
c     enough for us 
      PRINT*, 'convergencia en ', i, ' iteraciones'
      stop
      end program 

      function f(x)
      implicit none
      real*8 f,x
      f = (x**2 - 2)*EXP(-x**2)

      return
      end
