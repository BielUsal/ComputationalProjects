      program Bisection

      implicit none
      real*8 a, b, c, nc, error, tolerance, test, f, Secant
      integer i
      tolerance = 1 D-6
      a=0
      b=2
      c=0.
      nc =0.
      i=0
      error = 10.
      do while(error.GT.tolerance)
        nc = Secant(a,b)
        test = f(a)*f(nc)
        IF(test.GT.0) THEN
          a = nc
        ELSE
          b=nc
        END IF
        error = ABS(nc-c)
        c = nc
        i = i+1
      end do
      PRINT*, 'la raiz es ', c
      PRINT*, 'convergencia en ', i, ' iteraciones'
      stop
      end program 

      function Secant(a,b)
      implicit none
      real*8 Secant, a, f, b
      Secant = (f(b)*a-f(a)*b)/(f(b)-f(a))
      end
      
      function f(x)
      implicit none
      real*8 f,x
      f = (x**2 - 2)*EXP(-x**2)

      return
      end
