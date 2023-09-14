      program Secant
      implicit none
      real*8 a,b,error,tolerance,SecantUpdate,c,f,nc
      real*8 ainitial, binitial, test
      integer i, max, flag
      print*,'indica precision'
      read*,tolerance
      print*,'indica valor inicial 1'
      read*,ainitial
      print*,'indica valor inicial 1'
      read*,binitial
      a=ainitial
      b=binitial
      i=0
      flag=0
      max = 1000000
c     Initialization of some variables 
c      tolerance = 1. D-6
      error = 10.
      DO WHILE(error.GT.tolerance.AND.flag.EQ.0)
        IF(ABS(f(a)-f(b)).LT.1d-6)THEN
          PRINT*,'hay division por 0. Intentaremos regula falsi'
          flag=-1
        ELSE
        c = b
        b = SecantUpdate(a,b)
        a = c  
        END IF
        
        error = ABS(b-a)
        i = i+1
        IF(i.GT.max)THEN
          PRINT*,'no se ha encontrado raiz en 1000000 de iteraciones'
          flag = -1
        END IF
        END DO
c     If the badness flag is up, do regula falsi     
      IF(FLAG.EQ.0) THEN
        PRINT*, 'la raiz es ', b
        PRINT*, 'convergencia en ', i, ' iteraciones'
      ELSE
        a = ainitial
        b= binitial
        i=0
c     After initialization, regula falsi: Essentially combines the 
c     Bissection method with the secant and gives us one that is slower
c     but guarantees convergence
        do while(error.GT.tolerance)
          nc = SecantUpdate(a,b)
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
        END IF
      stop


      end program 

      function SecantUpdate(x,xold)
      implicit none
      real*8 SecantUpdate, x, f, xold
      SecantUpdate = x - f(x)*((x-xold)/(f(x)-f(xold)))
      end


      function f(x)
      implicit none
      real*8 f,x
      f = (x**2 - 2)*EXP(-x**2)
      return
      end
