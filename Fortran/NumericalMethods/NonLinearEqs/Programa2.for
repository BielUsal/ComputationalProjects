      program Newton
c     Simple implementation. It is possible to 'switch over' to an
c     implimentation with finite differences if the (analytic)
c     derivative is not known.
      implicit none
      real*8 root,rootold,error,tolerance,f,fprimeanal,fprimeaprox
      real*8 NRmethod
      integer i, max,flag
      print*,'indica precision'
      read*,tolerance
      print*,'indica valor inicial'
      read*,root
      max = 1000000
c      root=1
      flag=0
      i=0
      rootold=0
c      tolerance = 1. D-6
      error = 10.
c     Initialization of variables
      DO WHILE(error.GT.tolerance.AND.flag.EQ.0)
c     Check for badness. If horrible, stops. 
        IF(ABS(fprimeanal(root)).LT.1d-4)THEN
          print*, 'la derivada es muy cerca de 0 en algun punto.'
          print*, 'no converge el metodo.'
          flag =-1
        END IF
c       Updates by Newton-Rapson. All functions defined below. 
        root = NRmethod(root,f(root),fprimeanal(root))
        error = ABS(root - rootold)
        rootold=root
        i = i+1
        IF(i.GT.max)THEN
          PRINT*,'no se ha encontrado raiz en 1000000 de iteraciones'
          flag = -1
        END IF
          
      end do
c     Checks to see we aren't stuck in an infinite loop. 
      IF(FLAG.EQ.0) THEN
        PRINT*, 'la raiz es ', root
        PRINT*, 'convergencia en ', i, ' iteraciones'
      ENDIF
      stop
      end program 

      function NRmethod(x,efe,efeprima)
      implicit none
      real*8 NRmethod, x, efe, efeprima
      NRmethod= x - efe/efeprima
      end



      function fprimeanal(x)
      implicit none
      real*8 fprimeanal, x
      fprimeanal = -2*EXP(-x**2)*x*(x**2-3)
      return
      end

c      function fprimeaprox(x)
c      implicit none
c      real*8 f, fprimeaprox, x, h
c      h = 10 D-6
c      fprimeaprox=(f(x+h)-f(x))/h
c      return
c      end
c      I was having some problems with the convergence of this, so I
c      implemented a finite-difference derivative to check if the
c      program worked on some level

      function f(x)
      implicit none
      real*8 f,x
      f = (x**2 - 2)*EXP(-x**2)
      return
      end
