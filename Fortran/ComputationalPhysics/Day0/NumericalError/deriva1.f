c     Programa para derivar numéricamente la función f(x) en el intervalo [a,b] 

      program MAIN
      IMPLICIT none
      integer::N,i
      parameter (N=200)! número de celdas				
      real*8::a,b,h,func,fpa,fpc,fpp,epa,epc,epp,x,f
      a=0.0d0! límite inferior
      b=10.0d0! límite superior
      h=(b-a)/(N-1)! espaciado

      

      open (1,file='derivadas.dat')
      open (2,file='errores.dat')      

      do i=0,N
        x=a+i*h
        func = f(x)
        fpa=(f(x+h)-f(x))/h! f'(x) avanzada
        fpc=(f(x+h)-f(x-h))/(2.d0*h) ! f'(x) centrada
        fpp= (f(x+h)-2.d0*f(x)+f(x-h))/(h**2.d0)!f''(x)
        epa=abs((fpa-(-sin(x))))			! error absoluto en f'(x) avanzada
        epc=abs((fpc-(-sin(x))))! error absoluto en f'(x) centrada
        epp=abs((fpp-(-cos(x))))! error absoluto en f''(x)
        write(1,*) x,fpa,fpc,fpp			! escribimos las derivadas
        write(2,*) x,epa,epc,epp			! escribimos los errores absolutos en las derivadas
      end do

      close(1)
      close(2)      

      end program

      function f(x)
      implicit none
      real*8::f,x 
      f = cos(x)
      end function