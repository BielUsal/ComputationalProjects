c     Programa para derivar numéricamente la función f(x) en el punto x en función de h

      PROGRAM main
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      parameter (N=200)		! número de celdas
      parameter (x=1.0d0)		! punto en el que se evalúan las derivacas
      
      f0=cos(x)		! valor exacto de la función en el punto x
      fp0=-sin(x)		!  valor exacto de la derivada primera en el punto x
      fpp0=-cos(x)	!  valor exacto de la derivada segunda en el punto x
       
      open (1,FILE='error.dat')

      do i=1,N
      
        h=(1.0/3.0)**i			! h disminuye en un factor 3 cada vez

        fpa=(cos(x+h)-f0)/h						! f'(x), expresión adelantada
        fpc=(cos(x+h)-cos(x-h))/(2.d0*h)								! f'(x), expresión centrada
        fpp= (cos(x+h)-2.d0*cos(x)+cos(x-h))/(h**2)! f''(x)

        efpa=abs((fpa-fp0)/fp0)		! error relativo en f'(x), expresión adelantada
        efpb=abs((fpc-fp0)/fp0)! error relativo en f'(x), expresión centrada
        efpp=abs((fpp-fpp0)/fpp0)					! error relativo en f''(x)

        write(1,*) h,efpa,efpb,efpp	! escribimos los errores para cada valor del espaciado

      end do

      close(1)

      END
