      program Integrals
      implicit none
      real*8::a, b, f, integralexact
      real*8::intpm,inttrap,intsimp, intnc5
      real*8::intGL2,integralGL10,fgauss
      real*8:: epm, etrap, esimp, enc5,egl2,egL10
      a = 2.d0
      b=5.d0
      open(12,file='gl10.txt')
      integralGL10=0
      integralexact= 0.2224732644915359d0
c       open(12,file='results.txt')
      call Intgl10(a,b,integralGl10)
      epm = 100*abs(intpm(a,b)-integralexact)/integralexact
      etrap = 100*abs(inttrap(a,b)-integralexact)/integralexact
      esimp = 100*abs(intsimp(a,b)-integralexact)/integralexact
      enc5 = 100*abs(intnc5(a,b)-integralexact)/integralexact
      egl2 = 100*abs(intgl2(a,b)-integralexact)/integralexact
      egl10 = 100*abs(IntegralGL10-integralexact)/integralexact
      print*, intpm(a,b), epm  
      print*, inttrap(a,b), etrap
      print*, intsimp(a,b),esimp              
      print*,intnc5(a,b),enc5
      print*,intgl2(a,b),egl2
      print*,integralGl10,egl10
c        write(12,99) epm,
      end program 

      function intpm(a,b)
      implicit none
      real*8::intpm,a,b,f
      intpm = (b-a)*f((a+b)/2)
      return
      end function

      function inttrap(a,b)
      implicit none
      real*8::inttrap, a,b,f
      inttrap = (b-a)*(f(a)+f(b))/2
      return
      end function

      function intsimp(a,b)
      implicit none
      real*8::intsimp,a,b,f
      intsimp = (b-a)*(f(a)+4*f((a+b)/2)+f(b))/6
      return
      end function

      function intnc5(a,b)
      implicit none
      real*8::intnc5,a,b,f,h
      h = b-a
      intnc5 = (b-a)*(3*f(a)/90 + 16*f(a+h/4)/45 + 2*f(a+h/2)/15
     c+ 16*f(a+(3*h/4))/45+ 7*f(b)/90)
      return
      end function

      function intGL2(a,b)
      implicit none
      real*8::intGL2,a,b,fgauss
      intgl2 = 0.5*(b-a)*(1*fgauss(-0.577350269189626d0,a,b)
     c +1*fgauss(0.577350269189626d0,a,b))
      return
      end function
      
      subroutine intGL10(a,b,integral)
      implicit none
      real*8::integral,a,b,fgauss, alpha, omega, sum
      dimension::alpha(0:100), omega(0:100)
      integer::i
      do i=0,8
      read(12,*) alpha(i), omega(i)
      end do
      sum =0
      do i =0,8
        sum = sum + fgauss(alpha(i),a,b)*omega(i)
        end do
      integral = 0.5*(b-a)*sum
      return
      end

      function fgauss(x,a,b)
      implicit none
      real*8::f,fgauss,z,a,b,x
      z=0.5*(b+a+(b-a)*x)
      fgauss = f(z)
      return
      end function

      
      function f(x)
      implicit none
      real*8 f,x
      f = (1-cos(x))*exp(-x)
      return
      end
