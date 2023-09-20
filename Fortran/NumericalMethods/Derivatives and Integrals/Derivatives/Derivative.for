      program Derivatives
      implicit none
      real*8::a, h, f, fns, fs, fr, fexact, ens,es,er,pi
      integer::i
c      real::b, step, eens,ees,eer,
      fexact= 0.01173348873386704 d0
      a=1.d0
       open(12,file='results.txt')
       open(13,file='mathematica.txt')
      do i=0,9
        h = 0.1d0 -0.01d0*i
        ens = 100*abs(fns(a,h)-fexact)/fexact
        es = 100*abs(fs(a,h)-fexact)/fexact
        er = 100*abs(fr(a,h)-fexact)/fexact        
        write(12,99) h,fns(a,h),fs(a,h),fr(a,h),ens,es,er
        write(13,*) h, ens,es,er
      end do 
      print*, fs(a,0.1d0)
      
99    format(f5.3,2x,3(e12.6,1x),2x,3(e8.3,1x))
      end program 



      function fns(x,h)
      implicit none
      real*8::f, fns, x, h
      fns=(f(x+h)-f(x))/h
      return
      end

      function fs(x,h)
      implicit none
      real*8::f, fs, x, h
      fs=(f(x+h)-f(x-h))/(2*h)
      return
      end

      function fr(x,h)
      real*8::f, fr, x, h
      fr = -(f(x+h)-8*f(x+h/2)-f(x-h)+8*f(x-h/2))/(6*h)
      return 
      end


      
      function f(x)
      implicit none
      real*8 f,x, pi
      pi = 4.d0*datan(1.d0)
      f = dsin(pi*x/2)*(1-dexp(-2*pi*x))
      return
      end
