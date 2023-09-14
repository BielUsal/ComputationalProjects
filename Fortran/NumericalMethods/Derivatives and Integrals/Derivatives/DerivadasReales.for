      program Derivatives
      implicit none
      real*8::a, h, f, fns, fs, fr, fexact, ens,es,er,pi
      integer::i
      real::b, s, bens,bes,ber,bfs,bf,bfns,bfr
      fexact= 0.01173348873386704 d0
      a=1.d0
      b=1.e0
       open(12,file='results2.txt')
       open(13,file='hnosim.txt')
       open(14,file='hsim.txt')
       open(15,file='hri.txt')
       open(16,file='realhnosim.txt')
       open(17,file='realhsim.txt')
       open(18,file='realhri.txt')
      do i=0,9
        h = 0.1d0 -0.01d0*i
        s = 0.1e0 - 0.01e0*i
        ens = 100*abs(fns(a,h)-fexact)/fexact
        es = 100*abs(fs(a,h)-fexact)/fexact
        er = 100*abs(fr(a,h)-fexact)/fexact
        bens = 100*abs(bfns(b,s)-fexact)/fexact
        bes = 100*abs(bfs(b,s)-fexact)/fexact
        ber = 100*abs(bfr(b,s)-fexact)/fexact
                
        write(12,99) s,bfns(b,s),bfs(b,s),bfr(b,s),bens,bes,ber
        write(13,*) h,Log(ens) 
        write(14,*) h,Log(es)
        write(15,*) h,Log(er)
        write(16,*) h,Log(bens)
        write(17,*) h,Log(bes)
        write(18,*) h,Log(ber)
      end do 
      
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

      function bfns(x,h)
      implicit none
      real::bf, bfns, x, h
      bfns=(bf(x+h)-bf(x))/h
      return
      end

      function bfs(x,h)
      implicit none
      real::bf, bfs, x, h
      bfs=(bf(x+h)-bf(x-h))/(2*h)
      return
      end

      function bfr(x,h)
      real::bf, bfr, x, h
      bfr = -(bf(x+h)-8*bf(x+h/2)-bf(x-h)+8*bf(x-h/2))/(6*h)
      return 
      end


      
      function bf(x)
      implicit none
      real bf,x, pi
      pi = 4.e0*atan(1.e0)
      bf = sin(pi*x/2)*(1-exp(-2*pi*x))
      return
      end