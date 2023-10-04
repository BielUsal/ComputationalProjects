      program Sampler
      implicit none
      integer N
      parameter (N=100)
      real*8::a,b,h,x,f
      integer::i
      a = 0.d0
      b = 10.d0
      h = (b-a)/N
      open(21,file="sample.dat")
      do i=0,N
        x = a + h*i
        write(21,*) x, i, f(x)
      end do 

      end program
     
      function f(x)
      implicit none
      real*8 f,x
      f = dcos(x) 
      end function     

