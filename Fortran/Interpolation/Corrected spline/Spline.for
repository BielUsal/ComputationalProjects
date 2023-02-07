      program splines
      implicit none 
      
      real*8 :: x, f, d, a, b, c, xt,S,y
      real*8:: left,right
      integer :: n,i,flag, p,m,k
      dimension:: x(0:100), f(0:100), d(0:100), a(0:100)
      dimension:: b(0:100), c(0:100), xt(0:100)
      dimension::y(0:100)
      xt = 1
      flag = 0
      print*,"dar intervalo en forma de"
      print*, "a b"
      read*,left,right
      print*, "dar numero de datos"
      read*, n
      print*,"dar numero de puntos en interpolacion"
      read*,m
c     Reads the data into the arrays datax and datay 
      open(10, file = "data21.txt")
      do i=0,n
        read(10,*) x(i), f(i)
      end do
c    Calculate the spline coefficients !!!!!!  esto mejor en la subrutina     
      call Splinefy(x,f,n,a,b,c,d)
c      Uniform points.
      do k=0,m
        xt(k) = left + ((right-left)/m)*k
      end do
      do k=0,m
c       Locating where our point is
        i=0
        flag=0
        do while(flag.eq.0)
          if(xt(k).ge.x(i).and.xt(k).lt.x(i+1)) then
            p = i
            flag = -1
          end if
          i = i+1  
        end do
      s = a(p)*(xt(k)-x(p))**3 + b(p) * (xt(k)-x(p))**2
     c   + c(p) *(xt(k)-x(p))+ d(p)
      y(k) = s
      end do
      open(11,file="spline21.txt")
      do k=0,m
        write(11,*) xt(k), y(k)
      end do
      end program


      
      subroutine Splinefy(x,f,n,a,b,c,d)
      real*8:: x,f,a,b,c,d
      real*8::h,alpha,beta,gam,r
      dimension:: x(0:100), f(0:100), h(0:100), d(0:100), a(0:100)
      dimension:: b(0:100), c(0:100), alpha(0:100), beta(0:100)
      dimension:: gam(0:100), r(0:100)
      integer::n,i
      b(0) = 0
      b(n) = 0 
c      Define h
      do i=0,n-1
        h(i)=x(i+1)-x(i)
      end do
c     Define d      
      do i=0,n
        d(i)=f(i)
      end do
c     Define alpha, beta, gamma and r
      do i=2,n-1
        alpha(i)=h(i-1)
      end do
      do i=1,n-1
        beta(i)=2*(h(i)+h(i-1))
        r(i) = 3*(d(i+1)-d(i))/h(i) - 3*(d(i)-d(i-1))/h(i-1)
      end do
      do i=1,n-2
        gam(i)=h(i)
      end do

c    Gaussianly eliminate:
      do i=2,n-1
        beta(i) = beta(i)-gam(i-1)*alpha(i-1)/beta(i-1)
        r(i) = r(i)-r(i-1)*alpha(i)/beta(i-1)
      enddo
      b(n-1)=r(n-1)/beta(n-1)
c   Calculate a,b,c,d    
      do i=n-2,1,-1
        b(i)=(r(i)-gam(i)*b(i+1))/beta(i)
      end do
      do i=0,n-1
        a(i)=(b(i+1)-b(i))/(3*h(i))
      end do
      do i=0,n-1
        c(i)=(d(i+1)-d(i))/h(i) - (b(i+1)+2*b(i))*h(i)/3
      end do
c      do i=0,n-1
c        write(12,*) a(i),b(i),c(i),d(i)
c      end do   
      return 
      
      end subroutine

c      function EvalSpline(a,b,c,d,x,xlocal)
c      real*8::EvalSpline,a,b,c,d,x
c      EvalSpline = a*(x-xlocal)**3 + b * (x-xlocal)**2 + c*(x-xlocal)+ d
c      end function