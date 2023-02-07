      program DiffEq
      Implicit none
      real*8::x,y,h,yexact
      dimension::x(0:100),y(0:100),yexact(0:100)
      integer::n,i
      n=21
c     Setting up initial values
      x(0)=0.d0
      x(n)=2.0d0
      y(0)=1.0d0
c     Since all problems use the same x, here we just set it up
      h = (x(n)-x(0))/n
      do i=0,n
        x(i) = x(0)+h*i
        yexact(i)=dexp(-x(i)**2)
      end do

      
      call Euler(x,y,n,h,yexact)
      call Huen(x,y,n,h,yexact)
      call Midpoint(x,y,n,h,yexact)
      call RK4(x,y,n,h,yexact)
      call PredCor(x,y,n,h,yexact)
      end program

      subroutine Euler(x,y,n,h,yexact)
      implicit none
      real*8::x,y,h,f,yexact
      dimension::x(0:100),y(0:100),yexact(0:100)
      integer::i,n
      do i=0,n-1
        y(i+1)=y(i)+h*f(x(i),y(i))
      end do
c     Now we write them to a file for mathematica
      open(12,file="euler.txt")
      do i=0,n
        write(12,*) x(i),y(i)
      end do
      close(12)
      open(13,file="eulererror.txt")
      do i=1,n
        write(13,*) x(i), abs((y(i)-yexact(i)))/yexact(i)
      end do
      close(13)
      end subroutine

      subroutine Huen(x,y,n,h,yexact)
      implicit none
      real*8::x,y,h,f,yexact
      dimension::x(0:100),y(0:100),yexact(0:100)
      integer::i,n
      do i=0,n-1
        y(i+1)=y(i)+0.5*h*(f(x(i),y(i))+f(x(i+1),y(i)+h*f(x(i),y(i))))
      end do
c     Now we write them to a file for mathematica
      open(12,file="huen.txt")
      do i=0,n
        write(12,*) x(i),y(i)
      end do
      close(12)
      open(13,file="huenerror.txt")
      do i=1,n
        write(13,*) x(i), abs(y(i)-yexact(i))/yexact(i)
      end do
      close(13)
      end subroutine
      
      subroutine MidPoint(x,y,n,h,yexact)
      implicit none
      real*8::x,y,h,f,yexact
      dimension::x(0:100),y(0:100),yexact(0:100)
      integer::i,n
      do i=0,n-1
        y(i+1)=y(i)+h*f(x(i)+0.5*h,y(i)+0.5*h*f(x(i),y(i)))
      end do
c     Now we write them to a file for mathematica
      open(12,file="midpoint.txt")
      do i=0,n
        write(12,*) x(i),y(i)
      end do
      close(12)

      open(13,file="midpointerror.txt")
      do i=1,n
        write(13,*) x(i), abs(y(i)-yexact(i))/yexact(i)
      end do
      close(13)
      end subroutine
      

      subroutine RK4(x,y,n,h,yexact)
      implicit none
      real*8::x,y,h,f,a,b,c,d,yexact
      dimension::x(0:100),y(0:100),yexact(0:100)
      integer::n,i
      do i=0,n-1
      a=f(x(i),y(i))
      b=f(x(i)+0.5*h,y(i)+0.5*h*a)
      c=f(x(i)+0.5*h,y(i)+0.5*h*b)
      d=f(x(i)+h,y(i)+h*c)
        y(i+1)=y(i)+(h/6)*(a+2*b+2*c+d)
      end do
c     Now we write them to a file for mathematica
      open(12,file="rk4.txt")
      do i=0,n
        write(12,*) x(i),y(i)
      end do
      close(12)
      open(13,file="rk4error.txt")
      do i=1,n
        write(13,*) x(i), abs(y(i)-yexact(i))/yexact(i)
      end do
      close(13)
      end subroutine

      subroutine PredCor(x,y,n,h,yexact)
      implicit none
      real*8::x,y,h,f,yexact, ypred
      dimension::x(0:100),y(0:100),yexact(0:100)
      integer::n,i
      do i=0,2
        y(i)=Exp(-x(i)**2)
      end do
      do i=2,n-1
        ypred = y(i-2) + 0.75*h*(f(x(i-2),y(i-2))+3*f(x(i),y(i)))
        y(i+1)=y(i-2)+0.375*h*(f(x(i-2),y(i-2))+3*f(x(i-1),y(i-1))
     c  +3*f(x(i),y(i))+f(x(i+1),ypred))
      end do
c     Now we write them to a file for mathematica
      open(12,file="predcor.txt")
      do i=0,n
        write(12,*) x(i),y(i)
      end do
      close(12)
      open(13,file="predcorerror.txt")
      do i=3,n
        write(13,*) x(i), abs(y(i)-yexact(i))/yexact(i)
      end do
      close(13)
      
      end subroutine
      
      function f(x,y)
      implicit none
      real*8::x,y,f
      f = -2.d0 * x * y
      end function