      program SysDiffEq
      implicit none
      real*8::t,y,h
      dimension::t(0:1000),y(1:2,0:1000)
      integer:: i, n
      end program
      



      subroutine RK4(t,y,h,alpha,y0,v0,tf,n)
      implicit none
      real*8::t,y,h,f1,f2,K,alpha,y0,v0,tf
      dimension::t(0:1000),y(1:2,0:1000),K(1:2,1:4)
      integer::n,i,m
      n=tf/h
      do i=0,n
        t(i)= h * i
      end do
      y(1,0) = y0
      y(2,0)= v0
c     Note: \ddot{y} = -sin(y) - alpha*\dot{y}, i.e. \dot{y} = v and \dot{v} = -sin(y)- alpha *v
c    i.e, our diff equation will be: f1 = y2,f2 = -sin(y1) - alpha y2
      do i=0,n-1
        
      K(1,1)=f1(t(i),y(1,i),y(2,i),alpha)
      K(2,1)=f2(t(i),y(1,i),y(2,i),alpha)
      K(1,2) = f1(t(i)+0.5*h,y(1,i)+0.5*h*K(1,1),y(2,i)+0.5*h*K(2,1)
     c ,alpha)
      K(2,2) = f2(t(i)+0.5*h,y(1,i)+0.5*h*K(1,1),y(2,i)+0.5*h*K(2,1)
     c,alpha)
      K(1,3) = f1(t(i)+0.5*h,y(1,i)+0.5*h*K(1,2),y(2,i)+0.5*h*K(2,2)
      
     c,alpha)
      K(2,3) = f2(t(i)+0.5*h,y(1,i)+0.5*h*K(1,2),y(2,i)+0.5*h*K(2,2)
     c,alpha)
      K(1,4) = f1(t(i)+0.5*h,y(1,i)+h*K(1,3),y(2,i)+h*K(2,3),alpha)
      K(2,4) = f2(t(i)+0.5*h,y(1,i)+h*K(1,3),y(2,i)+h*K(2,3),alpha)
      y(1,i+1) = y(1,i) +(h/6)*(K(1,1)+2*K(1,2)+2*K(1,3)+K(1,4)) 
      y(2,i+1) = y(2,i) +(h/6)*(K(2,1)+2*K(2,2)+2*K(2,3)+K(2,4)) 
      end do
      end subroutine
      
      function f1(t,y,v,alpha)
      implicit none
      real*8::t,y,v,alpha,f1
      f1 = v
      end function
      
      function f2(t,y,v,alpha)
      implicit none
      real*8::t,y,v,f2,alpha
      f2 = -sin(y)- alpha*v
      end function