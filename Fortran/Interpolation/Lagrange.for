      program LagrangePoly
      implicit none
      real*8 x,y,datax,datay, l,prod,f,a,b
      integer i, n,j,k, m
      dimension datax(0:100), datay(0:100), l(0:100)
      dimension x(0:2000),y(0:2000)
      print*,"dar intervalo en forma de"
      print*, "a b"
      read*,a,b    
      print*, "dar numero de datos"
      read*, n
      print*,"dar numero de puntos en interpolacion"
      read*,m
c     Reads the data into the arrays datax and datay. Change to change 
c     file
      open(10, file = "data22.txt")
      do i=0,n
        read(10,*) datax(i), datay(i)
      end do
c   Uniform points between [0,4]
      do k=0,m
        x(k) = a + ((b-a)/m)*k
      end do
c     Iterate over x values
      do k=0,m
c     Calculates the 'l' coefficients of the polynomial. 
        do i=0,n
            prod = 1
            do j=0,n 
                if(j.EQ.i)THEN
                ELSE
                 prod = prod * (x(k)-datax(j))/(datax(i)-datax(j))
                END IF 
            END DO
            l(i) = prod
        END do
        y(k) = 0    
        do i=0,n
            y(k) = y(k) + l(i)*datay(i)
        end do 
      end do
c     Write our results to an output file. 
      open(11,file="poly22.txt")
      do k=0,m
        write(11,*) x(k), y(k)
      end do
        
      end program