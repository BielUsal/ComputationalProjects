      Program Gaussian
      Implicit none
      integer,parameter::n=4
      real*8,dimension(n,n)::A,Aorig
      real*8,dimension(n)::b,x
      open(12,file="5.2a.txt")
      call ReadMatrix(A,n)
      close(12)
      open(12,file="5.2b.txt")
      call ReadVector(b,n)
      close(12)
      Aorig = A
      call PrintMatrix(A,n)
      print*,""
      call PrintVector(b,n)
      print*,""
      print*,"---------------------------------------------------------"
      call GaussPivot(A,b,x,n)
      call PrintMatrix(A,n)
      print*,""
      call PrintVector(b,n)
      call solvetriang(A,b,x,n,-1)
      print*,""
      call PrintVector(x,n)
      print*,""
      call PrintVector(MatMul(Aorig,x),n)
      end program


c     Gaussian elimination. With pivoting. This makes your matrix into an
c     upper triangular matrix.
      subroutine GaussPivot(a,b,x,n)
      implicit none
      integer::i,j,k,n,piv
      real*8,dimension(n,n)::A
      real*8,dimension(n)::b,x
      real*8::temp,f
c     As suggested in the slides, we look for the biggest 'a(k,k)'' we 
c     can find.
      do k = 1, n-1
c     We initialize our pivot variable and run through all the values 
c     in the column to check if there is a bigger one            
        piv = k
        do i = k+1, n
            if (abs(a(i,k)).gt.abs(a(piv,k))) then 
                  piv = i
            end if
        end do
c     Now, piv has the biggest "a(k,k)", now we do the good'ol 
c     switchearoo   
        if (piv.ne.k) then
            do j = k, n
                temp = a(piv,j)
                a(piv,j) = a(k,j)
                a(k,j) = temp
            end do
            temp = b(piv)
            b(piv) = b(k)
            b(k) = temp
        end if
c     This is just normal gaussian elimination now.   
        do i = k+1, n
            f = a(i,k)/a(k,k)
            do j = k, n
                a(i,j) = a(i,j) -  a(k,j)*f
            enddo
            b(i) = b(i) - b(k)*f
          enddo
        enddo
      end subroutine

c     -------------------------------------------------------------------
c     Forward and backward substitution stuff
      
      subroutine solvetriang(a,b,x,n,direction)
      real*8,dimension(n,n)::A
      real*8,dimension(n)::b,x
      integer::n,direction
      select case(direction)
        case(1)
        call forsub(a,b,x,n)
        case(-1)
        call backsub(a,b,x,n)
        case default
        print*, "not gonna work"
      end select


      end subroutine
      subroutine backsub(a,b,x,n)
      implicit none
      real*8,dimension(n,n)::A
      real*8,dimension(n)::b,x
      integer::i,j,n

      x(n) = b(n) / a(n,n)
      do i = n-1, 1, -1
        x(i) = b(i)
        do j = i+1, n
            x(i) = x(i) - a(i,j) * x(j)
        end do
        x(i) = x(i) / a(i,i)
      end do

      end subroutine

      subroutine forsub(a,b,x,n)
      implicit none
      real*8,dimension(n,n)::A
      real*8,dimension(n)::b,x
      integer::i,j,n
      x(1)=b(1)/a(1,1)
      do i = 2, n
        x(i) = b(i)
        do j = 1, i-1
            x(i) = x(i) - a(i,j) * x(j)
        enddo
        x(i) = x(i) / a(i,i)
      enddo
      end subroutine

c     ------------------------------------------------------------------
c     Stuff to read/write matrices
      subroutine ReadVector(x,n)
      implicit none
      integer::i,n
      real*8,dimension(n)::x
      do i =1,n
            read(12,*) x(i)
      end do 
      end subroutine

      subroutine PrintVector(x,n)
      implicit none
      integer::i,n
      real*8,dimension(n)::x
      do i =1,n
            write(*,*) x(i)
      end do 
      end subroutine 

      subroutine PrintMatrix(A,n)
      implicit none
      integer::i,j,n
      real*8,dimension(n,n)::A
      do i =1,n
            write(*,*) (A(i,j),j=1,n)
      end do 
      end subroutine 

      subroutine ReadMatrix(A,n)
      implicit none
      integer::i,j
      integer::n
      real*8,dimension(n,n)::A
      do i =1,n
            read(12,*) (A(i,j),j=1,n)
      end do 
      end subroutine 

