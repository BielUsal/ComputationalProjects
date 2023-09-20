      Program BackAndForward
      Implicit none
      integer,parameter::n=4
      real*8,dimension(n,n)::A
      real*8,dimension(n)::b,x
      open(12,file="5.1a_inf.txt")
      call ReadMatrix(A,n)
      close(12)
      open(12,file="5.1b_inf.txt")
      call ReadVector(b,n)
      close(12)
      call solvetriang(a,b,x,n,1)
      print*, "Solucion del primero"
      call PrintVector(x,n)
      open(12,file="5.1a_sup.txt")
      call ReadMatrix(A,n)
      close(12)
      open(12,file="5.1b_sup.txt")
      call ReadVector(b,n)
      close(12)
      call solvetriang(a,b,x,n,-1)
      print*, "Solucion del segundo"
      call PrintVector(x,n)

      end program      

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

