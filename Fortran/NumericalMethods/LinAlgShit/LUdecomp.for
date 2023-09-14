      Program LUDecomp
      Implicit none
      integer,parameter::n=6
      real*8::det
      real*8,dimension(n,n)::A,L,U,test,Inv
      integer::i,j,k
      open(12,file="5.3a.txt")
      call ReadMatrix(A,n)
      close(12)
c     The above is initialization + reading the matrix from the file
c     this just saves some headache in defining the variable.
c      Instead of implementing it directly here, we made a subroutine 
c     for the LUer.
      call LUNoPivot(A,L,U,n)
      test = matmul(L,U)
      call PrintMatrix(test,n)
      print*, "L"
      call PrintMatrix(L,n)
      print*, "U"
      call PrintMatrix(U,n)
      call determinant(A,det,n)
      print*,"Determinant:",det
      call Inverse(A,Inv,n)
      call PrintMatrix(Inv,n)
      print*,"---------------------------------------------------------"
      call PrintMatrix(MatMul(A,Inv),n)

      end program
c    -------------------------------------------------------------------
c     Inverse
      subroutine Inverse(A,Inv,n)
      real*8,dimension(n,n)::A,L,U,Inv,Id, Z
      integer::i,j,k,n 
      call LUNoPivot(A,L,U,n)
      do i=1,n
        do j=1,n
          if(i.eq.j) then
            Id(i,j) = 1
          else 
            Id(i,j)=0
          end if
        end do
      end do
      do i=1,n
        call solvetriang(L,Id(n,:),Z(n,:),n,-1)
        call solvetriang(U,Z(n,:),Inv(n,:),n,1)
      end do
      call CleanMatrix(Inv,n)
      end subroutine
c    -------------------------------------------------------------------
c    LU decomp + determinant
      subroutine determinant(A,det,n)
      integer::i,n
      real*8,dimension(n,n)::A,L,U
      real*8::det
      call LUNoPivot(A,L,U,n)
      det =1
      do i=1,n
            det = det * U(i,i)
      end do
      end subroutine


      subroutine LUNoPivot(A,L,U,n)
      implicit none
      real*8::suma
      integer::i,j,k,n
      real*8,dimension(n,n)::A,L,U
      do i=1,n
            L(i,i) = 1
            U(1,i) = A(1,i)
      end do
      do k = 2, n
            do i = 1,k-1
                  suma = 0 
                  do j=1,i-1
                        suma = suma + L(k,j)*U(j,i)
                  end do
                  L(k,i)=(A(k,i)-suma)/(U(i,i))
            end do 
            do i =k,n
                  suma = 0 
                  do j=1,k-1
                        suma = suma + L(k,j)*U(j,i)
                  end do
            u(k,i)=a(k,i)-suma
      end do
      end do
      call CleanMatrix(L,n)
      call CleanMatrix(U,n)
      end subroutine

      subroutine CleanMatrix(A,n)
      implicit none
      integer::i,j,n
      real*8,dimension(n,n)::A
      do i =1, n
            do j=1,n 
                  if(abs(a(i,j)).lt.1.d-8) then
                        a(i,j)=0
                  end if
            end do
      end do
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

