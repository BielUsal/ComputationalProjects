      !-----------------------------------------------------------------
      !Typing and initializations
      Program random 
      Implicit none
      integer N, i, nbin !nbin is the number of bins
      Parameter (N=1E5)
      Parameter (nbin=100)
      real*8 vector_gen (N),x(N),lambda,P(nbin)
      real*8 a,b!,vector_read (N),a
      integer Iterations(N)
      call srand(1333)
      !-----------------------------------------------------------------
      !Open up the files
      open(15, file="serielambda1.txt")
      !open(16,file="serie.txt") !read vector example
      open(17,file="histogramalamb1.txt")
      !-----------------------------------------------------------------
      lambda = 1
      do i =1, N
      	Iterations(i) = i
      	vector_gen (i) = rand(0)
      	x(i) = -dlog(vector_gen(i))/lambda
      	write(15,*) Iterations(i), vector_gen(i), x(i)
      end do
      a = 0.00d0
      b = maxval(x)
      do i =1,nbin
      	P(i) = 0
      end do
      close(15)
      call Histogram(x,P,a,b)
      do i = 1,nbin
      	write(17,*) a+dble(i)*(b-a)/dble(nbin), P(i)/N
      end do

      !do i =1, N
      !	read(16,*) a, vector_read
      !end do
      !close(16)
      !write(*,*) vector_gen(15),vector_read(15)
      end

      subroutine Histogram(vector,bin,a,b)
      integer i,j, N, nbin
      Parameter (N=1E5)
      Parameter (nbin=100)
      real*8 vector(N),bin(nbin),a,b,h !V is the vector we want to histogram, bin represents the number in each bin
      h = (b-a)/DBLE(nbin)
      do i =1,N
      	j = int((vector(i)-a)/h)
      	bin(j) = bin(j)+1
      end do
      end subroutine