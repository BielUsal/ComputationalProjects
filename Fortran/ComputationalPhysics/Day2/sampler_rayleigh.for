      !-----------------------------------------------------------------
      !Typing and initializations
      Program random 
      Implicit none
      integer N, i, nbin !nbin is the number of bins
      Parameter (N=1E5)
      Parameter (nbin=100)
      real*8 vector_gen (N),x(N),lambda,P(nbin)
      real*8 a,b,s,sigma, average, variance, skew,pi!,vector_read (N),a
      Parameter (pi = 4*datan(1.d0))
      call srand(1333)
      !-----------------------------------------------------------------
      !Open up the files
      open(15, file="serielambda1.txt")
      !open(16,file="serie.txt") !read vector example
      open(17,file="histogramalamb1.txt")
      open(18,file="moments_rayleigh.txt")
      !-----------------------------------------------------------------
      !This is where we do the actual sampling
      
      s= 2 !
      sigma = sqrt(s/2)
      do i =1, N
        vector_gen (i) = rand(0)!This samples from a ~U(0,1)
        x(i) = dsqrt(-dlog(vector_gen(i))*s) !This transforms it into an ~Rayleigh(sigma)
        write(15,*) i, vector_gen(i), x(i) !This writes it to a file
      end do
      close(15)
      !-----------------------------------------------------------------
      !Now we do some histogramming. We have made a small little subroutine to do the dirty work so that we can easily use it later
      a = 0.0d0
      b = maxval(x) !Min and max value of our representation
      do i =1,nbin
        P(i) = 0
      end do
      call Histogram(x,P,a,b)
      close(17)
      call moments(x,average,variance,skew)
      write(18,*) average,variance,skew
      write(18,*) sigma, (4-pi)*0.5d0*sigma**2,(2*sqrt(pi)*(pi-3))
     &/(4-pi)**1.5d0
      end

      !-----------------------------------------------------------------
      subroutine Histogram(vector,bin,a,b)
      integer i,j, N, nbin
      Parameter (N=1E5)
      Parameter (nbin=100)
      real*8 vector(N),bin(nbin),a,b,h,x !'vector' is the vector we want to histogram, bin represents the number in each bin
      h = (b-a)/DBLE(nbin)!This uses a,b and nbin to calculate the length of each bin
      do i =1,N
        j = int((vector(i)-a)/h)!'locates' each of the vectors by rounding them up
        bin(j) = bin(j)+1!adds 1 to the bin where the vector was located
      end do
      do i = 1,nbin
        bnorm = b/(N*h)!Normalizes the bins
        x = a+dble(i)*(b-a)/dble(nbin)!Reparametrizes us to x
        write(17,*) x, bnorm!Writes to a file
      end do
      end subroutine

      subroutine moments(x,average,variance,skew)
      implicit none
      integer i, N
      Parameter (N=1E5)
      real*8 x(N), average, variance, skew
      average = 0
      variance=0
      skew =0
      do i =1,N
        average = average + x(i)
      end do
      average = average/N
      do i =1,N
        variance = variance + (x(i)-average)**2
        skew = skew + (x(i)-average)**3
      end do
      variance = variance/N 
      skew = skew/(variance**(1.5d0)*N)
      end subroutine