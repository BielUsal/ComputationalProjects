      !-----------------------------------------------------------------
      !Typing and initializations
      Program random 
      Implicit none
      integer N,Np, i, Ncor !Np is the number of points, 
      Parameter (Np=20001)
      Parameter (Ncor = 100)
      real*8 x(0:Np),t(0:Np)
      real*8 a,b,average,variance,skew,dt!,vector_read (Np),a
      real*8 R(-Ncor:Ncor)
      call srand(1333)
      !-----------------------------------------------------------------
      open(15,file="serie_temporal.txt")
      open(16,file="moments_autocor.txt")
      open(19,file="autocor.txt")
      !-----------------------------------------------------------------
      !This is where we do the actual sampling
      dt = 0.05d0
      do i =0, Np
        t(i) = -500+i*dt
      	x(i) = rand(0)!This samples from a ~U(0,1)
        write(15,*) t(i),x(i)
      end do
      close(15)
      call moments(x,average,variance,skew)
      write(16,*) average, variance, skew
      close(16)

      call Autocor(x,Ncor,average,R)
      do i = -Ncor,Ncor
        write(19,*) t(Np/2+i), R(i)
      end do 
      !-----------------------------------------------------------------
      !Now we do some histogramming. We have made a small little subroutine to do the dirty work so that we can easily use it later
      end
      !-----------------------------------------------------------------

      subroutine moments(x,average,variance,skew)
      implicit none
      integer i, Np
      Parameter (Np=20001)
      real*8 x(Np), average, variance, skew
      average = 0
      variance=0
      skew =0
      do i =1,Np
        average = average + x(i)
      end do
      average = average/Np
      do i =1,Np
        variance = variance + (x(i)-average)**2
        skew = skew + (x(i)-average)**3
      end do
      variance = variance/Np 
      skew = skew/(variance**(1.5d0)*Np)
      end subroutine

      subroutine Autocor(x,Ncor,average,R)
      integer Np, i,j, Ncor !nbin is the number of bins
      Parameter (Np=20001)
      real*8 x(Np), R(-Ncor:Ncor),average
      do j = -Ncor,Ncor
        R(j) = 0
        do i = Np/2-Ncor,Np/2+Ncor
          R(j) = R(j) + (x(i)-average)*(x(i+j)-average)
        end do
        R(j) = R(j)/((Np-2*Ncor)+1)
      end do
      end