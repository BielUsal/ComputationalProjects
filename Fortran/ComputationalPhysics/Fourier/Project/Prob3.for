      !-----------------------------------------------------------------
      !By 
      !Gabriel D'Andrade Furlanetto XDD204950
      !Alvaro Gamarra Ralero
      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      !Typing and initializations. Nothing very interesting other than parameters
      Program random 
      Implicit none
      integer N,i,Npf,Ncor 
      real*8 Tlim,dt,fmax,df,Tsamp
      Parameter(dt=10d-15)
      Parameter(N=1d6) 
      Parameter(Tlim = N*dt)
      Parameter(fmax = 5.d12)
      Parameter(df = 10.d9)
      Parameter(Npf = Int(fmax/df))
      Parameter(TSamp=5.0d-12)
      Parameter(Ncor=Int(TSamp/dt))
      real*8 v(N),t(N)
      real*8 Freq(Npf), R(Ncor),average, avg
      !-----------------------------------------------------------------
      open(15,file="Velocidades.dat")
      open(20,file="PSD_Vel6.txt")
      !-----------------------------------------------------------------
      !We now read the time-series we want to analyze
      do i =1,N
        read(15,*) v(i)
      enddo
      avg = average(v,Ncor)!We need the average to get the PSD, so we have a little function fetch that for us
      call Fourier(v,Freq,avg) !This subroutine does the Fourier transform and returns us 
c     the (squared) amplitude, which for this problem corresponds to the PSD
      do i = 1,Npf
        write(20,*) i*df, Freq(i) !This just writes the result
      enddo
      end
      !-----------------------------------------------------------------
      !------------------Functions and Subroutines----------------------
      !-----------------------------------------------------------------
      subroutine Fourier(x,Freq2,avg)
        implicit none
        integer N,i,Npf,j,Nsamp !For notation sake, i and Npf always refer to things in the frequency domain, N and j in the time domain   
        real*8 Tlim,dt,fmax,df,pi,avg
        Parameter(pi = 4*datan(1.0d0))
        Parameter(dt=10d-15)
        Parameter(N=1d6) ! Tlim/dt
        Parameter(Tlim = N*dt)
        Parameter(fmax = 5.d12)
        Parameter(df = 10.d9)
        Parameter(Npf = Int(fmax/df))
        real*8 x(N),t(N)
        real*8 Re(Npf),Im(Npf),Freq2(Npf)
        Nsamp=1D6 !Since we want to sample from both 1E4 and 1E6 points, this just selects that
        do i=1,Npf !This just initializes the variables to 0. Thought Fortran did this automatically, but I was having problems with this 
            Re(i) = 0
            Im(i) = 0
            Freq2(i) = 0
        enddo
        do i=1,Npf !Looping over the frequencies to define our vector
            do j =1,Nsamp !Integrating in time
                Re(i) = Re(i)+(x(j)-avg)*cos(2*pi*(i*df)*(j*dt))*dt !\int (x(t)-xbar)*cos(2\pi\omega t) dt
                Im(i) = Im(i)-(x(j)-avg)*sin(2*pi*(i*df)*(j*dt))*dt !-\int (x(t)-xbar)*sin(2\pi\omega t) dt
            enddo
            Freq2(i) = 2*(Re(i)**2 +Im(i)**2)/(Nsamp*dt) !This is |X(\omega)|^2
        enddo
      end

      function average(x,Np) !This just averages over the first Np*dt s, not that complicated
      implicit none
      integer i, N,Np
      Parameter (N=1d6)
      real*8 x(N), average
      average = 0
      do i =1,Np
        average = average + x(i)
      end do
      average = average/Np
      end 


      !-----------------------------------------------------------------
