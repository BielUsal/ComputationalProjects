      !-----------------------------------------------------------------
      !By 
      !Gabriel D'Andrade Furlanetto XDD204950
      !Alvaro Gamarra Ralero 04231602Q
      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      !Typing and initializations
      Program random 
      Implicit none
      integer N,i,Npf,Ncor 
      real*8 Tlim,dt,fmax,df,Tsamp
      Parameter(dt=10d-15)
      Parameter(N=1d6) ! Tlim/dt
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
      open(20,file="PSD_VelFFT.txt")
      !-----------------------------------------------------------------
      !We now read the time-series we want to analyze

      do i =1,N
        read(15,*) v(i)
      enddo
      call FastishFourier(v,Freq)!This subroutine does a Fourier transform and returns us 
c     the (squared) amplitude. We do this over many regions in our time-domain and average over them
      do i = 1,Npf
        write(20,*) i*df, Freq(i)
      enddo
      close(15)
      close(20)
      end
      !-----------------------------------------------------------------
      !------------------Functions and Subroutines----------------------
      !-----------------------------------------------------------------

      subroutine FastishFourier(x,Freq2)
        implicit none
        integer N,i,Npf,j,Nsamp,Nb,k,dN,Ni,Nf!i and Np are on the frequency domain, j and N on the time-domain 
c       and k and Nb are refering to subdivisions 
        real*8 Tlim,dt,fmax,df,pi,avg,average
        Parameter(pi = 4*datan(1.0d0))
        Parameter(dt=10d-15)
        Parameter(N=1d6) ! Tlim/dt
        Parameter(Tlim = N*dt)
        Parameter(fmax = 5.d12)
        Parameter(df = 10.d9)
        Parameter(Npf = Int(fmax/df))
        Parameter(Nb=100)
        Parameter(dN=N/Nb)
        real*8 x(N),t(N)
        real*8 Re(Npf),Im(Npf),Freq2(Npf),Freq(Nb,Npf) !Freq is a matrix, first index selects the 'bin' we are in, and second is just the frequency.
c       It can be thought of as a vector of squared-amplitude vectors from the original problem.
        !---------------------------------------------------------------
        do k=1,Nb
          do i=1,Npf
            Freq(k,i) = 0 !Initialize values to 0. Probably not necessary but I was paranoid at this 
c           point because I was having a strange divergence phenomenon related to not setting variables to 0.
          enddo
        enddo
        !---------------------------------------------------------------
        !The actual Fourier integrals and sums
        do k=1,Nb !Loop over all the bins
          do i=1,Npf !Reinitialize all "sum" variables to 0
            Re(i) = 0
            Im(i) = 0
            Freq2(i) = 0
          enddo
          do i=1,Npf !Loop over all the frequencies
            Ni = 1+(k-1)*dN !Where we start on the loop
            Nf=k*dN !Where we end
            avg = average(x,Ni,Nf) !We use our function to calculate the average in the region
              do j =Ni,Nf !Integrate over time
                  Re(i) = Re(i)+(x(j)-avg)*cos(2*pi*(i*df)*(j*dt))*dt!\int (x(t)-xbar)*cos(2\pi\omega t) dt
                  Im(i) = Im(i)-(x(j)-avg)*sin(2*pi*(i*df)*(j*dt))*dt!-\int (x(t)-xbar)*sin(2\pi\omega t) dt
              enddo
              Freq(k,i) = 2*(Re(i)**2 +Im(i)**2)/(dN*dt) !This is |X(\omega)|^2
          enddo
        enddo
        !----------------------------------------------------------------
        !Averaging over all the bins
        do i=1,Npf !Loop over the frequencies
          Freq2(i) = 0
          do k =1,Nb !Sum over the bins
            Freq2(i) =Freq2(i) +Freq(k,i)
          enddo
          Freq2(i) = Freq2(i)/Nb !Divide by Nb to get the average
        enddo 
      end
      !-----------------------------------------------------------------
      function average(x,Ni,Nf) !This function does the average over (Ni,Nf)
      implicit none
      integer i, N,Ni,Nf
      Parameter (N=1d6-1)
      real*8 x(N), average
      average = 0
      do i =Ni,Nf
        average = average + x(i)
      end do
      average = average/(Ni-Nf)
      end 


      !-----------------------------------------------------------------
