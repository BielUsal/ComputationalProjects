      !-----------------------------------------------------------------
      !Typing and initializations
      Program random 
      Implicit none
      integer N,i,Npf !Np is the number of points, 
      real*8 Tlim,dt,fmax,df!,vector_read (Np),a
      Parameter(Tlim = 10.0d0)
      Parameter(dt=5.0d-3)
      Parameter(N=Int(Tlim/dt)) ! Tlim/dt
      Parameter(fmax = 5.d0)
      Parameter(df = 2.d-3)
      Parameter(Npf = Int(fmax/df))
      real*8 x(-N:N),t(-N:N),g
      real*8 Re(-Npf:Npf),Im(-Npf:Npf),Freq(-Npf:Npf)
      !-----------------------------------------------------------------
      open(18,file="time_series.txt")
      open(19,file="fourier.txt")
      !-----------------------------------------------------------------
      !This is where we do the actual sampling
      do i =-N,N
        t(i) = i*dt 
        x(i) = g(t(i),2.0d0)
      write(18,*) t(i), x(i)
      enddo

      call Fourier(x,Re,Im,Freq)
      
      do i = -Npf,Npf
        write(19,*) i*df, Re(i), Im(i), Freq(i) 
      enddo
      !-----------------------------------------------------------------
      end

      function g(t,alpha)
        implicit none
        real*8 g, t,alpha
        g=t**2*dexp(-alpha*abs(t))
      end

      subroutine Fourier(x,Re,Im,Freq)
        implicit none
        integer N,i,Npf,j !Np is the number of points, 
        real*8 Tlim,dt,fmax,df,pi!,vector_read (Np),a
        Parameter(Tlim = 10.0d0)
        Parameter(pi = 4*datan(1.0d0))
        Parameter(dt=5.0d-3)
        Parameter(N=Int(Tlim/dt)) ! Tlim/dt
        Parameter(fmax = 5.d0)
        Parameter(df = 2.d-3)
        Parameter(Npf = Int(fmax/df))
        real*8 x(-N:N),t(-N:N)
        real*8 Re(-Npf:Npf),Im(-Npf:Npf),Freq(-Npf:Npf)
        do i=-Npf,Npf
            Re(i) = 0
            Im(i) = 0
            Freq(i) = 0
        enddo
        do i=-Npf,Npf !i*df represents our frequency
            do j = -N,N !j*dt represents our time variable
                Re(i) = Re(i)+x(j)*cos(2*pi*(i*df)*(j*dt))*dt
                Im(i) = Im(i)-x(j)*sin(2*pi*(i*df)*(j*dt))*dt
            enddo
            Freq(i) = sqrt(Re(i)**2 +Im(i)**2) !This is |X(\omega)|
        enddo
      end


      !-----------------------------------------------------------------
