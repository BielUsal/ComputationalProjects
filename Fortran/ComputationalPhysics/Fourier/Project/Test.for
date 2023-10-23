      !-----------------------------------------------------------------
      !Typing and initializations
      Program random 
      Implicit none
      integer N,i,Npf,Ncor 
      real*8 Tlim,dt,fmax,df,Tsamp
      Parameter(dt=10d-15)
      Parameter(N=1d6-1) ! Tlim/dt
      Parameter(Tlim = N*dt)
      Parameter(fmax = 5.d12)
      Parameter(df = 10.d9)
      Parameter(Npf = Int(fmax/df))
      Parameter(TSamp=5.0d-12)
      Parameter(Ncor=Int(TSamp/dt))
      real*8 v(0:N),t(0:N)
      real*8 Re(0:Npf), R(0:Ncor),average, avg
      !-----------------------------------------------------------------
      open(15,file="Velo_15_200k.dat")
      open(19,file="auto15200.txt")
      !-----------------------------------------------------------------
      !This is where we do the actual sampling
      do i =0,N
        read(15,*) v(i)
      enddo
      avg = average(v,N)
      write(17,*) avg
      call Autocor(v,N,Ncor,R,avg)
      do i =0,Ncor
        write(19,*) i*dt, R(i)
      enddo
      close(15)
      close(19)
      !-----------------------------------------------------------------
      end

      subroutine Fourier(x,Re)
        implicit none
        integer N,i,Npf,j,Ncor !Np is the number of points, 
        real*8 Tlim,dt,fmax,df,pi,Tsamp,test!,vector_read (Np),a
        Parameter(dt=10d-15)
        Parameter(N=1d6-1) ! Tlim/dt
        Parameter(df = 10.d9)
        Parameter(fmax = 5.d12)
        Parameter(Npf = Int(fmax/df))
        Parameter(TSamp=5.0d-12)
        Parameter(Ncor=Int(TSamp/dt))
        Parameter(pi=4*datan(1.d0))
        real*8 x(0:Ncor)
        real*8 Re(0:Npf)
        do i=0,Npf
            Re(i) = 0
        enddo
        do i=0,Npf !i*df represents our frequency
            do j = 0,Ncor !j*dt represents our time variable
                Re(i) = Re(i)+x(j)*dcos(2*pi*(i*df)*(j*dt))*dt
            enddo
            Re(i) = 4*Re(i)
        enddo
      end

      subroutine Autocor(v,N,Ncor,R,average)
      integer N, i,j, Ncor
      real*8 v(0:N), R(0:Ncor),average
      do j = 0,Ncor
        R(j) = 0
        do i = 0,N-Ncor
          R(j) = R(j) + (v(i)-average)*(v(i+j)-average)
        end do
        R(j) = R(j)/((N-Ncor)+1) !m^2/s^2. To get 10^14 cm^2/s^2, m^2/s^2 = 10^4 cm^2/s^2 = 10^-10 10^14 cm^2/s^2 
      end do
      end

      function average(x,Np)
      implicit none
      integer i, N,Np
      Parameter (N=1d6-1)
      real*8 x(N), average
      average = 0
      do i =1,Np
        average = average + x(i)
      end do
      average = average/Np
      end 
      !-----------------------------------------------------------------
