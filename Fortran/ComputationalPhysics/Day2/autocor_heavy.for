      !-----------------------------------------------------------------
      !Typing and initializations
      Program random 
      Implicit none
      integer N,Np, i, Ncor !Np is the number of points, 
      real*8 a,b,Ta,dt!,vector_read (Np),a
      Parameter(dt=0.1d-3)
      Parameter(Np=1d4)
      Parameter(Ncor=5d3)
      real*8 x(0:Np),t(0:Np)
      real*8 R(0:Ncor),test
      call srand(1333)
      !-----------------------------------------------------------------
      open(15,file="serie_temporal_heavy.txt")
      open(19,file="autocor_heavy.txt")
      !-----------------------------------------------------------------
      !This is where we do the actual sampling
      a = 5.0
      Ta = 0.2
      do i =0, Np
        t(i) = +i*dt
        if(t(i).lt.Ta) then
          x(i) = a
        else 
          x(i) = 0
        endif
        write(15,*) t(i),x(i)
      end do
      call Autocor(x,Ncor,R)
      do i = 0,Ncor
        write(19,*) t(i), R(i)
      end do 
      !-----------------------------------------------------------------
      !Now we do some histogramming. We have made a small little subroutine to do the dirty work so that we can easily use it later
      end
      !-----------------------------------------------------------------
      subroutine Autocor(x,Ncor,R)
      integer Np, i,j, Ncor !nbin is the number of bins
      Parameter(Np=1d4)
      real*8 x(0:Np), R(0:Ncor),average
      do j = 0,Ncor
        R(j) = 0
        do i = 0,Np-Ncor
          R(j) = R(j) + (x(i))*(x(i+j))
        end do
        R(j) = R(j)/((Np-Ncor)+1)
      end do
      end