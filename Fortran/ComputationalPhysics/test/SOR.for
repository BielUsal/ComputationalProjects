      !-----------------------------------------------------------------
      !By 
      !Gabriel D'Andrade Furlanetto 71080411S
      !Alvaro Gamarra Ralero 04231602Q
      !-----------------------------------------------------------------
      !Typing and initializations
      Program Poisson
      Implicit none
      integer N, i,j,k
      Parameter (N=25)
      real*8 xi,xf,yi,yf,h,error,omega
      real*8 rho(-2*N:2*N,-N:N),u(-2*N:2*N,-N:N)
      real*8 Ex(-2*N+1:2*N-1,-N+1:N-1),Ey(-2*N+1:2*N-1,-N+1:N-1)
      omega=1.9d0
      xi = -10.d0
      xf = 10.d0
      yi = -5.d0
      yf= 5.d0
      h = (xf-xi)/(4*N)
      open(20,file="V-SOR.dat")
      open(21,file="E-SOR.dat")
      open(22,file="error-SOR-1.8.dat")
      !initializing matrix to 0
      do i=-2*N,2*N
        do j=-N,N
          u(i,j)=0.d0
        enddo
      enddo
      !Boundary Conditions
      do i=-2*N,2*N
        u(i,N)=5.d0
        u(i,-N) = 0.d0
      enddo  
      do j=-N,N
        u(-2*N,j)=-0.d0
        u(2*N,j)=0.d0
      enddo
      !initializing rho
      do i=-2*N,2*N
        do j=-N,N
          rho(i,j)=0.d0
        enddo
      enddo
      k=1
      error=100
      do while (error>1d-5.and.k<100000)
        call PoissonSOR(u,rho,omega,h,N,error)
        write(22,*) k, error
        k=k+1
      enddo
      !Writing V(x,y)
      do i=-2*N,2*N
        do j=-N,N
          write(20,*) i*h,j*h,u(i,j)
        enddo
      enddo
      !Defining Ex(x,y)=-\partial_X V(x,y) and Ey(x,y)=-\partial_Y V(x,y)
      do i=-2*N+1,2*N-1
        do j=-N+1,N-1
          Ex(i,j) = -(u(i+1,j)-u(i-1,j))/(2.d0*h)
          Ey(i,j) = -(u(i,j+1)-u(i,j-1))/(2.d0*h)
        enddo
      enddo
      !write Ex(x,y) and Ey(x,y)
      do i=-2*N+1,2*N-1,2
        do j=-N+1,N-1,2
          write(21,*) i*h,j*h,Ex(i,j),Ey(i,j)
        enddo
      enddo

      end Program

      subroutine PoissonSOR(u,rho,omega,h,N,error)
      implicit none
      integer N, i,j
      real*8 rho(-2*N:2*N,-N:N),u(-2*N:2*N,-N:N),h,r(-2*N:2*N,-N:N)
      real*8 error,omega,uold(-2*N:2*N,-N:N)
      !Initialize uold
      do i=-2*N,2*N
        do j=-N,N
          uold(i,j)=u(i,j)
        enddo
      enddo
      !Over-relaxation
      do i=-2*N+1,2*N-1
        do j=-N+1,N-1
          u(i,j)=omega*0.25d0*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)
     &     +rho(i,j)*h**2) + (1.d0-omega)*u(i,j)
        enddo
      enddo
      !Enforcing v-N boundary conditions
c      do i=-2*N,2*N
c        u(i,-N) = u(i,-N+1)
c      enddo  
c      do j=-N,N
c        u(2*N,j)=u(2*N-1,j)
c      enddo
      !error calculation using the maximum norm
      error=0.d0
      do i=-2*N,2*N
        do j=-N,N
          error=max(error,abs(u(i,j)-uold(i,j)))
          uold(i,j)=u(i,j)
        enddo
      enddo


      end subroutine

