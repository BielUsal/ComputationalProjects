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
      real*8 xi,xf,yi,yf,h,error
      real*8 rho(-2*N:2*N,-N:N),u(-2*N:2*N,-N:N),unew(-2*N:2*N,-N:N)
      real*8 Ex(-2*N+1:2*N-1,-N+1:N-1),Ey(-2*N+1:2*N-1,-N+1:N-1)
      xi = -10.d0
      xf = 10.d0
      yi = -5.d0
      yf= 5.d0
      h = (xf-xi)/(2*N)
      open(20,file="V-GS.dat")
      open(21,file="E-GS.dat")
      open(22,file="error-GS.dat")
      !initializing matrix to 0
      do i=-2*N,2*N
        do j=-N,N
          u(i,j)=0.d0
        enddo
      enddo
      !Boundary Conditions
      do i=-2*N,2*N
        u(i,N)=5.d0
        u(i,-N) = u(i,-N+1)
      enddo  
      do j=-N,N
        u(-2*N,j)=-1.d0 
        u(2*N,j)=u(2*N-1,j)
      enddo
      !initializing rho
      do i=-2*N,2*N
        do j=-N,N
          rho(i,j)=0.d0
        enddo
      enddo
      !initializing unew
      do i=-2*N,2*N
        do j=-N,N
          unew(i,j)=u(i,j)
        enddo
      enddo
      i=1
      error=100
      do while (error>1d-5.and.k<10000)
        call PoissonGS(u,rho,h,N,error)
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
      do i=-2*N+1,2*N-1
        do j=-N+1,N-1
          write(21,*) i*h,j*h,Ex(i,j),Ey(i,j)
        enddo
      enddo

      end Program

      subroutine PoissonGS(u,rho,h,N,error)
      implicit none
      integer N, i,j
      real*8 rho(-2*N:2*N,-N:N),u(-2*N:2*N,-N:N),uold(-2*N:2*N,-N:N),h
      real*8 error
      !initializing uold
      do i=-2*N,2*N
        do j=-N,N
          uold(i,j)=u(i,j)
        enddo
      enddo
      !Gauss-Seidel evolution step
      do i=-2*N+1,2*N-1
        do j=-N+1,N-1
          u(i,j) = (u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)+h**2*rho(i,j))/4.d0 !\nabla^2 p(x,y) = \nabla \cdot C
        enddo
      enddo
      !Enforcing v-N boundary conditions
      do i=-2*N,2*N
        u(i,-N) = u(i,-N+1)
      enddo  
      do j=-N,N
        u(2*N,j)=u(2*N-1,j)
      enddo
      !Using maximum norm, check the error
      error=0.d0
      do i=-2*N,2*N
        do j=-N,N
          error=max(error,abs(u(i,j)-uold(i,j)))
          uold(i,j)=u(i,j)
        enddo
      enddo

      end subroutine

