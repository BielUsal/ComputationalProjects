      Program Coulomb
      Implicit none
      real*8:: k, sigma, hc, rho, rmin, rmax,mc,l,Eanal
      integer:: N, i, j, m,info,lwork,p
      Parameter (N=40)
      real*8:: H(N,N), B(N,N), rj(N),nu(N),eigvals(N),work(1+6*N+2*N**2)
      real*8::rsq(N,N),rms,rmsanal
      lwork = 1+6*N+2*N**2
      k = 0.52d0 !fm (?)
      sigma = 0.9255d0 !GeV fm^-1
      hc =0.19733d0 !GeV fm
      mc = 1.84d0 !GeV
      rmax = 110.0d0 !fm
      rmin = 0.1d0 !fm
      rho = (rmax/rmin)**(1.d0/(N-1))
      do j = 1,N
        rj(j) = rmin*(rho**(j-1))
        nu(j) = 1/(rj(j)**2)
      enddo
      write(*,*) "----------------Potencial Columbiano----------------"
      write(*,*) "----------------------Energias----------------------" 
      do j =0,2
        l = 1.d0*j
        call ConstructHamiltonian(nu,l,k,hc,0.0d0,mc,N,H,B,rsq)
        call DSYGV(1,'V','U',N,H,N,B,N,eigvals,work,lwork,info)
        write(*,*) " 			l =", l
        write(*,*) " Numerico			", "Analitico			", "Error"
        do i =1,5
          Eanal = -0.25d0*k**2 *mc/(i+l)**2
          write(*,*) eigvals(i), Eanal, abs(Eanal-eigvals(i))/abs(Eanal)
        enddo
      enddo
      !write(*,*) "---------------Radio Cuadratico Medio---------------"
      do m =0,2
        l = 1.d0*m
        call ConstructHamiltonian(nu,l,k,hc,0.0d0,mc,N,H,B,rsq)
        call DSYGV(1,'V','U',N,H,N,B,N,eigvals,work,lwork,info)
       ! write(*,*) " Numerico			", "Analitico			", "Error"
        rms = 0.0d0
      	do i = 1,40
      	  do j =1,40
            do p =1,40
               rms = rms + H(p,i)* rsq(i,j)*H(p,j)
           enddo
          enddo
        enddo
        rmsanal = 0.0d0
        write(*,*) rms
      enddo



      end program

      subroutine ConstructHamiltonian(nu,l,k,hc,sigma,mc,N,H,B,Rsq)
      Implicit none
      integer N,i, j
      real*8:: nu(N),l,k,hc,mc,sigma
      real*8::H(N,N),V(N,N), T(N,N), B(N,N), Rsq(N,N)

      do i = 1,N
      	do j = 1,N
          B(i,j) = ( 2*sqrt(nu(i)*nu(j))/( nu(i)+nu(j) ) )**(l+1.5d0)
          T(i,j) = (((2*l+3)*(hc)**2)/mc)*(2*nu(i)*nu(j))/(nu(i)+nu(j))
     &     *B(i,j)
          V(i,j) = -k*hc*Gamma(l+1.0d0)/Gamma(l+1.5d0)*B(i,j)
     &     *sqrt(nu(i)+nu(j))+sigma*(Gamma(l+2.0d0)/Gamma(l+1.5d0)
     &     *B(i,j)/sqrt(nu(i)+nu(j)))
          H(i,j) = T(i,j) + V(i,j)
          Rsq(i,j) = Gamma(l+2.5d0)/Gamma(l+1.5d0)*B(i,j)/(nu(i)+nu(j))!matrix element to calculate the rms radius    
        enddo
      enddo
      end subroutine 

