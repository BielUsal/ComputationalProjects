      Program Coulomb
      Implicit none
      real*8:: k, sigma, hc, rho, rmin, rmax,mc,l,Eanal,M0,pi,alpha
      integer:: N, i, j, m,info,lwork,p
      Parameter (N=50)
      real*8:: H(N,N), B(N,N), rj(N),nu(N),eigvals(N),work(1+6*N+2*N**2)
      real*8::rsq(N,N),rms(N),rmsanal,GammaFac(N)
      real*8:: Mn(N),C(N),psi0(N),Gam(N)
      Parameter(k = 0.52d0) !fm (?)
      Parameter(alpha=1.d0/137.036d0)
      Parameter(sigma = 0.9255d0) !GeV fm^-1
      Parameter(hc =0.19733d0) !GeV fm
      Parameter(mc = 1.84d0) !GeV
      Parameter(rmax = 70.0d0) !fm
      Parameter(rmin = 0.1d0) !fm
      Parameter(M0 =2.83747) !Gev
      Parameter(Pi = 4.d0*datan(1.d0))
      lwork = 1+6*N+2*N**2 !I found this in the documentation to be essentially _a_ minimum for it to work with eigvecs(i.e. it definitely works with this).
c     Correct implementation would imply using the search mode (lwork=-1) for our function and dynamic memory to reshape the work vector into it's optimal size.  
      rho = (rmax/rmin)**(1.d0/(N-1))
      do j = 1,N !These define the position parameters of our gaussians. In the end, nu(j) is the one we will actually use
        rj(j) = rmin*(rho**(j-1))
        nu(j) = 1/(rj(j)**2)
      enddo
      write(*,*) "-----------------------------------------------------"
      write(*,*) "----------------Potencial Columbiano----------------"
      write(*,*) "-----------------------------------------------------"
      write(*,*)""
      write(*,*) "----------------------Energias----------------------" 
      do j =0,2
        l = 1.d0*j
        call ConstructCornellHam(nu,l,k,hc,0.0d0,mc,N,H,B,rsq) !This subroutine takes in the parameters of the problem and outputs the Hamiltonian and an
c       auxiliary function to solve the root-mean squared problem
c       In this case, we take sigma to be 0 to have the Coulomb potential
        call DSYGV(1,'V','U',N,H,N,B,N,eigvals,work,lwork,info) !Lapack subroutine for generalized eigenvalue problems of symmetric matrices
c       eigvals is a vector of the eigenvalues and "H", the input matrix, is now a matrix of eigenvectors
        write(*,*) " 			l =", l
        write(*,*) " Numerico			", "Analitico			", "Error"
        do i =1,5
          Eanal = -0.25d0*k**2 *mc/(i+l)**2 !Analytic expression for the energy levels of the Coulomb potential
          write(*,*) eigvals(i),"GeV", Eanal,"GeV", 
     &     abs(Eanal-eigvals(i))/abs(Eanal)
        enddo
      enddo
      write(*,*) "---------------Radio Cuadratico Medio---------------"
      do m =0,2
        l = 1.d0*m
        write(*,*) " 			l =", l
        write(*,*) " Numerico			", "Analitico			", "Error"
        call ConstructCornellHam(nu,l,k,hc,0.0d0,mc,N,H,B,rsq)
        call DSYGV(1,'V','U',N,H,N,B,N,eigvals,work,lwork,info)
        call RootMeanSquared(H,Rsq,n,rms) !This subroutine just takes the eigenvectors and an auxiliary matrix and gives us a vector of the rms for 
c       the different energy levels
        do i =1,5
          rmsanal=(2*hc/(k*mc))**2*(i+l)**4*(1+1.5d0*(1-(l*(l+1)-1/3.d0)
     &      /(i+l)**2)) !Analytic expression for the root-mean squared radius for the Coulomb potential
          rmsanal = sqrt(rmsanal)
        write(*,*) rms(i),"fm",rmsanal,"fm", abs(rms(i)-rmsanal)/rmsanal
        enddo
      enddo
      write(*,*) "-----------------------------------------------------"!All will be a repeat but with sigma set to sigma now, until desintegration
      write(*,*) "-----------------Potencial Cornell-------------------"      
      write(*,*) "-----------------------------------------------------"
      write(*,*)""
      write(*,*) "----------------------Energias----------------------" 
      do j =0,2
        l = 1.d0*j
        call ConstructCornellHam(nu,l,k,hc,sigma,mc,N,H,B,rsq)
        call DSYGV(1,'V','U',N,H,N,B,N,eigvals,work,lwork,info)
        write(*,*) "l =", l
        do i =1,5
          write(*,*) eigvals(i)+M0 , "GeV"!The M0 adds the mass to fix the groundstate at the known mass of J/psi
        enddo
      enddo
      write(*,*) "---------------Radio Cuadratico Medio---------------"
      do m =0,2
        l = 1.d0*m
        write(*,*) "l =", l
        call ConstructCornellHam(nu,l,k,hc,sigma,mc,N,H,B,rsq)
        call DSYGV(1,'V','U',N,H,N,B,N,eigvals,work,lwork,info)
        call RootMeanSquared(H,Rsq,n,rms)
        do i =1,5
        write(*,*) rms(i), "fm"
        enddo
      enddo

      write(*,*) "-------------------Desintegracion-------------------" 
      call ConstructCornellHam(nu,0.0d0,k,hc,sigma,mc,N,H,B,rsq)
      call DSYGV(1,'V','U',N,H,N,B,N,eigvals,work,lwork,info)
      do j = 1,N
        C(j) = sqrt(2*(2*nu(j))**(1.5d0)/Gamma(1.5d0))
      enddo
      call CalcPsi(C,H,N,psi0) !Calculates psi(0) 
      do i = 1,5 !This calculates the decay width. It's a cumbersome formula but not that hard to implement
        Mn(i) = eigvals(i)+M0
        GammaFac(i) = 16*pi*(4/9.d0)*alpha**2*hc**3/Mn(i)**2
     &  *(1-4*k/pi)*1d6 !1d6 to get it in keV
        Gam(i) = GammaFac(i)*psi0(i)**2  
      enddo
      do i =1,5
        write(*,*) Gam(i),"keV"
      enddo
      end program

      subroutine ConstructCornellHam(nu,l,k,hc,sigma,mc,N,H,B,Rsq)
      Implicit none
      !This subroutine calculates the Hamiltonian and all associated matrices for the Cornell potential as linear combinations of Gaussian basis functions
      !Inputs: nu,l,k,hc,sigma,mc,N :: Parameters of the physical problem and dimension 
      !Outputs: H(N,N),B(N,N),Rsq(N,N) :: H is the Hamiltonian, B the generalized Eigval matrix, Rsq the Root-mean squared auxiliary matrix  
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
          Rsq(i,j) = Gamma(l+2.5d0)/Gamma(l+1.5d0)*B(i,j)/(nu(i)+nu(j))!matrix element to calculate the rms radius later
        enddo
      enddo
      end subroutine 

      subroutine RootMeanSquared(EigVec,Rsq,n,rms)
      !This subroutine calculates the root-mean square position given the eigenvectors of a hamiltonian. 
      !Inputs: Eigvec(N,N), Rsq(N,N), N :: Eigvec is the matrix of Eigenvectors, Rsq the Root-mean-squared auxiliary matrix, and N dimension
      !Output: rms(N) :: rms is the vector of the root mean square radius for each energy
      Integer N,i,j,k
      Real*8 EigVec(N,N), Rsq(N,N), Rms(N)
      do k = 1,N
        rms(k) = 0.0d0
      	do i = 1,N
      	  do j = 1,N
            rms(k) = rms(k) + EigVec(i,k) * Rsq(i,j) *EigVec(j,k)
          enddo
        enddo
        rms(k) = sqrt(rms(k))
      enddo
      end subroutine
      

      subroutine CalcPsi(C,EigVec,N,psi0)
      !This subroutine calculates the value of the wavefunction of l=0 at 0 given the eigenvectors and a vector related to the position parameters
      !which we will use to calculate the the decay-width of our potential.
      !Inputs: C(N), EigVec(N,N), N :: C is a vector that depends on the positions, EigVec is the matrix of eigenvectors of our system, N the dimension 
      !Outputs: psi0(N) :: Psi0 is the value of our wavefunction at 0 when l=0
      Integer N,i,j
      Real*8 EigVec(n,n), psi0(N),C(N),pi
      parameter (pi =4.d0*datan(1.d0))
      do i = 1,N
        psi0(i) = 0.d0
      	do j =1,N
          psi0(i) = psi0(i) + EigVec(j,i)*C(j)
          enddo
          psi0(i) = psi0(i)/sqrt(4*pi)
      enddo
      end subroutine