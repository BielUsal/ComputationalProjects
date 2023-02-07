	program oscilador_armonico
	implicit none
	integer  i,j,k,l,neq,nptos
	parameter(neq=2,nptos=100)
	real*8 tt(nptos),yy(neq,nptos)
	real*8 t_ini,t_fin,h,t,omega
	real*8 y(neq),yp(neq),yout(neq)
	common /parametros/ omega

	omega=1.d0
	t_ini=0.d0
	t_fin=30.d0
	h=(t_fin-t_ini)/(nptos-1.d0)
	do i=1,nptos
		tt(i)=t_ini+h*(i-1.0)
	end do

	!condicion inicial
	y(1)=1.0
	y(2)=0.d0
	yy(1,1)=y(1)	!guardamos las condiciones iniciales
	yy(2,1)=y(2)
	call derivs(t_ini,y,yp)
	do i=2,nptos
		t=tt(i)
		call rk4(y,yp,neq,t,h,yout,derivs)
		yy(1,i)=yout(1)
		yy(2,i)=yout(2)
		y(1) =yout(1)
		y(2) =yout(2)
		call derivs(t,y,yp)
	end do

	open(22,file='harmonicdata',status='unknown')
	do i=1,nptos
	write(22,*)tt(i),yy(1,i),yy(2,i)
	enddo
	close(22)
	end



      subroutine rk4(y,dydx,n,x,h,yout,derivs)
      integer n,nmax
      real*8 h,x,dydx(n),y(n),yout(n)
      external derivs
      parameter (nmax=50)
      integer i
      real*8 h6,hh,xh,dym(nmax),dyt(nmax),yt(nmax)
      hh=h*0.5d0
      h6=h/6.d0
      xh=x+hh
      do 11 i=1,n
        yt(i)=y(i)+hh*dydx(i)
11    continue
      call derivs(xh,yt,dyt)
      do 12 i=1,n
        yt(i)=y(i)+hh*dyt(i)
12    continue
      call derivs(xh,yt,dym)
      do 13 i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
13    continue
      call derivs(x+h,yt,dyt)
      do 14 i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.d0*dym(i))
14    continue
      return
      end




	subroutine derivs(t,y,yp)
	integer neq
	parameter(neq=2)
	real*8	y(neq),yp(neq),t,omega
	common /parametros/ omega

	yp(1)=y(2)
	yp(2)=-omega**2*y(1)
	return
	end


