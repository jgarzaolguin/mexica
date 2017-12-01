	subroutine becke_sp(mesh,r,x,f,rvx,xa)
	implicit double precision(a-h,o-z)
        include 'parameter'
	parameter(onethi=1.0/3.0)
	dimension r(ndim),x(ndim),rvx(2,ndim),xa(ndim),f(ndim)
	dimension drho(ndim),d2rho(ndim)
	pi=4.0*atan(1.0D00)
	c2=4.0*onethi
	fpi=4.0*pi
	h=r(2)
	do i=1,mesh
	  xa(i)=0.0
	end do
	nt=2
	do isp=1,nt
	  if (isp.eq.1) then
c	    CALL primder(mesh,h,x,drho)
c	    CALL seconder(mesh,h,x,d2rho) 
	    do i=1,mesh
		rho=x(i)
		if (dabs(rho).le.0.0) then
		  rvx(1,i)=0.0
		  term=0.0
		else
		  ri=r(i)
		  r2=r(2)
		  dr=drho(i)
		  d2r=d2rho(i)
		  rho43=rho**c2
		  vx=dabs(drho(i))/rho43
		  rcuad=ri*ri
		  term=2.0*fpi*rcuad*funcional_sp(rho43,vx)
		  rvx(1,i)=2.0*derfuncional_sp(ri,r2,rho,dr,d2r,vx)*ri
		end if
		xa(i)=xa(i)+term
	    end do
	  else
c            CALL primder(mesh,h,f,drho)
c            CALL seconder(mesh,h,f,d2rho)
	    do i=1,mesh
		rho=f(i)
		if (dabs(rho).le.0.0) then
		  rvx(2,i)=0.0
		  term=0.0
		else
		  ri=r(i)
		  r2=r(2)
		  rho43=rho**c2
		  dr=drho(i)
		  d2r=d2rho(i)
		  vx=dabs(drho(i))/rho43
		  rcuad=ri*ri
		  term=2.0*fpi*rcuad*funcional_sp(rho43,vx)
		  rvx(2,i)=2.0*derfuncional_sp(ri,r2,rho,dr,d2r,vx)*ri
		end if
		xa(i)=xa(i)+term
	    end do
	  end if
	end do
	return
	end
C	********************************************************
	double precision function funcional_sp(rho43,vx)
	implicit double precision(a-h,o-z)
	parameter(onethi=1.0/3.0)
	pi=4.0*atan(1.0D00)
	cx=1.5*(3.0/(4.0*pi))**onethi
	beta=0.0042
	funcional_sp=-rho43*(cx+beta*fbecke(vx))
	return
	end
C	****************************************************************
	double precision function derfuncional_sp(r,r2,rho,drho,d2rho,vx)
	implicit double precision(a-h,o-z)
	parameter(onethi=1.0/3.0)
	pi=4.0*atan(1.0D00)
	c2=4.0*onethi
	cx=1.5*(3.0/(4.0*pi))**onethi
	beta=0.0042
	rho13=rho**onethi
	rho43=rho**c2
	paux1=fbecke(vx)
	paux2=dfbecke(vx)
	paux3=d2fbecke(vx)
	if (r.lt.r2) then
	  paux4=1.0D35
	else
	  paux4=2.0/r
	end if
	pot1=-c2*cx*rho13
	pot2=-c2*beta*rho13*paux1
	pot3=beta*(-paux4+c2*rho13*vx)*paux2
	pot4=beta*(d2rho-c2*drho*drho/rho)*paux3/rho43
	derfuncional_sp=pot1+pot2+pot3+pot4
	return
	end
C	**********************************************************
	subroutine becke_sr(ri, r2, rho, drho, d2rho, auxe, term1)
	implicit double precision(a-h,o-z)
	parameter(onethi=1.0/3.0) 
	c2=4.0*onethi
	fac=2.0**onethi
        bound = 1d-20
	if (rho.le.bound) then
	  auxe = 0D00
	  term1 = 0D00
	else
	  rho43=rho**c2
	  dr=dabs(drho)
	  vx=fac*dr/rho43
	  auxe = funcional(rho43,vx)
	  term1 = derfuncional(ri,r2,rho,dr,d2rho,vx)
	end if
	return
	end
C     *****************************************************
	double precision function funcional(rho43,vx)
	implicit double precision(a-h,o-z)
	parameter(onethi=1.0/3.0)
	pi=4.0*atan(1.0D00)
	cd=0.75*(3.0/pi)**onethi
	fac=2.0**onethi
	beta=0.0042
	alfa=beta/fac
	funcional=-rho43*(cd+alfa*fbecke(vx))
	return
	end
C     *************************************************************
	double precision function derfuncional(r,r2,rho,drho,d2rho,vx)
	implicit double precision(a-h,o-z)
	parameter(onethi=1.0/3.0)
	pi=4.0*atan(1.0D00)
	c2=4.0*onethi
	cd=0.75*(3.0/pi)**onethi
	fac=2.0**onethi
	beta=0.0042
	alfa=beta/fac
	rho13=rho**onethi
	rho43=rho**c2
	paux1=fbecke(vx)
	paux2=dfbecke(vx)
	paux3=d2fbecke(vx)
	if (r.lt.r2) then
	  paux4=1.0D35
	else
	  paux4=2.0/r
	end if
	pot1=-c2*(cd+alfa*paux1)*rho13
	pot2=beta*(c2*drho/rho-paux4)*paux2
        if (dabs(d2rho).lt.1e-14.and.dabs(drho).lt.1e-14) then
          pot3 = 0d00
        else
          pot3=beta*fac*(d2rho-c2*drho*drho/rho)*paux3/rho43
        endif
	derfuncional=pot1+pot2+pot3
	return
	end
C	********************************************
	double precision function fbecke(x)
	implicit double precision(a-h,o-z)
	fbecke=x*x/denom(x)
	return
	end
C	****************************************************
	double precision function arcsinh(x)
	implicit double precision(a-h,o-z)
	arcsinh=dlog(x+raiz(x))
	return
	end
C	***********************************************
	double precision function dfbecke(x)
	implicit double precision(a-h,o-z)
	beta=0.0042
	t1=-3.0*beta*x*x
	t2=3.0*beta*x*raiz(x)*arcsinh(x)
	t3=denom(x)**2.0
	dfbecke=2.0*x*(t1+raiz(x)+t2)/(raiz(x)*t3)
	return
	end
C	******************************************************
	double precision function d2fbecke(x)
	implicit double precision(a-h,o-z)
	beta=0.0042
	cuad=x*x
	t1=-3.0*beta*cuad*(6.0+5.0*cuad-12.0*beta*cuad*raiz(x))
	t2=raiz(x)**3.0
	t3=-18.0*beta*beta*(x**3.0)*(2.0+x*x)*arcsinh(x)
	t4=denom(x)**3.0
	d2fbecke=2.0*(t1+t2+t3)/(t2*t4)
	return
	end
C	*****************************************************************
	double precision function raiz(x)
	implicit double precision(a-h,o-z)
	raiz=dsqrt(1.0+x*x)
	return
	end
C	**********************************************
	double precision function denom(x)
	implicit double precision(a-h,o-z)
	beta=0.0042
	denom=(1.0+6.0*beta*x*arcsinh(x))
	return
	end

