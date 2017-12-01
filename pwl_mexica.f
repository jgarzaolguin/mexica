	subroutine pw92sr(x,auxe,term1)
	implicit double precision(a-h,o-z)
	pi=4.0*atan(1.0D00)
	bound=1.0D-20
	pot=1.0/3.0
	cx=(3.0/(4.0*pi))**pot
	rhot=x
	if (dabs(rhot).lt.bound) then
	  auxe=0.0
	  term1=0.0
	else
	  rs=cx/(rhot**pot)
	  ec=ec0(rs)
	  cpc1=ec-pot*rs*decrs(rs,0.0d00)
	  term1=cpc1
	  auxe=rhot*ec
	end if
	return
	end
C	****************************************
	subroutine pw92_sp(mesh,r,x,f,rvx,xa)
	implicit double precision(a-h,o-z)
	include 'parameter'
C	parameter (mesh=441)
	dimension x(ndim),f(ndim)
	dimension xa(ndim),r(ndim),rvx(2,ndim)
	pi=4.0*atan(1.0D00)
	fpi=4.0*pi
	bound=1.0D-16
	d2fpw=1.709921
	pot=1.0/3.0
	cx=(3.0/(4.0*pi))**pot
	do i=1,mesh
	  rhot=x(i)+f(i)
	  if (dabs(rhot).lt.bound) then
	    auxe=0.0
	    term1=0.0
	    term2=0.0
	  else
	    psi=(x(i)-f(i))/rhot
	    rs=cx/(rhot**pot)
	    ec=ec0(rs)+alfac(rs)*fpw(psi)*(1.0-(psi**4.0))/d2fpw
     *         +(ec1(rs)-ec0(rs))*fpw(psi)*(psi**4.0)
	    cuad=r(i)*r(i)
	    cpc1=ec-pot*rs*decrs(rs,psi)
	    cpc2=decpsi(rs,psi)
	    term1=2.0*(cpc1-(psi-1.0)*cpc2)*r(i)
	    term2=2.0*(cpc1-(psi+1.0)*cpc2)*r(i)
	    auxe=2.0*fpi*cuad*rhot*ec
	  end if
	  xa(i)=xa(i)+auxe
	  rvx(1,i)=rvx(1,i)+term1
	  rvx(2,i)=rvx(2,i)+term2
	end do
	return
	end
C	****************************************************
	double precision function fpw(x)
	implicit double precision(a-h,o-z)
	pot=4.0/3.0
	t1=(1.0+x)**pot
	t2=(1.0-x)**pot
	t3=((2.0**pot)-2.0)
	fpw=(t1+t2-2.0)/t3
	return
	end
C	******************************************************
	double precision function gint(rs,a,a1,b1,b2,b3,b4,p)
	implicit double precision(a-h,o-z)
	t1=-2.0*a*(1.0+a1*rs)
	t2=b1*dsqrt(rs)
	t3=b2*rs
	t4=b3*(rs**1.5)
	t5=b4*(rs**(p+1))
	denom=2.0*a*(t2+t3+t4+t5)
	gint=t1*dlog(1.0+1.0/denom)
	return
	end
C	***********************************************
	double precision function ec0(rs)
	implicit double precision(a-h,o-z)
	a=0.031091
	a1=0.21370
	b1=7.5957
	b2=3.5876
	b3=1.6382
	b4=0.49294
	p=1.0
	ec0=gint(rs,a,a1,b1,b2,b3,b4,p)
	return
	end
C	****************************************************
	double precision function ec1(rs)
	implicit double precision(a-h,o-z)
	a=0.015545
	a1=0.20548
	b1=14.1189
	b2=6.1977
	b3=3.3662
	b4=0.62517
	p=1.0
	ec1=gint(rs,a,a1,b1,b2,b3,b4,p)
	return
	end
C	******************************************************
	double precision function alfac(rs)
	implicit double precision(a-h,o-z)
	a=0.016887
	a1=0.11125
	b1=10.357
	b2=3.6231
	b3=0.88026
	b4=0.49671
	p=1.0
	alfac=-gint(rs,a,a1,b1,b2,b3,b4,p)
	return
	end
C	***********************************************************
	double precision function dfpw(x)
	implicit double precision(a-h,o-z)
	pot=1.0/3.0
	cx=4.0*pot
	t1=(1.0+x)**pot
	t2=(1.0-x)**pot
	t3=((2.0**cx)-2.0)
	dfpw=cx*(t1-t2)/t3
	return
	end
C	**************************************************
	double precision function decpsi(rs,x)
	implicit double precision(a-h,o-z)
	x3=x**3.0
	x4=x3*x
	d2fpw=1.709921
	t1=4.0*x3*fpw(x)*(ec1(rs)-ec0(rs)-alfac(rs)/d2fpw)
	t2=alfac(rs)*(1.0-x4)/d2fpw+(ec1(rs)-ec0(rs))*x4
	decpsi=t1+dfpw(x)*t2
	return
	end
C	*********************************************************
	double precision function decrs(rs,x)
	implicit double precision(a-h,o-z)
	x4=x**4.0
	d2fpw=1.709921
	t1=dec0(rs)*(1.0-x4*fpw(x))
	t2=dec1(rs)*fpw(x)*x4
	t3=dalfac(rs)*fpw(x)*(1.0-x4)/d2fpw
	decrs=t1+t2+t3
	return
	end
C	*************************************************************
	double precision function dgint(rs,a,a1,b1,b2,b3,b4,p)
	implicit double precision(a-h,o-z)
	q0=-2.0*a*(1.0+a1*rs)
	t2=b1*dsqrt(rs)
	t3=b2*rs
	t4=b3*(rs**1.5)
	t5=b4*(rs**(p+1))
	q1=2.0*a*(t2+t3+t4+t5)
	dq1=a*(t2/rs+2.0*b2+3.0*t4/rs+2.0*(p+1)*t5/rs)
	dgint=-2.0*a*a1*dlog(1.0+1.0/q1)-q0*dq1/(q1*(q1+1.0))
	return
	end
C	************************************************************
	double precision function dec0(rs)
	implicit double precision(a-h,o-z)
	a=0.031091
	a1=0.21370
	b1=7.5957
	b2=3.5876
	b3=1.6382
	b4=0.49294
	p=1.0
	dec0=dgint(rs,a,a1,b1,b2,b3,b4,p)
	return
	end
C	***************************************************************
	double precision function dec1(rs)
	implicit double precision(a-h,o-z)
	a=0.015545
	a1=0.20548
	b1=14.1189
	b2=6.1977
	b3=3.3662
	b4=0.62517
	p=1.0
	dec1=dgint(rs,a,a1,b1,b2,b3,b4,p)
	return
	end
C	****************************************************************
	double precision function dalfac(rs)
	implicit double precision(a-h,o-z)
	a=0.016887
	a1=0.11125
	b1=10.357
	b2=3.6231
	b3=0.88026
	b4=0.49671
	p=1.0
	dalfac=-dgint(rs,a,a1,b1,b2,b3,b4,p)
	return
	end
