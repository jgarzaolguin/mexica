	subroutine lda_sp(mesh,r,x,f,rvx,xa)
	implicit double precision(a-h,o-z)
        include 'parameter'
	parameter(onethi=1.0/3.0)
	dimension r(ndim),x(ndim),rvx(2,ndim),xa(ndim)
	dimension f(ndim)
	pi=4.0*atan(1.0D00)
	fpi=4.0*pi
	bound=1.0D-20
	do i=1,mesh
	  xa(i)=0.0
	end do
	nt=2
	do isp=1,nt
	  do i=1,mesh
		if (isp.eq.1) then
		  if (dabs(x(i)).lt.bound) then
		    rvx(1,i)=0.0
		    term=0.0
		  else
		    rho=x(i)
		    cuad=r(i)*r(i)
		    term=2.0*fpi*cuad*funlda(rho)
		    rvx(1,i)=2.0*derfunlda(rho)*r(i)
		  end if
		else
		  if (dabs(f(i)).lt.bound) then
		    rvx(2,i)=0.0
		    term=0.0
		  else
		    rho=f(i)
		    cuad=r(i)*r(i)
		    term=2.0*fpi*cuad*funlda(rho)
		    rvx(2,i)=2.0*derfunlda(rho)*r(i)
		  end if
		end if
		xa(i)=xa(i)+term
	  end do
	end do
	return
	end
	double precision function funlda(rho)
	implicit double precision(a-h,o-z)
	parameter(onethi=1.0/3.0)
	pi=4.0*atan(1.0D00)
	c2=4.0*onethi
	cx=1.5*(3.0/(4.0*pi))**onethi
	rho43=rho**c2
	funlda=-cx*rho43
	return
	end
	double precision function derfunlda(rho)
	implicit double precision(a-h,o-z)
	parameter(onethi=1.0/3.0)
	pi=4.0*atan(1.0D00)
	c2=4.0*onethi
	cx=1.5*(3.0/(4.0*pi))**onethi
	rho13=rho**onethi
	derfunlda=-c2*cx*rho13
	return
	end

