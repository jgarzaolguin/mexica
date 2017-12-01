      subroutine x_slater_sp(rhoa, rhob, ex, vxa, vxb)
      implicit double precision(a-h,o-z)
      pi = 4d00*atan(1d00)
      fpi = 4d00*pi
      bound = 1d-20
      if (dabs(rhoa).lt.bound) then
        vxa = 0d00
        terma = 0d00
      else
        terma = funlda(rhoa)
        vxa = derfunlda(rhoa)
      end if
      if (dabs(rhob).lt.bound) then
        vxb = 0d00
        termb = 0d00
      else
        termb = funlda(rhob)
        vxa = derfunlda(rhob)
      end if
      ex = terma + termb
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

