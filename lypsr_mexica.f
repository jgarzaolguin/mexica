      subroutine lyp_sr(ri, rho, drho, d2rho, auxe, term1)
      implicit double precision(a-h,o-z)
      pi = 4d00*atan(1d00)
      a = 0.049d00
      b = 0.132d00
      c = 0.2533d00
      d = 0.349d00
      cf = (3d00*pi*pi)**(2d00/3d00)
      cf = 3d00*cf/10d00
      bound = 1d-20
      if (rho.le.bound) then
	enerc = 0.0d00
	potcor = 0.0d00
      else
        rho13 = rho**(1d00/3d00)
        rho23 = rho13**2d00
        rho53 = rho13**5d00
        exp1 = exp(c/rho13)
        if (exp1.gt.1d25) exp1 = 1d25
        drho13 = d + rho13
        term1 = -7d00*b*c*d - b*(7d00*c + 10d00*d)*rho13
     &          -3d00*b*rho23
        term1 = term1*drho**2d00 + 72d00*d*(b*cf + exp1)*rho**3d00
        term1 = term1 + 72d00*(b*cf + exp1)*rho13**10d00
        term1 = -a*term1/exp1
        enerc = term1/(72d00*rho53*drho13**2d00)
ccccccccccccccccccccccccccccccccccccccccccccccccc
c       Potential
ccc
        term2 = 7d00*b*c*c*d*d + b*c*(14d00*c - 25d00*d)*d*rho13
     &  + b*(7d00*c*c - 64d00*c*d - 40d00*d*d)*rho23
     &  - 3d00*b*(13d00*c + 23d00*d)*rho
        term2 = term2*drho*drho
        if (ri.eq.0d00) then
          rlaplac = 1d25
        else
          rlaplac = d2rho + 2d00*drho/ri
        endif
        term3 = (-15d00*b*drho*drho + 42d00*b*c*d*d*rlaplac)*rho13**4d00
        term3 = term3 + 12d00*b*d*(7d00*c + 5d00*d)*rlaplac*rho53
     &  + 6d00*b*(7d00*c + 13d00*d)*rlaplac*rho*rho 
     &  + 18d00*b*rlaplac*rho13**7d00
        term1 = term2 + term3 + 216d00*(b*cf + exp1)*rho**4 
     &  + 72d00*(b*cf*(c + 7d00*d) + 7d00*d*exp1)*rho13**11
     &  + 144d00*d*(b*cf*(c + 2d00*d) + 2d00*d*exp1)*rho13**10d00
     &  + 72d00*b*c*cf*d*d*rho**3d00
        potcor = -a*term1/(216d00*exp1*(drho13*rho)**3d00)
cccccccccccccccccccccccccccccccccccccccccccccccc
      end if
      auxe = enerc
      term1 = potcor
      return
      end
