c**********************
c Spin Polarized
c**********************        
        subroutine pbe_sp(mesh,r,x,f,rvx,xa)
        implicit double precision(a-h,o-z)
        include 'parameter'  
        dimension r(ndim),x(ndim),rvx(2,ndim),xa(ndim),f(ndim)
	dimension drho(ndim),d2rho(ndim)      
	parameter(onethi=1.d00/3.d00)
        parameter(fth=4.d00/3.d00)  
        pi=4.d00*atan(1.0D00)
        fpi=4.d00*pi
        cf=(rho*3.d00*pi**2.d00)**onethi
        ax=(-3.d00*(3.d00*pi**2.d00)**onethi)/fpi
	h=r(2)
	do i=1,mesh
	  xa(i)=0.0
	end do
	nt=2
	do isp=1,nt
	  if (isp.eq.1) then
c	    CALL primder(mesh,h,x,drho) Spin A
c	    CALL seconder(mesh,h,x,d2rho) 
	    do i=1,mesh
		rho=x(i)
		if (dabs(rho).le.0.0) then
		  rvx(1,i)=0.0
		  term=0.0
		else
		  ri=r(i)
	          r2=r(2)
               rho13=rho**onethi
               rho43=rho**fth 
                  dr=drho(i) 
                 d2r=d2rho(i)
                 sp=dabs(dr)/(2.d00*cf*rho)  
               rcuad=ri*ri
	        term=2.0*fpi*rcuad*funcPBE_sp(rho43,sp)
	    rvx(1,i)=2.0*derfuncPBE_sp(ri,r2,rho,dr,d2r,sp)*ri
	     
                end if
		xa(i)=xa(i)+term
	    end do
	  else
c CALL primder(mesh,h,f,drho) spin B
c CALL seconder(mesh,h,f,d2rho)
	    do  i=1,mesh
		rho=f(i)
		if (dabs(rho).le.0.0) then
              rvx(2,i)=0.0
		  term=0.0
		else
		  ri=r(i)
		  r2=r(2)
	       rho13=rho**onethi
               rho43=rho**fth 
	          dr=drho(i) 
                 d2r=d2rho(i)
                  cf=(rho*3.d00*pi**2.d00)**onethi
                  sp=dabs(dr)/(2.d00*cf*rho) 
               rcuad=ri*ri
            term=2.0*fpi*rcuad*funcPBE_sp(rho43,sp)
            rvx(2,i)=2.0*derfuncPBE_sp(ri,r2,rho,dr,d2r,sp)*ri
		end if
		xa(i)=xa(i)+term             
        end do
	  end if
	end do
	return
	end
C********************************
C Funcional-PBE-SP
C*****************************
      double precision function funcPBE_sp(rho43,sp)
      implicit double precision(a-h,o-z)
      parameter(onethi=1.d00/3.d00)
      parameter(fth=4.d00/3.d00)  
      pi=4.d00*atan(1.0D00)
      fpi=4.d00*pi
      ax=(-3.d00/fpi)*(3.d00*pi**2.d00)**onethi
      funcPBE_sp=ax*rho43*fx(sp)
      return
      end
c************************
c Der-Funcional PBE-SP 
c************************
       double precision function derfuncPBE_sp(r,r2,rho,drho,d2rho,sp)
       implicit double precision(a-h,o-z)
       parameter(onethi=1.d00/3.d00)
       parameter(fth=4.d00/3.d00) 
       pi=4.d00*atan(1.0D00)
       fpi=4.d00*pi
       rho13=rho**onethi
       dr=drho  
       d2r=d2rho
       alp=-(3.d00/pi)**onethi
       bet=-3.d00/(8.d00*pi)
c VXC=dExc/drho 
c Esquema 1  Sugerencia: Javier C.      
       p1=alp*rho13*(fx(sp)-sr*dfx(sp))
       p2=bet*dfx(sp)*(dr)/dabs(dr)
       derfuncPBE_sp=p1-p2
       return
       end
c*********************
c*********************
c******spin restricted
c*********************
      subroutine pbe_sr(rho, drho, auxe, term1, term2)
      implicit none
      double precision onethi, pi, bound, auxe, term1, term2,
     &                 rho, drho, rkf, varS, ener_PBE_x,
     &                 Fx, eps_x, der_PBE_sr_local,
     &                 der_PBE_sr_nolocal
      onethi=1.d00/3.d00
      pi=4.d00*datan(1.0D00)
      bound= 1d-20
      if (rho.le.bound) then
       auxe = 0D00
       term1 = 0D00
      else
        rkf = (3.d00*pi*pi*rho)**onethi
        varS=dabs(drho)/(2.d00*rkf*rho)  
        ener_PBE_x = rho*eps_x(rho)*Fx(varS)
        auxe=ener_PBE_x
        term1=der_PBE_sr_local(rho, varS)
        term2=der_PBE_sr_nolocal(rho, drho, varS)
      end if 
      return
      end
c*********************************************
c Derivada funcional para ser escrita sobre un grid
c Solamente para \'atomos
c************************************************
      subroutine pbegrid(rho, derho, varS, der_varS, pot_x)
      implicit none
      double precision rho, varS, der_varS, pot_x, derho
      double precision term1, term2, kappa, mu, cuadS, frac, pi
      double precision eps_x, Fx, der_Fx
      pi = 4d00*datan(1d00)
      kappa = 0.804d00
      mu = 0.21951d00
      term1 = 4d00*eps_x(rho)*(Fx(varS) - varS*der_Fx(varS))/3d00
      cuadS = varS*varS
      frac = mu/kappa
      term2 = kappa - 3d00*mu*cuadS
      term2 = term2/((1d00 + frac*cuadS)**3d00)
      term2 = 3d00*frac*term2/(4d00*pi)
      term2 = term2*der_varS
      pot_x = term1 + term2
      return
      end
c****************************************
c** Intercambio para el gas de electrones
c***************************************
      double precision function eps_x(rho)
      implicit none
      double precision onethi, pi, fpi, rkf, rho
      onethi = 1d00/3d00
      pi = 4d00*datan(1d00)
      fpi = 4d00*pi
      rkf = (3d00*pi*pi*rho)**onethi
      eps_x = -3d00*rkf/fpi
      return
      end
c***********************************************       
c** Definición de Función de exhacervamieto
c *********************************************
      double precision function Fx(varS)
      implicit none
      double precision pkap, pmu, ratio, den, varS
      pkap = 0.804d00
      pmu = 0.21951d00
      ratio = pmu/pkap 
      den = 1.d00+ratio*varS*varS
      Fx = 1.d00+pkap-pkap/den
      return
      end
c***********************************************
c** Derivada de Fx con respecto a la variable s
c***********************************************
      double precision function der_Fx(vars)
      implicit none
      double precision pkap, pmu, ratio, den, varS
      pkap = 0.804d00
      pmu = 0.21951d00
      ratio = pmu/pkap 
      den = 1.d00+ratio*varS*varS
      den = den*den
      der_Fx = 2d00*pmu*varS/den
      return
      end
c************************
c derivada funcional 
c************************
      double precision function der_PBE_sr_local(rho, varS)
      implicit none
      double precision rho, varS, term1
      double precision eps_x, Fx, der_Fx
      term1 = 4d00*eps_x(rho)*(Fx(varS) - varS*der_Fx(varS))/3d00
      der_PBE_sr_local = term1
      return
      end
c
      double precision function der_PBE_sr_nolocal(rho, drho, varS)
      implicit none
      double precision rho, drho, varS, term1, pi
      double precision der_Fx
      pi = 4d00*datan(1d00)
      term1 = -3d00*drho*der_Fx(varS)/(8d00*pi*dabs(drho))
      der_PBE_sr_nolocal = term1
      return
      end
c**********************************************
c**Derivadas de la funcion de exhacervamiento
c**********************************************
       double precision function dfx(varS)
       implicit double precision(a-h,o-z)
       pkap=0.804d00
       pmu=0.21951d00
       ratio=pmu/pkap 
       den=1.d00+ratio*varS**2.d00
       dfx=2.d00*pmu*varS/(den**2.d00)
       return
       end

       double precision function d2fx(varS)
       implicit double precision(a-h,o-z)
       pkap=0.804d00
       pmu=0.21951d00
       ratio=pmu/pkap 
       den=1.d00+ratio*varS**2.d00
       d2fx=(-8.d00*(pmu*varS**2.d00)**2.d00/(pkap*den**3.d00))
     &      + (2.d00*pmu/den**2.d00)
       return
       end
