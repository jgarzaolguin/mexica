/* Integrales monoelectronicas y bielectronicas
   para resolver las ecuaciones de Hartree-Fock
   en sistemas atomicos libres.
   Jorge Garza, Junio del 2006 */
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


 
int constant_normalization_free(int*    mu, 
				double* expo, 
				int*    np, 
				double* resultado)
{
 int argumento, interno_mu;
 double numerador, denominador, base, exponente, dentro;
 extern double factorial(int);
 interno_mu = *mu;
 base = (double)2.0*expo[interno_mu];
 exponente = (double)np[interno_mu] + (double)0.5;
 numerador = pow(base, exponente);
 argumento = 2*np[interno_mu];
 dentro = factorial(argumento);
 denominador = sqrt(dentro);
 *resultado = numerador/denominador;
 return 1;
 }

int constants_normalization_finite(int     mu, 
				   int    *np, 
				   int    *ang, 
				   double *zetas, 
				   double  Rc,
                     	           double *arreglo_factorial, 
				   double *arreglo_inv_factorial,
                                   char   *using_gamma,
                     		   double  gamma_couple,
                     		   double *const_n_minus, 
				   double *const_n_plus,
                     	           double *alphas)
{
 extern double intc(int a, double b, double r, double*, double*);
 extern double upper_incomplete_gamma(double Rc, int a, double b);
 double total, temp, temp_2, temp_3, exis, K_mu, alpha_mu, zeta_mu, factor1, factor2;

    if (strcmp(using_gamma,"YES") == 0){//label 1
  
	zeta_mu = zetas[mu];
	
	factor1 = gamma_couple/((1.f - gamma_couple)*Rc) + zeta_mu;
	factor2 = (double) (ang[mu] + np[mu])/Rc;
	
	if(factor1 < factor2)
	printf("\n Problems with condition for internal exponents in the basis set\n");

	temp       = (double) -1.f*(np[mu] + ang[mu]) + gamma_couple/(1.f - gamma_couple);
	temp_2     = (double) np[mu] + ang[mu];
	temp_3     = (double) (np[mu] + ang[mu]) + gamma_couple/(gamma_couple - 1.f);
	alpha_mu   = zeta_mu + temp/Rc;
        printf("alpha %lf\n", alpha_mu);
	*alphas    = alpha_mu;
	K_mu       = pow(Rc,-1.f*(temp_2 + 1.f))*exp(temp_3)/(1.f - gamma_couple);
	
	total      = Rc*Rc*
		     intc(2*np[mu], 2.f*zeta_mu, Rc, arreglo_factorial, arreglo_inv_factorial) 
                     -2.f*Rc*gamma_couple*
		     intc(2*np[mu] + 1, 2.f*zeta_mu, Rc, arreglo_factorial, arreglo_inv_factorial)   
                     +gamma_couple*gamma_couple*
		     intc(2*np[mu] + 2, 2.f*zeta_mu, Rc, arreglo_factorial, arreglo_inv_factorial);	
	total      = total*K_mu*K_mu;
	total      = upper_incomplete_gamma(Rc, 2*ang[mu], 2.f*alpha_mu) + total;

	total      = 1.f/sqrt(total);

   	*const_n_plus  = total;
   	*const_n_minus = K_mu*total;
      }//label 1
       else {//label 2
   
           zeta_mu = zetas[mu];
           
           temp_2 = (double) (np[mu] + ang[mu]);
           
           if(zeta_mu < temp_2/Rc)
               printf("\n Problems with condition for internal exponents in the basis set\n");;
           
           alpha_mu = (double) zeta_mu - temp_2/Rc;
          
           printf("alpha %lf\n", alpha_mu);
           
           *alphas = alpha_mu;
           
           K_mu = pow(Rc,-temp_2)*exp(temp_2);
           
           total = upper_incomplete_gamma(Rc, 2*ang[mu], 2.f*alpha_mu);
           total = K_mu*K_mu*intc(2*np[mu], 2.f*zeta_mu, Rc, arreglo_factorial, arreglo_inv_factorial) + total;
           
           total = 1.f/sqrt(total);
           *const_n_plus =  total;
           *const_n_minus = K_mu*total;
            
         
           }//label 2
 
 return 0;
 }

