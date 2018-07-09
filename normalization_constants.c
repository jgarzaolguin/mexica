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
                     	           double *alphas, int print_expo)
{ //esto es para stos
 extern double intc(int a, double b, double r, double*, double*);
 extern double upper_incomplete_gamma(double Rc, int a, double b);
 int result;
 double total, temp, temp_2, temp_3, exis, K_mu, alpha_mu, zeta_mu, factor1, factor2;

 result = 0;
    if (strcmp(using_gamma,"YES") == 0){//label 1
  
	zeta_mu = zetas[mu];
	
	factor1 = gamma_couple/((1.f - gamma_couple)*Rc) + zeta_mu;
	factor2 = (double) (ang[mu] + np[mu])/Rc;
	
	if(factor1 < factor2) {
	  printf("\n Problems with condition for internal exponents in the basis set\n");
          result = 1;
        }

	temp       = (double) -1.f*(np[mu] + ang[mu]) + gamma_couple/(1.f - gamma_couple);
	temp_2     = (double) np[mu] + ang[mu];
	temp_3     = (double) (np[mu] + ang[mu]) + gamma_couple/(gamma_couple - 1.f);
	alpha_mu   = zeta_mu + temp/Rc;
        if (print_expo == 1) printf("alpha %lf\n", alpha_mu);
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
           if(zeta_mu < temp_2/Rc) {
             printf("\n Problems with condition for internal exponents in the basis set\n");;
             result = 1;
             printf("Hola :( \n");
           }
           
           alpha_mu = (double) zeta_mu - temp_2/Rc;
           if (print_expo == 1) printf("alpha %lf\n", alpha_mu);
           *alphas = alpha_mu;
           K_mu = pow(Rc,-temp_2)*exp(temp_2);
           total = upper_incomplete_gamma(Rc, 2*ang[mu], 2.f*alpha_mu);
           total = K_mu*K_mu*intc(2*np[mu], 2.f*zeta_mu, Rc, arreglo_factorial, arreglo_inv_factorial) + total;
           
           total = 1.f/sqrt(total);
           *const_n_plus =  total;
           *const_n_minus = K_mu*total;
            
         
           }//label 2
 
 return result;
 } //esto es para stos

/* Mike */
/* ----------------------------------------------------------------------------------------------------------------- */
/* This is for GTO's */

 double constant_normalization_GTO(int nu, int n, double alpha, double *n_gto) {

    extern long long int factorial_mike(int );
    double result;
    long int factor;
    double pi;

    pi = ((double) 4)*atan(1.f);
    factor = ((double) pow(2, 2*n + 1));
    result = factor*((double) factorial_mike(n));
    result = result*pow(2.f*alpha, n + 0.5);
    result = result/(((double) factorial_mike(2*n))*sqrt(pi));
    result = sqrt(result);
    *n_gto = result;
    return 1;
 } /* utilizaré esta definición para los elementos de matriz */ 
 
 double cte_norm_gto(int *mu, int *np, double *expo, double *n_gto) {

    extern long long int factorial_mike(int );
    double result;
    long int factor;
    double pi;
    int index_mu;

    pi = ((double) 4)*atan(1.f);

    index_mu = *mu;
    factor = ((double) pow(2, 2*np[index_mu] + 1));
    result = factor*((double) factorial_mike(np[index_mu]));
    result = result*pow(2.f*expo[index_mu], np[index_mu] + 0.5);
    result = result/(((double) factorial_mike(2*np[index_mu]))*sqrt(pi));
    result = sqrt(result);
    *n_gto = result;
    return 1;
 } /* utilizaré esta definición para la construcción de la densidad radial */

  double constant_normalization_GTO_imp(int nu, int n, double alpha, double r0, double *n_gto_imp) {  /* únicamente para el caso de pared impenetrable */
    extern double lower_incomplete_gamma_function(int , double );
    double fact1, fact2, fact3, fact4, number1, number2, number3, y, result;
    
    number1 = ((double)n) + 0.5f;
    fact1 = sqrt(2.f*pow(2.f*alpha,number1));
    fact2 = 1.f + (2.f*((double)n) + 1.f)/(4.f*alpha*pow(r0,2.f));
    fact3 = 2.f/(r0*sqrt(2.f*alpha));
    number2 = ((double)n) - 0.5f;
    number3 = 2.f*((double)n) - 1.f;
    y = 2.f*alpha*pow(r0,2.f);
    fact4 = (pow(2.f*alpha,number2)*pow(r0,number3))/exp(y);

    result = fact2*lower_incomplete_gamma_function(2*n + 1, y) - fact3*lower_incomplete_gamma_function(2*(n + 1), y) - fact4;
    result = fact1/sqrt(result);

    *n_gto_imp = result;
    return 1;
    
 } /* utilizaré esta definición para los elementos de matriz */

   double cte_norm_gto_imp(int *nu, int *np, double *expo, double r0, double *n_gto_imp) {  /* únicamente para el caso de pared impenetrable */
    extern double lower_incomplete_gamma_function(int , double );
    double fact1, fact2, fact3, fact4, number1, number2, number3, y, result;
    int index_nu;

    index_nu = *nu;

    number1 = ((double)np[index_nu]) + 0.5f;
    fact1 = sqrt(2.f*pow(2.f*expo[index_nu],number1));
    fact2 = 1.f + (2.f*((double)np[index_nu]) + 1.f)/(4.f*expo[index_nu]*pow(r0,2.f));
    fact3 = 2.f/(r0*sqrt(2.f*expo[index_nu]));
    number2 = ((double)np[index_nu]) - 0.5f;
    number3 = 2.f*((double)np[index_nu]) - 1.f;
    y = 2.f*expo[index_nu]*pow(r0,2.f);
    fact4 = (pow(2.f*expo[index_nu],number2)*pow(r0,number3))/exp(y);

    result = fact2*lower_incomplete_gamma_function(2*np[index_nu] + 1, y) - fact3*lower_incomplete_gamma_function(2*(np[index_nu] + 1), y) - fact4;
    result = fact1/sqrt(result);

    *n_gto_imp = result;
    return 1;

 } /* utilizaré esta definición para la construcción de la densidad radial */

