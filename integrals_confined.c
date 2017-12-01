#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>



 double bielectronic_integral_confined(int     l, 
				       int     mu, 
				       int     nu, 
				       int     lam, 
				       int     sig, 
				       double  Rc, 
				       double  inv_Rc,
                	  	       double* expo, 
			               int    *np, 
				       double *arreglo_factorial,
                                       double *arreglo_inv_factorial)
 {
  double total;

  extern double intradc(int, int, int, int, int, double, double,  double*, int*, double*, double*);

  total = intradc(l, mu, nu, lam, sig, Rc, inv_Rc, expo, np, arreglo_factorial,
                  arreglo_inv_factorial);

  return(total);
 }

 double int1c(int k, int mu, int nu, double Rc, double* expo, int* np,
              double* arreglo_factorial, double* arreglo_inv_factorial)
 {
  int entero1;
  double doble1, term1, term2, term3, fact, total;

  extern double constc(int, double, double*, int*, double*, double*);
  extern double intc(int, double, double, double*, double*);

  fact = constc(mu, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial)*
         constc(nu, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial);
  entero1 = k + np[mu] + np[nu] - 2;
  doble1 = expo[mu] + expo[nu];
  term1 = intc(entero1, doble1, Rc, arreglo_factorial, arreglo_inv_factorial);
  entero1 = entero1 + 1;
  term2 = (double)2.0*intc(entero1, doble1, Rc, arreglo_factorial, arreglo_inv_factorial)/Rc;
  entero1 = entero1 + 1;
  term3 = intc(entero1, doble1, Rc, arreglo_factorial, arreglo_inv_factorial)/(Rc*Rc);
  total = fact*(term1 - term2 + term3);

  return(total);
 }

 double intc(int a, double b, double r, double *arreglo_factorial,
             double *arreglo_inv_factorial)
 {
  int i;
  double denom, suma, uno, total, doble1, doble2;

//jgo  extern double factorial(int);

  if(a <= 0)
  {
  if(a < 0)
  printf("\nEntero negativo incorrecto r^%d\n", a);
  else
  total = (1.f - exp(-r*b))/b;
  }
  else
  {
  doble1 = (double)(a + 1);
  denom = pow(b,doble1);
  doble2 = b*r;
  suma = (double)0.0;
  for (i = 1; i <= a; i++) {
    doble1 = (double)i;
    uno = pow(doble2, doble1);
//jgo    suma += uno/factorial(i);
    suma += uno*arreglo_inv_factorial[i];
  }
  suma += (double)1.0;
//jgo  total = factorial(a)*((double)1.0 - suma/exp(doble2))/denom;
  total = arreglo_factorial[a]*((double)1.0 - suma/exp(doble2))/denom;
  }
  return(total);
 }

 double constc(int mu, double Rc, double* expo, int* np, double *arreglo_factorial,
               double *arreglo_inv_factorial)
 {
  int entero1;
  double doble1, term1, term2, term3, total;

  extern double intc(int, double, double, double*, double*);

  entero1 = 2*np[mu];
  doble1 = (double)2.0*expo[mu];
  term1 = intc(entero1, doble1, Rc, arreglo_factorial, arreglo_inv_factorial);
  entero1 = entero1 + 1;
  term2 = (double)2.0*intc(entero1, doble1, Rc, arreglo_factorial, arreglo_inv_factorial)/Rc;
  entero1 = entero1 + 1;
  term3 = intc(entero1, doble1, Rc, arreglo_factorial, arreglo_inv_factorial)/(Rc*Rc);
  total = (double)1.0/sqrt(term1 - term2 + term3);

  return(total);
 }


 double int_p_conf(int b, int c, int q, int mu, int nu, int lamb,
                   int sig, int *np, double *expo,double Rc, double *arreglo_factorial,
                   double *arreglo_inv_factorial)
 {
  int n, entero1, entero2;
  double term1, term2, suma, doble1, doble2, los_dos, K3, total, Rc_cuad, temp1, temp2;

  extern double constc(int, double, double*, int*, double*, double*);
  extern double intc(int, double, double, double*, double*);
//jgo  extern double factorial(int);

  Rc_cuad = Rc*Rc;
  entero1 = c + np[lamb] + np[sig] - q;
//jgo  K3 = factorial(entero1);
  K3 = arreglo_factorial[entero1];
  doble1 = expo[lamb] + expo[sig];
  temp1 = (double) (entero1 + 1);
  K3 = K3/pow(doble1, temp1);

  doble2 = expo[mu] + expo[nu];
  entero2 = b + np[mu] + np[nu] - 2;

  term1 = intc(entero2, doble2, Rc, arreglo_factorial, arreglo_inv_factorial);
  entero2 = entero2 + 1;
  term1 = term1 - 2.f*intc(entero2, doble2, Rc, arreglo_factorial, arreglo_inv_factorial)/Rc;
  entero2 = entero2 + 1;
  term1 = term1 + intc(entero2, doble2, Rc, arreglo_factorial, arreglo_inv_factorial)/Rc_cuad;

  los_dos = doble1 + doble2;
  suma = 0.f;
  for (n = 0; n <= entero1; n++) {
    entero2 = b + n + np[mu] + np[nu] - 2;
    term2 = intc(entero2, los_dos, Rc, arreglo_factorial, arreglo_inv_factorial);
    entero2 = entero2 + 1;
    term2 = term2 - 2.f*intc(entero2, los_dos, Rc, arreglo_factorial, arreglo_inv_factorial)/Rc;
    entero2 = entero2 + 1;
    term2 = term2 + intc(entero2, los_dos, Rc, arreglo_factorial, arreglo_inv_factorial)/Rc_cuad;
//jgo    term2 = term2*pow(doble1, (double) n)/factorial(n);
    term2 = term2*pow(doble1, (double) n)*arreglo_inv_factorial[n];
    suma = suma + term2;
  }

  total = constc(mu, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial)*
          constc(nu, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial)*K3;
  total = total*(term1 - suma);
  return total;
  }

 double i_conf(int b, int c, int mu, int nu, int lamb, int sig, int *np, double *expo,
               double Rc, double* arreglo_factorial, double* arreglo_inv_factorial)
 {
  double term1, total;

  extern double int_p_conf(int b, int c, int q, int mu, int nu, int lamb,
                           int sig, int *np, double *expo, double Rc,
                           double* arreglo_factorial, double *arreglo_inv_factorial);
  extern double constc(int, double, double*, int*, double*, double*);

  term1 = int_p_conf(b, c, 2, mu, nu, lamb, sig, np, expo, Rc,
                     arreglo_factorial, arreglo_inv_factorial);
  term1 = term1 - 2.f*int_p_conf(b, c, 1, mu, nu, lamb, sig, np, expo, Rc,
                                 arreglo_factorial, arreglo_inv_factorial)/Rc;
  term1 = term1 + int_p_conf(b, c, 0, mu, nu, lamb, sig, np, expo, Rc,
                             arreglo_factorial, arreglo_inv_factorial)/(Rc*Rc);

  total = constc(lamb, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial)*
          constc(sig, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial)*term1;

  return(total);

 }

 double intradc(int l, int mu, int nu, int lam, int sig, double Rc, double inv_Rc,
                double* expo, int* np, double* arreglo_factorial,
                double* arreglo_inv_factorial)
 {
  int entero1, entero2, entero3;
  double doble1, doble2, doble3, total, temp1, temp2, exponente1;

//jgo  extern double factorial(int);
  extern double constc(int, double, double*, int*, double*, double*);
  extern double intc(int, double, double, double*, double*);
  extern double i_conf(int b, int c, int mu, int nu, int lamb,
                       int sig, int *np, double *expo, double Rc, double*, double*);

  entero1 = 1 - l;
  entero2 = 2 + l;

  doble1 = i_conf(entero1, entero2, mu, nu, lam, sig, np, expo, Rc, arreglo_factorial,
                  arreglo_inv_factorial);

  entero3 = np[mu] + np[nu] + l;
  exponente1 = expo[mu] + expo[nu];
  temp1 = intc(entero3, exponente1, Rc, arreglo_factorial, arreglo_inv_factorial);
  entero3 = np[mu] + np[nu] + l + 1;
//jgo  temp1 = temp1 - 2.f*intc(entero3, exponente1, Rc, arreglo_factorial, arreglo_inv_factorial)/Rc;
  temp1 = temp1 - 2.f*intc(entero3, exponente1, Rc, arreglo_factorial, arreglo_inv_factorial)*inv_Rc;
  entero3 = np[mu] + np[nu] + l + 2;
//jgo   temp1 = temp1 + intc(entero3, exponente1, Rc, arreglo_factorial, arreglo_inv_factorial)/(Rc*Rc);
  temp1 = temp1 + intc(entero3, exponente1, Rc, arreglo_factorial, arreglo_inv_factorial)*inv_Rc*inv_Rc;
  doble2 = constc(mu, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial)*
           constc(nu, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial)*temp1;

  entero3 = np[lam] + np[sig] - l - 1;
  exponente1 = expo[lam] + expo[sig];
  temp1 = intc(entero3, exponente1, Rc, arreglo_factorial, arreglo_inv_factorial);
  entero3 = np[lam] + np[sig] - l;
//jgo  temp1 = temp1 - 2.f*intc(entero3, exponente1, Rc, arreglo_factorial, arreglo_inv_factorial)/Rc;
  temp1 = temp1 - 2.f*intc(entero3, exponente1, Rc, arreglo_factorial, arreglo_inv_factorial)*inv_Rc;
  entero3 = np[lam] + np[sig] - l + 1;
//jgo  temp1 = temp1 + intc(entero3, exponente1, Rc, arreglo_factorial, arreglo_inv_factorial)/(Rc*Rc);
  temp1 = temp1 + intc(entero3, exponente1, Rc, arreglo_factorial, arreglo_inv_factorial)*inv_Rc*inv_Rc;
  doble2 = doble2*constc(lam, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial)*
           constc(sig, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial)*temp1;

  doble3 = i_conf(entero2, entero1, mu, nu, lam, sig, np, expo, Rc, arreglo_factorial,
                  arreglo_inv_factorial);

  total = doble1 + doble2 - doble3;

  return(total);
 }

