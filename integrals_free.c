// One- and two-electron integrals
// to solve HF equations
// for free atoms.
// Jorge Garza, June/2006
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>  // important for complex numbers



double bielectronic_integral_free(int     ang, 
				  int     mu, 
			   	  int     nu, 
				  int     lam, 
				  int     sig, 
			          int    *np, 
				  double *expo,
                                  double *arreglo_factorial, 
				  double *arreglo_inv_factorial)
{
 int i, ent1, ent2;
 double num1, num2, num3, total;

 extern double intl(int, int, int, double*, int*, double*);
 extern double int6l(int, int, int, int, int, int, int*, double*, double*, double*);

 ent1 = 1 - ang;
 ent2 = 2 + ang;
 num1 = int6l(mu, nu, lam, sig, ent1, ent2, np, expo, arreglo_factorial,
              arreglo_inv_factorial); 
 num2 = intl(mu, nu, ent2, expo, np, arreglo_factorial); 
 num2 = num2*intl(lam, sig, ent1, expo, np, arreglo_factorial); 
 num3 = int6l(mu, nu, lam, sig, ent2, ent1, np, expo, arreglo_factorial,
              arreglo_inv_factorial); 
 total = num1 + num2 - num3;
 return(total);
 }

double int6l(int mu, int nu, int lam, int sig, int a, int b, int *np, double *expo,
             double *arreglo_factorial, double *arreglo_inv_factorial)
{
 double num1, num3, num4, num5, num7, num8, suma, cont, coef1, coef2, total,
        result1, result2;
 int i, num2, num6, num9;
 extern double int1l(int, double, double*);
 extern int constant_normalization_free(int*, double*, int*, double*);
 extern double factorial(int);

 constant_normalization_free(&mu, expo, np, &result1);
 constant_normalization_free(&nu, expo, np, &result2);
 num1 = result1*result2;
 constant_normalization_free(&lam, expo, np, &result1);
 constant_normalization_free(&sig, expo, np, &result2);
 num1 = num1*result1*result2;
 num2 = b + np[lam] + np[sig] - 2;
 num3 = factorial(num2);
 num4 = (double)(num2 + 1);
 num5 = expo[lam] + expo[sig];
 coef1 = num1*num3/pow(num5, num4);
 num6 = np[mu] + np[nu] - 2 + a;
 num7 = expo[mu] + expo[nu];
 num8 = num5 + num7;
 suma = (double)0.0;
 for (i = 0; i<= num2; i++)
 {
  cont = (double)i;
// jgo  coef2 = pow(num5, cont)/factorial(i);
  coef2 = pow(num5, cont)*arreglo_inv_factorial[i];
  num9 = num6 + i;
  suma = suma + coef2*int1l(num9, num8, arreglo_factorial);
  }
 total = coef1*(int1l(num6, num7, arreglo_factorial) - suma);
 return(total);
 }

//jgo  double intl_Rc(int mu, int nu, int k, double* expo, int* np, double Rc)
//jgo  {
//jgo   int arg1;
//jgo   double arg2, final, result1, result2;
//jgo   extern double intc(int a, double b, double r, double* double*);
//jgo   extern int constant_normalization_free(int*, double*, int*, double*);
//jgo  
//jgo   arg1 = k + np[mu] + np[nu] - 2;
//jgo   arg2 = expo[mu] + expo[nu];
//jgo   constant_normalization_free(&mu, expo, np, &result1);
//jgo   constant_normalization_free(&nu, expo, np, &result2);
//jgo   final = result1*result2*intc(arg1, arg2, Rc);
//jgo  
//jgo   return(final);
//jgo   }

double intl(int mu, int nu, int k, double* expo, int* np, double *arreglo_factorial)
{
 int arg1;
 double arg2, final, result1, result2;
 extern double int1l(int, double, double*);
 extern int constant_normalization_free(int*, double*, int*, double*);

 arg1 = k + np[mu] + np[nu] - 2;
 arg2 = expo[mu] + expo[nu];
 constant_normalization_free(&mu, expo, np, &result1);
 constant_normalization_free(&nu, expo, np, &result2);
 final = result1*result2*int1l(arg1, arg2, arreglo_factorial);

 return(final);
 }

double int1l(int arg1, double arg2, double *arreglo_factorial)
{
//
// Evaluaci'on de la integral /int_0^{\infty} dr r^a e^{-b r}
//
 double numerador, denominador, arg3;
 extern double factorial(int);

 if (arg1 == 0)
   numerador = (double)1.0;
 else
//jgo   numerador = factorial(arg1);
   numerador = arreglo_factorial[arg1];
 arg3 = (double)(arg1 + 1);
 denominador = pow(arg2, arg3);

 return(numerador/denominador);
 }

double local_gamma(int a, double *arreglo_factorial)
{
 double total;
 extern double int1l(int arg1, double arg2, double *arreglo_factorial);

 if (a > 0) 
   total = int1l(a - 1, 1.f, arreglo_factorial);
 else {
   printf("ERROR ** Gamma function evaluated with negative exponent!!!\n");
   total = -1.f;
 }

 return total;
 }

//jgo double local_incomplete_gamma(int a, double r, double* arreglo_factorial)
//jgo {
//jgo  double total;
//jgo  extern double intc(int a, double b, double r);
//jgo  extern double local_gamma(int a, double* arreglo_factorial);
//jgo 
//jgo  if (a > 0)
//jgo    total = local_gamma(a, arreglo_factorial) - intc(a - 1, 1.f, r);
//jgo  else {
//jgo    printf("ERROR ** Gamma function evaluated with negative exponent!!!\n");
//jgo    total = -1.f;
//jgo  }
//jgo 
//jgo  return total;
//jgo  }
 


double factorial(int argumento)
{
 int contador;
 double total;
 total = (double)1.0;
 if (argumento != 0)
   for (contador = argumento; contador >= 1; contador--)
     total = ((double )contador)*total;
 return(total);
 }

double lineal(double pendiente, double ordenada, double x)
{
 double total;
 total = pendiente*x + ordenada;
 return(total);
}

double upper_incomplete_gamma(double Rc, int a, double b)
{
// Evaluation of 
//               \int_{R_c}^{\infty} r^{-a} exp(-b r)
// The method is based on
// the gauss quadrature with seven points in the interval [-1,1]
// and an adaptative procedure which is controled by 
// the parameter tol.
//
 int mu, espacios, elements, i, j;
 double  liminf, step, x1, x2, x3, x4, x5, x6, x7,
         w1, w2, w3, w4, w5, w6, w7, limsup,
         r1, r2, r3, r4, r5, r6, r7,
         f1, f2, f3, f4, f5, f6, f7,
         xinf, xsup, pendiente, ordenada, 
         integral, int1, tol1, dif, integralvieja, result,
         dif2, integraltotal, integralparcial, tol2, pi;
 extern double lineal(double, double, double);
 
if(a == 0)
{
integraltotal = exp(-b*Rc)/b;
}
else
{
 tol1 = 1.0e-15;
 tol2 = 1.0e-15;
 dif2 = 100.f;
 x1 = -0.9491079123427600f;
 w1 =  0.1294849661688884f;
 x2 = -0.7415311855993934f;
 w2 =  0.2797053914892702f;
 x3 = -0.4058451513773973f;
 w3 =  0.3818300505051192f;
 x4 =  0.0000000000000000f;
 w4 =  0.4179591836734694f;
 x5 = -x3;
 w5 =  w3;
 x6 = -x2;
 w6 =  w2;
 x7 = -x1;
 w7 =  w1;


 liminf = Rc;
 limsup = Rc + 80.f;
 integraltotal = 0.f;
 double decae = -(double) a;
 do {
   integralvieja = 0.f;
   i = 0;
   dif = 100.f;
   do {
     i = i + 1;
     espacios = 400 + 10*(i - 1);
     step = (limsup - liminf)/(double)(espacios);
     integral = 0.f;
     for (j = 1; j <= espacios; j++) {
       xinf = liminf + step*(double)(j - 1);
       xsup = liminf + step*(double)(j);
       pendiente = (xsup - xinf)/2.f;
       ordenada = (xsup + xinf)/2.f;
       r1 = lineal(pendiente, ordenada, x1);
       r2 = lineal(pendiente, ordenada, x2);
       r3 = lineal(pendiente, ordenada, x3);
       r4 = lineal(pendiente, ordenada, x4);
       r5 = lineal(pendiente, ordenada, x5);
       r6 = lineal(pendiente, ordenada, x6);
       r7 = lineal(pendiente, ordenada, x7);
       f1 = pow(r1,decae)*exp(-b*r1);
       f2 = pow(r2,decae)*exp(-b*r2);
       f3 = pow(r3,decae)*exp(-b*r3);
       f4 = pow(r4,decae)*exp(-b*r4);
       f5 = pow(r5,decae)*exp(-b*r5);
       f6 = pow(r6,decae)*exp(-b*r6);
       f7 = pow(r7,decae)*exp(-b*r7);

       int1 = pendiente*(w1*f1 + w2*f2 + w3*f3 + w4*f4 + w5*f5 + w6*f6 + w7*f7);
       integral = integral + int1;
     } // End For
     dif = fabs(integral - integralvieja);
     integralvieja = integral;
    } while (dif >= tol1 && i < 500); //End do 1
    integralparcial = integral;
    dif2 = fabs(integralparcial - integraltotal);
    limsup = limsup + 5.f;
    integraltotal = integralparcial;
 } while (dif2 > tol2 && limsup < 120.f); //End do 2    
// printf("Espacios en el grid  = %12d\n", espacios);
// printf("Infinito practico    = %12.5f\n", limsup);
// printf("integral             = %12.8f\n",integraltotal);
//printf("\nPROBAR_UPPER_INTEGRAL %1.18f\n", integraltotal);
} //End of else for a. 
 return(integraltotal);
 }


double use_upper_incomplete_gamma(double Rc, int a1, int a2, double b1, double b2)
{
// Evaluation of 
//               \int_{R_c}^{\infty} r^{-a1} exp(-b1 r) Int2(a2, b2, r1)
// The method is based on
// the gauss quadrature with seven points in the interval [-1,1]
// and an adaptative procedure which is controled by 
// the parameter tol.
//
 int mu, espacios, elements, i, j;
 double  liminf, step, x1, x2, x3, x4, x5, x6, x7,
         w1, w2, w3, w4, w5, w6, w7, limsup,
         r1, r2, r3, r4, r5, r6, r7,
         f1, f2, f3, f4, f5, f6, f7,
         xinf, xsup, pendiente, ordenada, 
         integral, int1, tol1, dif, integralvieja, result,
         dif2, integraltotal, integralparcial, tol2, pi;
 extern double lineal(double, double, double);
 extern double upper_incomplete_gamma(double Rc, int a, double b);

 if (a1 == 0) {
   integraltotal = exp(-b1*Rc)/b1;
 } else {
 tol1 = 1.0e-16;
 tol2 = 1.0e-16;
 dif2 = 100.f;
 x1 = -0.9491079123427600f;
 w1 =  0.1294849661688884f;
 x2 = -0.7415311855993934f;
 w2 =  0.2797053914892702f;
 x3 = -0.4058451513773973f;
 w3 =  0.3818300505051192f;
 x4 =  0.0000000000000000f;
 w4 =  0.4179591836734694f;
 x5 = -x3;
 w5 =  w3;
 x6 = -x2;
 w6 =  w2;
 x7 = -x1;
 w7 =  w1;
 

//printf("jgo: a1 = %d :: a2 = %d\n", a1, a2);
//printf("jgo: a1 = %f :: a2 = %f\n", (double) a1,(double) a2);


 liminf = Rc;
 limsup = Rc + 50.f;
 integraltotal = 0.f;
 double decae = -(double) a1;
 do {
   integralvieja = 0.f;
   i = 0;
   dif = 100.f;
   do {
     i = i + 1;
     espacios = 400 + 10*(i - 1);
     step = (limsup - liminf)/(double)(espacios);
     integral = 0.f;
     for (j = 1; j <= espacios; j++) {
       xinf = liminf + step*(double)(j - 1);
       xsup = liminf + step*(double)(j);
       pendiente = (xsup - xinf)/2.f;
       ordenada = (xsup + xinf)/2.f;
       r1 = lineal(pendiente, ordenada, x1);
       r2 = lineal(pendiente, ordenada, x2);
       r3 = lineal(pendiente, ordenada, x3);
       r4 = lineal(pendiente, ordenada, x4);
       r5 = lineal(pendiente, ordenada, x5);
       r6 = lineal(pendiente, ordenada, x6);
       r7 = lineal(pendiente, ordenada, x7);
       f1 = pow(r1,decae)*exp(-b1*r1)*upper_incomplete_gamma(r1, a2, b2);
       f2 = pow(r2,decae)*exp(-b1*r2)*upper_incomplete_gamma(r1, a2, b2);
       f3 = pow(r3,decae)*exp(-b1*r3)*upper_incomplete_gamma(r1, a2, b2);
       f4 = pow(r4,decae)*exp(-b1*r4)*upper_incomplete_gamma(r1, a2, b2);
       f5 = pow(r5,decae)*exp(-b1*r5)*upper_incomplete_gamma(r1, a2, b2);
       f6 = pow(r6,decae)*exp(-b1*r6)*upper_incomplete_gamma(r1, a2, b2);
       f7 = pow(r7,decae)*exp(-b1*r7)*upper_incomplete_gamma(r1, a2, b2);
       int1 = pendiente*(w1*f1 + w2*f2 + w3*f3 + w4*f4 + w5*f5 + w6*f6 + w7*f7);
       integral = integral + int1;
     } // End For
     dif = fabs(integral - integralvieja);
     integralvieja = integral;
    } while (dif >= tol1 && i < 500); //End do 1
    integralparcial = integral;
    dif2 = fabs(integralparcial - integraltotal);
    limsup = limsup + 5.f;
    integraltotal = integralparcial;
 } while (dif2 > tol2 && limsup < 120.f); //End do 2
// printf("Espacios en el grid  = %12d\n", espacios);
// printf("Infinito practico    = %12.5f\n", limsup);
// printf("integral  = %12.8f\n", integraltotal);
 }  // End of else for a1.
 return(integraltotal);
 }
/* ------------------------------------------------------------------------------------------ */
/* Mike */
/* Here I am going to start with the with the necessary integrals for the GTO basis functions */
 long long int factorial_mike(int n) {
        int i;
        long long int z;
        i = 1;
        z = 1;
        if(n == 0) {
           z = 1;
        }
        else{
           while(i <= n){
              z = z*i;
              i++;
           }
        }
        return(z);
 }
/* This function calculates the integral \int_{0}^{\infty} t^{(n/2)-1} e^{-t} dt = Gamma(n/2) = [(n/2) - 1]! for n = 1,2,3,4,5,... */
/* I have decided to insert this function due to its frequency with the matrix elements */
  double full_gamma_arg_2( int n1) {
        extern long long int factorial_mike(int );
        long long int factor2n;
        double gamma;
        double argument;
        int n2;
        double theta, pi;
        int n_0_1;

        pi = ((double) 4)*atan(1.f);
        argument = ((double)n1)/((double) 2);
        theta = argument*pi;
        n_0_1 = ((int) sin(theta));
        /* los únicos valores que puede tomar esta cantidad son "0" cuando (n_0_1/2)= entero; 1 ó -1 cuando (n_0_1/2) = 1/2, 3/2, etc. */
        if(n_0_1 == 0){
           n2 = ((int) n1/2);
           n2 = n2 - 1;
           gamma = (double) factorial_mike(n2);
        }
        else{
           n2 = n1 - 1;
           n2 = ((int) n2/2);
           factor2n = pow(2, n1 - 1);
           gamma = (double) factorial_mike(n1 - 1);
           gamma = gamma*sqrt(pi);
           gamma = gamma/((double) factor2n);
           gamma = gamma/((double) factorial_mike(n2));
        }
        return gamma;
 }
 /* This function calculates the integral \int_{0}^{y} t^{(n/2)-1} e^{-t} dt = \gamma(n/2,y) = lower incomplete gamma function n >= 2, n = 2,3,4,5,6,7,8,...  */
 double lower_incomplete_gamma_function(int n2, double y){
        int i, j, iglobal;
        double rho_j;
        double r0, r1, r[691], f_r[691];
        double factor1;
        int key;
        double pi, result;
        double a, b, c, d, f0, f1, f2, f0_2, f1_2, f2_2, x0, x1, x2, x0_2, x1_2, x2_2, c13, c12;
        factor1 = ((double) n2)/((double) 2);
        factor1 = factor1 - ((double) 1);
        key = 0;
        c12 = ((double) 1)/((double) 2);
        c13 = ((double) 1)/((double) 3);
        pi = ((double) 4)*atan(1.f);
 if(n2 == 1){ /* se tiene el caso \int_{0}^{y} t^{-1/2} e^{-t} dt = \gamma(1/2,y)  */
    result = sqrt(pi)*erf(sqrt(y));
 }
 else
    if(n2 == 2){ /* se tiene \int_{0}^y e^{-t} dt = \gamma(1,y); la solución es analítica y no hay necesidad de generar la malla */
      result = ((double) 1) - exp(-y);
    }
    else
       if(n2 == 3){   /* se tiene \int_{0}^y t^{1/2} e^{-t} dt = \gamma(3/2,y) */
          result = 0.5f*sqrt(pi)*erf(sqrt(y)) - sqrt(y)*exp(-y);
       }
       else
          if(n2 == 4){  /* se tiene \int_{0}^y t^{1} e^{-t} dt = \gamma(2,y) */
             result = ((double)1) - exp(-y)*(y + (double)1);
          }
          else
             if(n2 == 5){  /* se tiene \int_{0}^y t^{3/2} e^{-t} dt = \gamma(5/2,y) */
                result = 0.75f*sqrt(pi)*erf(sqrt(y)) - sqrt(y)*exp(-y)*(y + 1.5f);
             }   /* here we are ;) */
            else{ /* -------------------- Here begins the numerical integration -------------------- */
           for(j = 0; j <= 690; j++){ /* ---------- Here begins the construction of the grid ---------- */
              rho_j = ((double) j - 1);  //el nombre de esta variable no tiene nada que ver con la densidad electrónica
              rho_j = rho_j/((double) 30);
              rho_j = rho_j - ((double) 10);
              r0 = exp(rho_j);
              if(j == 0){
                 r[0] = ((double) 0);
                 f_r[0] = ((double) 0);
              }
              else{
                 r[j] = r0; //reasigno para crear el arreglo de la malla
                 f_r[j] = pow(r0,factor1)*exp(-r0); //para evaluar la función en cada punto de la malla
              }
              if(key == 0){
                 if(r0 > y){
                    if(j%2 != 0){ //para j "impar"
                       j = j + 1;
                       r[j] = y;
                       f_r[j] = pow(y,factor1)*exp(-y);
                       r1 = 0.5f*(y + r[j - 2]);
                       r[j - 1] = r1;
                       f_r[j - 1] = pow(r1,factor1)*exp(-r1);
                    }
                    else{   // para j "par"
                       r[j] = y;
                       f_r[j] = pow(y,factor1)*exp(-y);
                    }
                    key = 1; //reasignación
                 }
              }  //se controla que la malla pase exactamente por "y" y que el arreglo siempre termine siento un "j"---> par
              else{
                 f_r[j] = ((double) 0);  //no controlo hasta donde llega la malla, pero si que la función automáticamente sea cero más alla del límite de "y"
              }
           } /* ---------- Here ends the construction of the grid ---------- */
           result = ((double) 0);   //inicializamos result
           for(i = 0; i <= 344; i++){
              iglobal = 2*(i + 1);

              x0 = r[iglobal - 2];
              x1 = r[iglobal - 1];
              x2 = r[iglobal];

              x0_2 = pow(x0, 2.f);
              x1_2 = pow(x1, 2.f);
              x2_2 = pow(x2, 2.f);

              f0 = f_r[iglobal - 2];
              f1 = f_r[iglobal - 1];
              f2 = f_r[iglobal];

              f0_2 = pow(f0, 2.f);
              f1_2 = pow(f1, 2.f);
              f2_2 = pow(f2, 2.f);

              d = (x0 - x1)*(x0 - x2)*(x1 - x2);

              if(f2 == (double) 0){
                 a = ((double) 0);
                 b = ((double) 0);
                 c = ((double) 0);
              }
              else{
                 a = ((x1 - x2)*f0 + (x2 - x0)*f1 + (x0 - x1)*f2)/d;

                 b = ((f1 - f2)*x0_2 + (f2 - f0)*x1_2 + (f0 - f1)*x2_2)/d;

                 c = ((f2*x1 - f1*x2)*x0_2 + (f0*x2 - f2*x0)*x1_2 + (f1*x0 - f0*x1)*x2_2)/d;

                 result = result + c13*a*(x2*x2_2 - x0*x0_2) + c12*b*(x2_2 - x0_2) + c*(x2 - x0);
              }
//              printf("%20.15f \n", result);
           }
        } /* -------------------- Here ends the numerical integration -------------------- */
        return result;
 }

double upper_incomplete_gamma_5_points(double Rc, int a, double b, int n_points, int save_i, double *grid) {
// Evaluation of 
//               \int_{R_c}^{\infty} r^{-a} exp(-b r)
 int i, indice, blocks;
 double x0, x1, x2, x3, x4, f0, f1, f2, f3, f4;
 double decae, total_integral;

 extern double integral_five_points(double x_0, double x_1, double x_2, double x_3, double x_4,
                                    double f_0, double f_1, double f_2, double f_3, double f_4);

 if(a == 0)
   total_integral = exp(-b*Rc)/b;
 else {
   decae = -(double) a;
   blocks = (n_points - save_i + 1)/5;

   total_integral = 0.f;
   for (i = 0; i < blocks; i++) {
       indice = 4*i + save_i;
       x0 = grid[indice];     f0 = pow(x0,decae)*exp(-b*x0);
       x1 = grid[indice + 1]; f1 = pow(x1,decae)*exp(-b*x1);
       x2 = grid[indice + 2]; f2 = pow(x2,decae)*exp(-b*x2);
       x3 = grid[indice + 3]; f3 = pow(x3,decae)*exp(-b*x3);
       x4 = grid[indice + 4]; f4 = pow(x4,decae)*exp(-b*x4);
       total_integral = total_integral + integral_five_points(x0, x1, x2, x3, x4, f0, f1, f2, f3, f4);
   }
 }

 return(total_integral);
}

double use_upper_incomplete_gamma_5_points(double Rc, int a1, int a2, double b1, double b2, int n_points, int save_i, double *grid) {
// Evaluation of 
//               \int_{R_c}^{\infty} r^{-a1} exp(-b1 r) Int2(a2, b2, r1)
 int i, indice, blocks;
 double x0, x1, x2, x3, x4, f0, f1, f2, f3, f4;
 double decae, total_integral;

 extern double integral_five_points(double x_0, double x_1, double x_2, double x_3, double x_4,
                                    double f_0, double f_1, double f_2, double f_3, double f_4);
 extern double upper_incomplete_gamma_5_points(double Rc, int a, double b, int n_points, int save_i, double *grid);

 if (a1 == 0) {
   total_integral = exp(-b1*Rc)/b1;
 } else {
   decae = -(double) a1;
   blocks = (n_points - save_i + 1)/5;

   total_integral = 0.f;
   for (i = 0; i < blocks; i++) {
       indice = 4*i + save_i;
       x0 = grid[indice];     f0 = pow(x0,decae)*exp(-b1*x0)*upper_incomplete_gamma_5_points(x0, a2, b2, n_points, save_i, grid);
       x1 = grid[indice + 1]; f1 = pow(x1,decae)*exp(-b1*x1)*upper_incomplete_gamma_5_points(x1, a2, b2, n_points, save_i, grid);
       x2 = grid[indice + 2]; f2 = pow(x2,decae)*exp(-b1*x2)*upper_incomplete_gamma_5_points(x2, a2, b2, n_points, save_i, grid);
       x3 = grid[indice + 3]; f3 = pow(x3,decae)*exp(-b1*x3)*upper_incomplete_gamma_5_points(x3, a2, b2, n_points, save_i, grid);
       x4 = grid[indice + 4]; f4 = pow(x4,decae)*exp(-b1*x4)*upper_incomplete_gamma_5_points(x4, a2, b2, n_points, save_i, grid);
       total_integral = total_integral + integral_five_points(x0, x1, x2, x3, x4, f0, f1, f2, f3, f4);
   }
 }
 return(total_integral);
}

// Function to evaluate numerically \int_{x_0}^{x_4} f(x) dx with five points
//
double integral_five_points(double x_0, double x_1, double x_2, double x_3, double x_4,
                           double f_0, double f_1, double f_2, double f_3, double f_4) {
  double total, num1, num2, num3, num4, num5, denom1, denom2, denom3, denom4, denom5,
         term1, term2, term3, term4, term5;
  double dif_x4x0, dif_x4x1, dif_x4x3, dif_x3x2, dif_x3x1, dif_x3x0, dif_x4x2, dif_x2x1, dif_x2x0,
	 dif_x1x0;

  dif_x4x0 = x_4 - x_0;
  dif_x4x1 = x_4 - x_1;
  dif_x4x2 = x_4 - x_2;
  dif_x4x3 = x_4 - x_3;
  dif_x3x2 = x_3 - x_2;
  dif_x3x1 = x_3 - x_1;
  dif_x3x0 = x_3 - x_0;
  dif_x2x1 = x_2 - x_1;
  dif_x2x0 = x_2 - x_0;
  dif_x1x0 = x_1 - x_0;

  num1 = 3.f*x_0*x_0 + 10.f*x_1*x_2 - 5.f*x_0*(x_1 + x_2) + 4.f*x_0*x_4 - 5.f*(x_1 + x_2)*x_4 + 3.f*x_4*x_4;
  num1 = f_3*dif_x4x0*dif_x4x0*num1;
  denom1 = dif_x3x0*dif_x3x1*dif_x3x2*dif_x4x3;
  term1 = -num1/denom1;

  num2 = 3.f*x_0*x_0 + 10.f*x_1*x_3 - 5.f*x_0*(x_1 + x_3) + 4.f*x_0*x_4 - 5.f*(x_1 + x_3)*x_4 + 3.f*x_4*x_4;
  num2 = f_2*dif_x4x0*dif_x4x0*num2;
  denom2 = -dif_x2x0*dif_x2x1*dif_x3x2*dif_x4x2;
  term2 = -num2/denom2;

  num3 = 3.f*x_0*x_0 + 10.f*x_2*x_3 - 5.f*x_0*(x_2 + x_3) + 4.f*x_0*x_4 - 5.f*(x_2 + x_3)*x_4 + 3.f*x_4*x_4;
  num3 = f_1*dif_x4x0*dif_x4x0*num3;
  denom3 = dif_x1x0*dif_x2x1*dif_x3x1*dif_x4x1;
  term3 = -num3/denom3;

  num4 = 12.f*x_0*x_0*x_0 - 30.f*x_1*x_2*x_3 + 10.f*(x_2*x_3 + x_1*(x_2 + x_3))*x_4;
  num4 = num4 - 5.f*(x_1 + x_2 + x_3)*x_4*x_4 + 3.f*x_4*x_4*x_4 + 3.f*x_0*x_0*(-5.f*(x_1 + x_2 + x_3) + 3.f*x_4);
  num4 = num4 + 2.f*x_0*(10.f*(x_2*x_3 + x_1*(x_2 + x_3)) - 5.f*(x_1 + x_2 + x_3)*x_4 + 3.f*x_4*x_4);
  num4 = f_0*num4;
  denom4 = -dif_x1x0*dif_x2x0*dif_x3x0;
  term4 = -num4/denom4;

  num5 = 3.f*x_0*x_0*x_0 - 30.f*x_1*x_2*x_3 + 20.f*(x_2*x_3 + x_1*(x_2 + x_3))*x_4;
  num5 = num5 - 15.f*(x_1 + x_2 + x_3)*x_4*x_4 + 12.f*x_4*x_4*x_4 + x_0*x_0*(-5.f*(x_1 + x_2 + x_3) + 6.f*x_4);
  num5 = num5 + x_0*(10.f*(x_2*x_3 + x_1*(x_2 + x_3)) - 10.f*(x_1 + x_2 + x_3)*x_4 + 9.f*x_4*x_4);
  num5 = f_4*num5;
  denom5 = -dif_x4x1*dif_x4x2*dif_x4x3;
  term5 = num5/denom5;

  total = -dif_x4x0/60.f;
  total = total*(term1 + term2 + term3 + term4 + term5);
  return (total);
  }
/* ------------------------------------------------------------------------------------------ */
/* This is for the Debye plasma (only exponential without cos(r - r') ) */
/* The basis set functions for this case are only STOs */
/* Mike Nov 5, 2019 */
 double h_f1(int ang, int mu, int nu, int lam, int sig, int *np, double *expo, double alpha, int m, int k){
         extern double factorial(int );
         int s;
         double beta, b2, c0, c1, c2, n1, n2, e1, e2, hs;
         // parameters
         beta = (double) 1/alpha;
         b2 = (double) 2*beta;
         c0 = (double) -1;
         c1 = pow(c0, (double) k);
         c2 = pow(c0, (double) ang + 1);
         n1 = (double) np[mu] + np[nu];
         n2 = (double) np[sig] + np[lam];
         e1 = expo[mu] + expo[nu];
         e2 = e1 + expo[sig] + expo[lam];
         hs = (double) 0;
         for(s = 0; s <= np[mu] + np[nu] - k - 1; s++){
                 hs = hs + (factorial(np[sig] + np[lam] + s - m - 1)/factorial(s))*(c1/(pow(e1 - beta, (double) n1 - s - k)*pow(e2, (double) n2 + s - m)) +
                                 c2/(pow(e1 + beta, (double) n1 - s - k)*pow(e2 + b2, (double) n2 + s - m)));
         }
         return(hs);
  }

 double h_f2(int ang, int mu, int nu, int lam, int sig, int *np, double *expo, double alpha, int m){
         extern double h_f1(int , int , int , int , int , int *, double *, double , int , int );
         extern double factorial(int );
         int k;
         double beta, b2, c0, c1, c2, c3, n1, n2, e1, e2, hs;
         // parameters
         beta = (double) 1/alpha;
         b2 = (double) 2*beta;
         c0 = (double) -1;
         c2 = pow(c0, ang + 1);
         n1 = (double) np[mu] + np[nu];
         n2 = (double) np[sig] + np[lam];
         e1 = expo[mu] + expo[nu];
         e2 = expo[sig] + expo[lam];
         hs = (double) 0;
         for(k = 0; k <= ang; k++){
                 c1 = pow(c0, (double) k);
                 c3 = pow(b2, (double) k + 1);
                 hs = hs + ((factorial(ang + k)*factorial(np[mu] + np[nu] - k - 1))/(factorial(k)*factorial(ang - k)*c3))*
                         ((factorial(np[sig] + np[lam] - m - 1)/pow(e2 + beta, (double) n2 - m))*
                          (c1/pow(e1 - beta, (double) n1 - k) + c2/pow(e1 + beta, (double) n1 - k)) - h_f1(ang, mu, nu, lam, sig, np, expo, alpha, m, k));
         }
         return(hs);
 }

 double h_l(int ang, int mu, int nu, int lam, int sig, int *np, double *expo, double alpha){
         extern double h_f2(int , int , int , int , int , int *, double *, double , int );
         extern double factorial(int );
         int m;
         double beta, b2, b3, hs;
         // parameters
         beta = (double) 1/alpha;
         b2 = (double) 2*beta;
         hs = (double) 0;
         for(m = 0; m <= ang; m++){
                 b3 = pow(b2,m);
                 hs = hs + (factorial(ang + m)/(factorial(m)*factorial(ang - m)*b3*beta))*h_f2(ang, mu, nu, lam, sig, np, expo, alpha, m);
         }
         return(hs);
 } // This is ready 

 double g_f1(int ang, int mu, int nu, int lam, int sig, int *np, double *expo, double alpha, int m, int k){
         extern double factorial(int );
         int s;
         double beta, b2, c0, c1, c2, n2, e1, e2, hs;
         // parameters
         beta = (double) 1/alpha;
         b2 = (double) 2*beta;
         c0 = (double) -1;
         c1 = pow(c0, (double) m);
         c2 = pow(c0, (double) ang + 1);
         n2 = (double) np[sig] + np[lam];
         e1 = expo[mu] + expo[nu];
         e2 = e1 + expo[sig] + expo[lam];
         hs = (double) 0;
         for(s = 0; s <= np[mu] + np[nu] - k - 1; s++){
                 hs = hs + ((pow(e1 + beta, s)*factorial(np[sig] + np[lam] + s - m - 1))/factorial(s))*
                         (c1/pow(e2, n2 + s - m) + c2/pow(e2 + b2, n2 + s - m));
         }
         return(hs);
 }

 double g_f2(int ang, int mu, int nu, int lam, int sig, int *np, double *expo, double alpha, int m){
         extern double g_f1(int , int , int , int , int , int *, double *, double , int , int );
         extern double factorial(int );
         int k;
         double beta, b3, c3, n1, e1, hs;
         // parameters
         beta = (double) 1/alpha;
         n1 = (double) np[mu] + np[nu];
         e1 = expo[mu] + expo[nu];
         hs = (double) 0;
         for(k = 0; k <= ang; k++){
                 c3 = pow((double) 2,k);
                 b3 = pow(beta, (double) k + 1);
                 hs = hs + ((factorial(ang + k)*factorial(np[mu] + np[nu] - k - 1))/(factorial(k)*factorial(ang - k)*c3*b3*pow(e1 + beta, n1 - k)))*
                         g_f1(ang, mu, nu, lam, sig, np, expo, alpha, m, k);
         }
         return(hs);
 }

 double g_l(int ang, int mu, int nu, int lam, int sig, int *np, double *expo, double alpha){
         extern double g_f2(int , int , int , int , int , int *, double *, double , int );
         extern double factorial(int );
         int m;
         double beta, b2, c1, hs;
         // parameters
         beta = (double) 1/alpha;
         b2 = (double) 2*beta;
         hs = (double) 0;
         for(m = 0; m <= ang; m++){
                 c1 = pow(b2, (double) m + 1);
                 hs = hs + (factorial(ang + m)/(factorial(m)*factorial(ang - m)*c1))*g_f2(ang, mu, nu, lam, sig, np, expo, alpha, m);
         }
         return(hs);
 } // This is ready

 double bielectronic_integral_debye(int ang, int mu, int nu, int lam, int sig, int *np, double *expo, double alpha){
         extern double h_l(int , int , int , int , int , int *, double *, double );
         extern double g_l(int , int , int , int , int , int *, double *, double );
         extern int constant_normalization_free(int*, double*, int*, double*);
         double n_sto_mu, n_sto_nu, n_sto_lam, n_sto_sig, factor, rad;

         constant_normalization_free(&mu, expo, np, &n_sto_mu);
         constant_normalization_free(&nu, expo, np, &n_sto_nu);
         constant_normalization_free(&lam, expo, np, &n_sto_lam);
         constant_normalization_free(&sig, expo, np, &n_sto_sig);

         factor = ((double) 2*ang + 1)*n_sto_mu*n_sto_nu*n_sto_lam*n_sto_sig/alpha;
         rad = factor*(h_l(ang, mu, nu, lam, sig, np, expo, alpha) + g_l(ang, mu, nu, lam, sig, np, expo, alpha));

         return(rad);
 }
/* ------------------------------------------------------------------------------------------ */
/* This is for the Debye plasma with the factor cos(r - r') ) */
/* The basis set functions for this case are only STOs */
/* Mike Nov 22, 2019 */
  double complex h_f1_c(int ang, int mu, int nu, int lam, int sig, int *np, double *expo, double alpha, int m, int k, int phase){
         extern double factorial(int );
         int s;
         double beta, c0, c1, c2, n1, n2, e1, e2;
         double complex b, b2, hs;
         // convention
         beta = (double) 1/alpha;
         if(phase == 1)
                 b = beta - beta*I;
         else
                 b = beta + beta*I;
         // parameters
         b2 = (double) 2*b;
         c0 = (double) -1;
         c1 = pow(c0, (double) k);
         c2 = pow(c0, (double) ang + 1);
         n1 = (double) np[mu] + np[nu];
         n2 = (double) np[sig] + np[lam];
         e1 = expo[mu] + expo[nu];
         e2 = e1 + expo[sig] + expo[lam];
         hs = (double complex) 0;
         for(s = 0; s <= np[mu] + np[nu] - k - 1; s++){
                 hs = hs + (factorial(np[sig] + np[lam] + s - m - 1)/factorial(s))*(c1/(cpow(e1 - b, (double) n1 - s - k)*pow(e2, (double) n2 + s - m)) +
                                 c2/(cpow(e1 + b, (double) n1 - s - k)*cpow(e2 + b2, (double) n2 + s - m)));
         }
         return(hs);
  }
  double complex h_f2_c(int ang, int mu, int nu, int lam, int sig, int *np, double *expo, double alpha, int m, int phase){
         extern double complex h_f1_c(int , int , int , int , int , int *, double *, double , int , int , int );
         extern double factorial(int );
         int k;
         double beta, c0, c1, c2, n1, n2, e1, e2;
         double complex b, b2, c3, hs;
         // convention
         beta = (double) 1/alpha;
         if(phase == 1)
                 b = beta - beta*I;
         else
                 b = beta + beta*I;
         // parameters
         b2 = (double) 2*b;
         c0 = (double) -1;
         c2 = pow(c0, ang + 1);
         n1 = (double) np[mu] + np[nu];
         n2 = (double) np[sig] + np[lam];
         e1 = expo[mu] + expo[nu];
         e2 = expo[sig] + expo[lam];
         hs = (double complex) 0;
         for(k = 0; k <= ang; k++){
                 c1 = pow(c0, (double) k);
                 c3 = cpow(b2, (double) k + 1);
                 hs = hs + ((factorial(ang + k)*factorial(np[mu] + np[nu] - k - 1))/(factorial(k)*factorial(ang - k)*c3))*
                         ((factorial(np[sig] + np[lam] - m - 1)/cpow(e2 + b, (double) n2 - m))*
                          (c1/cpow(e1 - b, (double) n1 - k) + c2/cpow(e1 + b, (double) n1 - k)) - h_f1_c(ang, mu, nu, lam, sig, np, expo, alpha, m, k, phase));
         }
         return(hs);
 }

  double complex h_l_c(int ang, int mu, int nu, int lam, int sig, int *np, double *expo, double alpha, int phase){
         extern double complex h_f2_c(int , int , int , int , int , int *, double *, double , int , int );
         extern double factorial(int );
         int m;
         double beta;
         double complex b, b2, b3, hs;
         // convention
         beta = (double) 1/alpha;
         if(phase == 1)
                 b = beta - beta*I;
         else
                 b = beta + beta*I;
         // parameters
         b2 = (double) 2*b;
         hs = (double complex) 0;
         for(m = 0; m <= ang; m++){
                 b3 = cpow(b2,m);
                 hs = hs + (factorial(ang + m)/(factorial(m)*factorial(ang - m)*b3*b))*h_f2_c(ang, mu, nu, lam, sig, np, expo, alpha, m, phase);
         }
         return(hs);
 }

 double complex g_f1_c(int ang, int mu, int nu, int lam, int sig, int *np, double *expo, double alpha, int m, int k, int phase){
         extern double factorial(int );
         int s;
         double beta, c0, c1, c2, n2, e1, e2;
         double complex b, b2, hs;
         // convention
         beta = (double) 1/alpha;
         if(phase == 1)
                 b = beta - beta*I;
         else
                 b = beta + beta*I;
         // parameters
         b2 = (double) 2*b;
         c0 = (double) -1;
         c1 = pow(c0, (double) m);
         c2 = pow(c0, (double) ang + 1);
         n2 = (double) np[sig] + np[lam];
         e1 = expo[mu] + expo[nu];
         e2 = e1 + expo[sig] + expo[lam];
         hs = (double complex) 0;
         for(s = 0; s <= np[mu] + np[nu] - k - 1; s++){
                 hs = hs + ((cpow(e1 + b, s)*factorial(np[sig] + np[lam] + s - m - 1))/factorial(s))*
                         (c1/pow(e2, n2 + s - m) + c2/cpow(e2 + b2, n2 + s - m));
         }
         return(hs);
 }

 double complex g_f2_c(int ang, int mu, int nu, int lam, int sig, int *np, double *expo, double alpha, int m, int phase){
         extern double complex g_f1_c(int , int , int , int , int , int *, double *, double , int , int, int );
         extern double factorial(int );
         int k;
         double beta, c3, n1, e1;
         double complex b, b2, b3, hs;
         // convention
         beta = (double) 1/alpha;
         if(phase == 1)
                 b = beta - beta*I;
         else
                 b = beta + beta*I;
         // parameters
         n1 = (double) np[mu] + np[nu];
         e1 = expo[mu] + expo[nu];
         hs = (double) 0 + 0*I;  // ojo
         for(k = 0; k <= ang; k++){
                 c3 = pow((double) 2,k);
                 b3 = cpow(b, (double) k + 1);
                 hs = hs + ((factorial(ang + k)*factorial(np[mu] + np[nu] - k - 1))/(factorial(k)*factorial(ang - k)*c3*b3*cpow(e1 + b, n1 - k)))*
                         g_f1_c(ang, mu, nu, lam, sig, np, expo, alpha, m, k, phase);
         }
         return(hs);
 }

 double complex g_l_c(int ang, int mu, int nu, int lam, int sig, int *np, double *expo, double alpha, int phase){
         extern double complex g_f2_c(int , int , int , int , int , int *, double *, double , int, int );
         extern double factorial(int );
         int m;
         double beta;
         double complex c1, b, b2, hs;
         // convention
         beta = (double) 1/alpha;
         if(phase == 1)
                 b = beta - beta*I;
         else
                 b = beta + beta*I;
         // parameters
         b2 = (double) 2*b;
         hs = (double) 0 + 0*I;  // ojo
         for(m = 0; m <= ang; m++){
                 c1 = cpow(b2, (double) m + 1);
                 hs = hs + (factorial(ang + m)/(factorial(m)*factorial(ang - m)*c1))*g_f2_c(ang, mu, nu, lam, sig, np, expo, alpha, m, phase);
         }
         return(hs);
 }

 double complex bielectronic_integral_yukawa(int ang, int mu, int nu, int lam, int sig, int *np, double *expo, double alpha){
         extern double complex h_l_c(int , int , int , int , int , int *, double *, double, int );
         extern double complex g_l_c(int , int , int , int , int , int *, double *, double, int );
         extern int constant_normalization_free(int*, double*, int*, double*);
         double n_sto_mu, n_sto_nu, n_sto_lam, n_sto_sig, beta, factor, re;
         double complex b1, b2, rad, rad1, rad2;

         constant_normalization_free(&mu, expo, np, &n_sto_mu);
         constant_normalization_free(&nu, expo, np, &n_sto_nu);
         constant_normalization_free(&lam, expo, np, &n_sto_lam);
         constant_normalization_free(&sig, expo, np, &n_sto_sig);
         beta = (double) 1/alpha;
         b1 = beta - beta*I;
         b2 = beta + beta*I;
         factor = ((double) 2*ang + 1)*n_sto_mu*n_sto_nu*n_sto_lam*n_sto_sig/((double) 2);
         rad1 = b1*(h_l_c(ang, mu, nu, lam, sig, np, expo, alpha, 1) + g_l_c(ang, mu, nu, lam, sig, np, expo, alpha, 1));
         rad2 = b2*(h_l_c(ang, mu, nu, lam, sig, np, expo, alpha, 2) + g_l_c(ang, mu, nu, lam, sig, np, expo, alpha, 2));

         rad = factor*(rad1 + rad2);
//       printf("Re = %20.15f Im = %20.15f \n", creal(rad1), cimag(rad1));
//       printf("Im1 = %20.15f Im2 = %20.15f \n", cimag(rad1), cimag(rad2));
         return(rad);
 }


