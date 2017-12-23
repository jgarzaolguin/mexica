//Kinetic energy matrix
//
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
 #include <omp.h>
#else
 #define omp_get_thread_num() 0
#endif

void cinetica(char   *using_gamma,
              int     nt, 
	      double *matk, 
	      double *matkint, 
	      double *matkext, 
	      int    *np, 
	      int    *mang, 
	      int    *ncm, 
              double *expo, 
	      double  Rc,  
              double  gamma_couple,
              char   *bound,
              double *NC_minus, 
	      double *NC_plus, 
	      double *arreglo_factorial,
              double *arreglo_inv_factorial)
{
 int i, j, k;

 int index_i, index_j, total_elements;
 int ang_i, ang_j, ncm_i, ncm_j;
 double delta, delta1;
// double *array_temp, *array_temp_2;

 int entero1;
 double doble1, doble2, num1, num2, num3, num4, num5, num6, result1, result2;

 double total, total1, totaltemp;
 double alpha_nu, alpha_mu, zeta_mu, zeta_nu;
 double zetas, alphas;
 int n_mu, n_nu;
 int enes, eles;
 
 extern int indexes(int, int, int*, int*);
 extern double intl(int, int, int, double* , int*, double*);
 extern double int1c(int, int, int, double, double* , int*, double*, double*);
 extern double intc(int, double, double, double*, double*);
 extern double constc(int, double, double*, int*, double*, double*);
//mrb extern void memoria_double_uni(int, double**, char*);
 extern void memoria_double_uni(int      tamanio,
                                double **arreglo,
                                char    *titulo);

 extern int delta_kro_(int*, int*, double*);
 extern double upper_incomplete_gamma(double, int, double);
//jgo extern double intc(int a, double b, double r);

 total_elements = nt*nt;

#pragma omp parallel shared(total_elements, nt, matk, np, mang, ncm, expo, Rc, gamma_couple, NC_minus, NC_plus)

#pragma omp parallel private(i, j, k, index_i, index_j, ang_i, ang_j, ncm_i, ncm_j, delta, delta1, entero1, doble1, doble2, num1, num2, num3, num4, num5, num6, result1, result2, total, alpha_nu, alpha_mu, zetas, alphas, n_mu, n_nu, enes, eles, total1, zeta_nu, zeta_mu, totaltemp)

{
 int TID = omp_get_thread_num();
 if (strcmp(bound,"free") == 0) {
   #pragma omp for
   for (k = 0; k < total_elements; k++) {
     indexes(nt, k, &index_i, &index_j);
     index_i = index_i - 1;
     index_j = index_j - 1;
     i = index_i;
     j = index_j;
     ang_i = mang[i];
     ang_j = mang[j];
     delta_kro_(&ang_i, &ang_j, &delta);
     ncm_i = ncm[i];
     ncm_j = ncm[j];
     delta_kro_(&ncm_i, &ncm_j, &delta1);
     delta = delta*delta1;
     if (delta == (double)0) 
       matk[k] = (double)0;
     else {
        num1 = expo[j]*expo[j]*intl(i, j, 2, expo, np, arreglo_factorial);
        num2 = ((double)2)*expo[j]*np[j]*intl(i, j, 1, expo, np, arreglo_factorial);
        num3 = (double)(np[j]*(np[j] - 1) - mang[j]*(mang[j] + 1));
        num3 = num3*intl(i, j, 0, expo, np, arreglo_factorial);
        matk[k] = -((double)0.5)*(num1 - num2 + num3);
     }
   }
 }
 else if (strcmp(bound,"finite") == 0) {
   printf("Kinetic contribution to hamiltonian..\n");
   #pragma omp for
   for (k = 0; k <  total_elements; k++) {
     indexes(nt, k, &index_i, &index_j);
     index_i = index_i - 1;
     index_j = index_j - 1;
     ang_i = mang[index_i];
     ang_j = mang[index_j];
     //alpha_mu = expo[index_i];
     //alpha_nu = expo[index_j];
      zeta_mu = expo[index_i];
      zeta_nu = expo[index_j];
      n_mu =np[index_i];
      n_nu =np[index_j];
     delta_kro_(&ang_i, &ang_j, &delta);
     ncm_i = ncm[index_i];
     ncm_j = ncm[index_j];
     delta_kro_(&ncm_i, &ncm_j, &delta1);
     delta = delta*delta1;
     if (delta == (double)0)
       matk[k] = (double)0;
     else {
        
        if (strcmp(using_gamma,"YES") == 0){//label 1
         alpha_mu = zeta_mu  + (double) -1.f*(n_mu + ang_i)/Rc + gamma_couple/(Rc*(1.f - gamma_couple));
         alpha_nu = zeta_nu  + (double) -1.f*(n_nu + ang_j)/Rc + gamma_couple/(Rc*(1.f - gamma_couple));
         alphas = alpha_mu + alpha_nu;
         zetas = zeta_mu + zeta_nu;
         enes = n_mu + n_nu;
         eles = ang_i + ang_j;
         total1 = (double) (Rc*Rc*(n_nu - 1)*n_nu - Rc*Rc*ang_j*(ang_j + 1));
         total = total1*intc(enes - 2, zetas, Rc, arreglo_factorial, arreglo_inv_factorial);
         total1 = (double) (-2.f*Rc*n_nu*(gamma_couple*n_nu + Rc*zeta_nu) + 2.f*Rc*gamma_couple*ang_j*(ang_j + 1));
         total = total + total1*intc(enes - 1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial); 
         total1 = (double) (gamma_couple*gamma_couple*n_nu*n_nu + Rc*zeta_nu*(2.f*gamma_couple + Rc*zeta_nu) +
                            gamma_couple*n_nu*(gamma_couple + 4.f*Rc*zeta_nu) - gamma_couple*gamma_couple*ang_j*(ang_j + 1));
         total = total + total1*intc(enes, zetas, Rc, arreglo_factorial, arreglo_inv_factorial);
         total1 = (double) (-2.f*gamma_couple*zeta_nu*(gamma_couple + gamma_couple*n_nu + Rc*zeta_nu));
         total = total + total1*intc(enes + 1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial);
         total1 = gamma_couple*gamma_couple*zeta_nu*zeta_nu;
         total = total + total1*intc(enes + 2, zetas, Rc, arreglo_factorial, arreglo_inv_factorial);
         total = total*NC_minus[index_i]*NC_minus[index_j];
         total1 =(double) (ang_j*(1 + ang_j) - ang_j*(1 + ang_j));
         totaltemp = total1*upper_incomplete_gamma(Rc, eles + 2, alphas);
         total1 = (double) 2.f*ang_j*alpha_nu;
         totaltemp = totaltemp + total1*upper_incomplete_gamma(Rc, eles + 1, alphas);
         total1 = alpha_nu*alpha_nu;
         totaltemp = totaltemp + total1*upper_incomplete_gamma(Rc, eles, alphas);
         totaltemp = NC_plus[index_i]*NC_plus[index_j]*totaltemp;
         matkint[k]=-0.5*total; 
         matkext[k]=-0.5*totaltemp;
         total = -0.5*(total + totaltemp); 
         matk[k] = total;
        }//label 1
         else {//label 2
             alpha_mu = zeta_mu  - (double) (n_mu + ang_i)/Rc;
             alpha_nu = zeta_nu  - (double) (n_nu + ang_j)/Rc;
             alphas = alpha_mu + alpha_nu;
             zetas = zeta_mu + zeta_nu;
             enes = n_mu + n_nu;
             eles = ang_i + ang_j;
             total1 = (double) (n_nu - 1)*n_nu - ang_j*(ang_j + 1);
             total = total1*intc(enes - 2, zetas, Rc, arreglo_factorial, arreglo_inv_factorial);
             total1 = (double) -2.f*n_nu*zeta_nu;
             total = total + total1*intc(enes - 1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial);
             total1 = (double) zeta_nu*zeta_nu;
             total = total + total1*intc(enes, zetas, Rc, arreglo_factorial, arreglo_inv_factorial);
             total = total*NC_minus[index_i]*NC_minus[index_j];
             total1 =(double) ang_j*(1 + ang_j) - ang_j*(1 + ang_j);
             totaltemp = total1*upper_incomplete_gamma(Rc, eles + 2, alphas);
             total1 = (double) 2.f*ang_j*alpha_nu;
             totaltemp = totaltemp + total1*upper_incomplete_gamma(Rc, eles + 1, alphas);
             total1 = alpha_nu*alpha_nu;
             totaltemp = totaltemp + total1*upper_incomplete_gamma(Rc, eles, alphas);
             totaltemp = totaltemp*NC_plus[index_i]*NC_plus[index_j];
             matkint[k] = -0.5*total;
             matkext[k] = -0.5*totaltemp;
             total = -0.5*(total + totaltemp);
             matk[k] = total;
             }//label 2
     }
   }
 } else {
   #pragma omp for
   for (k = 0; k < total_elements; k++) {
     indexes(nt, k, &index_i, &index_j);
     index_i = index_i - 1;
     index_j = index_j - 1;
     i = index_i;
     j = index_j;
     ang_i = mang[i];
     ang_j = mang[j];
     delta_kro_(&ang_i, &ang_j, &delta);
     ncm_i = ncm[i];
     ncm_j = ncm[j];
     delta_kro_(&ncm_i, &ncm_j, &delta1);
     delta = delta*delta1;
     if (delta == (double)0)
       matk[k] = (double)0;
     else {
       doble1 = expo[i] + expo[j];
       result1 = constc(i, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial);
       result2 = constc(j, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial);
       doble2 = (double)2.0*result1*result2/Rc;
       entero1 = np[i] + np[j] - 1;
       num1 = doble2*intc(entero1, doble1, Rc, arreglo_factorial, arreglo_inv_factorial);
       num1 = -num1*(double)np[j];
       entero1++;
       num2 = doble2*intc(entero1, doble1, Rc, arreglo_factorial, arreglo_inv_factorial);
       num2 = num2*(Rc*expo[j] + (double)np[j] )/Rc;
       entero1++;
       num3 = doble2*intc(entero1, doble1, Rc, arreglo_factorial, arreglo_inv_factorial);
       num3 = -num3*expo[j]/Rc;
       entero1 = np[j]*( np[j] - 1 ) - mang[j]*( mang[j] + 1 );
       num4 = (double)entero1*int1c(0, i, j, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial);
       entero1 = 2*np[j];
       num5 = -(double)entero1*expo[j]*int1c(1, i, j, Rc, expo, np, arreglo_factorial,
                                             arreglo_inv_factorial);
       num6 = expo[j]*expo[j]*int1c(2, i, j, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial);

       matk[k] = -(double)0.5*(num1 + num2 + num3 + num4 + num5 + num6);
     }
   }
 }
 }//Termina omp

 }

