// Overlap matrix
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

void traslape(char   *using_gamma,
              int     nt, 
	      double *mats, 
	      int    *np, 
	      int    *mang,
              int    *ncm, 
	      double *expo, 
	      double  Rc,
              double  gamma_couple,
	      char   *bound, 
	      double  U_0,
              double *NC_minus, 
	      double *NC_plus, 
	      double *arreglo_factorial,
              double *arreglo_inv_factorial)
{
 int i, j, k, index_i, index_j, total_elements;
 int ang_i, ang_j, ncm_i, ncm_j;
 double delta, delta1, part_real;
 
 double total;
 int enes, eles;
 double alphas, zetas;

 extern int indexes(int, int, int*, int*);
 extern double intl(int, int, int, double* , int*, double* );
 extern double int1c(int, int, int, double, double* , int*, double*, double*);
 extern double upper_incomplete_gamma(double, int, double);
//mrb extern void memoria_double_uni(int, double**, char*);
 extern void memoria_double_uni(int      tamanio,
                                double **arreglo,
                                char    *titulo);

 extern int delta_kro_(int*, int*, double*);
 extern double intc(int, double, double, double*, double*);

 total_elements = nt*nt;

#pragma omp parallel shared(total_elements, nt, mats, np, mang, ncm, expo, Rc, gamma_couple, bound, U_0, NC_minus, NC_plus, arreglo_factorial, arreglo_inv_factorial) 

#pragma omp parallel private(i, j, k, index_i, index_j, ang_i, ang_j, ncm_i, ncm_j, delta, delta1, total, enes, eles, alphas, zetas)

{
 int TID = omp_get_thread_num();
 if (strcmp(bound,"free") == 0) {
   #pragma omp for
   for (k = 0; k < total_elements; k++) {
     indexes(nt, k, &index_i, &index_j);
     index_i = index_i - 1;
     index_j = index_j - 1;
//     index_i = k/nt;
//     part_real = fmod(k,nt);
//     index_j = (int) part_real;
     ang_i = mang[index_i];
     ang_j = mang[index_j];
     delta_kro_(&ang_i, &ang_j, &delta);
     ncm_i = ncm[index_i];
     ncm_j = ncm[index_j];
     delta_kro_(&ncm_i, &ncm_j, &delta1);
     delta = delta*delta1;
     if (delta == (double)0) 
       mats[k] = (double)0;
     else 
       mats[k] = intl(index_i, index_j, 2, expo, np, arreglo_factorial);
   }
 }
 else if (strcmp(bound,"finite") == 0) {
   #pragma omp for
   for (k = 0; k <  total_elements; k++) { //label for omp
     indexes(nt, k, &index_i, &index_j);
     index_i = index_i - 1;
     index_j = index_j - 1;
    ang_i = mang[index_i];
    ang_j = mang[index_j];
     delta_kro_(&ang_i, &ang_j, &delta);
     ncm_i = ncm[index_i];
     ncm_j = ncm[index_j];
     delta_kro_(&ncm_i, &ncm_j, &delta1);
     delta = delta*delta1;
     if (delta == (double)0)
       mats[k] = (double)0;
     else {

        if (strcmp(using_gamma,"YES") == 0) {//label 1

           //  alphas = expo[index_i] + expo[index_j];
           zetas = expo[index_i] + expo[index_j];
           enes = np[index_i] + np[index_j];
           eles = ang_i + ang_j;	

           //zetas = alphas + (double) (enes + eles)/Rc  + 2.f*gamma_couple/(Rc*(gamma_couple - 1.f));
           alphas = zetas + (double) -1.f*(enes + eles)/Rc  + 2.f*gamma_couple/(Rc*(1.f - gamma_couple));

           total = Rc*Rc*intc(enes, zetas, Rc, arreglo_factorial, arreglo_inv_factorial);
          
           total = total -2.f*Rc*gamma_couple*intc(enes + 1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial); 
           
           total = total + gamma_couple*gamma_couple*intc(enes + 2, zetas, Rc, arreglo_factorial, arreglo_inv_factorial);
          
           total = total*NC_minus[index_i]*NC_minus[index_j];
          
           total = total + NC_plus[index_i]*NC_plus[index_j]*upper_incomplete_gamma(Rc, eles, alphas);
          
           mats[k] = total;
          }//label 1
           else {//label 2
              

               //alphas = expo[index_i] + expo[index_j];
               zetas = expo[index_i] + expo[index_j];
               enes = np[index_i] + np[index_j];
               eles = ang_i + ang_j;

               //zetas = alphas + (double) (enes + eles)/Rc;
               alphas = zetas - (double) (enes + eles)/Rc;

               total = NC_minus[index_i]*NC_minus[index_j]*intc(enes, zetas, Rc, arreglo_factorial, arreglo_inv_factorial);
               total = total + NC_plus[index_i]*NC_plus[index_j]*upper_incomplete_gamma(Rc, eles, alphas);


               mats[k] = total;
               
               }//label 2


     }
   }//label for omp
 } else {
   #pragma omp for
   for (k = 0; k <  total_elements; k++) {
     indexes(nt, k, &index_i, &index_j);
     index_i = index_i - 1;
     index_j = index_j - 1;
     ang_i = mang[index_i];
     ang_j = mang[index_j];
     delta_kro_(&ang_i, &ang_j, &delta);
     ncm_i = ncm[index_i];
     ncm_j = ncm[index_j];
     delta_kro_(&ncm_i, &ncm_j, &delta1);
     delta = delta*delta1;
     if (delta == (double)0)
       mats[k] = (double)0;
     else
       mats[k] = int1c(2, index_i, index_j, Rc, expo, np, arreglo_factorial, 
                       arreglo_inv_factorial);
   }
 }
 }//Termina omp
}
