//External potential matrix
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
void expected_value_r(int     nt, 
		      int     r_exp, 
		      int    *np, 
		      int    *mang,
                      int    *ncm, 
		      double *expo, 
		      double  Rc, 
		      char   *bound,
                      double *vectsfin, 
		      int     elecalfa, 
                      double  gamma_couple,
                      double *NC_minus,
                      double *NC_plus,
		      double *arreglo_factorial,
                      double *arreglo_inv_factorial,
                      char   *using_gamma,
                      double *zeta,
                      double *grid)
{
 int i, j, k;

 int index_i, index_j, element_k, total_elements, proc_come;
 int bloque, restante, do_ini, do_fin, total_size, sizebi;
 int ang_i, ang_j, ncm_i, ncm_j, orbital;
 int enes, eles;
 int compara, points;
 double alphas, zetas;
 double delta, delta1, part_real, suma, mat_r_i_j;
 double *array_temp, *array_temp_2;
 double value_integral_ext,
        value_integral_int;
 double partial;                       //new 
 double r;

 extern int indexes(int, int, int*, int*);
 extern double intl(int, int, int, double*, int*, double*);
 extern double int1c(int, int, int, double, double*, int*, double*, double*);
//mrb extern void memoria_double_uni(int, double**, char*);
 extern void memoria_double_uni(int      tamanio,
                                double **arreglo,
                                char    *titulo);
 extern double delta_kro_(int*, int*, double*);

 extern double intc(int a, double b, double r, double*, double*);
 extern double upper_incomplete_gamma(double, int, double);

 extern double finite_orbital_int(char   *using_gamma,
                                   int     nt,
                                   int     orbital,
                                   double  r,
                                   double  Rc,
                                   double *zeta,
                                   int    *np,
                                   double *vectors,
                                   double  gamma_couple,
                                   double *N_minus);
 extern double finite_orbital_ext(int     nt,
                                   int     orbital,
                                   double  r,
                                   double  Rc,
                                   double *alfa,
                                   int    *mang,
                                   double *vectors,
                                   double  gamma_couple,
                                   double *N_plus);

 extern double numerical_int(double   *grid,
                             double   *grid_fun,
                             int       points);

 points = 600;
 array_temp = (double *) malloc (points * sizeof (double));

 total_elements = nt*nt;

//#pragma omp parallel shared(nt, Rc, total_elements, mang, ncm, expo, np, z, matv) private(index_i, index_j, i, j, ang_i, ang_j, ncm_i, ncm_j, k, delta, delta1)
//{
// int TID = omp_get_thread_num();
 if (strcmp(bound,"free") == 0) {
//   #pragma omp for
   for (orbital = 0; orbital < elecalfa; orbital++) {
   suma = 0.f;
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
       mat_r_i_j = (double)0;
     else 
       mat_r_i_j = intl(i, j, 2 + r_exp, expo, np, arreglo_factorial);
//     suma = suma + matp[element2]*mat_r_i_j;
     suma = suma + vectsfin[i*nt + orbital]*vectsfin[j*nt + orbital]*mat_r_i_j;
   }
   printf("Orbital %d, <r^%d>=%10.5f\n", orbital, r_exp, suma); 
   }
 }
 else
 if(strcmp(bound,"confined") == 0) {
//   #pragma omp for
   for (orbital = 0; orbital < elecalfa; orbital++) {
   suma = 0.f;
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
       mat_r_i_j = (double)0;
     else 
     mat_r_i_j = int1c(i, j, 2 + r_exp, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial);
     suma = suma + vectsfin[i*nt + orbital]*vectsfin[j*nt + orbital]*mat_r_i_j;
   }
   printf("Orbital %d, <r^%d>=%10.5f\n", orbital, r_exp, suma); 
   }
 }
 else
 {
//   #pragma omp for
  double pi;
  pi = 4.f*atan(1.f);
   for (i = 0; i < points; i++) array_temp[i] = 0.f;
   for (orbital = 0; orbital < elecalfa; orbital++) {
      suma = 0.f;
      for (i = 0; i < points; i++) {
        r = grid[i];
        if (r < Rc )
          partial = finite_orbital_int(using_gamma, nt, orbital, r, Rc, expo, np,
                                       vectsfin, gamma_couple, NC_minus);
        else
          partial = finite_orbital_ext(nt, orbital, r, Rc, zeta, mang, vectsfin,
                                       gamma_couple, NC_plus);

        if (i == 0) 
          array_temp[i] = 0.f;
        else
          array_temp[i] = pow(r,(double) r_exp)*partial*partial; //This array will be multiplied by r^2 in
                                                                 //numerical_int function
   }
   suma = numerical_int(grid, array_temp, points);
   printf("Orbital %d, <r^%d>=%10.5f\n", orbital, r_exp, suma);
   }

 }
 free(array_temp);
//} Cierra openmp
}


void potencial(char   *using_gamma,
               int     nt, 
	       double  z, 
	       double *matv, 
	       int    *np, 
	       int    *mang,
               int    *ncm, 
	       double *expo, 
	       double  Rc, 
               double  gamma_couple,
               char   *bound, 
               int     iter_pol,
               double  charge_int,
	       double  U_0,
               double *NC_minus, 
	       double *NC_plus, 
	       double *arreglo_factorial,
               double *arreglo_inv_factorial,
               double elec,
               int plasma)
{
 int i, j, k;

 int index_i, index_j, total_elements;
 double delta, delta1;

 double total, adicional;
 int ang_i, ang_j, ncm_i, ncm_j;
 int enes, eles;	
 double alphas, zetas;


// double *array_temp, *array_temp_2;

 extern int indexes(int, int, int*, int*);
 extern double intl(int, int, int, double*, int*, double*);
 extern double int1c(int, int, int, double, double* , int*, double*, double*);
//mrb extern void memoria_double_uni(int, double**, char*);
 extern void memoria_double_uni(int      tamanio,
                                double **arreglo,
                                char    *titulo);
 extern double delta_kro_(int*, int*, double*);
 extern double intc(int a, double b, double r, double*, double*);
 extern double upper_incomplete_gamma(double, int, double);
 extern int constant_normalization_free(int*, double*, int*, double*);


 total_elements = nt*nt;

#pragma omp parallel shared(total_elements, nt, z, matv, np, mang, ncm, expo, Rc, gamma_couple, bound, U_0, NC_minus, NC_plus, arreglo_factorial, arreglo_inv_factorial, iter_pol, charge_int) private(i, j, k, index_i, index_j, delta, delta1, total, ang_i, ang_j, ncm_i, ncm_j, enes, eles, alphas, zetas, adicional)



{
 int TID = omp_get_thread_num();
 if (strcmp(bound,"free") == 0) {
 double result1, result2;
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
       matv[k] = (double)0;
     else
	{ 
          matv[k] = -z*intl(i, j, 1, expo, np, arreglo_factorial);  
       }
   }
 }
 else 

 if (strcmp(bound,"dielectricnc") == 0) {
 double result1, result2;
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
       matv[k] = (double)0;
     else
        {
          constant_normalization_free(&i, expo, np, &result1);
          constant_normalization_free(&j, expo, np, &result2);
          matv[k] = -z*(result1*result2*intc(np[i] + np[j] - 1, expo[i] + expo[j], Rc, arreglo_factorial, arreglo_inv_factorial) + (intl(i, j, 1, expo, np, arreglo_factorial) - result1*result2*intc(np[i] + np[j] - 1, expo[i] + expo[j], Rc, arreglo_factorial, arreglo_inv_factorial))/U_0);

       }
   }
 }
 else 

if (strcmp(bound,"finite") == 0) {
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
       matv[k] = (double)0;
     else {

     if (strcmp(using_gamma,"YES") == 0){//label 1
        //alphas = expo[index_i] + expo[index_j];
        zetas = expo[index_i] + expo[index_j];
        enes = np[index_i] + np[index_j];
        eles = ang_i + ang_j;

        //zetas = alphas + (double) (enes + eles)/Rc + 2.f*gamma_couple/(Rc*(gamma_couple - 1.f));
        alphas = zetas + (double) -1.f*(enes + eles)/Rc + 2.f*gamma_couple/(Rc*(1.f - gamma_couple));

        total = Rc*Rc*intc(enes - 1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial);
        total = total - 2.f*Rc*gamma_couple*intc(enes, zetas, Rc, arreglo_factorial, arreglo_inv_factorial);
        total = total +  gamma_couple*gamma_couple*intc(enes + 1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial);
        total = -z*total*NC_minus[index_i]*NC_minus[index_j];
        total = total + 0.5*U_0*NC_plus[index_i]*NC_plus[index_j]*upper_incomplete_gamma(Rc, eles, alphas);
 
        matv[k] = total;
       }//label 1
        else {//label 2

        //alphas = expo[index_i] + expo[index_j];
        //zetas = alphas + (double) (enes + eles)/Rc;
        
        zetas = expo[index_i] + expo[index_j];
        enes = np[index_i] + np[index_j];
        eles = ang_i + ang_j;

        alphas  = zetas - (double) (enes + eles)/Rc;


        total = -z*NC_minus[index_i]*NC_minus[index_j]*intc(enes - 1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial);
        total = total + 0.5*U_0*NC_plus[index_i]*NC_plus[index_j]*upper_incomplete_gamma(Rc, eles, alphas);

        matv[k] = total;


            }//label 2


     }
   }
 } 
 else

 if (strcmp(bound,"dielectricc") == 0 && strcmp(using_gamma,"NO") == 0) {
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
       matv[k] = (double)0;
     else {

        zetas = expo[index_i] + expo[index_j];
        enes = np[index_i] + np[index_j];
        eles = ang_i + ang_j;

        alphas  = zetas - (double) (enes + eles)/Rc;

        total = -z*NC_minus[index_i]*NC_minus[index_j]*intc(enes - 1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial);
        total = total - z*(NC_plus[index_i]*NC_plus[index_j]*upper_incomplete_gamma(Rc, eles + 1, alphas))/U_0;



        matv[k] = total;


     }
   }
 } 
 else
  if (strcmp(bound,"polarization") == 0 && strcmp(using_gamma,"NO") == 0) {
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
       matv[k] = (double)0;
     else {

        zetas = expo[index_i] + expo[index_j];
        enes = np[index_i] + np[index_j];
        eles = ang_i + ang_j;

        alphas  = zetas - (double) (enes + eles)/Rc;

        total = -z*NC_minus[index_i]*NC_minus[index_j]*intc(enes - 1, zetas, Rc, arreglo_factorial, arreglo_inv_factorial);
        //mrb total = total - z*(NC_plus[index_i]*NC_plus[index_j]*upper_incomplete_gamma(Rc, eles + 1, alphas))/U_0 + (1.f - 1.f/U_0)*charge_int*NC_plus[index_i]*NC_plus[index_j]*upper_incomplete_gamma(Rc, eles + 1, alphas);
        if (iter_pol > 1)
        total = total + NC_plus[index_i]*NC_plus[index_j]*upper_incomplete_gamma(Rc, eles + 1, alphas)*(-z/U_0 + (1.f - 1.f/U_0)*charge_int);

        //mrb total = total + NC_plus[index_i]*NC_plus[index_j]*upper_incomplete_gamma(Rc, eles + 1, alphas)*(-z/U_0 + (1.f - 1.f/U_0)*charge_int);
        matv[k] = total;


     }
   }
  }
 else {
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
       matv[k] = (double)0;
     else {
      matv[k] = -z*int1c(1, i, j, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial);
      adicional = 0e00;
      if (plasma == 1) {
        adicional = 3e00*int1c(2, index_i, index_j, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial);
        adicional = adicional - int1c(4, index_i, index_j, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial)/(Rc*Rc);
        adicional = (z - elec)*adicional/(2e00*Rc);
      }
      matv[k] = matv[k] + adicional;
     }
   }
 }

 } //termina omp
}

