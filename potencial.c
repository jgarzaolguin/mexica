//External potential matrix
//
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

//#ifdef _OPENMP
// #include <omp.h>
//#else
// #define omp_get_thread_num() 0
//#endif
void expected_value_r(int     nt, 
		      int     r_exp, 
		      int    *np, 
		      int    *mang,
                      int    *ncm, 
		      double *expo, 
		      double  Rc, 
		      char   *bound,
                      char   *basis,  // new
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
 int n_mu, n_nu;                // mike
 double n_gto_mu, n_gto_nu, n_gto_imp_mu, n_gto_imp_nu, y;     // mike
 double g0, g1, g2, g3, g4, g5, g6;                            // mike
 double alphas, zetas;
 double delta, delta1, part_real, suma, mat_r_i_j;
 double *array_temp, *array_temp_2;
 double value_integral_ext,
        value_integral_int;
 double num1, pi;
 double partial;                        
 double r;

 extern int indexes(int, int, int*, int*);
 extern double intl(int, int, int, double*, int*, double*);
 extern double int1c(int, int, int, double, double*, int*, double*, double*);
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
 extern double full_gamma_arg_2(int );
 extern double constant_normalization_GTO(int , int , double , double *);
 extern double constant_normalization_GTO_imp(int , int , double , double , double *);  /* mike */
 extern double lower_incomplete_gamma_function(int, double ); 

 points = 600;
 array_temp = (double *) malloc (points * sizeof (double));
 total_elements = nt*nt;

 //#pragma omp parallel shared(nt, Rc, total_elements, mang, ncm, expo, np, z, matv) private(index_i, index_j, i, j, ang_i, ang_j, ncm_i, ncm_j, k, delta, delta1)
 //{
 // int TID = omp_get_thread_num();
 if(strcmp(basis,"gtos") == 0){     // begins gtos
    for(orbital = 0; orbital < elecalfa; orbital++){
       suma = ((double) 0);
       for(k = 0; k < total_elements; k++){
          indexes(nt, k, &index_i, &index_j);
          index_i = index_i - 1;
          index_j = index_j - 1;
          i = index_i;
          j = index_j;
          n_mu =np[i];
          n_nu =np[j];
          ang_i = mang[i];
          ang_j = mang[j];
          delta_kro_(&ang_i, &ang_j, &delta);
          ncm_i = ncm[i];
          ncm_j = ncm[j];
          delta_kro_(&ncm_i, &ncm_j, &delta1);
          delta = delta*delta1; 
          if(delta == ((double) 0)){
             mat_r_i_j = (double)0;
          }
          else {  /* begins principal else */
             if(strcmp(bound,"confined") == 0){
//                printf("¡Aun estoy trabajando en eso! \n");
                y = (expo[i] + expo[j])*pow(Rc,2.f);
                /* estas dos integrales son las fundamentales, el resto pueden ser puestas en términos de estas */
                g0 = lower_incomplete_gamma_function(n_mu + n_nu - 1, y);
                g1 = lower_incomplete_gamma_function(n_mu + n_nu, y);
                /* aprovecharemos estas definiciones y utilizaré una reducción de orden */
                g2 = 0.5f*((double) n_mu + n_nu - 1.f)*g0 - pow(y,0.5f*((double) n_mu + n_nu - 1.f))*exp(-y);
                g3 = 0.5f*((double) n_mu + n_nu)*g1 - pow(y,0.5f*((double) n_mu + n_nu))*exp(-y);
                g4 = 0.25f*((double) n_mu + n_nu + 1.f)*((double) n_mu + n_nu - 1.f)*g0 - 
                     0.5f*pow(y,0.5f*((double) n_mu + n_nu - 1.f))*exp(-y)*((double) n_mu + n_nu + 1.f + 2.f*y);
                g5 = 0.25f*((double) n_mu + n_nu)*((double) n_mu + n_nu + 2.f)*g1 - 
                     0.5f*pow(y,0.5f*((double) n_mu + n_nu))*exp(-y)*((double) n_mu + n_nu + 2.f + 2.f*y);
                g6 = 0.5f*((double) n_mu + n_nu + 3.f)*g4 - pow(y,0.5f*((double) n_mu + n_nu + 3.f))*exp(-y);

                constant_normalization_GTO_imp(i, n_mu, expo[i], Rc, &n_gto_imp_mu);
                constant_normalization_GTO_imp(j, n_nu, expo[j], Rc, &n_gto_imp_nu);
                if(r_exp == -1){
                   mat_r_i_j = (0.5f*n_gto_imp_mu*n_gto_imp_nu)/pow(expo[i] + expo[j],0.5f*((double) n_mu + n_nu));
                   mat_r_i_j = mat_r_i_j*( g3/((expo[i] + expo[j])*pow(Rc,2.f)) - (2.f*g2)/(sqrt(expo[i] + expo[j])*Rc) + g1 ); 
                   suma = suma + vectsfin[i*nt + orbital]*vectsfin[j*nt + orbital]*mat_r_i_j;
                }
                else{
                   if(r_exp == 1){
                      mat_r_i_j = (0.5f*n_gto_imp_mu*n_gto_imp_nu)/pow(expo[i] + expo[j],0.5f*((double) n_mu + n_nu + 2.f));
                      mat_r_i_j = mat_r_i_j*( g5/((expo[i] + expo[j])*pow(Rc,2.f)) - (2.f*g4)/(sqrt(expo[i] + expo[j])*Rc) + g3 ); 
                      suma = suma + vectsfin[i*nt + orbital]*vectsfin[j*nt + orbital]*mat_r_i_j;
                   }
                   else{
                      if(r_exp == 2){
                         mat_r_i_j = (0.5f*n_gto_imp_mu*n_gto_imp_nu)/pow(expo[i] + expo[j],0.5f*((double) n_mu + n_nu + 3.f));
                         mat_r_i_j = mat_r_i_j*( g6/((expo[i] + expo[j])*pow(Rc,2.f)) - (2.f*g5)/(sqrt(expo[i] + expo[j])*Rc) + g4 ); 
                         suma = suma + vectsfin[i*nt + orbital]*vectsfin[j*nt + orbital]*mat_r_i_j;
                      }
                   }
                } 
             }
             else{     /* begins the rest of the cases: finite, dielec, parabolic and free */
                constant_normalization_GTO(i, n_mu, expo[i], &n_gto_mu);
                constant_normalization_GTO(j, n_nu, expo[j], &n_gto_nu);
                if(r_exp == -1){
                   num1 = ((double) n_mu + n_nu)/((double) 2);
                   mat_r_i_j = (n_gto_mu*n_gto_nu*full_gamma_arg_2(n_mu + n_nu))/(((double) 2)*pow(expo[i] + expo[j], num1));
                   suma = suma + vectsfin[i*nt + orbital]*vectsfin[j*nt + orbital]*mat_r_i_j;
                }
                else{
                   if(r_exp == 1){
                      num1 = ((double) n_mu + n_nu + 2)/((double) 2);
                      mat_r_i_j = (n_gto_mu*n_gto_nu*full_gamma_arg_2(n_mu + n_nu + 2))/(((double) 2)*pow(expo[i] + expo[j], num1));
                      suma = suma + vectsfin[i*nt + orbital]*vectsfin[j*nt + orbital]*mat_r_i_j;
                   }
                   else{
                      if(r_exp == 2){
                         num1 = ((double) n_mu + n_nu + 3)/((double) 2);
                         mat_r_i_j = (n_gto_mu*n_gto_nu*full_gamma_arg_2(n_mu + n_nu + 3))/(((double) 2)*pow(expo[i] + expo[j], num1));
                         suma = suma + vectsfin[i*nt + orbital]*vectsfin[j*nt + orbital]*mat_r_i_j;
                      }                   
                   }
                }
             }         /* ends the rest of the cases: finite, dielec, parabolic and free */
          }       /* ends principal else */
       }
       printf("Orbital %d, <r^%d>=%10.5f\n", orbital, r_exp, suma);
    }
 }  // ends gtos
 else{     
//    if(strcmp(basis,"stos") == 0){  // begins stos
       if(strcmp(bound,"free") == 0 || strcmp(bound,"debye") == 0 || strcmp(bound,"yukawa") == 0) {
///   #pragma omp for
          for(orbital = 0; orbital < elecalfa; orbital++) {
             suma = 0.f;
             for(k = 0; k < total_elements; k++) {
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
                if(delta == (double)0) 
                   mat_r_i_j = (double)0;
                else 
                   mat_r_i_j = intl(i, j, 2 + r_exp, expo, np, arreglo_factorial);
//		suma = suma + matp[element2]*mat_r_i_j;
                suma = suma + vectsfin[i*nt + orbital]*vectsfin[j*nt + orbital]*mat_r_i_j;
             }
             printf("Orbital %d, <r^%d>=%10.5f\n", orbital, r_exp, suma); 
          }
       }
       else
          if(strcmp(bound,"confined") == 0) {
//   #pragma omp for
             for(orbital = 0; orbital < elecalfa; orbital++) {
                suma = 0.f;
                for(k = 0; k < total_elements; k++) {
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
                   if(delta == (double)0)
                      mat_r_i_j = (double)0;
                   else 
                      mat_r_i_j = int1c(i, j, 2 + r_exp, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial);
                   suma = suma + vectsfin[i*nt + orbital]*vectsfin[j*nt + orbital]*mat_r_i_j;
                }
                printf("Orbital %d, <r^%d>=%10.5f\n", orbital, r_exp, suma); 
             }
          }
          else {
//   #pragma omp for
             double pi;
             pi = 4.f*atan(1.f);
             for(i = 0; i < points; i++) array_temp[i] = 0.f;
                for(orbital = 0; orbital < elecalfa; orbital++) {
                   suma = 0.f;
                   for(i = 0; i < points; i++) {
                      r = grid[i];
                      if(r < Rc )
                         partial = finite_orbital_int(using_gamma, nt, orbital, r, Rc, expo, np,
                                                      vectsfin, gamma_couple, NC_minus);
                      else
                         partial = finite_orbital_ext(nt, orbital, r, Rc, zeta, mang, vectsfin,
                                                   gamma_couple, NC_plus);
                      if(i == 0) 
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
//    }   // ends stos
 }  
//} Cierra openmp
} //aquí cierra expected value

 double complex v_ne_yukawa(int n_mu, int n_nu, double expo_mu, double expo_nu, double alpha){
 // This is for the debye potential with the factor cos(\beta * r)
         extern double factorial(int );
         double beta, e1, n1, c1;
         double complex b1, b2, partial;

         beta = (double) 1/alpha;
         b1 = beta - beta*I;
         b2 = beta + beta*I;
         n1 = (double) n_mu + n_nu;
         e1 = expo_mu + expo_nu;
         c1 = (double) 1;
         partial = factorial(n_mu + n_nu - 1)*(c1/cpow(e1 + b1,n1) + c1/cpow(e1 + b2,n1));
//       printf("Re = %20.15f, Im = %20.15f \n", creal(partial), cimag(partial));
//       return(creal(partial));
         return(partial);
 }

void potencial(char   *using_gamma,
               int     nt, 
	       double  z, 
	       double *matv, 
	       int    *np, 
	       int    *mang,
               int    *ncm, 
	       double *expo,
               char   *basis, 
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
 int n_mu, n_nu;                                                   // mike
 double n_gto_mu, n_gto_nu, n_gto_imp_mu, n_gto_imp_nu, limit;     // mike
 double num2, num3, num4, num5, num6, y, n_sto_mu, n_sto_nu;       // mike
 double num0, num1, pi;                                            // mike
 double alphas, zetas;
 double result1, result2;
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
// mike--------------------------------------------------------------------
 extern double full_gamma_arg_2(int );
 extern double lower_incomplete_gamma_function(int, double );
 extern double constant_normalization_GTO(int , int , double , double *);
 extern double constant_normalization_GTO_imp(int , int , double , double , double *);
 extern long long int factorial_mike(int );
 extern double constant_normalization_sto(int , int , double , double *);
 extern double complex v_ne_yukawa(int , int , double , double , double );

 total_elements = nt*nt;
 pi = 4.f*atan(1.f);

 if (strcmp(basis,"gtos") == 0) { //aquí empiezan los cálculos para los gtos
     for (k = 0; k < total_elements; k++) {
          indexes(nt, k, &index_i, &index_j);
          index_i = index_i - 1;
          index_j = index_j - 1;
          i = index_i;
          j = index_j;
          n_mu =np[i];
          n_nu =np[j];
          ang_i = mang[i];
          ang_j = mang[j];
          delta_kro_(&ang_i, &ang_j, &delta);
          ncm_i = ncm[i];
          ncm_j = ncm[j];
          delta_kro_(&ncm_i, &ncm_j, &delta1);
          delta = delta*delta1;
          if (delta == (double)0)
            matv[k] = (double)0;
          else { //begins principal else
             num1 = ((double) n_mu + n_nu);
             constant_normalization_GTO(i, n_mu, expo[i], &n_gto_mu);
             constant_normalization_GTO(j, n_nu, expo[j], &n_gto_nu);
             if(strcmp(bound,"free") == 0){
                matv[k] = -((double) z)*n_gto_mu*n_gto_nu/(((double) 2)*pow(expo[i] + expo[j], num1/2.f));

                matv[k] = matv[k]*((double) full_gamma_arg_2(n_mu + n_nu));
             } 
             else{
                if(strcmp(bound,"dielectricc") == 0) {
                   matv[k] = -((double) z)*n_gto_mu*n_gto_nu/(((double) 2)*pow(expo[i] + expo[j], num1/2.f));

                   matv[k] = matv[k]*((U_0 - ((double) 1))*lower_incomplete_gamma_function(n_mu + n_nu, (expo[i] + expo[j])*pow(Rc,2.f))/U_0 +
                                      ((double) full_gamma_arg_2(n_mu + n_nu))/U_0 );
                } 
                else{
                   if(strcmp(bound,"finite") == 0){
                     matv[k] = -n_gto_mu*n_gto_nu/(((double) 2)*pow(expo[i] + expo[j], num1/2.f));

                     matv[k] = matv[k]*(
                               ((double) z)*lower_incomplete_gamma_function(n_mu + n_nu, (expo[i] + expo[j])*pow(Rc,2.f)) - 
                               ( U_0/pow(expo[i] + expo[j],0.5f) )*
                               ( ((double) full_gamma_arg_2(n_mu + n_nu + 1)) -
                                 lower_incomplete_gamma_function(n_mu + n_nu + 1, (expo[i] + expo[j])*pow(Rc,2.f)) )
                               );
//                    printf("V_mu,nu[%d] = %f \n", k, matv[k]);
                   }
                   else{
                      if(strcmp(bound,"parabolic") == 0){
                         matv[k] = -n_gto_mu*n_gto_nu/(((double) 2)*pow(expo[i] + expo[j], num1/2.f));

                         matv[k] = matv[k]*(
                                   ((double) z)*lower_incomplete_gamma_function(n_mu + n_nu, (expo[i] + expo[j])*pow(Rc,2.f)) - 
                                     pow(U_0,2.f)/(((double) 2)*pow(expo[i] + expo[j],1.5f))*
                                     ( ((double) full_gamma_arg_2(n_mu + n_nu + 3)) - 
                                       lower_incomplete_gamma_function(n_mu + n_nu + 3, (expo[i] + expo[j])*pow(Rc,2.f)) )
                                   );
                      }
                      else{
                         if(strcmp(bound,"confined") == 0){
                            limit = (expo[i] + expo[j])*pow(Rc,2.f); 
                            constant_normalization_GTO_imp(i, n_mu, expo[i], Rc, &n_gto_imp_mu);
                            constant_normalization_GTO_imp(j, n_nu, expo[j], Rc, &n_gto_imp_nu);
 
                            matv[k] = -(((double)z)*n_gto_imp_mu*n_gto_imp_nu)/( ((double)2)*pow(expo[i] + expo[j],num1/2.f));

                            matv[k] = matv[k]*(
                                      lower_incomplete_gamma_function(n_mu + n_nu, limit) -
                                      (((double)2)*lower_incomplete_gamma_function(n_mu + n_nu + 1, limit))/(Rc*sqrt(expo[i] + expo[j])) +
                                      lower_incomplete_gamma_function(n_mu + n_nu + 2, limit)/(pow(Rc,2.f)*(expo[i] + expo[j]))
                            );
 
//                            printf("V_mu,nu[%d] = %10.15f \n", k, matv[k]);  /* Ya quedo mike */
                         }
                      }
                   }
                }
             }
          } // ends principal else
     }
          /* Atención, al parecer hay una especie de reasignación pues epsilon ya está definida en data.c sin embargo, para evitar poner mas argumentos en la función de potencial lo que se hizo fue reemplazar U_0 por la constante dieléctrica */
 } 
 else {//aquí empieza el else para los stos

//#pragma omp parallel shared(total_elements, nt, z, matv, np, mang, ncm, expo, Rc, gamma_couple, bound, U_0, NC_minus, NC_plus, arreglo_factorial, arreglo_inv_factorial, iter_pol, charge_int) private(i, j, k, index_i, index_j, delta, delta1, total, ang_i, ang_j, ncm_i, ncm_j, enes, eles, alphas, zetas, adicional)

//{ // begins omp
// int TID = omp_get_thread_num();
 if (strcmp(bound,"free") == 0) {
//   #pragma omp for
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
//	  printf("V_free[%d] = %f \n", k, matv[k]); // mike
       }
   }
 }
 else
      if(strcmp(bound,"debye") == 0){ // begins Debye
//         printf("lambda = %f \n", U_0);
//         printf("1/lambda = %f \n", 1.f/U_0);
         for(k = 0; k < total_elements; k++) {
            indexes(nt, k, &index_i, &index_j);
            index_i = index_i - 1;
            index_j = index_j - 1;
            i = index_i;
            j = index_j;
            n_mu =np[i];
            n_nu =np[j];
            ang_i = mang[i];
            ang_j = mang[j];
            delta_kro_(&ang_i, &ang_j, &delta);
            ncm_i = ncm[i];
            ncm_j = ncm[j];
            delta_kro_(&ncm_i, &ncm_j, &delta1);
            delta = delta*delta1;
            if(delta == (double)0)
               matv[k] = (double)0;
            else{
               constant_normalization_sto(i, n_mu, expo[i], &n_sto_mu);
               constant_normalization_sto(j, n_nu, expo[j], &n_sto_nu);

               num0 = ((double) 1)/U_0;  // 1/lambda_d
               num1 = n_mu + n_nu;
               num2 = expo[i] + expo[j] + num0;
               num3 = pow(num2, (double) num1);

               matv[k] = -z*n_sto_mu*n_sto_nu*factorial_mike(n_mu + n_nu - 1)/num3;

//               printf("V_debye[%d,%d] = %f \n", i,j, matv[k]); // mike
            }
         }
      } // ends Debye
      else
              if(strcmp(bound,"yukawa") == 0){
//                    printf("lambda = %f \n", U_0);
                      for(k = 0; k < total_elements; k++) {
                              indexes(nt, k, &index_i, &index_j);
                              index_i = index_i - 1;
                              index_j = index_j - 1;
                              i = index_i;
                              j = index_j;
                              n_mu =np[i];
                              n_nu =np[j];
                              ang_i = mang[i];
                              ang_j = mang[j];
                              delta_kro_(&ang_i, &ang_j, &delta);
                              ncm_i = ncm[i];
                              ncm_j = ncm[j];
                              delta_kro_(&ncm_i, &ncm_j, &delta1);
                              delta = delta*delta1;
                              if(delta == (double) 0)
                                      matv[k] = (double)0;
                              else{
                                      constant_normalization_sto(i, n_mu, expo[i], &n_sto_mu);
                                      constant_normalization_sto(j, n_nu, expo[j], &n_sto_nu);
                                      matv[k] = -z*n_sto_mu*n_sto_nu*v_ne_yukawa(n_mu, n_nu, expo[i], expo[j], U_0)/((double) 2);
//                                    printf("V_yukawa[%d,%d] = %f \n", i,j, matv[k]); // mike
                              }
                      }
              }  // ends yukawa
              else	      
 if (strcmp(bound,"dielectricnc") == 0) {
//   #pragma omp for
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
//   #pragma omp for
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
//	printf("V_finite[%d] = %f \n", k, matv[k]); // mike


            }//label 2


     }
   }
 } 
 else

 if (strcmp(bound,"dielectricc") == 0 && strcmp(using_gamma,"NO") == 0) {
//   #pragma omp for
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
//   #pragma omp for
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
//   #pragma omp for
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

// } //termina omp
}// aquí termina el else para las stos
}

