//Kinetic energy matrix
//
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

//#ifdef _OPENMP
// #include <omp.h>
//#else
// #define omp_get_thread_num() 0
//#endif

void cinetica(char   *using_gamma,
              int     nt, 
	      double *matk, 
	      double *matkint, 
	      double *matkext, 
              char   *basis,
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
 double n_gto_mu, n_gto_nu, n_gto_imp_mu, n_gto_imp_nu, limit;    // mike
 double a0, a1, a2, a3, a4, a5, a6, g0, g1, g2, g3, g4, g5, g6;   // mike
 
 extern int indexes(int, int, int*, int*);
 extern double intl(int, int, int, double* , int*, double*);
 extern double int1c(int, int, int, double, double* , int*, double*, double*);
 extern double intc(int, double, double, double*, double*);
 extern double constc(int, double, double*, int*, double*, double*);
 extern void memoria_double_uni(int      tamanio,
                                double **arreglo,
                                char    *titulo);

 extern int delta_kro_(int*, int*, double*);
 extern double upper_incomplete_gamma(double, int, double);

 extern double lower_incomplete_gamma_function(int, double ); 
 extern double full_gamma_arg_2(int );
 extern double constant_normalization_GTO(int , int , double , double *);
 extern double constant_normalization_GTO_imp(int , int , double , double , double *);

 double pi;
 pi = ((double) 4)*atan(1.f);

// printf("Kinetic contribution to hamiltonian..\n");

 total_elements = nt*nt;

 if(strcmp(basis,"gtos") == 0) { /* Aquí empieza la maquinaria para gtos */
       for(k = 0; k < total_elements; k++) {
          indexes(nt, k, &index_i, &index_j);
          index_i = index_i - 1;
          index_j = index_j - 1;
          i = index_i;
          j = index_j;
          n_mu =np[i];
          n_nu =np[j];
          ang_i = mang[i];   /* número cuántico de momento ángular */
          ang_j = mang[j];
          delta_kro_(&ang_i, &ang_j, &delta);
          ncm_i = ncm[i];
          ncm_j = ncm[j];
          delta_kro_(&ncm_i, &ncm_j, &delta1);
          delta = delta*delta1;
          if(delta == (double)0){
            matk[k] = (double)0;
          }
          else {   // begins principal else
             num1 = ((double) n_mu + n_nu - 1.f);
             num2 = ((double) n_mu + n_nu + 1.f);
             num3 = ((double) n_mu + n_nu + 3.f);

             constant_normalization_GTO(i, n_mu, expo[i], &n_gto_mu);
             constant_normalization_GTO(j, n_nu, expo[j], &n_gto_nu);

             num4 = (-0.25f)*n_gto_mu*n_gto_nu*(
                                                (((double) n_nu*(n_nu - 1)) - ((double) ang_j*(ang_j + 1)))/pow(expo[i] + expo[j],num1/2.f) -
                                                (expo[j]*((double) 2.f*n_nu + 1.f)*num1)/pow(expo[i] + expo[j],num2/2.f) +
                                                (pow(expo[j],2.f)*num2*num1)/pow(expo[i] + expo[j],num3/2.f)
                                               );
             if(strcmp(bound,"free") == 0 || strcmp(bound,"FREE") == 0){
                matk[k] = num4*full_gamma_arg_2(n_mu + n_nu - 1);
             }   // ends free
             else{
                if(strcmp(bound,"confined") == 0){  /* Confinement by impenetrable walls */ 
                   a6 = 4.f*pow(expo[j],2.f)/pow(Rc,2.f);
                   a5 = 8.f*pow(expo[j],2.f)/Rc;
                   a4 = 4.f*pow(expo[j],2.f) - expo[j]*((double)4*n_nu + (double)6)/pow(Rc,2.f);                   
                   a3 = 8.f*expo[j]*((double)n_nu + 1)/Rc;
                   a2 = ((double)n_nu*n_nu + n_nu - ang_j*ang_j - ang_j)/pow(Rc,2.f) - 2.f*expo[j]*((double)2*n_nu + 1);
                   a1 = 2.f*((double) ang_j*ang_j + ang_j - n_nu*n_nu)/Rc;
                   a0 = ((double) n_nu*n_nu - n_nu - ang_j*ang_j - ang_j); 

                   limit = (expo[i] + expo[j])*pow(Rc,2.f);

                   g0 = lower_incomplete_gamma_function(n_mu + n_nu - 1, limit);
                   g1 = lower_incomplete_gamma_function(n_mu + n_nu, limit);
                   g2 = 0.5f*((double) n_mu + n_nu - 1)*g0 - pow(limit,num1/2.f)*exp(-limit);
                   g3 = 0.5f*((double) n_mu + n_nu)*g1 - pow(limit,((double) n_mu + n_nu)/2.f)*exp(-limit);
                   g4 = 0.25f*((double) n_mu + n_nu + 1)*((double) n_mu + n_nu - 1)*g0 - 
                        0.5f*pow(limit,num1/2.f)*exp(-limit)*(2.f*limit + ((double)n_mu + n_nu + 1));
                   g5 = 0.25f*((double)n_mu + n_nu)*((double)n_mu + n_nu + 2)*g1 - 
                        0.5f*pow(limit,((double) n_mu + n_nu)/2.f)*exp(-limit)*(2.f*limit + ((double)n_mu + n_nu + 2));


                   g6 = 0.125f*((double)n_mu + n_nu + 3)*((double)n_mu + n_nu + 1)*((double)n_mu + n_nu - 1)*g0 - 
                        0.25*pow(limit,num1/2.f)*exp(-limit)*(
                        4.f*pow(limit,2.f) + 2.f*((double)n_mu + n_nu + 3)*limit + ((double)n_mu + n_nu + 3)*((double)n_mu + n_nu + 1)
                        );

                   constant_normalization_GTO_imp(i, n_mu, expo[i], Rc, &n_gto_imp_mu);
                   constant_normalization_GTO_imp(j, n_nu, expo[j], Rc, &n_gto_imp_nu);
                   
                   num5 = -0.25f*n_gto_imp_mu*n_gto_imp_nu;
                   num5 = num5/pow(expo[i] + expo[j],num1/2.f);
                   
                   matkint[k] = a6*g6/pow(expo[i] + expo[j],3.0f) -
                                a5*g5/pow(expo[i] + expo[j],2.5f) +  
                                a4*g4/pow(expo[i] + expo[j],2.0f) +  
                                a3*g3/pow(expo[i] + expo[j],1.5f) +  
                                a2*g2/(expo[i] + expo[j]) +  
                                a1*g1/sqrt(expo[i] + expo[j]) +  
                                a0*g0;  
                   matkint[k] = num5*matkint[k];
                   matkext[k] = ((double) 0); 
                   matk[k] = matkint[k];     /* ya quedo */
//                   printf("matk_[%d] = %10.10f \n", k, matk[k]);
                }   // ends confined by impenetrable walls
                else{
                      matkint[k] = num4*lower_incomplete_gamma_function(n_mu + n_nu - 1, (expo[i] + expo[j])*pow(Rc,2.f));
                      matkext[k] = num4*(full_gamma_arg_2(n_mu + n_nu - 1) - lower_incomplete_gamma_function(n_mu + n_nu - 1, (expo[i] + expo[j])*pow(Rc,2.f)));
                      matk[k] = num4*full_gamma_arg_2(n_mu + n_nu - 1);
//                      printf("T_mu,nu[%d] = %f \n", k, matk[k]);
                }   // ends the rest of the cases dielectric, parabolic and finite
             }
          }   // ends principal else 
       }   // ends for
    }   // ends GTOs
 else { /* Aquí empieza toda la maquinaria para los stos */
//#pragma omp parallel shared(total_elements, nt, matk, np, mang, ncm, expo, Rc, gamma_couple, NC_minus, NC_plus)

//#pragma omp parallel private(i, j, k, index_i, index_j, ang_i, ang_j, ncm_i, ncm_j, delta, delta1, entero1, doble1, doble2, num1, num2, num3, num4, num5, num6, result1, result2, total, alpha_nu, alpha_mu, zetas, alphas, n_mu, n_nu, enes, eles, total1, zeta_nu, zeta_mu, totaltemp)

//{ // begins omp
// int TID = omp_get_thread_num();
 if (strcmp(bound,"free") == 0 || strcmp(bound,"debye") == 0 || strcmp(bound,"yukawa") == 0 || strcmp(bound,"baimbetov") == 0) {
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
       matk[k] = (double)0;
     else {
        num1 = expo[j]*expo[j]*intl(i, j, 2, expo, np, arreglo_factorial);
        num2 = ((double)2)*expo[j]*np[j]*intl(i, j, 1, expo, np, arreglo_factorial);
        num3 = (double)(np[j]*(np[j] - 1) - mang[j]*(mang[j] + 1));
        num3 = num3*intl(i, j, 0, expo, np, arreglo_factorial);
        matk[k] = -((double)0.5)*(num1 - num2 + num3);
//	printf("K_free[%d] = %f \n", k, matk[k]); // mike
     }
   }
 }
 else if (strcmp(bound,"finite") == 0) {
//   #pragma omp for
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
//	     printf("K_finite[%d] = %f \n", k, matk[k]); // mike
             }//label 2
     }
   }
 } else {
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
// }//Termina omp
 }

 }

