#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

int wf_closed_shell(int z, char *using_gamma, int compara, char *bound, int nt, int elecalfa, int elecbeta,
                    double Rc, double *expo, int *np, double *zetas, int *mang, int *ncm, double *vectsfinalfa,
                    double *vectsfinbeta, char *tipo, double *NC_minus, double *NC_plus, double gamma_couple,
                    double *grid, double *grid_rho, double *grid_rho_beta, double *grid_der, double *grid_der_beta,
                    int n_points, 
                    double *arreglo_factorial, double *arreglo_inv_factorial,
                    double *matp, double *mats, int *iter, int save_i, int print_vectors, double *cusp_kato) {

 extern int grid_rhorad(int z, char *using_gamma, int compara, char *bound, int nt, int elecalfa, int elecbeta,
                        double Rc, double *expo, int *np, double *zetas, int *mang, double *vectsfinalfa,
                        double *vectsfinbeta, char *tipo, double *NC_minus, double *NC_plus, double gamma_couple,
                        double *grid, double *grid_rho, double *grid_rho_beta, int n_points, int n_boundary,
                        double *arreglo_factorial, double *arreglo_inv_factorial);

 extern int grid_rhorad_orbital(int z, char *using_gamma, int compara, char *bound, int nt, int elecalfa, int elecbeta,
                        double Rc, double *expo, int *np, double *zetas, int *mang, double *vectsfinalfa,
                        double *vectsfinbeta, char *tipo, double *NC_minus, double *NC_plus, double gamma_couple,
                        double *grid, double *grid_rho, double *grid_rho_beta, int n_points, int n_boundary,
                        double *arreglo_factorial, double *arreglo_inv_factorial, int selected_orb);

 extern int grid_derrad(int z, char *using_gamma, int compara, char *bound, int nt, int elecalfa, int elecbeta,
                        double Rc, double *expo, int *np, double *zetas, int *mang, double *vectsfinalfa,
                        double *vectsfinbeta, char *tipo, double *NC_minus, double *NC_plus, double gamma_couple,
                        double *grid, double *grid_rho, double *grid_rho_beta, int n_points, int n_boundary,
                        double *arreglo_factorial, double *arreglo_inv_factorial);

 extern void memoria_double_uni(int tamanio, double **arreglo, char *titulo);

 extern double numerical_int(double *grid, double *grid_fun, int points);

 extern void expected_value_r(int nt, int r_exp, int *np, int *mang, int *ncm, double *expo, 
                              double  Rc, char *bound, double *vectsfin, int elecalfa, double gamma_couple,
                              double *NC_minus, double *NC_plus, double *arreglo_factorial,
                              double *arreglo_inv_factorial, char *using_gamma, double *zeta, double *grid);

 grid_rhorad(z, using_gamma, compara, bound, nt, elecalfa, 0, Rc, expo, np, zetas, mang, vectsfinalfa, NULL,
             tipo, NC_minus, NC_plus, gamma_couple, grid, grid_rho, grid_rho_beta, n_points, save_i, 
             arreglo_factorial, arreglo_inv_factorial);

 grid_derrad(z, using_gamma, compara, bound, nt, elecalfa, 0, Rc, expo, np, zetas, mang, vectsfinalfa, NULL,
             tipo, NC_minus, NC_plus, gamma_couple, grid, grid_der, grid_der_beta, n_points, save_i,
             arreglo_factorial, arreglo_inv_factorial);

 int h, mu, nu, elemento1, elemento2, selected_orb;
 double sum, rho_0, drho_0, SHAN, Pi;
 double *array_i;

 array_i = NULL;
 memoria_double_uni(n_points*(sizeof(double)), &array_i, "Array_i");

 Pi = 4.f*atan(1.f);

 if (print_vectors == 1) {
   sum = 0.f;
   for (mu = 0; mu < nt; mu++) {
     for (nu = 0; nu < nt; nu++) {
       elemento1 = mu*nt + nu;
       elemento2 = nu*nt + mu;
       sum = sum + matp[elemento1]*mats[elemento2];
     }
   }
   printf("---------------------------\n");
   printf("Number of electrons: = %f\n", sum);
   printf("---------------------------\n");
   expected_value_r(nt, -1, np, mang, ncm, expo, Rc, bound, vectsfinalfa, elecalfa, gamma_couple,
                    NC_minus, NC_plus, arreglo_factorial, arreglo_inv_factorial, using_gamma,
                    zetas, grid);
   printf("---------------------------\n");
   expected_value_r(nt, 1,np, mang, ncm, expo, Rc, bound, vectsfinalfa, elecalfa, gamma_couple,
                    NC_minus, NC_plus, arreglo_factorial, arreglo_inv_factorial, using_gamma,
                    zetas, grid);
   printf("---------------------------\n");
   expected_value_r(nt, 2, np, mang, ncm, expo, Rc, bound, vectsfinalfa, elecalfa, gamma_couple,
                    NC_minus, NC_plus, arreglo_factorial, arreglo_inv_factorial, using_gamma,
                    zetas, grid);
   printf("---------------------------\n");
   for(h = 0; h < n_points; h++) array_i[h] = 0.f;
   for(h = 1; h < n_points - 1; h++) {
     if(grid_rho[h] > 1e-20)
        array_i[h] = -1.f*grid_rho[h]*log(grid_rho[h]);
   }

   SHAN = 4.f*Pi*numerical_int(grid, array_i, n_points);

   if(strcmp(bound,"free") == 0 && Rc == 0.f)
      printf("\nShannon entropy(%s,%s)= %5.4lf\n", tipo, bound, SHAN);
   else
      printf("\nShannon entropy(%s,%s,Rc: %3.3lf)= %5.4lf\n", tipo, bound, Rc, SHAN);

   printf("\nElectrons from numerical integration= %5.8lf\n\n", 4.f*Pi*numerical_int(grid, grid_rho, n_points));

   FILE *workout;
   char nameout[200];

   if (strcmp(bound, "free") == 0 && Rc == 0.f)
     sprintf(nameout, "%s_%s_rho_drho_+drho_rdf_divdrhorho", bound, tipo);
   else
     sprintf(nameout, "%3.3f_%s_%s_rho_drho_+drho_rdf_divdrhorho", Rc, bound, tipo);
   workout = fopen(nameout, "w");

   for (h = 0; h < n_points; h++)
      if(grid_rho[h] > 1e-16)
        fprintf(workout,"%5.4f  %24.14f  %24.14f  %24.14f  %24.14f %24.14f %24.14f\n", grid[h], grid_rho[h],
                grid_der[h], -1.f*grid_der[h], 4.f*Pi*grid[h]*grid[h]*grid_rho[h], grid_der[h]/grid_rho[h],
                -1.f*grid_der[h]/(2.f*z*grid_rho[h]));

    fclose(workout);

 } //imprime full

 printf("rho(0)                = %.4f\n", grid_rho[0]);
 printf("rho'(0)               = %.4f\n", grid_der[0]);
 printf("KATO CUSP = %4.4f\n", -grid_der[0]/(2.f*z*grid_rho[0]));
 *cusp_kato = -grid_der[0]/(2.f*z*grid_rho[0]);
 if (-grid_der[0]/(2.f*z*grid_rho[0])  < 0.98 || -grid_der[0]/(2.f*z*grid_rho[0]) > 1.20)
    *iter = 1e7;

 free(array_i);
 array_i = 0;

 return 0;
}
