#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

int grid_rhorad(int           z,
                char         *using_gamma,
                int           compara,
                char         *bound,
                char         *basis, 
                int           nt,
                int           elecalfa,
                int           elecbeta,
                double        Rc,
                double       *expo,
                int          *np,
                double       *zetas,
                int          *mang,
                double       *vectsfinalfa,
                double       *vectsfinbeta,
                char         *tipo,
                double       *NC_minus,
                double       *NC_plus,
                double        gamma_couple,
                double       *grid,
                double       *grid_rho,
                double       *grid_rho_beta,
                int           n_points,
                int           n_boundary,
                double       *arreglo_factorial,
                double       *arreglo_inv_factorial)
 
 {//function

 extern double rho_radial_confined(int     nt,
                                   int     elec,
                                   double  r,
                                   double  Rc,
                                   double *expo,
                                   int    *np,
                                   double *vectors,
                                   char   *tipo,
                                   double *arreglo_factorial,
                                   double *arreglo_inv_factorial);

 extern double rho_radial_finite_int(char   *using_gamma,
                                     int     nt,
                                     int     elec,
                                     double  r,
                                     double  Rc,
                                     double *zeta,
                                     int    *np,
                                     double *vectors,
                                     char   *tipo,
                                     double  gamma_couple,
                                     double *N_minus);

 extern double rho_radial_finite_ext(int     nt,
                                     int     elec,
                                     double  r,
                                     double  Rc,
                                     double *alfa,
                                     int    *mang,
                                     double *vectors,
                                     char   *tipo,
                                     double  gamma_couple,
                                     double *N_plus);

 extern double rho_radial_free(int     nt,
                               int     elec,
                               double  r,
                               double *expo,
                               int    *np,
                               double *vectors,
                               char   *tipo);

 extern double rho_radial_gto(int nt, int elec, double r, double* expo, int* np, double* vectors, char* tipo);

 extern double rho_radial_gto_imp(int nt, int elec, double r, double r0, double* expo, int* np, double* vectors, char* tipo);

 int unsigned i;
 double p1;

 if(strcmp(bound,"free") == 0){
    if(strcmp(basis,"STOs") == 0){
       for(i = 0; i < n_points; i++) { 
         p1 = grid[i];
         if(compara == 0) {
            grid_rho[i] = rho_radial_free(nt, elecalfa, p1, expo, 
                                          np, vectsfinalfa, tipo);
         } 
         else {
            grid_rho[i] = rho_radial_free(nt, elecalfa, p1, expo,
                                          np, vectsfinalfa, tipo);
            
            grid_rho_beta[i] = rho_radial_free(nt, elecbeta, p1, expo,
                                               np, vectsfinbeta, tipo);
         }
       }
    }  /* ends stos */
    else {
       if(strcmp(basis,"GTOs") == 0){
          for(i = 0; i < n_points; i++) { 
             p1 = grid[i];
             if(compara == 0) {
                grid_rho[i] = rho_radial_gto(nt, elecalfa, p1, expo, np, vectsfinalfa, tipo);
             } 
             else {
                grid_rho[i] = rho_radial_gto(nt, elecalfa, p1, expo, np, vectsfinalfa, tipo);
                grid_rho_beta[i] = rho_radial_gto(nt, elecbeta, p1, expo, np, vectsfinbeta, tipo);
             }
//          printf("rho[%d] = %f \n", i, grid_rho[i]);
          }
       }  /* ends gtos */
    }
 } /* ends Free */
 else {
    if(strcmp(bound,"finite") == 0){
       if(strcmp(basis,"STOs") == 0){
          for(i = 0; i < n_boundary; i++) { //for save_temp
             p1 = grid[i];
             if(compara == 0) {
                grid_rho[i] = rho_radial_finite_int(using_gamma, nt, elecalfa, p1, Rc, expo, np,
                                                    vectsfinalfa, tipo, gamma_couple, NC_minus);
                } 
                else {
                   grid_rho[i] = rho_radial_finite_int(using_gamma, nt, elecalfa, p1, Rc, expo, np,
                                                       vectsfinalfa, tipo, gamma_couple, NC_minus);  
                   grid_rho_beta[i] = rho_radial_finite_int(using_gamma, nt, elecbeta, p1, Rc, expo, np,
                                                            vectsfinbeta, tipo, gamma_couple, NC_minus);
                }
          } //for save_temp 
          for(i = n_boundary; i < n_points; i++) { //for save_temp
             p1 = grid[i];
             if(compara == 0) {
                grid_rho[i] = rho_radial_finite_ext(nt, elecalfa, p1, Rc, zetas, mang, 
                                                    vectsfinalfa, tipo, gamma_couple, NC_plus);
             } 
             else {
                grid_rho[i] = rho_radial_finite_ext(nt, elecalfa, p1, Rc, zetas, mang,
                                                    vectsfinalfa, tipo, gamma_couple, NC_plus);
                grid_rho_beta[i] =  rho_radial_finite_ext(nt, elecbeta, p1, Rc, zetas, mang,
                                                          vectsfinbeta, tipo, gamma_couple, NC_plus);
             }
          }
       }  /* ends stos*/
       else{
          if(strcmp(basis,"GTOs") == 0){
             for(i = 0; i < n_points; i++) {
                p1 = grid[i];
                if(compara == 0) {
                   grid_rho[i] = rho_radial_gto(nt, elecalfa, p1, expo, np, vectsfinalfa, tipo);
                }
                else {
                   grid_rho[i] = rho_radial_gto(nt, elecalfa, p1, expo, np, vectsfinalfa, tipo);
                   grid_rho_beta[i] = rho_radial_gto(nt, elecbeta, p1, expo, np, vectsfinbeta, tipo);
                }
             }
          }  /* ends gtos */
       }
    }  /*ends finite */
    else {
       if(strcmp(bound,"parabolic") == 0){
          if(strcmp(basis,"GTOs") == 0){
             for(i = 0; i < n_points; i++) {
                p1 = grid[i];
                if(compara == 0) {
                   grid_rho[i] = rho_radial_gto(nt, elecalfa, p1, expo, np, vectsfinalfa, tipo);
                }
                else {
                   grid_rho[i] = rho_radial_gto(nt, elecalfa, p1, expo, np, vectsfinalfa, tipo);
                   grid_rho_beta[i] = rho_radial_gto(nt, elecbeta, p1, expo, np, vectsfinbeta, tipo);
                }
//             printf("rho[%d] = %f \n", i, grid_rho[i]);
             }
          }  /* ends gtos*/
       } /* ends parabolic */
       else {
          if(strcmp(bound,"confined") == 0){
             if(strcmp(basis,"STOs") == 0){
                for(i = 0; i < n_points; i++) { //for save_temp
                   p1 = grid[i];
                   if(compara == 0) { 
                      grid_rho[i] = rho_radial_confined(nt, elecalfa, p1, Rc, expo, np, vectsfinalfa,
                                                        tipo, arreglo_factorial, arreglo_inv_factorial);
                   }
                   else {
                      grid_rho[i] = rho_radial_confined(nt, elecalfa, p1, Rc, expo, np, vectsfinalfa,
                                                        tipo, arreglo_factorial, arreglo_inv_factorial);
                      grid_rho_beta[i] = rho_radial_confined(nt, elecbeta, p1, Rc, expo, np, vectsfinbeta,
                                                             tipo, arreglo_factorial, arreglo_inv_factorial);
                   }
                }
             }  /* ends stos */
             else{
                if(strcmp(basis,"GTOs") == 0){
                   for(i = 0; i < n_points; i++) { //for save_temp
                      p1 = grid[i];
                      if(compara == 0) {
                         grid_rho[i] = rho_radial_gto_imp(nt, elecalfa, p1, Rc, expo, np, vectsfinalfa, tipo);
                      }
                      else {
                         grid_rho[i] = rho_radial_gto_imp(nt, elecalfa, p1, Rc, expo, np, vectsfinalfa, tipo);
                         grid_rho_beta[i] = rho_radial_gto_imp(nt, elecbeta, p1, Rc, expo, np, vectsfinbeta, tipo);
                      }
                   }
                }
             } /* ends gtos */
          } /* ends confined */
          else {
             if(strcmp(bound,"dielectricc") == 0){
                if(strcmp(basis,"GTOs") == 0){
                   for(i = 0; i < n_points; i++) {
                      p1 = grid[i];
                      if(compara == 0) {
                         grid_rho[i] = rho_radial_gto(nt, elecalfa, p1, expo, np, vectsfinalfa, tipo);
                      }
                      else {
                         grid_rho[i] = rho_radial_gto(nt, elecalfa, p1, expo, np, vectsfinalfa, tipo);
                         grid_rho_beta[i] = rho_radial_gto(nt, elecbeta, p1, expo, np, vectsfinbeta, tipo);
                      }
                   }
                }  /* ends gtos */
             } /* ends dielec */
          }
       }
    }
 }

  return 0;
}//function

int grid_rhorad_orbital(int           z,
                char         *using_gamma,
                int           compara,
                char         *bound,
                int           nt,
                int           elecalfa,
                int           elecbeta,
                double        Rc,
                double       *expo,
                int          *np,
                double       *zetas,
                int          *mang,
                double       *vectsfinalfa,
                double       *vectsfinbeta,
                char         *tipo,
                double       *NC_minus,
                double       *NC_plus,
                double        gamma_couple,
                double       *grid,
                double       *grid_rho,
                double       *grid_rho_beta,
                int           n_points,
                int           n_boundary,
                double       *arreglo_factorial,
                double       *arreglo_inv_factorial,
                int           selected_orb)
 
 {//function

 extern double rho_orbital_radial_finite_int(char   *using_gamma,
                                     int     nt,
                                     int     elec,
                                     double  r,
                                     double  Rc,
                                     double *zeta,
                                     int    *np,
                                     double *vectors,
                                     char   *tipo,
                                     double  gamma_couple,
                                     double *N_minus,
                                     int     selected_orb);

 extern double rho_orbital_radial_finite_ext(int     nt,
                                     int     elec,
                                     double  r,
                                     double  Rc,
                                     double *alfa,
                                     int    *mang,
                                     double *vectors,
                                     char   *tipo,
                                     double  gamma_couple,
                                     double *N_plus,
                                     int     selected_orb);

 int unsigned i;
 double p1;


 for (i = 0; i < n_boundary; i++) { //for save_temp
   p1 = grid[i];
   if (compara == 0) {
     grid_rho[i] = rho_orbital_radial_finite_int(using_gamma, nt, elecalfa, p1, Rc, expo, np,
                                    vectsfinalfa, tipo, gamma_couple, NC_minus, selected_orb);
   } else {
         grid_rho[i] = rho_orbital_radial_finite_int(using_gamma, nt, elecalfa, p1, Rc, expo, np,
                                        vectsfinalfa, tipo, gamma_couple, NC_minus, selected_orb);
         grid_rho_beta[i] = rho_orbital_radial_finite_int(using_gamma, nt, elecbeta, p1, Rc, expo, np,
                                             vectsfinbeta, tipo, gamma_couple, NC_minus, selected_orb);
   }
 } //for save_temp
         
 for (i = n_boundary; i < n_points; i++) { //for save_temp
   p1 = grid[i];
   if (compara == 0) {
     grid_rho[i] = rho_orbital_radial_finite_ext(nt, elecalfa, p1, Rc, zetas, mang, 
                           vectsfinalfa, tipo, gamma_couple, NC_plus, selected_orb);
   } else {
         grid_rho[i] = rho_orbital_radial_finite_ext(nt, elecalfa, p1, Rc, zetas, mang,
                               vectsfinalfa, tipo, gamma_couple, NC_plus, selected_orb);
         
         grid_rho_beta[i] =  rho_orbital_radial_finite_ext(nt, elecbeta, p1, Rc, zetas, mang,
                                     vectsfinbeta, tipo, gamma_couple, NC_plus, selected_orb);
     }
 }
 return 0;
}//function

