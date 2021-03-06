#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

 int grid_secderrad(int           z,                                    //Derivada radial "derrad"
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
                    double       *grid_secder,
                    double       *grid_secder_beta,
                    int           n_points,
                    int           n_boundary,
                    double       *arreglo_factorial,
                    double       *arreglo_inv_factorial)
 
 {//function

 extern double der_rho_radial_confined(int     nt,
                                       int     elec,
                                       double  r,
                                       double  Rc,
                                       double *expo,
                                       int    *np,
                                       double *vectors,
                                       char   *tipo,
                                       double *arreglo_factorial,
                                       double *arreglo_inv_factorial);

 extern double der_rho_radial_finite_int(char   *using_gamma,
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

 extern  double der_rho_radial_finite_ext(int     nt,
                                          int     elec,
                                          double  r,
                                          double  Rc,
                                          double *alfa,
                                          int    *mang,
                                          double *vectors,
                                          char   *tipo,
                                          double  gamma_couple,
                                          double *N_plus);

 extern double der_rho_radial_free(int     nt, 
                                   int     elec, 
                                   double  r, 
                                   double *expo,
                                   int    *np, 
                                   double *vectors, 
                                   char   *tipo);

 extern double sec_der_rho_radial_free(int     nt, 
                                       int     elec, 
                                       double  r, 
                                       double *expo,
                                       int    *np, 
                                       double *vectors, 
                                       char   *tipo);

 extern double sec_der_rho_radial_finite_int(char   *using_gamma,
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

 extern double sec_der_rho_radial_finite_ext(int     nt,
                                             int     elec,
                                             double  r,
                                             double  Rc,
                                             double *alfa,
                                             int    *mang,
                                             double *vectors,
                                             char   *tipo,
                                             double  gamma_couple,
                                             double *N_plus);

 extern double sec_der_rho_radial_confined(int nt, int elec, double r, double Rc, double* expo,
                                           int* np, double* vectors, char* tipo, double *arreglo_factorial,
                                           double *arreglo_inv_factorial);

 extern double sec_der_rho_radial_gto(int nt, int elec, double r, double* expo, int* np, double* vectors, char* tipo);
 extern double sec_der_rho_radial_gto_imp(int nt, int elec, double r, double r0, double* expo, int* np, double* vectors, char* tipo);




 int unsigned i;
 double p1;

 if(strcmp(basis,"stos") == 0){     /* begins stos */
   if (strcmp(bound,"free") == 0 || strcmp(bound,"debye") == 0 || strcmp(bound,"yukawa") == 0 || strcmp(bound,"baimbetov") == 0) {
     for (i = 0; i < n_points; i++) { //for save_temp
        p1 = grid[i];
        if (compara == 0) {
          grid_secder[i] = sec_der_rho_radial_free(nt, elecalfa, p1, expo, np,                  //Capa cerrada 
                                                   vectsfinalfa, tipo);

        } else {
          grid_secder[i] = sec_der_rho_radial_free(nt, elecalfa, p1, expo, np,                  //Capa abierta
                                                   vectsfinalfa, tipo);

          grid_secder_beta[i] = sec_der_rho_radial_free(nt, elecbeta, p1, expo, np,       
                                                        vectsfinbeta, tipo);
          }
     }

  } else
    if (strcmp(bound,"confined") == 0) { //if confined
       for (i = 0; i < n_boundary; i++) { //for save_temp
          p1 = grid[i];
          if (compara == 0) {
            grid_secder[i] = sec_der_rho_radial_confined(nt, elecalfa, p1, Rc, expo, np, vectsfinalfa, 
                                                         tipo, arreglo_factorial, arreglo_inv_factorial);
          } else {
            grid_secder[i] = sec_der_rho_radial_confined(nt, elecalfa, p1, Rc, expo, np, vectsfinalfa,
                                                         tipo, arreglo_factorial, arreglo_inv_factorial);
   
            grid_secder_beta[i] = sec_der_rho_radial_confined(nt, elecbeta, p1, Rc, expo, np, vectsfinbeta,
                                                              tipo, arreglo_factorial, arreglo_inv_factorial);
            }
       } //for save_temp
    } else 
      if (strcmp(bound,"finite") == 0) { //if finite
         for (i = 0; i < n_boundary; i++) { //for save_temp
            p1 = grid[i];
            if (compara == 0) {
              grid_secder[i] =  sec_der_rho_radial_finite_int(using_gamma, nt, elecalfa, p1, Rc, expo, np, 
                                                              vectsfinalfa, tipo, gamma_couple, NC_minus);
            } else {
              grid_secder[i] =  sec_der_rho_radial_finite_int(using_gamma, nt, elecalfa, p1, Rc, expo, np,
                                                              vectsfinalfa, tipo, gamma_couple, NC_minus);

              grid_secder_beta[i] = sec_der_rho_radial_finite_int(using_gamma, nt, elecbeta, p1, Rc, expo, np,
                                                                  vectsfinbeta, tipo, gamma_couple, NC_minus);
              }
         } //for save_temp
         
         for (i = n_boundary; i < n_points; i++) { //for save_temp
            p1 = grid[i];
            if (compara == 0) {
              grid_secder[i] = sec_der_rho_radial_finite_ext(nt, elecalfa, p1, Rc, zetas, mang, vectsfinalfa,
                                                             tipo, gamma_couple, NC_plus);
            } else {
              grid_secder[i] = sec_der_rho_radial_finite_ext(nt, elecalfa, p1, Rc, zetas, mang, vectsfinalfa,
                                                             tipo, gamma_couple, NC_plus);

              grid_secder_beta[i] = sec_der_rho_radial_finite_ext(nt, elecbeta, p1, Rc, zetas, mang, vectsfinbeta,
                                                                  tipo, gamma_couple, NC_plus);
              }
            }
                
      } //if finite
  } /* ends stos */
  else{
     if(strcmp(basis,"gtos") == 0){
        if(strcmp(bound,"confined") == 0){
           for(i = 0; i < n_boundary; i++) { //for save_temp
              p1 = grid[i];
              if(compara == 0) {
                 grid_secder[i] = sec_der_rho_radial_gto_imp(nt, elecalfa, p1, Rc, expo, np, vectsfinalfa, tipo);
              } 
              else {
                 grid_secder[i] = sec_der_rho_radial_gto_imp(nt, elecalfa, p1, Rc, expo, np, vectsfinalfa, tipo);

                 grid_secder_beta[i] = sec_der_rho_radial_gto_imp(nt, elecbeta, p1, Rc, expo, np, vectsfinbeta, tipo);

              }
           } //for save_temp
        }  /* ends confined */
        else {     /* begins the rest of the cases: free, finite, dielec, parabolic */
           for(i = 0; i < n_points; i++) { //for save_temp
              p1 = grid[i];
              if(compara == 0) {
                 grid_secder[i] = sec_der_rho_radial_gto(nt, elecalfa, p1, expo, np, vectsfinalfa, tipo);
              } 
              else {
                 grid_secder[i] = sec_der_rho_radial_gto(nt, elecalfa, p1, expo, np, vectsfinalfa, tipo);

                 grid_secder_beta[i] = sec_der_rho_radial_gto(nt, elecbeta, p1, expo, np, vectsfinbeta, tipo);
              }
           } 
        }  /* ends the rest of the cases */
     }     /* ends gtos */
  }

//////////////////////
  return 0;
 }//function

