#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

int print_out_array(int points, double *grid, double *array, char *name_file) {
  int i;
  FILE *target;

  target = fopen(name_file,"w");

  for (i = 0; i < points; i++)
    fprintf(target,"%20.8f  %20.8f\n", grid[i], array[i]);

  fclose(target);

  return 0;
}

int x_slater(double  rho, 
             double *energy,
             double *potential) 
{
 double coef13, coef43, coefx;
 
 coef13 = 1.f/3.f;
 coef43 = 4.f*coef13;
 coefx = -0.7386f;
 *energy = coefx*pow(rho,coef43);
 *potential = coef43*coefx*pow(rho,coef13);
 return 0;
}

double xc_energy(double *correlationc,
                 double *exchangex,
                 int     points, 
                 int     compara,
                 char   **save_dft,
                 double *weight_dft,
                 int     flag_dft,
                 double *energy_array,
                 double *grid, 
                 double *rho_alpha, 
                 double *rho_beta, 
                 double *derho_alpha, 
                 double *derho_beta,
                 double *d2rho_alpha, 
                 double *d2rho_beta) 
{

 extern double numerical_int(double   *grid,
                             double   *grid_fun,
                             int       points);


 extern int x_slater(double rho, double *energy, double *potential);
 extern void x_slater_sp(double *rhoa, double *rhob, double *energy, double *vxa, double *vxb);
 extern void pw92sr(double *rho_local, double *local_energy, double *potential);
 extern void becke_sr(double *ri, double *r2, double *rho_local, double *drho_local, double *d2rho_local,
                      double *local_energy, double *potential);
 extern void pbe_sr(double *rho_local,double *drho_local,
                    double *local_energy,double *potential, double *nolocalpot);
 extern void lyp_sr(double *ri, double *rho, double *drho, double *d2rho, double *local_energy, double *potential);
//////////////
  int cont;
  double ri, r2;
  double rho_local, local_energy, local_energy_corr, local_potential, nolocal_pot, total_e_xc, Pi,
         drho_local, d2rho_local;
  //Open-Shell
  double rho_local_a, rho_local_b, local_potential_a, local_potential_b;

  Pi = 4.f*atan(1.f);
  int dft;
  double weight;
  total_e_xc = 0.f;

  double exchange; 
  exchange = 0.f;
 
  double correlation;
  correlation = 0.f;

  if (compara == 0) {
     for (dft = 1; dft < flag_dft; dft++)
      {
      weight = weight_dft[dft];
     
      if(strcmp(save_dft[dft], "slater") == 0) {
     
            for (cont = 0; cont < points; cont++) {
              rho_local = rho_alpha[cont];
              x_slater(rho_local, &local_energy, &local_potential);
              energy_array[cont] = local_energy;
            }
            exchange   = 4.f*Pi*numerical_int(grid, energy_array, points)*weight;
            exchangex[dft] = exchange;
            total_e_xc = total_e_xc + exchange;
       } else
       if(strcmp(save_dft[dft], "pbe") == 0) {
         for (cont = 0; cont < points; cont++) {
           rho_local = rho_alpha[cont];
           drho_local = derho_alpha[cont];
           pbe_sr_(&rho_local, &drho_local, &local_energy, &local_potential, &nolocal_pot);
           energy_array[cont] = local_energy;
          }
        exchange = 4.f*Pi*numerical_int(grid, energy_array, points)*weight;
        exchangex[dft] = exchange;
        total_e_xc = total_e_xc + exchange;
       } else
       if(strcmp(save_dft[dft], "becke") == 0) {
         for (cont = 0; cont < points; cont++) {
           rho_local = rho_alpha[cont];
           drho_local = derho_alpha[cont];
           d2rho_local = d2rho_alpha[cont];
           ri = grid[cont];
           r2 = grid[1];
           becke_sr_(&ri, &r2, &rho_local, &drho_local, &d2rho_local, &local_energy, &local_potential);
           energy_array[cont] = local_energy;
         }
         exchange   = 4.f*Pi*numerical_int(grid, energy_array, points)*weight;
         exchangex[dft] = exchange;
               		total_e_xc = total_e_xc + exchange;
                 
                 }  else
                        if(strcmp(save_dft[dft], "pwl") == 0) {
                          for (cont = 0; cont < points; cont++) {
                          rho_local = rho_alpha[cont];
                     	  pw92sr_(&rho_local, &local_energy, &local_potential);
                          energy_array[cont] = local_energy;
                         }
                         correlation = 4.f*Pi*numerical_int(grid, energy_array, points)*weight;
                         correlationc[dft] = correlation;
                         total_e_xc  = total_e_xc + correlation;
                 }  else 
                               if(strcmp(save_dft[dft], "lyp") == 0) {
                                  for (cont = 0; cont < points; cont++) {
                                  rho_local = rho_alpha[cont];
                                  drho_local = derho_alpha[cont];
                         	  d2rho_local = d2rho_alpha[cont];
                         	  ri = grid[cont];
                                  lyp_sr_(&ri, &rho_local, &drho_local, &d2rho_local, &local_energy, &local_potential);
                                  energy_array[cont] = local_energy;
                                 }
                      		 correlation = 4.f*Pi*numerical_int(grid, energy_array, points)*weight;
                       		 correlationc[dft] = correlation;
                       	         total_e_xc  = total_e_xc + correlation;
                          }  else
                                  if(strcmp(save_dft[dft],"hf") != 0)
                                  printf("\n>>>> XC FUCTIONAL NOT IMPLEMENTED YET %s\n", save_dft[dft]);
     
                }
  } else {//Open-Shell
      for (dft = 1; dft < flag_dft; dft++)
      {
      weight = weight_dft[dft];

      if(strcmp(save_dft[dft], "slater") == 0) {

            for (cont = 0; cont < points; cont++) {
              rho_local_a = rho_alpha[cont];
              rho_local_b = rho_beta[cont];
              x_slater_sp_(&rho_local_a, &rho_local_b, &local_energy, &local_potential_a, &local_potential_b);
              energy_array[cont] = local_energy;
            }
            total_e_xc = total_e_xc + 4.f*Pi*numerical_int(grid, energy_array, points)*weight;
      } else
 //remove or apply changes here       if(strcmp(save_dft[dft], "becke") == 0) {
 //remove or apply changes here
 //remove or apply changes here             for (cont = 0; cont < points; cont++) {
 //remove or apply changes here               rho_local = rho_alpha[cont];
 //remove or apply changes here               drho_local = derho_alpha[cont];
 //remove or apply changes here               d2rho_local = d2rho_alpha[cont];
 //remove or apply changes here               ri = grid[cont];
 //remove or apply changes here               r2 = grid[1];
 //remove or apply changes here               becke_sr_(&ri, &r2, &rho_local, &drho_local, &d2rho_local, &local_energy, &local_potential);
 //remove or apply changes here               energy_array[cont] = local_energy;
 //remove or apply changes here             }
 //remove or apply changes here             total_e_xc = total_e_xc + 4.f*Pi*numerical_int(grid, energy_array, points)*weight;
 //remove or apply changes here       }  else
 //remove or apply changes here           if(strcmp(save_dft[dft], "perdew") == 0) {
 //remove or apply changes here
 //remove or apply changes here                 for (cont = 0; cont < points; cont++) {
 //remove or apply changes here                   rho_local = rho_alpha[cont];
 //remove or apply changes here                   pw92sr_(&rho_local, &local_energy, &local_potential);
 //remove or apply changes here                   energy_array[cont] = local_energy;
 //remove or apply changes here                 }
 //remove or apply changes here              total_e_xc = total_e_xc + 4.f*Pi*numerical_int(grid, energy_array, points)*weight;
 //remove or apply changes here           }  else
                 if(strcmp(save_dft[dft],"hf") != 0)
                   printf("\n>>>> NOT HAVE THIS FUCTIONAL %s\n", save_dft[dft]);

      }
    }


  return (total_e_xc);
}


 void xc_potential (int           points,
                    int           compara,
                    int           point_interm,
                    char        **save_dft,
                    double       *weight_dft,
                    int           flag_dft,
                    char         *using_gamma,
                    char         *bound,
                    double        Rc,
                    int           nt,
                    double       *expo,
                    int          *np,
                    double       *zetas,
                    int          *mang,
                    int          *ncm,
                    double       *NC_minus,
                    double       *NC_plus,
                    double        gamma_couple,
                    double       *pot_array,
                    double       *pot_array_beta,
                    double       *local_int,
                    double       *local_int_beta,
                    double       *grid,
                    double       *rho_alpha,
                    double       *rho_beta,
                    double       *derho_alpha, 
                    double       *derho_beta,
                    double       *d2rho_alpha, 
                    double       *d2rho_beta, 
                    double       *mat_ks,
                    double       *mat_ks_beta,
                    double       *arreglo_factorial,
                    double       *arreglo_inv_factorial,
                    double       *nolocal_pot_array)
 {
 extern double numerical_int(double   *grid,
                             double   *grid_fun,
                             int       points);
// Definition of different basis set fucntions
/*free*/
 extern double sto(int mu,double r,double *expo,int *np);
/*impenetrable*/
 extern double msto(int mu,double  r,double  Rc,double *expo,int *np,double *arreglo_factorial,double *arreglo_inv_factorial);
/*penetrable*/
 extern double msto_finite_int(char *using_gamma,int mu,double r,double Rc,int *np,double gamma_couple,double *zeta,double *N_minus);
 extern double msto_finite_ext(int mu,double  r,double  Rc,int *mang,double  gamma_couple,double *alfa,double *N_plus);

// Derivatives
/*free*/
 extern double der_sto_r(int mu, double r, double* expo, int* np);
/*impenetrable*/ 
 extern double der_msto_r(int mu, double r, double Rc, double* expo, int* np, double *arreglo_factorial, double *arreglo_inv_factorial);
/*penetrable(intern,extern)*/
 extern double der_msto_finite_int(char *using_gamma,int mu,double r,double Rc,int *np,double gamma_couple,double *zeta,double *N_minus);
 extern double der_msto_finite_ext(int mu,double r,double Rc,int *mang,double gamma_couple,double *alfa,double *N_plus);

 extern int indexes(int, int, int*, int*);
 extern int delta_kro_(int*, int*, double*);
 extern void pw92sr(double *rho_local, double *local_energy, double *potential);
 extern void becke_sr(double *ri, double *r2, double *rho_local, double *drho_local, double *d2rho_local,
                      double *local_energy, double *potential);
 extern void lyp_sr(double *ri, double *rho, double *drho, double *d2rho, double *local_energy, double *potential);

 extern void pbe_sr(double *rho_local,double *drho_local,double *local_energy,double *potential, double *nolocalpot);

 extern int x_slater(double rho, double *energy, double *potential);
 extern void x_slater_sp(double *rhoa, double *rhob, double *energy, double *vxa, double *vxb);
 double Pi, rho_local, local_energy, local_potential, nolocal_pot, local_potential_corr,drho_local, d2rho_local;
 double ri, r2;
 int cont, k, i, otro_k;
 int total_elements, index_i, index_j, ang_i, ang_j, ncm_i, ncm_j;
 double p1,temp_potnol;;
 

 otro_k = 0;

 Pi = 4.f*atan(1.f);

  double delta, delta1;
  total_elements = nt*nt;

  double temp_add;
  double temp_add_beta;
  double weight; 
  int     dft;
  temp_add = 0.f;
  temp_add_beta = 0.f;
  dft      = 0;
  //Open-Shell
  double rho_local_a, rho_local_b, local_potential_a, local_potential_b;


  for (cont = 0; cont < points; cont++) {
    pot_array[cont] = 0.f;
    nolocal_pot_array[cont] = 0.f;
  }

  if(compara == 0) {
    for (dft = 1; dft < flag_dft; dft++)
     {
     weight = weight_dft[dft];
     if(strcmp(save_dft[dft], "slater") == 0) {
           for (cont = 0; cont < points; cont++) {
             rho_local = rho_alpha[cont];
             if (rho_local > 1e-16) 
              x_slater(rho_local, &local_energy, &local_potential);
             else 
              local_potential = 0.f;
             temp_add = pot_array[cont];
             pot_array[cont] = temp_add + local_potential*weight;
           }
     } else
        
        if(strcmp(save_dft[dft], "becke") == 0) {
              for (cont = 0; cont < points; cont++) {
                rho_local = rho_alpha[cont];
                drho_local = derho_alpha[cont];
                d2rho_local = d2rho_alpha[cont];
                ri = grid[cont];
                r2 = grid[1];
                if (rho_local > 1e-16) 
                  becke_sr_(&ri, &r2, &rho_local, &drho_local, &d2rho_local, &local_energy, &local_potential);
                else
                  local_potential = 0.f;             
                temp_add = pot_array[cont];
                pot_array[cont] = temp_add + local_potential*weight;
              }
        }  else
            if(strcmp(save_dft[dft], "pwl") == 0) {
                  
                  for (cont = 0; cont < points; cont++) {
                    rho_local = rho_alpha[cont];
                    if (rho_local > 1e-16) 
                     pw92sr_(&rho_local, &local_energy, &local_potential);
                    else 
                     local_potential = 0.f;
                    temp_add = pot_array[cont];
                    pot_array[cont] = temp_add + local_potential*weight;
                  }
            }  else 
                 if(strcmp(save_dft[dft], "lyp") == 0) {
                       for (cont = 0; cont < points; cont++) {
                         rho_local = rho_alpha[cont];
                         drho_local = derho_alpha[cont];
                         d2rho_local = d2rho_alpha[cont];
                         ri = grid[cont];
                         if (rho_local > 1e-16) 
                           lyp_sr_(&ri, &rho_local, &drho_local, &d2rho_local, &local_energy, &local_potential);
                         else
                           local_potential = 0.f;             
                         temp_add = pot_array[cont];
                         pot_array[cont] = temp_add + local_potential*weight;
                       }
                  } else
                  if(strcmp(save_dft[dft], "pbe") == 0) {
                    for (cont = 0; cont < points; cont++) {
                      rho_local = rho_alpha[cont];
                      drho_local = derho_alpha[cont];
                      if (rho_local > 1e-16) 
                        pbe_sr_(&rho_local, &drho_local, &local_energy, &local_potential, &nolocal_pot);
                      else {
                        local_potential = 0.f;             
                        nolocal_pot = 0.f;             
                      }
                      temp_add = pot_array[cont];
                      pot_array[cont] = temp_add + local_potential*weight;
                      nolocal_pot_array[cont] = nolocal_pot*weight;
                    }
                  }  else
                  if(strcmp(save_dft[dft],"hf") != 0)
                              printf("\n>>>> NOT HAVE THIS FUCTIONAL %s\n", save_dft[dft]);
     
     }



     for (k = 0; k < total_elements; k++) {
       indexes(nt, k, &index_i, &index_j);
       index_i = index_i - 1;
       index_j = index_j - 1;
       if (index_j >= index_i ) {
         ang_i = mang[index_i];
         ang_j = mang[index_j];
         delta_kro_(&ang_i, &ang_j, &delta);
         ncm_i = ncm[index_i];
         ncm_j = ncm[index_j];
         delta_kro_(&ncm_i, &ncm_j, &delta1);
         delta = delta*delta1;
         if (delta == (double)0)
           mat_ks[k] = (double) 0.f;
         else {//else
//CONF. PENETRABLE    
           if(strcmp(bound,"finite") == 0) {
             for (i = 0; i < points; i++) {
               p1 = grid[i];
                 if (i < point_interm){ 
                    local_int[i]  = pot_array[i]*msto_finite_int(using_gamma, index_i, p1, Rc, np, gamma_couple, expo, NC_minus)*msto_finite_int(using_gamma, index_j, p1, Rc, np, gamma_couple, expo, NC_minus);
                  
                   temp_potnol=msto_finite_int(using_gamma, index_i, p1, Rc, np, gamma_couple, expo, NC_minus)*der_msto_finite_int(using_gamma,index_j, p1, Rc, np, gamma_couple,expo, NC_minus);

                  temp_potnol = temp_potnol +msto_finite_int(using_gamma, index_j, p1, Rc, np, gamma_couple, expo, NC_minus)*der_msto_finite_int(using_gamma, index_i, p1, Rc, np, gamma_couple, expo, NC_minus);
                   local_int[i] = local_int[i] + nolocal_pot_array[i]*temp_potnol;
	       
                } else {
                 local_int[i]  = pot_array[i]*msto_finite_ext(index_i, p1, Rc, mang, gamma_couple, zetas, NC_plus)*msto_finite_ext(index_j, p1, Rc, mang, gamma_couple, zetas, NC_plus);
                    
                   temp_potnol= msto_finite_ext(index_i, p1, Rc, mang, gamma_couple, zetas, NC_plus)*der_msto_finite_ext(index_j, p1, Rc, mang, gamma_couple, zetas, NC_plus);
                   
                  temp_potnol=temp_potnol+msto_finite_ext(index_j, p1, Rc, mang, gamma_couple, zetas, NC_plus)*der_msto_finite_ext(index_i, p1, Rc, mang, gamma_couple, zetas, NC_plus);
                   local_int[i]=local_int[i] + nolocal_pot_array[i]*temp_potnol;
                   
             }
            }
          } else
//LIBRE    
           if(strcmp(bound,"free") == 0) {
              for (i = 0; i < points; i++) {
                            p1 = grid[i];
                 local_int[i]  = pot_array[i]*sto(index_i, p1, expo, np)*sto(index_j, p1, expo, np);
               temp_potnol = sto(index_i, p1, expo, np)*der_sto_r(index_j, p1, expo, np);
               temp_potnol = temp_potnol + sto(index_j, p1, expo, np)*der_sto_r(index_i, p1, expo, np);
               local_int[i]  = local_int[i] + nolocal_pot_array[i]*temp_potnol;
             }
           } else
     
///CONF. IMPENETRABLE

           if(strcmp(bound,"confined") == 0) {
             for (i = 0; i < points; i++) {
             p1 = grid[i];
             local_int[i]  = pot_array[i]*msto(index_i, p1, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial)*msto(index_j, p1, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial);
             temp_potnol = msto(index_i, p1, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial)*der_msto_r(index_j, p1, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial);
             temp_potnol = temp_potnol + msto(index_j, p1, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial)*der_msto_r(index_i, p1, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial);
             local_int[i]  = local_int[i] + nolocal_pot_array[i]*temp_potnol;
             }
           } else 
              printf("\nI do not have correct flag\n");
           mat_ks[k] = numerical_int(grid, local_int, points);
         }
       }
       else {
         otro_k = index_j*nt + index_i;
         mat_ks[k] = mat_ks[otro_k];
       }
     }

  } else {//Open-Shell
        for (cont = 0; cont < points; cont++) {
           pot_array_beta[cont] = 0.f;
        }
        for (dft = 1; dft < flag_dft; dft++)
         {
         weight = weight_dft[dft];
         if(strcmp(save_dft[dft], "slater") == 0) {
      
               for (cont = 0; cont < points; cont++) {
                 rho_local_a = rho_alpha[cont];
                 rho_local_b = rho_beta[cont];
                 x_slater_sp_(&rho_local_a, &rho_local_b, &local_energy, &local_potential_a, &local_potential_b);

                 temp_add = pot_array[cont];
                 pot_array[cont] = temp_add + local_potential_a*weight;

                 temp_add_beta = pot_array_beta[cont];
                 pot_array_beta[cont] = temp_add_beta + local_potential_b*weight;
               }
         } else
      
//mrb remove or apply changes           if(strcmp(save_dft[dft], "becke") == 0) {
//mrb remove or apply changes     
//mrb remove or apply changes                 for (cont = 0; cont < points; cont++) {
//mrb remove or apply changes                   rho_local = rho_alpha[cont];
//mrb remove or apply changes                   drho_local = derho_alpha[cont];
//mrb remove or apply changes                   d2rho_local = d2rho_alpha[cont];
//mrb remove or apply changes                   ri = grid[cont];
//mrb remove or apply changes                   r2 = grid[1];
//mrb remove or apply changes                   becke_sr_(&ri, &r2, &rho_local, &drho_local, &d2rho_local, &local_energy, &local_potential);
//mrb remove or apply changes                   temp_add = pot_array[cont];
//mrb remove or apply changes                   pot_array[cont] = temp_add + local_potential*weight;
//mrb remove or apply changes                 }
//mrb remove or apply changes           }  else
//mrb remove or apply changes               if(strcmp(save_dft[dft], "perdew") == 0) {
//mrb remove or apply changes     
//mrb remove or apply changes                     for (cont = 0; cont < points; cont++) {
//mrb remove or apply changes                       rho_local = rho_alpha[cont];
//mrb remove or apply changes                       pw92sr_(&rho_local, &local_energy, &local_potential);
//mrb remove or apply changes                       temp_add = pot_array[cont];
//mrb remove or apply changes                       pot_array[cont] = temp_add + local_potential*weight;
//mrb remove or apply changes                     }
//mrb remove or apply changes               }  else
                    if(strcmp(save_dft[dft],"hf") != 0)
                      printf("\n>>>> NOT HAVE THIS FUCTIONAL %s\n", save_dft[dft]);
      
         }
         double msto_all;
         for (k = 0; k < total_elements; k++) {
           indexes(nt, k, &index_i, &index_j);
           index_i = index_i - 1;
           index_j = index_j - 1;
           if (index_j >= index_i ) {
             ang_i = mang[index_i];
             ang_j = mang[index_j];
             delta_kro_(&ang_i, &ang_j, &delta);
             ncm_i = ncm[index_i];
             ncm_j = ncm[index_j];
             delta_kro_(&ncm_i, &ncm_j, &delta1);
             delta = delta*delta1;
             if (delta == (double)0) {
               mat_ks[k] = (double) 0.f;
               mat_ks_beta[k] = (double) 0.f;
             }
             else {//else
         
               if(strcmp(bound,"finite") == 0) {
                 for (i = 0; i < points; i++) {
                   p1 = grid[i];
                   if (i < point_interm) {
                    msto_all = msto_finite_int(using_gamma, index_i, p1, Rc, np, gamma_couple, expo, NC_minus)*msto_finite_int(using_gamma, index_j, p1, Rc, np, gamma_couple, expo, NC_minus);
                    local_int[i]  = pot_array[i]*msto_all;
                    local_int_beta[i]  = pot_array_beta[i]*msto_all;
                   }
                   else { 
                    msto_all = msto_finite_ext(index_i, p1, Rc, mang, gamma_couple, zetas, NC_plus)*msto_finite_ext(index_j, p1, Rc, mang, gamma_couple, zetas, NC_plus);
                    local_int[i]  = pot_array[i]*msto_all;
                    local_int_beta[i]  = pot_array_beta[i]*msto_all;
                   }
                 }
                } else
         
               if(strcmp(bound,"free") == 0) {
                 for (i = 0; i < points; i++) {
                   p1 = grid[i];
                   msto_all = sto(index_i, p1, expo, np)*sto(index_j, p1, expo, np);
                   local_int[i]  = pot_array[i]*msto_all;
                   local_int_beta[i]  = pot_array_beta[i]*msto_all;
                 }
               } else
         
               if(strcmp(bound,"confined") == 0) {
                 for (i = 0; i < points; i++) {
                 p1 = grid[i];
                 msto_all = msto(index_i, p1, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial)*msto(index_j, p1, Rc, expo, np, arreglo_factorial, arreglo_inv_factorial);
                 local_int[i]  = pot_array[i]*msto_all;
                 local_int_beta[i]  = pot_array_beta[i]*msto_all;
                 }
               } else 
                  printf("\nI do not have correct flag\n");
               mat_ks[k] = numerical_int(grid, local_int, points);
               mat_ks_beta[k] = numerical_int(grid, local_int_beta, points);
             }
           }
           else {
             otro_k = index_j*nt + index_i;
             mat_ks[k] = mat_ks[otro_k];
             mat_ks_beta[k] = mat_ks_beta[otro_k];
           }
         }
    }
 
}

