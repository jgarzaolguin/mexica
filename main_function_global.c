#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


int  main_finite_global(double   z,
                               int      elecalfa,
                               int      elecbeta,
                               char    *tipo,
                               char    *using_gamma,
                               char    *correlation,
                               char    *propagador,
                               char   **save_dft,
                               double *weight_dft,
                               int     flag_dft,
                               char    *bound,
                               char    *espin,
                               double   Rc,
                               double   tol,
                               double   mezcla,
                               int      orbital,
                               int      maxiter,
                               int      nt,
                               int     *np,
                               int     *mang,
                               int     *ncm,
                               double   gamma_couple,
                               double  *expo,
                               char    *opt,
                               double   step,
                               int     *opt_flag,
                               char    *nombre,
                               double   epsilon,
                               char    *kind_of_cal,
                               double   count_temp,
                               double   count_final,
                               int      steps,
                               double  *energia,
                               int      plasma,
                               char    *properties)

{//label 1

extern int optimiza_main(int     nt,
                         char   *opt,
                         int     elecalfa,
                         int     elecbeta,
                         double  z,
                         int     orbital,
                         double  tol,
                         char   *using_gamma,
                         char   *correlation,
                         char   *propagador,
                         char   **save_dft,
                         double *weight_dft,
                         int     flag_dft,
                         double  mezcla,
                         double *expo,
                         int    *np,
                         int    *mang,
                         int    *ncm,
                         double  gamma_couple,
                         char   *tipo,
                         int     maxiter,
                         double  Rc,
                         char   *bound,
                         char   *espin,
                         double  step,
                         int    *opt_flag,
                         double  epsilon,
                         char   *config,
                         char   *nombre,
                         int     imprime, 
                         int     plasma,
                         char   *properties);

  extern int           scf(int     nt,
                           int     elecalfa,
                           int     elecbeta,
                           double  z,
                           int     orbital,
                           double  tol,
                           char   *using_gamma,
                           char   *correlation,
                           char   *propagador,
                           char   **save_dft,
                           double *weight_dft,
                           int     flag_dft,
                           double  mezcla,
                           double *expo,
                           int    *np,
                           int    *mang,
                           int    *ncm,
                           double  gamma_couple,
                           char   *tipo,
                           int     maxiter,
                           double  Rc,
                           char   *bound,
                           char   *espin,
                           double *total_energy,
                           int     print_vectors,
                           double  epsilon,
                           int     imprime,
                           int     plasma,
                           double *cusp_kato,
                           char   *properties);



  int i, j;
  int prue1, prue2, prue3, prue4;

  double x0, x1, x2, energia_0, energia_1, energia_tmp, temp_gam;

  double energy_array[500], gamma_array[500];
  
  char config[4];

  double add_tol;

  double exponente_interno_temporal,
         factor1,
         factor2;
  
  int cont;
  FILE* archivo_opt;
  double cusp_kato;



 if(strcmp(using_gamma,"YES") == 0){//label 4

    if(strcmp(kind_of_cal,"YES") != 0){//label 2
       count_temp = gamma_couple;
      
      gamma_array[50] = count_temp;
      gamma_array[51] = count_temp;
      cont = 0;
      
      gamma_array[cont] = count_temp;
      
      
//mrb     for (i = 0; i < nt; i++)
//mrb            {
//mrb            if (opt_flag[i] != 0)
//mrb             {
//mrb              exponente_interno_temporal = expo[i];
//mrb              factor2 = (double) (mang[i] + np[i])/Rc;
//mrb              factor1 = factor2 - count_temp/((1.f - count_temp)*Rc) - 2.f*step;
//mrb     
//mrb               if(exponente_interno_temporal <= factor1)
//mrb                 {
//mrb                   expo[i] = fabs(factor1) + 4.f;
//mrb                 }
//mrb              }
//mrb            }
      
      printf("checking_gamma %1.10f\n", count_temp);
      
      optimiza_main(nt,
                    opt,
                    elecalfa,
                    elecbeta,
                    z,
                    orbital,
                    tol,
                    using_gamma,
                    correlation,
                    propagador,
                    save_dft,
                    weight_dft,
                    flag_dft,
                    mezcla,
                    expo,
                    np,
                    mang,
                    ncm,
                    count_temp,
                    tipo,
                    maxiter,
                    Rc,
                    bound,
                    espin,
                    step,
                    opt_flag,
                    epsilon,
                    config,
                    nombre,
                    0,
                    plasma,
                    properties);
      
      
      scf(nt,
          elecalfa,
          elecbeta,
          z,
          orbital,
          tol,
          using_gamma,
          correlation,
          propagador,
          save_dft,
          weight_dft,
          flag_dft,
          mezcla,
          expo,
          np,
          mang,
          ncm,
          count_temp,
          tipo,
          maxiter,
          Rc,
          bound,
          espin,
          &energia_0,
          1,
          epsilon,
          1,
          plasma,
         &cusp_kato,
          properties);
      
      
      energy_array[cont] = energia_0;
      *energia = energia_0;
      
      printf("------------------------------------\n");
      printf("------------------------------------\n");
      
     printf("\n ¿opt %1.9f  %15.4f (%15.4f Ryd) Gamma and Energy\n", gamma_array[cont],
                                                                     energy_array[cont],
                                                                     2.f*energy_array[cont]);
      
      printf("------------------------------------\n");
      printf("------------------------------------\n");
      
      }//label 2
      else{ //label 3
          double expo_2[800];
          double expo_save_0[800];
          double expo_save_1[800];
          double expo_save_2[800];
          for (i = 0; i < nt; i++) {
             expo_save_0[i]    = 0.f;
             expo_save_1[i]    = 0.f;
             expo_save_2[i]    = 0.f;
             }

        printf("jgo: optimizo gamma\n");
        cont = steps;
        add_tol = (double) (count_final - count_temp)/steps;
        for(j = 0; j < cont; j++) { //label cont
          energia_0 = 0.f;
          energia_1 = 0.f;
          for (i = 0; i < nt; i++)
             expo_2[i]    = 0.f;
          gamma_array[250 + 2*j] = count_temp;
          x0 = count_temp;

//mrb         for (i = 0; i < nt; i++) {
//mrb            if (opt_flag[i] != 0) {
//mrb              exponente_interno_temporal = expo[i];
//mrb              factor1 = x0/((1.f - x0)*Rc)
//mrb                      + exponente_interno_temporal;
//mrb              factor2 = (double) (mang[i] + np[i])/Rc + 2.f*step;
//mrb               if(factor1 < factor2) 
//mrb                   expo[i] = factor1 + factor2 + 4.f;
//mrb            }
//mrb         }
          printf("checking_gamma1 %1.10f\n", x0);
          optimiza_main(nt,
                        opt,
                           elecalfa,
                           elecbeta,
                           z,
                           orbital,
                           tol,
                           using_gamma,
                           correlation,
                           propagador,
                           save_dft,
                           weight_dft,
                           flag_dft,
                           mezcla,
                           expo,
                           np,
                           mang,
                           ncm,
                           x0,
                           tipo,
                           maxiter,
                           Rc,
                           bound,
                           espin,
                           step,
                           opt_flag,
                           epsilon,
                           config,
                           nombre,
                           0, 
                           plasma,
                           properties);
          prue1 = scf(nt,
                         elecalfa,
                         elecbeta,
                         z,
                         orbital,
                         tol,
                         using_gamma,
                         correlation,
                         propagador,
                         save_dft,
                         weight_dft,
                         flag_dft,
                         mezcla,
                         expo,
                         np,
                         mang,
                         ncm,
                         x0,
                         tipo,
                         maxiter,
                         Rc,
                         bound,
                         espin,
                         &energia_0,
                         1,
                         epsilon,
                         0, 
                         plasma,
                         &cusp_kato,
                         properties);
          count_temp = count_temp + add_tol;
          gamma_array[250 + 2*j + 1] = count_temp;
          x1 = count_temp;
//mrb         for (i = 0; i < nt; i++)
//mrb            {
//mrb              expo_2[i]=expo[i];
//mrb            if (opt_flag[i] != 0)
//mrb             {
//mrb              exponente_interno_temporal = expo_2[i];
//mrb     
//mrb            factor1 = x1/((1.f - x1)*Rc)
//mrb                      + exponente_interno_temporal;
//mrb     
//mrb            factor2 = (double) (mang[i] + np[i])/Rc + 2.f*step;
//mrb     
//mrb               if(factor1 < factor2)
//mrb                 {
//mrb                   expo_2[i] = factor1 + factor2 + 4.f;
//mrb                 }
//mrb              }
//mrb            }
      
      
      printf("checking_gamma2 %1.10f\n", x1);
             optimiza_main(nt,
                           opt,
                           elecalfa,
                           elecbeta,
                           z,
                           orbital,
                           tol,
                           using_gamma,
                           correlation,
                           propagador,
                           save_dft,
                           weight_dft,
                           flag_dft,
                           mezcla,
                           expo_2,
                           np,
                           mang,
                           ncm,
                           x1,
                           tipo,
                           maxiter,
                           Rc,
                           bound,
                           espin,
                           step,
                           opt_flag,
                           epsilon,
                           config,
                           nombre,
                           0, 
                           plasma,
                           properties);
      
       prue2 =       scf(nt,
                         elecalfa,
                         elecbeta,
                         z,
                         orbital,
                         tol,
                         using_gamma,
                         correlation,
                         propagador,
                         save_dft,
                         weight_dft,
                         flag_dft,
                         mezcla,
                         expo_2,
                         np,
                         mang,
                         ncm,
                         x1,
                         tipo,
                         maxiter,
                         Rc,
                         bound,
                         espin,
                         &energia_1,
                         1,
                         epsilon,
                         0,
                         plasma,
                         &cusp_kato,
                         properties);
      
      
      if (prue1 != 0) energia_0 = 1000.f;
      if (prue2 != 0) energia_1 = 1001.f;
      if (energia_0 < -2000.f) energia_0 = 1000.f;
      if (energia_1 < -2000.f) energia_1 = 1001.f;
      
      
             do {
               x2 =  (x0 + x1)/2e00;
      
      if (energia_0 == 1000.f && energia_1 == 1001.f) energia_1 = 1001.f;
      
              if (energia_0 < energia_1) {
                  x1 = x2;
      printf("checking_gamma3 %1.10f\n", x1);

//mrb         for (i = 0; i < nt; i++)
//mrb            {
//mrb            if (opt_flag[i] != 0)
//mrb             {
//mrb              exponente_interno_temporal = expo_2[i];
//mrb
//mrb            factor1 = x1/((1.f - x1)*Rc)
//mrb                      + exponente_interno_temporal;
//mrb
//mrb            factor2 = (double) (mang[i] + np[i])/Rc + 2.f*step;
//mrb
//mrb               if(factor1 < factor2)
//mrb                 {
//mrb                   expo_2[i] = factor1 + factor2 + 4.f;
//mrb                 }
//mrb              }
//mrb            }
      
               optimiza_main(nt,
                             opt,
                             elecalfa,
                             elecbeta,
                             z,
                             orbital,
                             tol,
                             using_gamma,
                             correlation,
                             propagador,
                             save_dft,
                             weight_dft,
                             flag_dft,
                             mezcla,
                             expo_2,
                             np,
                             mang,
                             ncm,
                             x1,
                             tipo,
                             maxiter,
                             Rc,
                             bound,
                             espin,
                             step,
                             opt_flag,
                             epsilon,
                             config,
                             nombre,
                             0, 
                             plasma,
                             properties);
      
       prue3 =        scf(nt,
                          elecalfa,
                          elecbeta,
                          z,
                          orbital,
                          tol,
                          using_gamma,
                          correlation,
                          propagador,
                          save_dft,
                          weight_dft,
                          flag_dft,
                          mezcla,
                          expo_2,
                          np,
                          mang,
                          ncm,
                          x1,
                          tipo,
                          maxiter,
                          Rc,
                          bound,
                          espin,
                          &energia_1,
                          1,
                          epsilon,
                          0,
                          plasma,
                          &cusp_kato,
                          properties);
      
      
      if(prue3 != 0) energia_1 = 1001.f;
      if (energia_1 < -2000.f) energia_1 = 1001.f;
      
               } else {
                       x0 = x2;
                       printf("checking_gamma4 %1.10f\n", x0);
      
//mrb                     for (i = 0; i < nt; i++)
//mrb                        {
//mrb                        if (opt_flag[i] != 0)
//mrb                         {
//mrb                          exponente_interno_temporal = expo[i];
//mrb            
//mrb                        factor1 = x0/((1.f - x0)*Rc)
//mrb                                  + exponente_interno_temporal;
//mrb            
//mrb                        factor2 = (double) (mang[i] + np[i])/Rc + 2.f*step;
//mrb            
//mrb                           if(factor1 < factor2)
//mrb                             {
//mrb                               expo[i] = factor1 + factor2 + 4.f;
//mrb                             }
//mrb                          }
//mrb                        }

      
                       optimiza_main(nt,
                                     opt,
                                     elecalfa,
                                     elecbeta,
                                     z,
                                     orbital,
                                     tol,
                                     using_gamma,
                                     correlation,
                                     propagador,
                                     save_dft,
                                     weight_dft,
                                     flag_dft,
                                     mezcla,
                                     expo,
                                     np,
                                     mang,
                                     ncm,
                                     x0,
                                     tipo,
                                     maxiter,
                                     Rc,
                                     bound,
                                     espin,
                                     step,
                                     opt_flag,
                                     epsilon,
                                     config,
                                     nombre,
                                     0,
                                     plasma,
                                     properties);
              
               prue4 =        scf(nt,
                                  elecalfa,
                                  elecbeta,
                                  z,
                                  orbital,
                                  tol,
                                  using_gamma,
                                  correlation,
                                  propagador,
                                  save_dft,
                                  weight_dft,
                                  flag_dft,
                                  mezcla,
                                  expo,
                                  np,
                                  mang,
                                  ncm,
                                  x0,
                                  tipo,
                                  maxiter,
                                  Rc,
                                  bound,
                                  espin,
                                  &energia_0,
                                  1,
                                  epsilon,
                                  0,  
                                  plasma,
                                  &cusp_kato,
                                  properties);
              
              if(prue4 != 0) energia_0 = 1000.f;
              if (energia_0 < -2000.f) energia_0=1000.f;
      
                }
      
      printf("¡Testing Hartree energy by gamma: %15.4f, %15.4f y gamma:  %2.9f   %2.9f y diff %f\n", energia_0, energia_1, x0, x1, fabs(x1 - x0));
      printf("¡Testing Ryd energy by gamma: %15.4f, %15.4f y gamma:  %2.9f   %2.9f y diff %f\n", 2.f*energia_0, 2.f*energia_1, x0, x1, fabs(x1 - x0));
      
     } while(fabs(x1 - x0) > 1e-4);

      if (j >= 2) {
         for (i = 0; i < nt; i++)
            expo_save_0[i] =  expo_save_1[i];
        }
      
      
      printf("------------------------------------\n");
              if(energia_0 < energia_1){
              gamma_array[j]=x0;
      
      
          printf("Reading initial basis set from %s\n", nombre);
          archivo_opt = fopen(nombre,"r");
          printf(":)  Final exponents 0:\n");
          double temp_dble;
          int temp_int;
          i = -1;
      
          while (fscanf(archivo_opt,"%s %lf %d", config, &temp_dble, &temp_int) != EOF)
          { 
            if (config[1] == 'S') i = i + 1;
            if (config[1] == 'P') i = i + 3;
            if (config[1] == 'D') i = i + 5;
            if (config[1] == 'F') i = i + 7;
            if (config[1] == 'G') i = i + 9;
            if (config[1] == 'H') i = i + 11;
            if (config[1] == 'I') i = i + 13;
            printf(":)  %s  %8.14lf %d \n", config, expo[i], temp_int);
          }
          fclose(archivo_opt);
          printf("** SCF after optimization **\n");

          for (i = 0; i < nt; i++)
              expo_save_2[i] =  expo[i];
          

          energia_0=0e00;
          prue3 = scf(nt,
                      elecalfa,
                      elecbeta,
                      z,
                      orbital,
                      tol,
                      using_gamma,
                      correlation,
                      propagador,
                      save_dft,
                      weight_dft,
                      flag_dft,
                      mezcla,
                      expo,
                      np,
                      mang,
                      ncm,
                      x0,
                      tipo,
                      maxiter,
                      Rc,
                      bound,
                      espin,
                      &energia_0,
                      1,
                      epsilon ,
                      0,
                      plasma,
                      &cusp_kato,
                      properties);
      
                            if(prue3 != 0){ energia_0 = 1000.f;
                                  energy_array[j]=energia_0;
                                          }
                                       else {

                                         energy_array[j]=energia_0;
      printf("\n%d ¿ Gamma %1.14f Energy %15.4f (%15.4f Ryd) Range_of_gamma_from_%1.9f_to_%1.9f\n", j, gamma_array[j],
                                                                                                 energy_array[j],
                                                                                                 2.f*energy_array[j],
                                                                                                 gamma_array[250 + 2*j],
                                                                                                 gamma_array[250 +2*j+1]);
                                            }
      
                                        }
                              else{
                              gamma_array[j]=x1;
      

                           printf("Reading initial basis set from %s\n", nombre);
                           archivo_opt = fopen(nombre,"r");
                           printf(":)  Final exponents 1:\n");
                           double temp_dble;
                           int temp_int;
                           i = -1;
                       
                           while (fscanf(archivo_opt,"%s %lf %d", config, &temp_dble, &temp_int) != EOF)
                           { 
                             if (config[1] == 'S') i = i + 1;
                             if (config[1] == 'P') i = i + 3;
                             if (config[1] == 'D') i = i + 5;
                             if (config[1] == 'F') i = i + 7;
                             if (config[1] == 'G') i = i + 9;
                             if (config[1] == 'H') i = i + 11;
                             if (config[1] == 'I') i = i + 13;
                             printf(":)  %s  %8.14lf %d \n", config, expo_2[i], temp_int);
                           }
                           fclose(archivo_opt);
                           printf("** SCF after optimization **\n");

                           for (i = 0; i < nt; i++)
                               expo_save_2[i] =  expo_2[i];
      
                           energia_1 = 0e00;
                           prue4 =   scf(nt,
                                         elecalfa,
                                         elecbeta,
                                         z,
                                         orbital,
                                         tol,
                                         using_gamma,
                                         correlation,
                                         propagador,
                                         save_dft,
                                         weight_dft,
                                         flag_dft,
                                         mezcla,
                                         expo_2,
                                         np,
                                         mang,
                                         ncm,
                                         x1,
                                         tipo,
                                         maxiter,
                                         Rc,
                                         bound,
                                         espin,
                                         &energia_1,
                                         1,
                                         epsilon,
                                         0,
                                         plasma,
                                         &cusp_kato,
                                         properties);
      
                            if(prue4 != 0) { energia_1 = 1001.f;
                                                 energy_array[j]=energia_1;
                                           }
                                 else {

                                          energy_array[j]=energia_1;
      
                                      printf("\n%d ¿ Gamma %1.14f Energy %15.4f (%15.4f Ryd) Range_of_gamma_from_%1.9f_to_%1.9f\n", j, 
                                                                                                                                    gamma_array[j],
                                                                                                                                    energy_array[j],
                                                                                                                                    2.f*energy_array[j],
                                                                                                                                    gamma_array[250 + 2*j],
                                                                                                                                    gamma_array[250 +2*j+1]);
                                       }
      
                                 } //Ending for for range of gamma

               if(j >= 2) {
                  printf("\n Checking energy\n");
                  int minimum_j;
                  if (energy_array[j-2] > energy_array[j - 1] && energy_array[j - 1] < energy_array[j]){

                      minimum_j = j - 1; 

                           printf("Reading initial basis set from %s\n", nombre);
                           archivo_opt = fopen(nombre,"r");
                           printf("Fopt:  Final exponents 0:\n");
                           double temp_dble;
                           int temp_int;
                           i = -1;

                           while (fscanf(archivo_opt,"%s %lf %d", config, &temp_dble, &temp_int) != EOF)
                           {
                             if (config[1] == 'S') i = i + 1;
                             if (config[1] == 'P') i = i + 3;
                             if (config[1] == 'D') i = i + 5;
                             if (config[1] == 'F') i = i + 7;
                             if (config[1] == 'G') i = i + 9;
                             if (config[1] == 'H') i = i + 11;
                             if (config[1] == 'I') i = i + 13;
                             printf("Fopt:  %s  %8.14lf %d \n", config, expo_save_1[i], temp_int);
                           }
                           fclose(archivo_opt);

                           printf("** SCF after optimization **\n");

                           scf(nt,
                               elecalfa,
                               elecbeta,
                               z,
                               orbital,
                               tol,
                               using_gamma,
                               correlation,
                               propagador,
                               save_dft,
                               weight_dft,
                               flag_dft,
                               mezcla,
                               expo_save_1,
                               np,
                               mang,
                               ncm,
                               gamma_array[minimum_j],
                               tipo,
                               maxiter,
                               Rc,
                               bound,
                               espin,
                               &energia_tmp,
                               1,
                               epsilon,
                               1, 
                               plasma,
                               &cusp_kato,
                               properties);

                      printf("\n0 Fopt: Minimum %d ¿ Gamma %1.14f  Energy %15.4f (%15.4f Ryd) Range_of_gamma_from_%1.9f_to_%1.9f\n", minimum_j,
                                                                                                                                     gamma_array[minimum_j],
                                                                                                                                     energy_array[minimum_j],
                                                                                                                                     2.f*energy_array[minimum_j],
                                                                                                                                     gamma_array[250 + 2*minimum_j],
                                                                                                                                     gamma_array[250 +2*minimum_j+1]);
                      *energia = energia_tmp;
                      break;
                    }


 
                  if (energy_array[j-2] < energy_array[j-1] && energy_array[j-1] < energy_array[j]){

                      minimum_j = j - 2; 

                           printf("Reading initial basis set from %s\n", nombre);
                           archivo_opt = fopen(nombre,"r");
                           printf("Fopt:  Final exponents 1:\n");
                           double temp_dble;
                           int temp_int;
                           i = -1;

                           while (fscanf(archivo_opt,"%s %lf %d", config, &temp_dble, &temp_int) != EOF)
                           {
                             if (config[1] == 'S') i = i + 1;
                             if (config[1] == 'P') i = i + 3;
                             if (config[1] == 'D') i = i + 5;
                             if (config[1] == 'F') i = i + 7;
                             if (config[1] == 'G') i = i + 9;
                             if (config[1] == 'H') i = i + 11;
                             if (config[1] == 'I') i = i + 13;
                             printf("Fopt:  %s  %8.14lf %d \n", config, expo_save_0[i], temp_int);
                           }
                           fclose(archivo_opt);

                           printf("** SCF after optimization **\n");

                           scf(nt,
                               elecalfa, 
                               elecbeta, 
                               z,
                               orbital,  
                               tol,
                               using_gamma,
                               correlation,
                               propagador,
                               save_dft,
                               weight_dft,
                               flag_dft,
                               mezcla,   
                               expo_save_0,
                               np,
                               mang,
                               ncm,
                               gamma_array[minimum_j],
                               tipo,
                               maxiter,
                               Rc,
                               bound,
                               espin,
                               &energia_tmp,
                               1,
                               epsilon,
                               1,
                               plasma,
                               &cusp_kato,
                               properties);


                      printf("\n1 Fopt: Minimum %d ¿ Gamma %1.14f Energy %15.4f (%15.4f Ryd) Range_of_gamma_from_%1.9f_to_%1.9f\n", minimum_j, 
                                                                                                                                    gamma_array[minimum_j],
                                                                                                                                    energy_array[minimum_j],
                                                                                                                                    2.f*energy_array[minimum_j],
                                                                                                                                    gamma_array[250 + 2*minimum_j],
                                                                                                                                    gamma_array[250 +2*minimum_j+1]);
                      *energia = energia_tmp;
                      break;
                   }
                 }
               if (j == cont - 1){

                           printf("Reading initial basis set from %s\n", nombre);
                           archivo_opt = fopen(nombre,"r");
                           printf("Fopt: Final exponents 2:\n");
                           double temp_dble;
                           int temp_int;
                           i = -1;

                           while (fscanf(archivo_opt,"%s %lf %d", config, &temp_dble, &temp_int) != EOF)
                           {
                             if (config[1] == 'S') i = i + 1;
                             if (config[1] == 'P') i = i + 3;
                             if (config[1] == 'D') i = i + 5;
                             if (config[1] == 'F') i = i + 7;
                             if (config[1] == 'G') i = i + 9;
                             if (config[1] == 'H') i = i + 11;
                             if (config[1] == 'I') i = i + 13;
                             printf("Fopt:  %s  %8.14lf %d \n", config, expo_save_2[i], temp_int);
                           }
                           fclose(archivo_opt);

                           printf("** SCF after optimization **\n");

                           scf(nt,
                               elecalfa, 
                               elecbeta, 
                               z,
                               orbital,  
                               tol,
                               using_gamma,
                               correlation,
                               propagador,
                               save_dft,
                               weight_dft,
                               flag_dft,
                               mezcla,   
                               expo_save_2,
                               np,
                               mang,
                               ncm,
                               gamma_array[j],
                               tipo,
                               maxiter,
                               Rc,
                               bound,
                               espin,
                               &energia_tmp,
                               1,
                               epsilon,
                               1,
                               plasma,
                               &cusp_kato,
                               properties);



                      printf("\n2 Fopt: Minimum %d ¿ Gamma %1.14f Energy %15.4f (%15.4f Ryd) Range_of_gamma_from_%1.9f_to_%1.9f\n", j,           
                                                                                                                                    gamma_array[j],
                                                                                                                                    energy_array[j],
                                                                                                                                    2.f*energy_array[j],
                                                                                                                                    gamma_array[250 + 2*j],
                                                                                                                                    gamma_array[250 +2*j+1]);
                      *energia = energia_tmp;
                 }


               for (i = 0; i < nt; i++) {
                   expo_save_1[i] =  expo_save_2[i];
                   expo_save_2[i] =  0.f;
                  }

            }//label cont    
  
      }//label 3 
      
  }//label 4
  else{ //label 5


      count_temp = gamma_couple;
      
      gamma_array[50] = count_temp;
      gamma_array[51] = count_temp;
      cont = 0;
      
      gamma_array[cont] = count_temp;
      
      
//mrb     for (i = 0; i < nt; i++)
//mrb            {
//mrb            if (opt_flag[i] != 0)
//mrb             {
//mrb     
//mrb            factor1 =  expo[i];
//mrb     
//mrb            factor2 = (double) (mang[i] + np[i])/Rc + 2.f*step;
//mrb     
//mrb               if(factor1 < factor2)
//mrb                 {
//mrb                   expo[i] = factor1 + factor2 + 4.f;
//mrb                 }
//mrb              }
//mrb            }
      
      
      printf("MRB_gamma %1.10f\n", count_temp);
      optimiza_main(nt,
                   opt,
                   elecalfa,
                   elecbeta,
                   z,
                   orbital,
                   tol,
                   using_gamma,
                   correlation,
                   propagador,
                   save_dft,
                   weight_dft,
                   flag_dft,
                   mezcla,
                   expo,
                   np,
                   mang,
                   ncm,
                   count_temp,
                   tipo,
                   maxiter,
                   Rc,
                   bound,
                   espin,
                   step,
                   opt_flag,
                   epsilon,
                   config,
                   nombre,
                   0,
                   plasma,
                   properties);
      
      scf(nt,
         elecalfa,
         elecbeta,
         z,
         orbital,
         tol,
         using_gamma,
         correlation,
         propagador,
         save_dft,
         weight_dft,
         flag_dft,
         mezcla,
         expo,
         np,
         mang,
         ncm,
         count_temp,
         tipo,
         maxiter,
         Rc,
         bound, 
         espin,
         &energia_0,
         1,
         epsilon,
         1,
         plasma,
         &cusp_kato,
         properties);
      
      
      energy_array[cont]=energia_0;
      *energia = energia_0;
//jgo      printf("------------------------------------\n");
//jgo      printf("------------------------------------\n");
//jgo      printf("\nCalculation done \n");
      
//jgo      printf("\n:) Energy %15.4f (%15.4f Ryd)\n", energy_array[cont], 2.f*energy_array[cont]);
      
//jgo      printf("------------------------------------\n");
//jgo      printf("------------------------------------\n");

      }//label 5

 return 0;

 }//label 1
