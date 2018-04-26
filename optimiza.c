#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>





void optimiza(int     nt, 
              char   *opt,
	      int     elecalfa, 
	      int     elecbeta,
              double  z, 
	      int     orbital, 
	      double  tol, 
              char   *using_gamma,
              char   *correlation,
              char   *propagador,
              char  **save_dft,
              double *weight_dft,
              int     flag_dft,
	      double  mezcla, 
              char   *basis,
	      double *expo, 
	      int    *np,
              int    *mang, 
	      int    *ncm, 
	      double  gamma_couple,
              char   *tipo, 
	      int     maxiter, 
	      double  Rc, 
              char   *nombre,
	      char   *bound,
              char   *espin, 
	      double  step,
              int    *opt_flag, 
	      double  epsilon, 
	      int     imprime,
              int     plasma,
              char   *properties)
{
 int i, j, k;
 int index, todos, ang_test, deg, test_scf;
 double step_work, energia; 
 double energia_0;
 int bt;
 double cusp_kato;
 int unsigned ident;
 
 int expo_no_correcto,
     checking_expo;

     checking_expo    = 0;
     expo_no_correcto = 0;

 int q, expo_bound, bandera;
 double x, x1, x0, x2, expo_diff[500];
 double expo_opt[800];

 double ener_array[3];
 double maximo, dif_1, dif_2, dif_3, dif_max, tol_opt;
 double dif_1_expo,
        dif_2_expo,
        dif_3_expo,
        dif_max_expo,
        tol_opt_expo;

 for (i = 0; i < 3; i++) 
    ener_array[i] = 0.f;

 double maxdiff;


 extern	int check_expotents_finite(char   *using_gamma,
                                   double  gamma_couple,
	                           double  Rc,
	                           int    *np,
	                           int    *mang,
	                           double *expo,
	                           int     index);

 extern int scf(int      nt, 
                int      elecalfa, 
                int      elecbeta,
                double   z, 
                int      orbital,
                double   tol, 
                char    *using_gamma,
                char    *correlation,
                char    *propagador,
                char   **save_dft,
                double *weight_dft,
                int     flag_dft,
                double   mezcla,
                char    *basis,
                double  *expo, 
                int     *np, 
                int     *mang, 
                int     *ncm, 
                double  gamma_couple,
                char   *tipo, 
                int     maxiter, 
                double  Rc,
                char   *bound, 
                char   *espin, 
                double *total_energy, 
                int     print_vects, 
                double  epsilon, 
                int     imprime,
                int     plasma,
                double *cusp_kato,
                char   *properties);


 extern int checking_linear_dependece(int       nt,
                                      int      *bt,
                                      double   *expo,
                                      int      *mang,
                                      char     *nombre);



 int optimiza_new(int     nt,
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
                  char   *nombre,
                  char   *bound,
                  char   *espin,
                  double  step,
                  int    *opt_flag,
                  double  epsilon,
                  int     imprime,
                  int     plasma);


  if(strcmp(opt,"optfull") == 0)
     optimiza_new(nt,
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
                  gamma_couple,
                  tipo,
                  maxiter,
                  Rc,
                  nombre,
                  bound,
                  espin,
                  step,
                  opt_flag,
                  epsilon,
                  imprime,
                  plasma);

 int unsigned todos_antes;
 FILE* file_basis; 
 FILE* archivo_opt;
 int  unsigned file;
 int  unsigned file_ii;
 char nombre_ii[80];
 char config[4];
 file    = 0;
 file_ii = 1;

 double temp_dble;
 int temp_int;

 energia_0 = 1e9;  // energy_0 asignation
    
 sprintf(nombre_ii, "%s_restart", nombre);
    do {
       step_work    = step;
       while (step_work >= (double)0.0001) {
        todos_antes  = 0;
        todos        = 0;
        do {
            expo_bound = 0;
            bandera    = 0;
            deg        = 1;
            ang_test   = mang[todos];
            switch(ang_test) {
              case(1) : deg = 3; break;
              case(2) : deg = 5; break;
              case(3) : deg = 7; break;
              case(4) : deg = 9; break;
              case(5) : deg = 11; break;
              case(6) : deg = 13; break;
              case(7) : deg = 15; break;
              default : deg =  1; break;
            }
               if (opt_flag[todos] == 0) {
                 bandera = 1;
                   for (i = 1; i <= deg; i++) {
                     index = todos + i - 1;
                     printf("0 Expo %d final value = %f\n", index,
                                                            expo[index]);
                   }
                   todos = index + 1;
               } 
               else {  //else1
                  test_scf = ((int) 0);      // test_scf asignation
                  x  = expo[todos];  // se asigna cada exponente a x
 //***************************************************************************************
                  x1 = x;            // valor de referencia  x1, el exponente permanece igual 
                  for(i = 1; i <= deg; i++) {
                     index = todos + i - 1;
                     expo[index] = x1;
                  }
                  test_scf = scf(nt,
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
                                 basis,
                                 expo,
                                 np,
                                 mang,
                                 ncm,
                                 gamma_couple,
                                 tipo,
                                 maxiter,
                                 Rc,
                                 bound,
                                 espin,
                                 &energia,
                                 0,
                                 epsilon,
                                 0,
                                 plasma,
                                 &cusp_kato,
                                 properties);

                  if(test_scf == 0)
                     ener_array[1] = energia;
                  else
                     ener_array[1] = 1e11; 

 //***************************************************************************************                    
                  x0 = x + step_work;   //  x0, al exponente se le agrega el tamaño del paso             
                  for(i = 1; i <= deg; i++) {
                     index = todos + i - 1;
                     expo[index] = x0;
                  }              
                  test_scf = scf(nt,
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
                                 basis,
                                 expo,
                                 np,
                                 mang,
                                 ncm,
                                 gamma_couple,
                                 tipo,
                                 maxiter,
                                 Rc,
                                 bound,
                                 espin,
                                 &energia,
                                 0,
                                 epsilon,
                                 0,
                                 plasma,
                                 &cusp_kato,
                                 properties);

                  if(test_scf == 0)
                     ener_array[0] = energia;
                  else 
                     ener_array[0] = 1e12;

 //***************************************************************************************                    
                  x2 = x - step_work;   //  x2, al exponente se le resta el tamaño del paso
//                  if(x2 > 0.1) {  //if x2
                     for(i = 1; i <= deg; i++) {
                        index = todos + i - 1;
                        expo[index] = x2;
                     }
//                     if(strcmp(basis, "STOs") || strcmp(basis, "stos") == 0 ) {   
//                        if(strcmp(bound, "finite") || strcmp(bound, "dielectricc") == 0 || strcmp(bound,"polarization") == 0) {
//                           checking_expo = check_expotents_finite(using_gamma,
//                                                                  gamma_couple,
//                                                                  Rc,
//                                                                  np,
//                                                                  mang,
//                                                                  expo,
//                                                                  todos);
//                        }
//                        if(checking_expo == 0) {
                           test_scf = scf(nt,
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
                                          basis,
                                          expo,
                                          np,
                                          mang,
                                          ncm,
                                          gamma_couple,
                                          tipo,
                                          maxiter,
                                          Rc,
                                          bound,
                                          espin,
                                          &energia,
                                          0,
                                          epsilon,
                                          0,
                                          plasma,
                                          &cusp_kato,
                                          properties);
//                        }
//                        else {
//                           test_scf = 1;
//                        }
//                     } //  if basis
//                  }    // if x2
//                  else
//                     test_scf = 1;

                  if(test_scf == 0)
                     ener_array[2] = energia;
                  else 
                     ener_array[2] = 1e12;

 //***************************************************************************************                   
 // We have the next cases... 
                  if(ener_array[0] < ener_array[1]) {      // first case: E0 < E1
                     for(i = 1; i <= deg; i++) {
                        index = todos + i - 1;
                        expo[index] = x0;
                     } 
                     energia = ener_array[0];
                     printf("\n case_1 E0 < E1     %4.10lf  %4.25lf \n", ener_array[0], x0);
                  } 
                  else 
                     if(ener_array[2] < ener_array[1]) {   // second case: E2 < E1
                        for(i = 1; i <= deg; i++) {
                           index = todos + i - 1;
                           expo[index] = x2;
                        } 
                        energia = ener_array[2];
                        printf("\n case_2 E2 < E1     %4.10lf  %4.25lf \n", ener_array[2], x2);
                     } 
                     else{                                 // third case: E1 is still smaller than E2 and E0
                        for(i = 1; i <= deg; i++) {
                           index = todos + i - 1;
                           expo[index] = x1;
                        }
                        energia = ener_array[1];
                        printf("\n case_3 E1 < E2, E0 %4.10lf  %4.25lf \n", ener_array[1], x1);
                     }
                     todos_antes = todos;
                     todos        = index + 1;
                 }//else1
        } while (todos < nt);

        step_work = step_work/10.f;
       }
    
       if(energia < energia_0) {
          energia_0 = energia;
          ident = 0;
          for(i = 0; i < nt; i++)
             expo_opt[i] = expo[i];
          printf("\nAQUIVOY_I\n");

          file++;
          if(file == 2*file_ii) {
             file_ii += file_ii;
             file_basis = fopen(nombre_ii,"w");
             archivo_opt = fopen(nombre,"r");
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
                fprintf(file_basis,"%s  %8.25lf %d \n", config, expo_opt[i], temp_int);
             }
             fclose(archivo_opt);
             fclose(file_basis);
          }
       }
       else 
         if(energia >= energia_0) {
            printf("\nAQUIVOY\n");
            ident = 1;
            for (i = 0; i < nt; i++){
//                expo[i] = expo_opt[i]; //aquí únicamente estoy pasando un cero (original)
                expo_opt[i] = expo[i]; 
//                printf(" expo %lf \n", expo[i]);         //mike :(
//                printf(" expo_opt %lf \n", expo_opt[i]); //mike :(
            }
         }

      } while (ident == 0);

      remove(nombre_ii);

            if(strcmp(bound,"free") == 0)
             sprintf(nombre_ii, "basis_opt_%s", bound);
            else
             sprintf(nombre_ii, "basis_opt_%s_Rc_%3.4lf", bound, Rc);

            file_basis = fopen(nombre_ii,"w");
            archivo_opt = fopen(nombre,"r");
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
               fprintf(file_basis,"%s  %8.25lf %d \n", config, expo_opt[i], temp_int);
             }
          fclose(archivo_opt);
          fclose(file_basis);

 
  

// // do {
//   maxdiff = 0.f;
//   q = q + 1;
//   printf("Step %d in the optimization process\n", q);
// 
//   maximo = (double)1.e10;
//   expo_bound = 0;
//   // tol_opt = (double)1.e-8;
//   tol_opt = (double)1.e-7;
//   tol_opt_expo = (double)1.e-12;
// //     step_work = step;
//   step_work = 0.01;
// 
//   while (step_work >= (double)0.0001) {
// 
//    printf("Step for the optimization = %f\n", step_work);
//    todos = 0;
// 
//    do {
// 
//        expo_bound = 0;
//        bandera    = 0;
//        deg        = 1;
//        ang_test   = mang[todos];
//        switch(ang_test) {
//          case(1) : deg = 3; break;
//          case(2) : deg = 5; break;
//          case(3) : deg = 7; break;
//          case(4) : deg = 9; break;
//          case(5) : deg = 11; break;
//          case(6) : deg = 13; break;
//          case(7) : deg = 15; break;
//          default : deg =  1; break;
//        }
// 
// 
//        if (opt_flag[todos] == 0) {
//          bandera = 1;
//            for (i = 1; i <= deg; i++) {
//              index = todos + i - 1;
//              printf("0 Expo %d final value = %f\n", index,
//                                                     expo[index]);
//            }
//            todos = index + 1;
//        } else {//label 1
//                x  = expo[todos];
//                x1 = x;
//                for (i = 1; i <= deg; i++) {
//                  index       = todos + i - 1;
//                  expo[index] = x1;
//                }
// 
//                if (strcmp(bound,"finite") == 0) {
//                       expo_no_correcto =  check_expotents_finite(using_gamma,
//                                                                  gamma_couple,
//                                                                  Rc,
//                                                                  np,
//                                                                  mang,
//                                                                  expo,
//                                                                  todos);
//                } else
//                      expo_no_correcto = 0;
// 
//                checking_expo = checking_linear_dependece(nt,
//                                                          &bt,
//                                                          expo,
//                                                          mang,
//                                                          nombre);
// 
// 
//                if (expo_no_correcto == 1 || checking_expo == 1) {
//                  test_scf = 1;
//                  energia  = maximo;
//                } else
//                      test_scf = scf(nt,
//                                     elecalfa,
//                                     elecbeta,
//                                     z,
//                                     orbital,
//                                     tol,
//                                     using_gamma,
//                                     correlation,
//                                     propagador,
//                                     mezcla,
//                                     expo,
//                                     np,
//                                     mang,
//                                     ncm,
//                                     gamma_couple,
//                                     tipo,
//                                     maxiter,
//                                     Rc,
//                                     bound,
//                                     espin,
//                                     &energia,
//                                     0,
//                                     epsilon,
//                                     0,
//                                     plasma,
// 
//                if (test_scf == 0) {
//                  ener_array[1] = energia;
//                  printf("1 Expo %d, value=%f, energ = %f\n", todos,
//                                                              x1,
//                                                              energia);
//                } else {
//                        printf("Convergence problems \n");
//                        ener_array[1] = maximo;
//                        for (i = 1; i <= deg; i++) {
//                          index       = todos + i - 1;
//                          expo[index] = x;
//                        }
//                  }
//                x0 = x - step_work;
//                if (x0 <= 0.f) {
//                  x0            = x;
//                  checking_expo = 1;
//                } else {
//                        for (i = 1; i <= deg; i++) {
//                           index       = todos + i - 1;
//                           expo[index] = x0;
//                        }
// 
//                        if (strcmp(bound,"finite") == 0) {
//                               expo_no_correcto =  check_expotents_finite(using_gamma,
//                                                                          gamma_couple,
//                                                                          Rc,
//                                                                          np,
//                                                                          mang,
//                                                                          expo,
//                                                                          todos);
//                        } else
//                              expo_no_correcto = 0;
// 
//                        checking_expo = checking_linear_dependece(nt,
//                                                            &bt,
//                                                            expo,
//                                                            mang,
//                                                            nombre);
//                  }
// 
//                if (expo_no_correcto == 1 || checking_expo == 1) {
//                  test_scf = 1;
//                  x0       = x;
//                  energia  = maximo;
//                } else
//                      test_scf = scf(nt,
//                                     elecalfa,
//                                     elecbeta,
//                                     z,
//                                     orbital,
//                                     tol,
//                                     using_gamma,
//                                     correlation,
//                                     propagador,
//                                     mezcla,
//                                     expo,
//                                     np,
//                                     mang,
//                                     ncm,
//                                     gamma_couple,
//                                     tipo,
//                                     maxiter,
//                                     Rc,
//                                     bound,
//                                     espin,
//                                     &energia,
//                                     0,
//                                     epsilon,
//                                     0,
//                                     plasma);
// 
//                if (test_scf == 0) {
//                  ener_array[0] = energia;
//                  printf("2 Expo %d, value=%f, energ = %f\n", todos,
//                                                              x0,
//                                                              energia);
//                } else {
//                        printf("Convergence problems \n");
//                        ener_array[0] = maximo;
//                        for (i = 1; i <= deg; i++) {
//                          index       = todos + i - 1;
//                          expo[index] = x;
//                        }
//                  }
// 
//                x2 = x + step_work;
//                for (i = 1; i <= deg; i++) {
//                  index       = todos + i - 1;
//                  expo[index] = x2;
//                }
//                if (strcmp(bound,"finite") == 0) {
//                       expo_no_correcto =  check_expotents_finite(using_gamma,
//                                                                  gamma_couple,
//                                                                  Rc,
//                                                                  np,
//                                                                  mang,
//                                                                  expo,
//                                                                  todos);
//                } else
//                      expo_no_correcto = 0;
// 
//                checking_expo = checking_linear_dependece(nt,
//                                                          &bt,
//                                                          expo,
//                                                          mang,
//                                                          nombre);
// 
// 
//                if (expo_no_correcto == 1 || checking_expo == 1) {
//                  test_scf = 1;
//                  x2       = x;
//                  energia  = maximo;
//                } else
//                      test_scf = scf(nt,
//                                     elecalfa,
//                                     elecbeta,
//                                     z,
//                                     orbital,
//                                     tol,
//                                     using_gamma,
//                                     correlation,
//                                     propagador,
//                                     mezcla,
//                                     expo,
//                                     np,
//                                     mang,
//                                     ncm,
//                                     gamma_couple,
//                                     tipo,
//                                     maxiter,
//                                     Rc,
//                                     bound,
//                                     espin,
//                                     &energia,
//                                     0,
//                                     epsilon,
//                                     0,
//                                     plasma);
// 
//                if (test_scf == 0) {
//                  ener_array[2] = energia;
//                  printf("3 Expo %d, value=%f, energ = %f\n", todos,
//                                                              x2,
//                                                              energia);
//                } else {
//                        printf("Convergence problems \n");
//                        ener_array[2] = maximo;
//                        for (i = 1; i <= deg; i++) {
//                          index       = todos + i - 1;
//                          expo[index] = x;
//                    }
//                  }
// 
//                dif_1        = fabs(ener_array[0] - ener_array[1]);
//                dif_2        = fabs(ener_array[0] - ener_array[2]);
//                dif_3        = fabs(ener_array[2] - ener_array[1]);
//                dif_max      = (dif_1 > dif_2 ? dif_1 : dif_2);
//                dif_max      = (dif_max > dif_3 ? dif_max : dif_3);
// 
//                dif_1_expo   = fabs(x0 - x1);
//                dif_2_expo   = fabs(x0 - x2);
//                dif_3_expo   = fabs(x2 - x1);
//                dif_max_expo = (dif_1_expo > dif_2_expo ? dif_1_expo : dif_2_expo);
//                dif_max_expo = (dif_max_expo > dif_3_expo ? dif_max_expo : dif_3_expo);
// 
//                printf("x0 = %f, e0 = %18.10f\n", x0,
//                                                  ener_array[0]);
//                printf("x1 = %f, e1 = %18.10f\n", x1,
//                                                  ener_array[1]);
//                printf("x2 = %f, e2 = %18.10f\n", x2,
//                                                  ener_array[2]);
//                printf("dif max = %14.10f\n", dif_max);
//                printf("dif max expo = %14.10f\n", dif_max_expo);
//                if (dif_max <= tol_opt || dif_max_expo <= tol_opt_expo) {
//                    bandera = 1;
//                    for (i = 1; i <= deg; i++) {
//                      index       = todos + i - 1;
//                      expo[index] = x1;
//                      printf("4 Expo %d final value = %f\n", index,
//                                                             expo[index]);
//                    }
//                    todos = index + 1;
//                } else {//label 2
//                        do {//label 3
//                            if (ener_array[1] < ener_array[0]) {//labelcomp1
//                              if (ener_array[1] < ener_array[2]) {
//                                 bandera = 1;
//                                 for (i = 1; i <= deg; i++) {
//                                   index       = todos + i - 1;
//                                   expo[index] = x1;
//                                   printf("5 Expo %d final value = %f\n", index,
//                                                                          expo[index]);
//                                 }
//                                 todos = index + 1;
//                              } else {//labelcomp2
//                                      x0            = x1;
//                                      ener_array[0] = ener_array[1];
//                                      x1            = x2;
//                                      ener_array[1] = ener_array[2];
//                                      x2            = x1 + step_work;
//                                      for (i = 1; i <= deg; i++) {
//                                         index       = todos + i - 1;
//                                         expo[index] = x2;
//                                      }
//                                      if (strcmp(bound,"finite") == 0) {
//                                             expo_no_correcto =  check_expotents_finite(using_gamma,
//                                                                                        gamma_couple,
//                                                                                        Rc,
//                                                                                        np,
//                                                                                        mang,
//                                                                                        expo,
//                                                                                        todos);
//                                      } else
//                                            expo_no_correcto = 0;
// 
//                                      checking_expo = checking_linear_dependece(nt,
//                                                                                &bt,
//                                                                                expo,
//                                                                                mang,
//                                                                                nombre);
// 
// 
//                                      if (expo_no_correcto == 1 || checking_expo == 1) {
//                                        test_scf = 1;
//                                        energia  = maximo;
//                                      } else
//                                            test_scf = scf(nt,
//                                                           elecalfa,
//                                                           elecbeta,
//                                                           z,
//                                                           orbital,
//                                                           tol,
//                                                           using_gamma,
//                                                           correlation,
//                                                           propagador,
//                                                           mezcla,
//                                                           expo,
//                                                           np,
//                                                           mang,
//                                                           ncm,
//                                                           gamma_couple,
//                                                           tipo,
//                                                           maxiter,
//                                                           Rc,
//                                                           bound,
//                                                           espin,
//                                                           &energia,
//                                                           0,
//                                                           epsilon,
//                                                           0,
//                                                           plasma);
//                                      if (test_scf == 0) {
//                                        ener_array[2] = energia;
//                                        printf("6 Expo %d, value=%f, energ = %f\n",todos,
//                                                                                   x2,
//                                                                                   energia);
//                                      } else {
//                                              printf("Convergence problems \n");
//                                              ener_array[2] = maximo;
//                                      }
//                              }//labelcomp2
//             /*labelcomp1*/ } else {//labelcomp7
//                                    if (ener_array[1] < ener_array[2]) {//labelcomp3
//                                      x2            = x1;
//                                      ener_array[2] = ener_array[1];
//                                      x1            = x0;
//                                      ener_array[1] = ener_array[0];
//                                      x0            = x1 - step_work;
//                                      if (x0 <= 0.f) {
//                                        x0         = 0.01f;
//                                        expo_bound = 1;
//                                      }
//                                      for (i = 1; i <= deg; i++) {
//                                        index       = todos + i - 1;
//                                        expo[index] = x0;
//                                      }
//                                      if (strcmp(bound,"finite") == 0) {
//                                             expo_no_correcto =  check_expotents_finite(using_gamma,
//                                                                                        gamma_couple,
//                                                                                        Rc,
//                                                                                        np,
//                                                                                        mang,
//                                                                                        expo,
//                                                                                        todos);
//                                      } else
//                                            expo_no_correcto = 0;
// 
//                                      checking_expo = checking_linear_dependece(nt,
//                                                                                &bt,
//                                                                                expo,
//                                                                                mang,
//                                                                                nombre);
// 
// 
//                                      if(expo_no_correcto == 1 || checking_expo == 1) {
//                                      test_scf = 1;
//                                      energia  = maximo;
//                                      } else
//                                            test_scf = scf(nt,
//                                                           elecalfa,
//                                                           elecbeta,
//                                                           z,
//                                                           orbital,
//                                                           tol,
//                                                           using_gamma,
//                                                           correlation,
//                                                           propagador,
//                                                           mezcla,
//                                                           expo,
//                                                           np,
//                                                           mang,
//                                                           ncm,
//                                                           gamma_couple,
//                                                           tipo,
//                                                           maxiter,
//                                                           Rc,
//                                                           bound,
//                                                           espin,
//                                                           &energia,
//                                                           0,
//                                                           epsilon,
//                                                           0,
//                                                           plasma);
//                                      if (test_scf == 0) {
//                                        ener_array[0] = energia;
//                                        printf("7 Expo %d, value=%f, energ = %f\n", todos,
//                                                                                    x0,
//                                                                                    energia);
//                                      } else {
//                                          printf("Convergence problems \n");
//                                          ener_array[0] = maximo;
//                                        }
//                    /*labelcomp3*/  } else {//labelcomp6
//                                            if (ener_array[2] < ener_array[0]) {//labelcomp4
//                                              x0            = x1;
//                                              ener_array[0] = ener_array[1];
//                                              x1            = x2;
//                                              ener_array[1] = ener_array[2];
//                                              x2            = x1 + step_work;
//                                              for (i = 1; i <= deg; i++) {
//                                                index       = todos + i - 1;
//                                                expo[index] = x0;
//                                              }
//                                              if (strcmp(bound,"finite") == 0) {
//                                                     expo_no_correcto =  check_expotents_finite(using_gamma,
//                                                                                                gamma_couple,
//                                                                                                Rc,
//                                                                                                np,
//                                                                                                mang,
//                                                                                                expo,
//                                                                                                todos);
//                                              } else
//                                                    expo_no_correcto = 0;
// 
//                                              checking_expo = checking_linear_dependece(nt,
//                                                                                        &bt,
//                                                                                        expo,
//                                                                                        mang,
//                                                                                        nombre);
// 
// 
//                                              if (expo_no_correcto == 1 || checking_expo == 1) {
//                                                test_scf = 1;
//                                                energia  = maximo;
//                                              } else
//                                                    test_scf = scf(nt,
//                                                                   elecalfa,
//                                                                   elecbeta,
//                                                                   z,
//                                                                   orbital,
//                                                                   tol,
//                                                                   using_gamma,
//                                                                   correlation,
//                                                                   propagador,
//                                                                   mezcla,
//                                                                   expo,
//                                                                   np,
//                                                                   mang,
//                                                                   ncm,
//                                                                   gamma_couple,
//                                                                   tipo,
//                                                                   maxiter,
//                                                                   Rc,
//                                                                   bound,
//                                                                   espin,
//                                                                   &energia,
//                                                                   0,
//                                                                   epsilon,
//                                                                   0,
//                                                                   plasma);
//                                              if (test_scf == 0) {
//                                                ener_array[2] = energia;
//                                                printf("8 Expo %d, value=%f, energ = %f\n", todos,
//                                                                                            x2,
//                                                                                            energia);
//                                              } else {
//                                                      printf("Convergence problems \n");
//                                                      ener_array[2] = maximo;
//                                                }
//                        /*labelcomp4*/      } else {//labelcomp5
//                                                    x2            = x1;
//                                                    ener_array[2] = ener_array[1];
//                                                    x1            = x0;
//                                                    ener_array[1] = ener_array[0];
//                                                    x0            = x1 - step_work;
//                                                    if (x0 <= 0.f) {
//                                                       x0         = 0.01f;
//                                                       expo_bound = 1;
//                                                    }
//                                                    for (i = 1; i <= deg; i++) {
//                                                       index       = todos + i - 1;
//                                                       expo[index] = x0;
//                                                    }
//                                                    if (strcmp(bound,"finite") == 0) {
//                                                           expo_no_correcto =  check_expotents_finite(using_gamma,
//                                                                                                      gamma_couple,
//                                                                                                      Rc,
//                                                                                                      np,
//                                                                                                      mang,
//                                                                                                      expo,
//                                                                                                      todos);
//                                                    } else
//                                                          expo_no_correcto = 0;
// 
//                                                    checking_expo = checking_linear_dependece(nt,
//                                                                                              &bt,
//                                                                                              expo,
//                                                                                              mang,
//                                                                                              nombre);
// 
//                                                    if (expo_no_correcto == 1 || checking_expo == 1) {
//                                                      test_scf = 1;
//                                                      energia  = maximo;
//                                                    } else
//                                                          test_scf = scf(nt,
//                                                                         elecalfa,
//                                                                         elecbeta,
//                                                                         z,
//                                                                         orbital,
//                                                                         tol,
//                                                                         using_gamma,
//                                                                         correlation,
//                                                                         propagador,
//                                                                         mezcla,
//                                                                         expo,
//                                                                         np,
//                                                                         mang,
//                                                                         ncm,
//                                                                         gamma_couple,
//                                                                         tipo,
//                                                                         maxiter,
//                                                                         Rc,
//                                                                         bound,
//                                                                         espin,
//                                                                         &energia,
//                                                                         0,
//                                                                         epsilon,
//                                                                         0,
//                                                                         plasma);
//                                                        if (test_scf == 0) {
//                                                          ener_array[0] = energia;
//                                                          printf("9 Expo %d, value=%f, energ = %f\n", todos,
//                                                                                                      x0,
//                                                                                                      energia);
//                                                        } else {
//                                                                printf("Convergence problems \n");
//                                                                ener_array[0] = maximo;
//                                                        }
// 
//                                              }//labelcomp5
//                                      }//labelcomp6
//                              }//labelcomp7
//                            dif_1        = fabs(ener_array[0] - ener_array[1]);
//                            dif_2        = fabs(ener_array[0] - ener_array[2]);
//                            dif_3        = fabs(ener_array[2] - ener_array[1]);
//                            dif_max      = (dif_1 > dif_2 ? dif_1 : dif_2);
//                            dif_max      = (dif_max > dif_3 ? dif_max : dif_3);
// 
//                            dif_1_expo   = fabs(x0 - x1);
//                            dif_2_expo   = fabs(x0 - x2);
//                            dif_3_expo   = fabs(x2 - x1);
//                            dif_max_expo = (dif_1_expo > dif_2_expo ? dif_1_expo : dif_2_expo);
//                            dif_max_expo = (dif_max_expo > dif_3_expo ? dif_max_expo : dif_3_expo);
// 
//                            if (expo_bound == 1) {
//                                dif_max         = 0.f;
//                                opt_flag[todos] = 0;
//                            }
//                            if (dif_max <= tol_opt || dif_max_expo <= tol_opt_expo) {
//                               bandera = 1;
//                               for (i = 1; i <= deg; i++) {
//                                  index       = todos + i - 1;
//                                  expo[index] = x1;
//                                  printf("10 Expo %d final value = %f\n", index, expo[index]);
//                               }
//                            todos = index + 1;
//                            }
//                            printf("dif max = %14.10f\n", dif_max);
// 
//                        } while (bandera == 0);//label 3
//                  }//label 2
//              }//label 1
//      } while(todos < nt);
// 
//      step_work = step_work/(double)10.;
//     }
//     printf("After partial optimization\n");
//     for (i = 0; i < nt; i++)
//        printf("@ Partial opt Expo %d : %20.12f\n", i, expo[i]);
// 
//   for (j = 0; j < nt; j++)
//      {
//       if (fabs(expo[j] - expo_diff[j]) > maxdiff)
//         maxdiff = fabs(expo[j] - expo_diff[j]);
//       expo_diff[j] = expo[j];
//   }
// // } while (maxdiff > 1.e-5);


 }
