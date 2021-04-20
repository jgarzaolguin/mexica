#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


void open_range(int     n, 
                int     b, 
                double *p) {

       int  i, j;
       double adicion;

       double rango_menor;
       double rango_mayor;
       double suma;
       int    pasos;
       rango_menor = 0.1f;
       rango_mayor = 40.f;
       pasos       = n;

       adicion     = (double) (rango_mayor - rango_menor)/pasos;
       suma        = rango_menor;


       for (i = 0; i < b; i++){
          for (j = 0; j < n; j++){

           p[j + i*n] = suma;
         // if (i == 0)
         // p[j + i*n] = suma;
         // else
         // p[j + i*n] = pow(suma, (double) i);
           
           printf("%8.4lf  ", p[j + i*(n + 1)]); 
           suma = suma + adicion + (double) 0.15*(i + 1);
       }
          suma = rango_menor;
          //printf("\n");
      }
          //printf("\n");

     }
                                                                       


void optimiza(int     nt, 
	      int     elecalfa, 
	      int     elecbeta,
              double  z, 
	      int     orbital, 
	      double  tol, 
              char   *using_gamma,
              char   *correlation,
              char   *propagador,
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
	      double  gamma_nicp, 
	      int     imprime,
              int plasma)
{
 int i, j, k, sizeint, sizedouble;
 int bandera, index, blanco, todos, ang_test, deg, test_scf, expo_bound;
 double x, x0, x1, x2, step_work;
 int bt;
 
 double energia, ener_array[3];
 double maximo, dif_1, dif_2, dif_3, dif_max, tol_opt;
 double dif_1_expo,
        dif_2_expo,
        dif_3_expo,
        dif_max_expo,
        tol_opt_expo;

// double exponente_interno_temporal,
//	factor1,
//	factor2;

 int expo_no_correcto,
     checking_expo;

     checking_expo    = 0;
     expo_no_correcto = 0;

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
                double   mezcla,
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
                double  gamma_nicp, 
                int     imprime,
                int plasma);

 extern int checking_linear_dependece(int       nt,
                                      int      *bt,
                                      double   *expo,
                                      int      *mang,
                                      char     *nombre);

 void open_range(int     n,
                 int     b,
                 double *p);


 int a, 
     total;
 double *p;

 double gamma_nicp; // mike, nicp --> nonideal classical plasma

 a = 500;

 checking_linear_dependece(nt,
                           &bt,
                           expo,
                           mang,
                           nombre);
 total = a*bt;
 
 p     = (double*)malloc(total*sizeof(double));
 if (p == NULL) {
   free(p);
   p = 0;
   printf("Error en arreglo p, insuficiente RAM\n");
   exit(1);
 }
 
 open_range(a, 
            bt, 
            p);

//mrb for (i = 0; i < bt; i++) {
//mrb    for (j = 0; j < a; j++){
//mrb       printf("%8.4lf(%d)  ", p[j + i*a], j + i*a);
//mrb       }
//mrb     printf("\n");
//mrb }


 srand(time(NULL));
 double ppast, diff;
 int s, q;
 double p_temp[bt];
 double expo_opt[500];
 double expo_diff[500];
 double energia_0;
 int    pa;
 
 for (j = 0; j < nt; j++) {
    expo_diff[j] = 0.f;
    expo_opt[j] = 0.f;
 }

 energia_0 = 1.e9;

 for (q = 0; q < 10; q++) {//labelmrb2
    pa = 100*(1 + q);
    for (j = 0; j < 200; j++) {//labelmrb1
       for (i = 0; i < bt; i++) {
    
          if (i < 1) {
            s = rand() % pa;
            p_temp[i] = p[s];
            ppast = p[s];
          } else {
    
                 s = rand() % (pa*(i + 1) - pa*i) + pa*i;
    
                 diff = fabs(ppast - p[s]);
                 if (diff < 0.13) {
                   do {
                        s    = rand() % (pa*(i + 1) - pa*i) + pa*i;
                        diff = fabs(ppast - p[s]);
                   } while (ppast == p[s] || diff < 0.13);
                 }
                   p_temp[i] = p[s];
                   ppast = p[s];
            }
    //mrb      printf("\n");
    //mrb      printf("%8.4lf(%d)  ", p_temp[i], i);
       }
    
    todos = 0;
    k     = 0;
    do {
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
    
        for (i = 1; i <= deg; i++) { 
        index = todos + i - 1;
        expo[index] = p_temp[k];
        printf("\n e %lf t  %lf\n", expo[index], p_temp[k]);
        }
        k++;
        todos = index + 1;
    }while(todos < nt); 
    
    test_scf = scf(nt,
                   elecalfa,
                   elecbeta,
                   z,
                   orbital,
                   tol,
                   using_gamma,
                   correlation,
                   propagador,
                   mezcla,
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
                   gamma_nicp,
                   0,
                   plasma);
    
     if (test_scf != 0)
       energia = 1.e10;
    
     if (energia < energia_0){
          energia_0 = energia;
          for (i = 0; i < nt; i++)
             expo_opt[i] = expo[i];
          printf("\n HOLA2 energia %lf", energia);
       }
    
        printf("\n");
    }//labelmrb1
 }//labelmrb2

 printf("\nHOLA_1\n");
 for (i = 0; i < nt; i++)
    expo[i] = expo_opt[i];

 double maxdiff;
 do {
     maxdiff = 0.f;
     q = q + 1;
     printf("Step %d in the optimization process\n", q);

     maximo = (double)1.e10;
     expo_bound = 0;
     // tol_opt = (double)1.e-8;
     tol_opt = (double)1.e-7;
     tol_opt_expo = (double)1.e-12;
     step_work = step;
     
     while (step_work >= (double)0.0001) {

      printf("Step for the optimization = %f\n", step_work);
      todos = 0;
     
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
          } else {//label 1
     
                  /*Testing_each_other_expo_test_linear_dependence x0, x1 and x2*/
                  
                //mrb x  = expo[todos];
                //mrb x0 = x - step_work;
                //mrb for (i = 1; i <= deg; i++) {
                //mrb   index       = todos + i - 1;
                //mrb   expo[index] = x0;
                //mrb }
                //mrb checking_expo = checking_linear_dependece(nt,
                //mrb                                           &bt,
                //mrb                                           expo,
                //mrb                                           mang,
                //mrb                                           nombre);
                //mrb if (checking_expo == 0) {
                //mrb   x1 = x;
                //mrb   for (i = 1; i <= deg; i++) {
                //mrb     index       = todos + i - 1;
                //mrb     expo[index] = x1;
                //mrb   }
                //mrb   checking_expo = checking_linear_dependece(nt,
                //mrb                                             &bt,
                //mrb                                             expo,
                //mrb                                             mang,
                //mrb                                             nombre);
                //mrb   if (checking_expo == 0) {
                //mrb     x2 = x + step_work;
                //mrb     for (i = 1; i <= deg; i++) {
                //mrb       index       = todos + i - 1;
                //mrb       expo[index] = x2;
                //mrb     }
                //mrb     checking_expo = checking_linear_dependece(nt,
                //mrb                                               &bt,
                //mrb                                               expo,
                //mrb                                               mang,
                //mrb                                               nombre);
                //mrb   }
                //mrb }
                //mrb if (checking_expo == 1) {
                //mrb   for (i = 1; i <= deg; i++) {
                //mrb      index       = todos + i - 1;
                //mrb      expo[index] = x + 0.15;
                //mrb   }                
                //mrb }
                //mrb /*Testing_each_other_expo_test_finite_restriction x0, x1 and x2*/
                //mrb if (strcmp(bound,"finite") == 0) {
                //mrb }
     
                  x  = expo[todos];
                  x1 = x;
                  for (i = 1; i <= deg; i++) {
                    index       = todos + i - 1;
                    expo[index] = x1;
                  }
           
                  if (strcmp(bound,"finite") == 0) {
                    i = 1;
                     do
                        {
                         index            = todos + i - 1;
                         expo_no_correcto =  check_expotents_finite(using_gamma,
                                                                    gamma_couple,
                                                                    Rc,
                                                                    np,
                                                                    mang,
                                                                    expo,
                                                                    index);
                         if(expo_no_correcto == 1)
                         break;
                         i++;
                        } while(i <= deg);
                  } else
                        expo_no_correcto = 0;
           
                  checking_expo = checking_linear_dependece(nt,
                                                            &bt,
                                                            expo,
                                                            mang,
                                                            nombre);
                  
           
                  if (expo_no_correcto == 1 || checking_expo == 1) {
                    test_scf = 1;
                    energia  = maximo;
                  } else
                        test_scf = scf(nt, 
                                       elecalfa, 
                                       elecbeta, 
                                       z, 
                                       orbital, 
                                       tol, 
                                       using_gamma,
                                       correlation,
                                       propagador,
                                       mezcla,
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
                                       gamma_nicp, 
                                       0,
                                       plasma);
           
                  if (test_scf == 0) {
                    ener_array[1] = energia;
                    printf("1 Expo %d, value=%f, energ = %f\n", todos, 
                                                                x1, 
                                                                energia);
                  } else {
                          printf("Convergence problems \n");
                          ener_array[1] = maximo;
                          for (i = 1; i <= deg; i++) {
                            index       = todos + i - 1;
                            expo[index] = x;
                          }
                    }
           
                  x0 = x - step_work;
                  if (x0 <= 0.f) { 
                    x0            = x;
                    checking_expo = 1;
                  } else {
                          for (i = 1; i <= deg; i++) {
                             index       = todos + i - 1;
                             expo[index] = x0;
                          }
           
                          if (strcmp(bound,"finite") == 0) {
                            i = 1;
                             do
                                {
                                 index            = todos + i - 1;
                                 expo_no_correcto =  check_expotents_finite(using_gamma,
                                                                            gamma_couple,
                                                                            Rc,
                                                                            np,
                                                                            mang,
                                                                            expo,
                                                                            index);
                                 break;
                                 
                                 i++;
                                } while(i <= deg);
                          } else
                                expo_no_correcto = 0;
           
                          checking_expo = checking_linear_dependece(nt,
                                                              &bt,
                                                              expo,
                                                              mang,
                                                              nombre);
                    }
           
                  if (expo_no_correcto == 1 || checking_expo == 1) {
                    test_scf = 1;
                    x0       = x;
                    energia  = maximo;
                  } else
                        test_scf = scf(nt, 
                                       elecalfa, 
                                       elecbeta, 
                                       z, 
                                       orbital, 
                                       tol, 
                                       using_gamma,
                                       correlation,
                                       propagador,
                                       mezcla,
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
                                       gamma_nicp, 
                                       0,
                                       plasma);
           
                  if (test_scf == 0) {
                    ener_array[0] = energia;
                    printf("2 Expo %d, value=%f, energ = %f\n", todos, 
                                                                x0, 
                                                                energia);
                  } else {
                          printf("Convergence problems \n");
                          ener_array[0] = maximo;
                          for (i = 1; i <= deg; i++) {
                            index       = todos + i - 1;
                            expo[index] = x;
                          }
                    }
           
                  x2 = x + step_work;
                  for (i = 1; i <= deg; i++) {
                    index       = todos + i - 1;
                    expo[index] = x2;
                  }
           
                  if (strcmp(bound,"finite") == 0) {
                    i = 1;
                     do
                        {
                         index            = todos + i - 1;
                         expo_no_correcto =  check_expotents_finite(using_gamma,
                                                                    gamma_couple,
                                                                    Rc,
                                                                    np,
                                                                    mang,
                                                                    expo,
                                                                    index);
                         if(expo_no_correcto == 1)
                         break;
                         i++;
                        } while(i <= deg);
                  } else
                        expo_no_correcto = 0;
           
                  checking_expo = checking_linear_dependece(nt,
                                                            &bt,
                                                            expo,
                                                            mang,
                                                            nombre);
           
           
                  if (expo_no_correcto == 1 || checking_expo == 1) {
                    test_scf = 1;
                    x2       = x;
                    energia  = maximo;
                  } else 
                        test_scf = scf(nt, 
                                       elecalfa, 
                                       elecbeta, 
                                       z, 
                                       orbital, 
                                       tol, 
                                       using_gamma,
                                       correlation,
                                       propagador,
                                       mezcla,
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
                                       gamma_nicp, 
                                       0,
                                       plasma);
           
                  if (test_scf == 0) {
                    ener_array[2] = energia;
                    printf("3 Expo %d, value=%f, energ = %f\n", todos, 
                                                                x2, 
                                                                energia);
                  } else {
                          printf("Convergence problems \n");
                          ener_array[2] = maximo;
                          for (i = 1; i <= deg; i++) {
                            index       = todos + i - 1;
                            expo[index] = x;
                      }
                    }
             
                  dif_1        = fabs(ener_array[0] - ener_array[1]);
                  dif_2        = fabs(ener_array[0] - ener_array[2]);
                  dif_3        = fabs(ener_array[2] - ener_array[1]);
                  dif_max      = (dif_1 > dif_2 ? dif_1 : dif_2); 
                  dif_max      = (dif_max > dif_3 ? dif_max : dif_3);
           
                  dif_1_expo   = fabs(x0 - x1);
                  dif_2_expo   = fabs(x0 - x2);
                  dif_3_expo   = fabs(x2 - x1);
                  dif_max_expo = (dif_1_expo > dif_2_expo ? dif_1_expo : dif_2_expo); 
                  dif_max_expo = (dif_max_expo > dif_3_expo ? dif_max_expo : dif_3_expo);
           
                  printf("x0 = %f, e0 = %18.10f\n", x0, 
                                                    ener_array[0]);
                  printf("x1 = %f, e1 = %18.10f\n", x1, 
                                                    ener_array[1]);
                  printf("x2 = %f, e2 = %18.10f\n", x2, 
                                                    ener_array[2]);
                  printf("dif max = %14.10f\n", dif_max);
                  printf("dif max expo = %14.10f\n", dif_max_expo);
           
                  if (dif_max <= tol_opt || dif_max_expo <= tol_opt_expo) {
                      bandera = 1;
                      for (i = 1; i <= deg; i++) {
                        index       = todos + i - 1;
                        expo[index] = x1;
                        printf("4 Expo %d final value = %f\n", index, 
                                                               expo[index]);
                      }
                      todos = index + 1;
                  } else {//label 2
                          do {//label 3
                              if (ener_array[1] < ener_array[0]) {//labelcomp1
                                if (ener_array[1] < ener_array[2]) {
                                   bandera = 1;
                                   for (i = 1; i <= deg; i++) {
                                     index       = todos + i - 1;
                                     expo[index] = x1;
                                     printf("5 Expo %d final value = %f\n", index, 
                                                                            expo[index]);
                                   }
                                   todos = index + 1;
                                } else {//labelcomp2
                                        x0            = x1;
                                        ener_array[0] = ener_array[1];
                                        x1            = x2;
                                        ener_array[1] = ener_array[2];
                                        x2            = x1 + step_work;
                                        for (i = 1; i <= deg; i++) {
                                           index       = todos + i - 1;
                                           expo[index] = x2;
                                        }
                                        if (strcmp(bound,"finite") == 0) {
                                          i = 1;
                                           do
                                              {
                                               index            = todos + i - 1;
                                               expo_no_correcto =  check_expotents_finite(using_gamma,
                                                                                          gamma_couple,
                                                                                          Rc,
                                                                                          np,
                                                                                          mang,
                                                                                          expo,
                                                                                          index);
                                               if(expo_no_correcto == 1)
                                               break;
                                               i++;
                                              } while(i <= deg);
                                        } else
                                              expo_no_correcto = 0;
           
                                        checking_expo = checking_linear_dependece(nt,
                                                                                  &bt,
                                                                                  expo,
                                                                                  mang,
                                                                                  nombre);
           
           
                                        if (expo_no_correcto == 1 || checking_expo == 1) {
                                          test_scf = 1;
                                          energia  = maximo;
                                        } else 
                                              test_scf = scf(nt, 
                                                             elecalfa, 
                                                             elecbeta, 
                                                             z,
                                                             orbital, 
                                                             tol, 
                                                             using_gamma,
                                                             correlation,
                                                             propagador,
                                                             mezcla,
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
                                                             gamma_nicp, 
                                                             0,
                                                             plasma);
                                        if (test_scf == 0) {
                                          ener_array[2] = energia;
                                          printf("6 Expo %d, value=%f, energ = %f\n",todos, 
                                                                                     x2, 
                                                                                     energia);
                                        } else {
                                                printf("Convergence problems \n");
                                                ener_array[2] = maximo;
                                        }
                                }//labelcomp2
               /*labelcomp1*/ } else {//labelcomp7
                                      if (ener_array[1] < ener_array[2]) {//labelcomp3
                                        x2            = x1;
                                        ener_array[2] = ener_array[1];
                                        x1            = x0;
                                        ener_array[1] = ener_array[0];
                                        x0            = x1 - step_work;
                                        if (x0 <= 0.f) {
                                          x0         = 0.01f;
                                          expo_bound = 1;
                                        }
                                        for (i = 1; i <= deg; i++) {
                                          index       = todos + i - 1;
                                          expo[index] = x0;
                                        }
                                        if (strcmp(bound,"finite") == 0) {
                                          i = 1;
                                           do
                                              {
                                               index            = todos + i - 1;
                                               expo_no_correcto =  check_expotents_finite(using_gamma,
                                                                                          gamma_couple,
                                                                                          Rc,
                                                                                          np,
                                                                                          mang,
                                                                                          expo,
                                                                                          index);
                                               if(expo_no_correcto == 1)
                                               break;
                                               i++;
                                              } while(i <= deg);
                                        } else
                                              expo_no_correcto = 0;
           
                                        checking_expo = checking_linear_dependece(nt,
                                                                                  &bt,
                                                                                  expo,
                                                                                  mang,
                                                                                  nombre);
           
           
                                        if(expo_no_correcto == 1 || checking_expo == 1) {
                                        test_scf = 1;
                                        energia  = maximo;
                                        } else
                                              test_scf = scf(nt, 
                                                             elecalfa, 
                                                             elecbeta, 
                                                             z,
                                                             orbital, 
                                                             tol, 
                                                             using_gamma,
                                                             correlation,
                                                             propagador,
                                                             mezcla,
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
                                                             gamma_nicp, 
                                                             0,
                                                             plasma);
                                        if (test_scf == 0) {
                                          ener_array[0] = energia;
                                          printf("7 Expo %d, value=%f, energ = %f\n", todos, 
                                                                                      x0, 
                                                                                      energia);
                                        } else {
                                            printf("Convergence problems \n");
                                            ener_array[0] = maximo;
                                          }
                      /*labelcomp3*/  } else {//labelcomp6
                                              if (ener_array[2] < ener_array[0]) {//labelcomp4
                                                x0            = x1;
                                                ener_array[0] = ener_array[1];
                                                x1            = x2;
                                                ener_array[1] = ener_array[2];
                                                x2            = x1 + step_work;
                                                for (i = 1; i <= deg; i++) {
                                                  index       = todos + i - 1;
                                                  expo[index] = x0;
                                                }
                                                if (strcmp(bound,"finite") == 0) {
                                                  i = 1;
                                                   do
                                                      {
                                                       index            = todos + i - 1;
                                                       expo_no_correcto =  check_expotents_finite(using_gamma,
                                                                                                  gamma_couple,
                                                                                                  Rc,
                                                                                                  np,
                                                                                                  mang,
                                                                                                  expo,
                                                                                                  index);
                                                       if(expo_no_correcto == 1)
                                                       break;
                                                       i++;
                                                      }
                                                   while(i <= deg);
                                                } else
                                                      expo_no_correcto = 0;
           
                                                checking_expo = checking_linear_dependece(nt,
                                                                                          &bt,
                                                                                          expo,
                                                                                          mang,
                                                                                          nombre);
           
           
                                                if (expo_no_correcto == 1 || checking_expo == 1) {
                                                  test_scf = 1;
                                                  energia  = maximo;
                                                } else
                                                      test_scf = scf(nt, 
                                                                     elecalfa, 
                                                                     elecbeta, 
                                                                     z,
                                                                     orbital, 
                                                                     tol, 
                                                                     using_gamma,
                                                                     correlation,
                                                                     propagador,
                                                                     mezcla,
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
                                                                     gamma_nicp, 
                                                                     0,
                                                                     plasma);
                                                if (test_scf == 0) {
                                                  ener_array[2] = energia;
                                                  printf("8 Expo %d, value=%f, energ = %f\n", todos, 
                                                                                              x2, 
                                                                                              energia);
                                                } else {
                                                        printf("Convergence problems \n");
                                                        ener_array[2] = maximo;
                                                  }
                          /*labelcomp4*/      } else {//labelcomp5
                                                      x2            = x1;
                                                      ener_array[2] = ener_array[1];
                                                      x1            = x0;
                                                      ener_array[1] = ener_array[0];
                                                      x0            = x1 - step_work;
                                                      if (x0 <= 0.f) {
                                                         x0         = 0.01f;
                                                         expo_bound = 1;
                                                      }
                                                      for (i = 1; i <= deg; i++) {
                                                         index       = todos + i - 1;
                                                         expo[index] = x0;
                                                      }
                                                      if (strcmp(bound,"finite") == 0) {
                                                        i = 1;
                                                         do
                                                            {
                                                             index            = todos + i - 1;
                                                             expo_no_correcto =  check_expotents_finite(using_gamma,
                                                                                                        gamma_couple,
                                                                                                        Rc,
                                                                                                        np,
                                                                                                        mang,
                                                                                                        expo,
                                                                                                        index);
                                                             if(expo_no_correcto == 1)
                                                             break;
                                                             i++;
                                                          } while(i <= deg);
                                                      } else
                                                            expo_no_correcto = 0;
           
                                                      checking_expo = checking_linear_dependece(nt,
                                                                                                &bt,
                                                                                                expo,
                                                                                                mang,
                                                                                                nombre);
           
                                                      if (expo_no_correcto == 1 || checking_expo == 1) {
                                                        test_scf = 1;
                                                        energia  = maximo;
                                                      } else
                                                            test_scf = scf(nt, 
                                                                           elecalfa, 
                                                                           elecbeta, 
                                                                           z,
                                                                           orbital, 
                                                                           tol, 
                                                                           using_gamma,
                                                                           correlation,
                                                                           propagador,
                                                                           mezcla,
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
                                                                           gamma_nicp, 
                                                                           0,
                                                                           plasma);
                                                          if (test_scf == 0) {
                                                            ener_array[0] = energia;
                                                            printf("9 Expo %d, value=%f, energ = %f\n", todos, 
                                                                                                        x0, 
                                                                                                        energia);
                                                          } else {
                                                                  printf("Convergence problems \n");
                                                                  ener_array[0] = maximo;
                                                          }
           
                                                }//labelcomp5
                                        }//labelcomp6
                                }//labelcomp7
           
                              dif_1        = fabs(ener_array[0] - ener_array[1]);
                              dif_2        = fabs(ener_array[0] - ener_array[2]);
                              dif_3        = fabs(ener_array[2] - ener_array[1]);
                              dif_max      = (dif_1 > dif_2 ? dif_1 : dif_2); 
                              dif_max      = (dif_max > dif_3 ? dif_max : dif_3);
           
                              dif_1_expo   = fabs(x0 - x1);
                              dif_2_expo   = fabs(x0 - x2);
                              dif_3_expo   = fabs(x2 - x1);
                              dif_max_expo = (dif_1_expo > dif_2_expo ? dif_1_expo : dif_2_expo); 
                              dif_max_expo = (dif_max_expo > dif_3_expo ? dif_max_expo : dif_3_expo);
           
                              if (expo_bound == 1) {
                                  dif_max         = 0.f;
                                  opt_flag[todos] = 0;
                              }
                              if (dif_max <= tol_opt || dif_max_expo <= tol_opt_expo) {
                                 bandera = 1;
                                 for (i = 1; i <= deg; i++) {
                                    index       = todos + i - 1;
                                    expo[index] = x1;
                                    printf("10 Expo %d final value = %f\n", index, expo[index]);
                                 }
                              todos = index + 1;
                              }
                              printf("dif max = %14.10f\n", dif_max);
           
                          } while (bandera == 0);//label 3
                    }//label 2
                }//label 1
        } while(todos < nt);
     
        step_work = step_work/(double)10.;
       }
       printf("After partial optimization\n");
       for (i = 0; i < nt; i++) 
          printf("@ Partial opt Expo %d : %20.12f\n", i, expo[i]);

     for (j = 0; j < nt; j++)
        {
         if (fabs(expo[j] - expo_diff[j]) > maxdiff)
           maxdiff = fabs(expo[j] - expo_diff[j]);
         expo_diff[j] = expo[j];
     }
 } while (maxdiff > 1.e-5);
 free(p);
 p = 0;
 }
