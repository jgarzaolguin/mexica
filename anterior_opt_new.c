#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

 double alpha_value(int i, double incr, double valpha){
        double value;
        value =(double) valpha + incr*i;
        return (value);
        }

 double beta_value(int i, double incr, double vbeta){
        double value;
        value =(double) vbeta + incr*i;
        return (value);
        }

 int optimiza_new(int     nt,
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
                  int     imprime,
                  int     plasma) {

 
 extern double alpha_value(int i, double incr, double valpha);
 extern double beta_value(int i, double incr, double vbeta);
 
 extern int check_expotents_finite(char   *using_gamma,
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
                int     imprime,
                int     plasma,
                double *cusp_kato);

 int i, j, k;
 int index, todos, ang_test, deg, test_scf;
 double step_work, energia;
 int bt;
 double cusp_kato;

 int expo_no_correcto,
     checking_expo;

     checking_expo    = 0;
     expo_no_correcto = 0;

 double  factor1, factor2;
 double  expalpha_opt[36], expbeta_opt[36];
 double  expalpha[36], expbeta[36];
 double  rango_a[36], rango_b[36];
 double  valpha_opt[36];
 double  vbeta_opt[36];
 int     k1[36];
 int     total_a[36], total_b[36];
 int     pasos;
 double  incr;
 double  array_cusp_kato[100];
 double  array_energy[100];
 double  *array_alpha, *array_beta;
 double  energia_0;
 double  kato_0;
 int     index_t;
 int unsigned q, w;
 int    ident;
 int    ident_test;
 int    rand_1;
 int    rand_2;
 int    iteraciones;
 int    identify;
 double diff_kato_0;
 double diff_kato_1;

  array_alpha     = (double*)malloc(5000*sizeof(double));
 if (array_alpha == NULL) {
   free(array_alpha);
   array_alpha = 0;
   printf("Error en arreglo alpha, insuficiente RAM\n");
   exit(1);
 }

  array_beta     = (double*)malloc(5000*sizeof(double));
 if (array_beta == NULL) {
   free(array_beta);
   array_beta = 0;
   printf("Error en arreglo beta, insuficiente RAM\n");
   exit(1);
 }

 double alpha_add;
 double beta_add;
 alpha_add = 0.f;
 beta_add  = 0.f;

 if (z >= 2) {
   alpha_add = 0.2;
   beta_add  = 1.2;
 } else {
      alpha_add = 0.01;
      beta_add  = 1.1;
   }
 

 pasos = 100000;
 
 for (j = 0; j < 36; j++) {
  rango_a[j]          = 3.f;
  rango_b[j]          = 3.3;
  total_a[j]          = (int) pasos*rango_a[j];
  total_b[j]          = (int) pasos*rango_b[j];
  expalpha_opt[j]     = 0.0;
  expbeta_opt[j]      = 0.0;
  valpha_opt[j]       = alpha_add;
  vbeta_opt[j]        = beta_add;
  expalpha[j]         = 0.f;
  expbeta[j]          = 0.f;
  k1[j]               = 0;
 }

 incr   = (double) rango_a[0]/total_a[0];

//mrb1 for (j = 0; j < 36; j++) 
//mrb1   printf("\n KONDA %lf  %lf  %d %f\n", alpha_value(total_a[j] - 1, incr, valpha_opt[j]), beta_value(total_b[j] - 1, incr, vbeta_opt[j]), total_b[j], incr);

 srand((unsigned) time(NULL));


 energia_0  = 1e9;
 kato_0     = 1.1;
 ident      = 0;
 ident_test = 0;

 iteraciones = 11000;

 for (w = 0; w <= 2; w++) {
      for (q = 0; q < iteraciones; q++) {
          index_t   = 0;
          k         = 0;
          factor1   = 0.f;
          factor2   = 0.f;
          printf("\n STEP (%d) %d\n", w, q);

           for (j = 0; j < 36; j++) {
              rand_1     = rand()%(total_a[j]);
              rand_2     = rand()%(total_b[j]);
              expalpha[j]  = alpha_value(rand_1, incr, valpha_opt[j]);
              expbeta[j]   = beta_value(rand_2, incr, vbeta_opt[j]);
           }

           for (j = 0; j < 36; j++)
           k1[j]    = 0;
           todos    = 0;
           index    = 0;

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

               if (np[todos] == 1)
               switch(deg) {
                 case(1)  : k = k1[0]; index_t = 0; break;
               }
               if (np[todos] == 2)
               switch(deg) {
                 case(1)  : k = k1[1]; index_t = 1; break;
                 case(3)  : k = k1[2]; index_t = 2; break;
               }
               if (np[todos] == 3)
               switch(deg) {
                 case(1)  : k = k1[3]; index_t = 3; break;
                 case(3)  : k = k1[4]; index_t = 4; break;
                 case(5)  : k = k1[5]; index_t = 5; break;
               }
               if (np[todos] == 4)
               switch(deg) {
                 case(1)  : k = k1[6]; index_t = 6; break;
                 case(3)  : k = k1[7]; index_t = 7; break;
                 case(5)  : k = k1[8]; index_t = 8; break;
                 case(7)  : k = k1[9]; index_t = 9; break;
               }
               if (np[todos] == 5)
               switch(deg) {
                 case(1)  : k = k1[10]; index_t = 10; break;
                 case(3)  : k = k1[11]; index_t = 11; break;
                 case(5)  : k = k1[12]; index_t = 12; break;
                 case(7)  : k = k1[13]; index_t = 13; break;
                 case(9)  : k = k1[14]; index_t = 14; break;
               }
               if (np[todos] == 6)
               switch(deg) {
                 case(1)  : k = k1[15]; index_t = 15; break;
                 case(3)  : k = k1[16]; index_t = 16; break;
                 case(5)  : k = k1[17]; index_t = 17; break;
                 case(7)  : k = k1[18]; index_t = 18; break;
                 case(9)  : k = k1[19]; index_t = 19; break;
                 case(11) : k = k1[20]; index_t = 20; break;
               }
               if (np[todos] == 7)
               switch(deg) {
                 case(1)  : k = k1[21]; index_t = 21; break;
                 case(3)  : k = k1[22]; index_t = 22; break;
                 case(5)  : k = k1[23]; index_t = 23; break;
                 case(7)  : k = k1[24]; index_t = 24; break;
                 case(9)  : k = k1[25]; index_t = 25; break;
                 case(11) : k = k1[26]; index_t = 26; break;
                 case(13) : k = k1[27]; index_t = 27; break;
               }
               if (np[todos] == 8)
               switch(deg) {
                 case(1)  : k = k1[28]; index_t = 28; break;
                 case(3)  : k = k1[29]; index_t = 29; break;
                 case(5)  : k = k1[30]; index_t = 30; break;
                 case(7)  : k = k1[31]; index_t = 31; break;
                 case(9)  : k = k1[32]; index_t = 32; break;
                 case(11) : k = k1[33]; index_t = 33; break;
                 case(13) : k = k1[34]; index_t = 34; break;
                 case(15) : k = k1[35]; index_t = 35; break;
               }
                 if (k == 0)
                   factor1 =(double) expalpha[index_t];
                   else
                       factor1 =(double) expalpha[index_t]*pow(expbeta[index_t],(double) k);

               printf("expoalpha %lf expobeta %lf total %lf (k %d|||deg %d) expo_opt %lf %lf\n",(double) expalpha[index_t],
                                                                                                   (double) expbeta[index_t],
                                                                                                   factor1, 
                                                                                                   k, 
                                                                                                   deg, 
                                                                                                   expalpha_opt[index_t], 
                                                                                                   expbeta_opt[index_t]);
               for (i = 1; i <= deg; i++) {
                  index = todos + i - 1;
                  expo[index] = factor1;
               }


               if (strcmp(bound,"finite") == 0) {
                      expo_no_correcto =  check_expotents_finite(using_gamma,
                                                                 gamma_couple,
                                                                 Rc,
                                                                 np,
                                                                 mang,
                                                                 expo,
                                                                 todos);
               }

               if (np[todos] == 1)
               switch(deg) {
                 case(1)  : k1[0]++; break;
               }
               if (np[todos] == 2)
               switch(deg) {
                 case(1)  : k1[1]++; break;
                 case(3)  : k1[2]++; break;
               }
               if (np[todos] == 3)
               switch(deg) {
                 case(1)  : k1[3]++; break;
                 case(3)  : k1[4]++; break;
                 case(5)  : k1[5]++; break;
               }
               if (np[todos] == 4)
               switch(deg) {
                 case(1)  : k1[6]++; break;
                 case(3)  : k1[7]++; break;
                 case(5)  : k1[8]++; break;
                 case(7)  : k1[9]++; break;
               }
               if (np[todos] == 5)
               switch(deg) {
                 case(1)  : k1[10]++; break;
                 case(3)  : k1[11]++; break;
                 case(5)  : k1[12]++; break;
                 case(7)  : k1[13]++; break;
                 case(9)  : k1[14]++; break;
               }
               if (np[todos] == 6)
               switch(deg) {
                 case(1)  : k1[15]++; break;
                 case(3)  : k1[16]++; break;
                 case(5)  : k1[17]++; break;
                 case(7)  : k1[18]++; break;
                 case(9)  : k1[19]++; break;
                 case(11) : k1[20]++; break;
               }
               if (np[todos] == 7)
               switch(deg) {
                 case(1)  : k1[21]++; break;
                 case(3)  : k1[22]++; break;
                 case(5)  : k1[23]++; break;
                 case(7)  : k1[24]++; break;
                 case(9)  : k1[25]++; break;
                 case(11) : k1[26]++; break;
                 case(13) : k1[27]++; break;
               }
               if (np[todos] == 8)
               switch(deg) {
                 case(1)  : k1[28]++; break;
                 case(3)  : k1[29]++; break;
                 case(5)  : k1[30]++; break;
                 case(7)  : k1[31]++; break;
                 case(9)  : k1[32]++; break;
                 case(11) : k1[33]++; break;
                 case(13) : k1[34]++; break;
                 case(15) : k1[35]++; break;
               }


               if (expo_no_correcto == 1)
                 todos = nt;
                 else
                   todos = index + 1;

           } while(todos < nt);


           if (expo_no_correcto == 0)
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
                            0,
                            plasma,
                           &cusp_kato);
              else
                  test_scf = 1;

            if (test_scf != 0) {
              energia   = 1.e10;
              //cusp_kato = 1000.f;
              diff_kato_0 = 1000.f;
            } else {             
                  diff_kato_0 = fabs(cusp_kato - 1.f);
              }

            if (energia < energia_0 && diff_kato_0 <  kato_0) {
               energia_0 = energia;
               kato_0 = diff_kato_0;
               array_cusp_kato[ident] = cusp_kato;
               array_energy[ident]    = energia;
               printf("\n%d HOLA1 energy %25.12lf Kato %lf\n", ident, energia, array_cusp_kato[ident]);
               for (j = 0; j < 36; j++) {
                    array_alpha[ident*35 + j] = expalpha[j];
                    array_beta[ident*35 + j]  =  expbeta[j];
                    //mrb1 printf("\nEXPO_OPT a %f b %f || a %f b %f\n", expalpha_opt[j], expbeta_opt[j], array_alpha[ident*35 + j], array_beta[ident*35 + j]);
               }
              ident++;
            }
   } //primer for q
 if (w == 2) break;

//identifica la mejor condicion de kato cercana a 1 y mejor energia
 int unsigned cut;

 if (w == 0 || w == 1) {
    if (w == 1) {
     pasos = pasos*10;
     cut = w;
    }
    identify  = ident - 1;
    energia_0 = 1e9;
    for (i = 0; i <= identify; i++) {
       if (array_energy[i] < energia_0) {
         energia_0 = array_energy[i];
         ident_test = i;
       }
    }
    energia_0 = array_energy[ident_test];
    printf("\n %d 1Good energy cusp kato test %f energy %lf\n", cut, array_cusp_kato[ident_test], array_energy[ident_test]);
   }

//mrb2 if (w == 1) {
//mrb2   identify = ident - 1;
//mrb2   diff_kato_0 = 0.f;
//mrb2   diff_kato_1 = 10000.f;

//mrb2   for (i = 0; i <= identify; i++) {
//mrb2//mrb   printf("All cusp kato test %f\n", array_cusp_kato[i]);
//mrb2      diff_kato_0 = fabs(array_cusp_kato[i] - 1.f);
//mrb2      if (diff_kato_0 < diff_kato_1){
//mrb2        diff_kato_1 = diff_kato_0;
//mrb2        ident_test = i;
//mrb2      }
//mrb2   }
//mrb2  energia_0 = array_energy[ident_test];
//mrb2//mrb  kato_0 = array_cusp_kato[ident_test];
//mrb2  printf("\n %d 2Good cusp kato test %f energy %lf\n", cut, array_cusp_kato[ident_test], array_energy[ident_test]);
//mrb2 }

   for (j = 0; j < 36; j++) {
      expalpha_opt[j] = array_alpha[ident_test*35 + j];
      expbeta_opt[j]  = array_beta[ident_test*35 + j];
     
  if (expalpha_opt[j] - (double) 1.f/(1 + cut) > 0.f)
    valpha_opt[j] = expalpha_opt[j] - (double) 1.f/(1 + cut);
     else
  if (expalpha_opt[j] - (double) 0.9/(1 + cut) > 0.f)
    valpha_opt[j] = expalpha_opt[j] - (double) 0.9/(1 + cut);
    else
  if (expalpha_opt[j] - (double) 0.8/(1 + cut) > 0.f)
    valpha_opt[j] = expalpha_opt[j] - (double) 0.8/(1 + cut);
    else
  if (expalpha_opt[j] - (double) 0.7/(1 + cut) > 0.f)
    valpha_opt[j] = expalpha_opt[j] - (double) 0.7/(1 + cut);
    else
  if (expalpha_opt[j] - (double) 0.6/(1 + cut) > 0.f)
    valpha_opt[j] = expalpha_opt[j] - (double) 0.6/(1 + cut);
    else
   if (expalpha_opt[j] - (double) 0.5/(1 + cut) > 0.f)
     valpha_opt[j] = expalpha_opt[j] - (double) 0.5/(1 + cut);
     else
   if (expalpha_opt[j] - (double) 0.4/(1 + cut) > 0.f)
     valpha_opt[j] = expalpha_opt[j] - (double) 0.4/(1 + cut);
     else
   if (expalpha_opt[j] - (double) 0.3/(1 + cut) > 0.f)
     valpha_opt[j] = expalpha_opt[j] - (double) 0.3/(1 + cut);
     else
   if (expalpha_opt[j] - (double) 0.2/(1 + cut) > 0.f)
     valpha_opt[j] = expalpha_opt[j] - (double) 0.2/(1 + cut);
     else
   if (expalpha_opt[j] - (double) 0.1/(1 + cut) > 0.f)
     valpha_opt[j] = expalpha_opt[j] - (double) 0.1/(1 + cut);
     else
     valpha_opt[j] = expalpha_opt[j];
   
   
  if (expbeta_opt[j]  - (double) 1.f/(1 + cut) > 1.1)
     vbeta_opt[j]  = expbeta_opt[j]  - (double) 1.f/(1 + cut);
      else
  if (expbeta_opt[j]  - (double) 0.9/(1 + cut) > 1.1)
     vbeta_opt[j]  = expbeta_opt[j]  - (double) 0.9/(1 + cut);
      else
  if (expbeta_opt[j]  - (double) 0.8/(1 + cut) > 1.1)
     vbeta_opt[j]  = expbeta_opt[j]  - (double) 0.8/(1 + cut);
      else
  if (expbeta_opt[j]  - (double) 0.7/(1 + cut) > 1.1)
     vbeta_opt[j]  = expbeta_opt[j]  - (double) 0.7/(1 + cut);
      else
  if (expbeta_opt[j]  - (double) 0.6/(1 + cut) > 1.1)
     vbeta_opt[j]  = expbeta_opt[j]  - (double) 0.6/(1 + cut);
      else
   if (expbeta_opt[j]  - (double) 0.5/(1 + cut) > 1.1)
      vbeta_opt[j]  = expbeta_opt[j]  - (double) 0.5/(1 + cut);
       else
   if (expbeta_opt[j]  - (double) 0.4/(1 + cut) > 1.1)
      vbeta_opt[j]  = expbeta_opt[j]  - (double) 0.4/(1 + cut);
       else
   if (expbeta_opt[j]  - (double) 0.3/(1 + cut) > 1.1)
      vbeta_opt[j]  = expbeta_opt[j]  - (double) 0.3/(1 + cut);
       else
   if (expbeta_opt[j]  - (double) 0.2/(1 + cut) > 1.1)
      vbeta_opt[j]  = expbeta_opt[j]  - (double) 0.2/(1 + cut);
       else
   if (expbeta_opt[j]  - (double) 0.1/(1 + cut) > 1.1)
      vbeta_opt[j]  = expbeta_opt[j]  - (double) 0.1/(1 + cut);
       else
      vbeta_opt[j]  = expbeta_opt[j];
   
    //mrb rango_a[j] = expalpha_opt[j] + (double) 1.0f/(1 + cut);
    //mrb rango_b[j] = expbeta_opt[j] + (double) 1.0f/(1 + cut);
     rango_a[j] = expalpha_opt[j] + (double) 0.5f/(1 + cut);
     rango_b[j] = expbeta_opt[j] + (double) 0.5f/(1 + cut);
   
     total_a[j] = (int) pasos*rango_a[j];
     total_b[j] = (int) pasos*rango_b[j];
   
     incr   = (double) rango_a[0]/total_a[0];
 
   //mrb1 printf("\nHOLA2 %f %f %f %f %f %f %d %d %f\n", expalpha_opt[j], expbeta_opt[j], valpha_opt[j], vbeta_opt[j], rango_a[j], rango_b[j], total_a[j], total_b[j], incr);
  }
  iteraciones = 5000; //iteraciones/4 + iteraciones/2;
 

 } //segundo for w


 identify = ident - 1;
 //mrb1 *identify = ident - 1;
                           
 //mrb1 for (j = 0; j <= identify; j++)
 //mrb1  printf("KATOO %lf\n", array_cusp_kato[j]);
 //mrb1  printf("\n");



 ident = 0;
 diff_kato_0 = 0.f;
 diff_kato_1 = 10000.f;
 for (i = 0; i <= identify; i++) {
    printf("All cusp kato test %f\n", array_cusp_kato[i]);
    diff_kato_0 = fabs(array_cusp_kato[i] - 1.f);
    if (diff_kato_0 < diff_kato_1){
      diff_kato_1 = diff_kato_0;
      ident = i;
    }
 }
 printf("\nCUSP KATO FINAL %f\n", array_cusp_kato[ident]);

 for (i = 0; i < 36; i++) {
 //mrb1 printf("alpha %lf  beta %lf", array_alpha[ident*35 + i], array_beta[ident*35 + i]);
    expalpha_opt[i] = array_alpha[ident*35 + i];
    expbeta_opt[i]  = array_beta[ident*35 + i];
 }
 printf("\n"); 
 for (j = 0; j < 36; j++) {
   if (k1[j] != 0) {
     printf(" %d FINAL alpha %lf  beta %lf\n", j, expalpha_opt[j], expbeta_opt[j]);
   }
 }
 printf("\n"); 

 ////////////////
 ////////////////Crea exponentes

 index_t   = 0;
 k         = 0;
 factor1   = 0.f;
 factor2   = 0.f;

 for (j = 0; j < 36; j++)
    k1[j]    = 0;
 todos    = 0;
 index    = 0;

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

     if (np[todos] == 1)
     switch(deg) {
       case(1)  : k = k1[0]; index_t = 0; break;
     }
     if (np[todos] == 2)
     switch(deg) {
       case(1)  : k = k1[1]; index_t = 1; break;
       case(3)  : k = k1[2]; index_t = 2; break;
     }
     if (np[todos] == 3)
     switch(deg) {
       case(1)  : k = k1[3]; index_t = 3; break;
       case(3)  : k = k1[4]; index_t = 4; break;
       case(5)  : k = k1[5]; index_t = 5; break;
     }
     if (np[todos] == 4)
     switch(deg) {
       case(1)  : k = k1[6]; index_t = 6; break;
       case(3)  : k = k1[7]; index_t = 7; break;
       case(5)  : k = k1[8]; index_t = 8; break;
       case(7)  : k = k1[9]; index_t = 9; break;
     }
     if (np[todos] == 5)
     switch(deg) {
       case(1)  : k = k1[10]; index_t = 10; break;
       case(3)  : k = k1[11]; index_t = 11; break;
       case(5)  : k = k1[12]; index_t = 12; break;
       case(7)  : k = k1[13]; index_t = 13; break;
       case(9)  : k = k1[14]; index_t = 14; break;
     }
     if (np[todos] == 6)
     switch(deg) {
       case(1)  : k = k1[15]; index_t = 15; break;
       case(3)  : k = k1[16]; index_t = 16; break;
       case(5)  : k = k1[17]; index_t = 17; break;
       case(7)  : k = k1[18]; index_t = 18; break;
       case(9)  : k = k1[19]; index_t = 19; break;
       case(11) : k = k1[20]; index_t = 20; break;
     }
     if (np[todos] == 7)
     switch(deg) {
       case(1)  : k = k1[21]; index_t = 21; break;
       case(3)  : k = k1[22]; index_t = 22; break;
       case(5)  : k = k1[23]; index_t = 23; break;
       case(7)  : k = k1[24]; index_t = 24; break;
       case(9)  : k = k1[25]; index_t = 25; break;
       case(11) : k = k1[26]; index_t = 26; break;
       case(13) : k = k1[27]; index_t = 27; break;
     }
     if (np[todos] == 8)
     switch(deg) {
       case(1)  : k = k1[28]; index_t = 28; break;
       case(3)  : k = k1[29]; index_t = 29; break;
       case(5)  : k = k1[30]; index_t = 30; break;
       case(7)  : k = k1[31]; index_t = 31; break;
       case(9)  : k = k1[32]; index_t = 32; break;
       case(11) : k = k1[33]; index_t = 33; break;
       case(13) : k = k1[34]; index_t = 34; break;
       case(15) : k = k1[35]; index_t = 35; break;
     }
       if (k == 0)
         factor1 =(double) expalpha_opt[index_t];
         else
             factor1 =(double) expalpha_opt[index_t]*pow(expbeta_opt[index_t],(double) k);

     //mrb1              printf("expoalpha %lf expobeta %lf total %lf (k %d|||deg %d)\n",(double) expalpha_opt[index_t],
     //mrb1                                                                                         (double) expbeta_opt[index_t],
     //mrb1                                                                                         factor1,
     //mrb1                                                                                         k,
     //mrb1                                                                                         deg);
     for (i = 1; i <= deg; i++) {
       index = todos + i - 1;
       expo[index] = factor1;
     }


     if (strcmp(bound,"finite") == 0) {
            expo_no_correcto =  check_expotents_finite(using_gamma,
                                                       gamma_couple,
                                                       Rc,
                                                       np,
                                                       mang,
                                                       expo,
                                                       todos);
     }
     if (np[todos] == 1)
     switch(deg) {
       case(1)  : k1[0]++; break;
     }
     if (np[todos] == 2)
     switch(deg) {
       case(1)  : k1[1]++; break;
       case(3)  : k1[2]++; break;
     }
     if (np[todos] == 3)
     switch(deg) {
       case(1)  : k1[3]++; break;
       case(3)  : k1[4]++; break;
       case(5)  : k1[5]++; break;
     }
     if (np[todos] == 4)
     switch(deg) {
       case(1)  : k1[6]++; break;
       case(3)  : k1[7]++; break;
       case(5)  : k1[8]++; break;
       case(7)  : k1[9]++; break;
     }
     if (np[todos] == 5)
     switch(deg) {
       case(1)  : k1[10]++; break;
       case(3)  : k1[11]++; break;
       case(5)  : k1[12]++; break;
       case(7)  : k1[13]++; break;
       case(9)  : k1[14]++; break;
     }
     if (np[todos] == 6)
     switch(deg) {
       case(1)  : k1[15]++; break;
       case(3)  : k1[16]++; break;
       case(5)  : k1[17]++; break;
       case(7)  : k1[18]++; break;
       case(9)  : k1[19]++; break;
       case(11) : k1[20]++; break;
     }
     if (np[todos] == 7)
     switch(deg) {
       case(1)  : k1[21]++; break;
       case(3)  : k1[22]++; break;
       case(5)  : k1[23]++; break;
       case(7)  : k1[24]++; break;
       case(9)  : k1[25]++; break;
       case(11) : k1[26]++; break;
       case(13) : k1[27]++; break;
     }
     if (np[todos] == 8)
     switch(deg) {
       case(1)  : k1[28]++; break;
       case(3)  : k1[29]++; break;
       case(5)  : k1[30]++; break;
       case(7)  : k1[31]++; break;
       case(9)  : k1[32]++; break;
       case(11) : k1[33]++; break;
       case(13) : k1[34]++; break;
       case(15) : k1[35]++; break;
     }


     if (expo_no_correcto == 1)
        exit(1);
        else
          todos = index + 1;

     }while(todos < nt);

 free(array_alpha);
 array_alpha = 0;
 free(array_beta);
 array_beta = 0;

 return (1);
}
