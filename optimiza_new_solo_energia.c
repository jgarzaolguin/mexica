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
 double  *array_alpha, *array_beta;
 double  energia_0;
 int     index_t;
 int unsigned q, w;
 int    ident;
 int    rand_1;
 int    rand_2;
 int    iteraciones;
 int    identify;
 double diff_kato_0;
 double diff_kato_1;

  array_alpha     = (double*)malloc(2000*sizeof(double));
 if (array_alpha == NULL) {
   free(array_alpha);
   array_alpha = 0;
   printf("Error en arreglo alpha, insuficiente RAM\n");
   exit(1);
 }

  array_beta     = (double*)malloc(2000*sizeof(double));
 if (array_beta == NULL) {
   free(array_beta);
   array_beta = 0;
   printf("Error en arreglo beta, insuficiente RAM\n");
   exit(1);
 }



 pasos = 100000;
 
 for (j = 0; j < 36; j++) {
  rango_a[j]          = 5.f;
  rango_b[j]          = 5.f;
  total_a[j]          = (int) pasos*rango_a[j];
  total_b[j]          = (int) pasos*rango_b[j];
  expalpha_opt[j]     = 0.0;
  expbeta_opt[j]      = 0.0;
  valpha_opt[j]       = 0.01;
  vbeta_opt[j]        = 1.1;
  expalpha[j]         = 0.f;
  expbeta[j]          = 0.f;
  k1[j]               = 0;
 }

 incr   = (double) rango_a[0]/total_a[0];

//mrb1 for (j = 0; j < 36; j++) 
//mrb1   printf("\n KONDA %lf  %lf  %d %f\n", alpha_value(total_a[j] - 1, incr, valpha_opt[j]), beta_value(total_b[j] - 1, incr, vbeta_opt[j]), total_b[j], incr);

 srand((unsigned) time(NULL));


 energia_0 = 1e9;
 ident     = 0;

 iteraciones = 10000;

 for (w = 0; w < 5; w++) {
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
                           &cusp_kato);
              else
                  test_scf = 1;

            if (test_scf != 0)
              energia = 1.e10;

            if (energia < energia_0) {
               energia_0 = energia;
               array_cusp_kato[ident] = cusp_kato;
               printf("\n HOLA1 energy %25.12lf %d %lf\n", energia, ident, array_cusp_kato[ident]);
               for (j = 0; j < 36; j++) {
                    array_alpha[ident*35 + j] = expalpha[j];
                    array_beta[ident*35 + j]  =  expbeta[j];
                    expalpha_opt[j]           = expalpha[j];
                    expbeta_opt[j]            = expbeta[j];
                    //mrb1 printf("\nEXPO_OPT a %f b %f || a %f b %f\n", expalpha_opt[j], expbeta_opt[j], array_alpha[ident*35 + j], array_beta[ident*35 + j]);
               }
              ident++;
            }
   } //primer for q

   pasos = pasos*10;
   
   for (j = 0; j < 36; j++) {
     
//mrb if (expalpha_opt[j] - (double) 1.f/(1 + w) > 0.f)
//mrb   valpha_opt[j] = expalpha_opt[j] - (double) 1.f/(1 + w);
//mrb    else
//mrb if (expalpha_opt[j] - (double) 0.9/(1 + w) > 0.f)
//mrb   valpha_opt[j] = expalpha_opt[j] - (double) 0.9/(1 + w);
//mrb   else
//mrb if (expalpha_opt[j] - (double) 0.8/(1 + w) > 0.f)
//mrb   valpha_opt[j] = expalpha_opt[j] - (double) 0.8/(1 + w);
//mrb   else
//mrb if (expalpha_opt[j] - (double) 0.7/(1 + w) > 0.f)
//mrb   valpha_opt[j] = expalpha_opt[j] - (double) 0.7/(1 + w);
//mrb   else
//mrb if (expalpha_opt[j] - (double) 0.6/(1 + w) > 0.f)
//mrb   valpha_opt[j] = expalpha_opt[j] - (double) 0.6/(1 + w);
//mrb   else
   if (expalpha_opt[j] - (double) 0.5/(1 + w) > 0.f)
     valpha_opt[j] = expalpha_opt[j] - (double) 0.5/(1 + w);
     else
   if (expalpha_opt[j] - (double) 0.4/(1 + w) > 0.f)
     valpha_opt[j] = expalpha_opt[j] - (double) 0.4/(1 + w);
     else
   if (expalpha_opt[j] - (double) 0.3/(1 + w) > 0.f)
     valpha_opt[j] = expalpha_opt[j] - (double) 0.3/(1 + w);
     else
   if (expalpha_opt[j] - (double) 0.2/(1 + w) > 0.f)
     valpha_opt[j] = expalpha_opt[j] - (double) 0.2/(1 + w);
     else
   if (expalpha_opt[j] - (double) 0.1/(1 + w) > 0.f)
     valpha_opt[j] = expalpha_opt[j] - (double) 0.1/(1 + w);
     else
     valpha_opt[j] = expalpha_opt[j];
   
   
//mrb if (expbeta_opt[j]  - (double) 1.f/(1 + w) > 1.1)
//mrb    vbeta_opt[j]  = expbeta_opt[j]  - (double) 1.f/(1 + w);
//mrb     else
//mrb if (expbeta_opt[j]  - (double) 0.9/(1 + w) > 1.1)
//mrb    vbeta_opt[j]  = expbeta_opt[j]  - (double) 0.9/(1 + w);
//mrb     else
//mrb if (expbeta_opt[j]  - (double) 0.8/(1 + w) > 1.1)
//mrb    vbeta_opt[j]  = expbeta_opt[j]  - (double) 0.8/(1 + w);
//mrb     else
//mrb if (expbeta_opt[j]  - (double) 0.7/(1 + w) > 1.1)
//mrb    vbeta_opt[j]  = expbeta_opt[j]  - (double) 0.7/(1 + w);
//mrb     else
//mrb if (expbeta_opt[j]  - (double) 0.6/(1 + w) > 1.1)
//mrb    vbeta_opt[j]  = expbeta_opt[j]  - (double) 0.6/(1 + w);
//mrb     else
   if (expbeta_opt[j]  - (double) 0.5/(1 + w) > 1.1)
      vbeta_opt[j]  = expbeta_opt[j]  - (double) 0.5/(1 + w);
       else
   if (expbeta_opt[j]  - (double) 0.4/(1 + w) > 1.1)
      vbeta_opt[j]  = expbeta_opt[j]  - (double) 0.4/(1 + w);
       else
   if (expbeta_opt[j]  - (double) 0.3/(1 + w) > 1.1)
      vbeta_opt[j]  = expbeta_opt[j]  - (double) 0.3/(1 + w);
       else
   if (expbeta_opt[j]  - (double) 0.2/(1 + w) > 1.1)
      vbeta_opt[j]  = expbeta_opt[j]  - (double) 0.2/(1 + w);
       else
   if (expbeta_opt[j]  - (double) 0.1/(1 + w) > 1.1)
      vbeta_opt[j]  = expbeta_opt[j]  - (double) 0.1/(1 + w);
       else
      vbeta_opt[j]  = expbeta_opt[j];
   
//mrb     rango_a[j] = expalpha_opt[j] + (double) 1.f/(1 + w);
//mrb     rango_b[j] = expbeta_opt[j] + (double) 1.f/(1 + w);
      rango_a[j] = expalpha_opt[j] + (double) 0.5/(1 + w);
      rango_b[j] = expbeta_opt[j] + (double) 0.5/(1 + w);
   
     total_a[j] = (int) pasos*rango_a[j];
     total_b[j] = (int) pasos*rango_b[j];
   
     incr   = (double) rango_a[0]/total_a[0];
 
   //mrb1 printf("\nHOLA2 %f %f %f %f %f %f %d %d %f\n", expalpha_opt[j], expbeta_opt[j], valpha_opt[j], vbeta_opt[j], rango_a[j], rango_b[j], total_a[j], total_b[j], incr);
  }
  iteraciones = 5000;//iteraciones/4 + iteraciones/2;
 

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
