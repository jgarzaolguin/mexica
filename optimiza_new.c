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
                  char   **save_dft,
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

 extern int checking_linear_dependece(int       nt,
                                      int      *bt,
                                      double   *expo,
                                      int      *mang,
                                      char     *nombre);


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
 double  expo_opt[800];
 double  expalpha_opt[36], expbeta_opt[36];
 double  expalpha[36], expbeta[36];
 double  rango_a[36], rango_b[36];
 double  valpha_opt[36];
 double  vbeta_opt[36];
 double  restric_a[36];
 int     k1[36];
 int     total_a[36], total_b[36];
 double  incr_a;
 double  incr_b;
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
 int unsigned cut;
 double range_inter_a;
 double range_inter_b;

 double alpha_add;
 double beta_add;
 alpha_add = 0.f;
 beta_add  = 0.f;

 if(z >= 2) {
    alpha_add = 0.2;
    beta_add  = 1.2;
 } 
 else {
    alpha_add = 0.01;
    beta_add  = 1.1;
 }

 int unsigned  pasos_a;
 pasos_a = 5000;
 int unsigned  pasos_b;
 pasos_b = 5000;

 for(j = 0; j < 36; j++) {
    rango_a[j]          = 3.5;
    rango_b[j]          = 3.2;
    total_a[j]          = (int) pasos_a*rango_a[j];
    total_b[j]          = (int) pasos_b*rango_b[j];
    expalpha_opt[j]     = 0.0;
    expbeta_opt[j]      = 0.0;
    valpha_opt[j]       = alpha_add;
    vbeta_opt[j]        = beta_add;
    expalpha[j]         = 0.f;
    expbeta[j]          = 0.f;
    k1[j]               = 0;
    restric_a[j]        = 0.f;
 }
 
 for(j = 0; j < nt; j++){
    expo_opt[j] = 0.f;
 }      
 char bound_pol[50];     //mike
 if(strcmp(basis,"GTOs") == 0) {
    sprintf(bound_pol,"%s", bound);
    strcpy(bound,"free");
 }

 if(strcmp(bound,"finite") == 0 || strcmp(bound,"dielectricc") == 0 || strcmp(bound,"polarization") == 0 || strcmp(bound,"parabolic") == 0  ||
    strcmp(bound,"confined") == 0) {     /* begins "if" all cases with gtos */

    todos = 0;
    index = 0;
    do{
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
       if(np[todos] == 1)
          switch(deg) {
                 case(1)  : index_t = 0; break;
          }
       if(np[todos] == 2)
          switch(deg) {
                 case(1)  : index_t = 1; break;
                 case(3)  : index_t = 2; break;
          }
       if(np[todos] == 3)
          switch(deg) {
                 case(1)  : index_t = 3; break;
                 case(3)  : index_t = 4; break;
                 case(5)  : index_t = 5; break;
          }
       if(np[todos] == 4)
          switch(deg) {
                 case(1)  : index_t = 6; break;
                 case(3)  : index_t = 7; break;
                 case(5)  : index_t = 8; break;
                 case(7)  : index_t = 9; break;
          }
       if(np[todos] == 5)
          switch(deg) {
                 case(1)  : index_t = 10; break;
                 case(3)  : index_t = 11; break;
                 case(5)  : index_t = 12; break;
                 case(7)  : index_t = 13; break;
                 case(9)  : index_t = 14; break;
          }
       if(np[todos] == 6)
           switch(deg) {
                 case(1)  : index_t = 15; break;
                 case(3)  : index_t = 16; break;
                 case(5)  : index_t = 17; break;
                 case(7)  : index_t = 18; break;
                 case(9)  : index_t = 19; break;
                 case(11) : index_t = 20; break;
           }
       if(np[todos] == 7)
           switch(deg) {
                 case(1)  : index_t = 21; break;
                 case(3)  : index_t = 22; break;
                 case(5)  : index_t = 23; break;
                 case(7)  : index_t = 24; break;
                 case(9)  : index_t = 25; break;
                 case(11) : index_t = 26; break;
                 case(13) : index_t = 27; break;
           }
       if(np[todos] == 8)
           switch(deg) {
                 case(1)  : index_t = 28; break;
                 case(3)  : index_t = 29; break;
                 case(5)  : index_t = 30; break;
                 case(7)  : index_t = 31; break;
                 case(9)  : index_t = 32; break;
                 case(11) : index_t = 33; break;
                 case(13) : index_t = 34; break;
                 case(15) : index_t = 35; break;
           }
      
       if(strcmp(using_gamma,"YES") == 0){
          factor1  = gamma_couple/((1.f - gamma_couple)*Rc);
          factor2  = (double) (mang[todos] + np[todos])/Rc;
          valpha_opt[index_t]  = 0.01 + factor1 + factor2; 
          //restric_a[index_t]   = valpha_opt[index_t];
          rango_a[index_t] = valpha_opt[index_t] + 3.5;
       } 
       else {
          factor2 = (double) (np[todos] + mang[todos]);
          valpha_opt[index_t]  = 0.01 + factor2/Rc; 
          //restric_a[index_t]   = valpha_opt[index_t];
          rango_a[index_t] = valpha_opt[index_t] + 3.5;
       }

       for(i = 1; i <= deg; i++) {
          index = todos + i - 1;
       }
       todos = index + 1;
    } while(todos < nt);
    for(j = 0; j < 36; j++) {
       rango_b[j]          = 3.2;
       total_a[j]          = (int) pasos_a*rango_a[j];
       total_b[j]          = (int) pasos_b*rango_b[j];
       expalpha_opt[j]     = 0.0;
       expbeta_opt[j]      = 0.0;
       expalpha[j]         = 0.f;
       expbeta[j]          = 0.f;
       k1[j]               = 0;
    }
 }     /* ends "if" for all cases with gtos */ 
 else 

 incr_a   = (double) rango_a[0]/total_a[0];
 incr_b   = (double) rango_b[0]/total_b[0];

//mrb1 for (j = 0; j < 36; j++) 
//mrb1   printf("\n KONDA %lf  %lf  %d %f\n", alpha_value(total_a[j] - 1, incr, valpha_opt[j]), beta_value(total_b[j] - 1, incr, vbeta_opt[j]), total_b[j], incr);

 srand((unsigned) time(NULL));

 energia_0  = 1e9;
 kato_0     = 1.1;
 ident      = 0;
 ident_test = 0;

 iteraciones = 1000;     // mike


 w = 0;
// for (w = 0; w < 1; w++) {
      for (q = 0; q < iteraciones; q++) {
          index_t   = 0;
          k         = 0;
          factor1   = 0.f;
          factor2   = 0.f;
          printf("\n STEP (%d) %d\n", w, q);

           for (j = 0; j < 36; j++) {
              rand_1     = rand()%(total_a[j]);
              rand_2     = rand()%(total_b[j]);
              expalpha[j]  = alpha_value(rand_1, incr_a, valpha_opt[j]);
              expbeta[j]   = beta_value(rand_2, incr_b, vbeta_opt[j]);
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


               if(strcmp(bound,"finite") == 0 || strcmp(bound,"dielectricc") == 0 || strcmp(bound,"polarization") == 0 || strcmp(bound,"parabolic") == 0 ||
                          strcmp(bound,"confined") == 0) {
                      expo_no_correcto =  check_expotents_finite(using_gamma,
                                                                 gamma_couple,
                                                                 Rc,
                                                                 np,
                                                                 mang,
                                                                 expo,
                                                                 todos);
               }
//mrb1              if (expo_no_correcto == 1) {
//mrb1                for (j = 0; j < 36; j++) {
//mrb1                 rango_a[j]          = 3.3;
//mrb1                 rango_b[j]          = 3.3;
//mrb1                 total_a[j]          = (int) pasos*rango_a[j];
//mrb1                 total_b[j]          = (int) pasos*rango_b[j];
//mrb1                 valpha_opt[j]       = alpha_add;
//mrb1                 vbeta_opt[j]        = beta_add;
//mrb1                }
//mrb1                incr_a   = (double) rango_a[0]/total_a[0];
//mrb1                incr_b   = (double) rango_b[0]/total_b[0];
//mrb1              }

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


           expo_no_correcto = checking_linear_dependece(nt,
                                                        &bt,
                                                        expo,
                                                        mang,
                                                        nombre);
           expo_no_correcto = 0;


           if (expo_no_correcto == 0){           //mike
            if (strcmp(basis,"GTOs") == 0) {
              sprintf(bound,"%s", bound_pol);
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
                          &cusp_kato);
           } 
           else
              test_scf = 1;

           if (test_scf != 0) {
              energia   = 1.e10;
              diff_kato_0 = 1000.f;
           } 
           else {  /* mike */
              if(strcmp(basis, "STOs") == 0){           
                 diff_kato_0 = fabs(cusp_kato - 1.f);
              }
              else{
                 diff_kato_0 = 0.f;
              }
           }

           if(energia < energia_0 && diff_kato_0 < 1.e-02 /*kato_0*/) {
              energia_0 = energia;
               //kato_0 = diff_kato_0;
              printf("\nHOLA Energy %6.12lf Kato %lf\n", energia, cusp_kato);
               for (j = 0; j < 36; j++) {
                   expalpha_opt[j] = expalpha[j];
                   expbeta_opt[j]  = expbeta[j];
               }
              for (j = 0; j < nt; j++)
                 expo_opt[j] = expo[j];
            }
     if (strcmp(basis,"GTOs") == 0) {
      sprintf(bound_pol,"%s", bound);
      strcpy(bound,"free");
     }

   } //primer for q

   
//mrb if (w == 1) break;

//identifica la mejor condicion de kato cercana a 1 y mejor energia
//mrb if (w == 1) { 
//mrb    pasos = pasos*2;
//mrb    cut = w;
//mrb }

//mrb   for (j = 0; j < 36; j++) {
     
//mrb if (expalpha_opt[j] - (double) 1.f/(1 + cut) > restric_a[j]) {
//mrb   range_inter_a = (double) 1.f/(1 + cut);
//mrb   valpha_opt[j] = expalpha_opt[j] - range_inter_a;
//mrb }
//mrb    else
//mrb if (expalpha_opt[j] - (double) 0.9/(1 + cut) > restric_a[j]) {
//mrb   range_inter_a = (double) 0.9/(1 + cut);
//mrb   valpha_opt[j] = expalpha_opt[j] - range_inter_a;
//mrb }
//mrb   else
//mrb if (expalpha_opt[j] - (double) 0.8/(1 + cut) > restric_a[j]) {
//mrb   range_inter_a = (double) 0.8/(1 + cut);
//mrb   valpha_opt[j] = expalpha_opt[j] - range_inter_a;
//mrb }
//mrb   else
//mrb if (expalpha_opt[j] - (double) 0.7/(1 + cut) > restric_a[j]) {
//mrb   range_inter_a = (double) 0.7/(1 + cut);
//mrb   valpha_opt[j] = expalpha_opt[j] - range_inter_a;
//mrb }
//mrb   else
//mrb if (expalpha_opt[j] - (double) 0.6/(1 + cut) > restric_a[j]) {
//mrb   range_inter_a = (double) 0.6/(1 + cut);
//mrb   valpha_opt[j] = expalpha_opt[j] - range_inter_a;
//mrb }
//mrb   else
//mrb  if (expalpha_opt[j] - (double) 0.5/(1 + cut) > restric_a[j]) {
//mrb   range_inter_a = (double) 0.5/(1 + cut);
//mrb    valpha_opt[j] = expalpha_opt[j] - range_inter_a;
//mrb }
//mrb    else
//mrb  if (expalpha_opt[j] - (double) 0.4/(1 + cut) > restric_a[j]) {
//mrb   range_inter_a = (double) 0.4/(1 + cut);
//mrb    valpha_opt[j] = expalpha_opt[j] - range_inter_a;
//mrb }
//mrb    else
//mrb  if (expalpha_opt[j] - (double) 0.3/(1 + cut) > restric_a[j]) {
//mrb   range_inter_a = (double) 0.3/(1 + cut);
//mrb    valpha_opt[j] = expalpha_opt[j] - range_inter_a;
//mrb }
//mrb    else
//mrb  if (expalpha_opt[j] - (double) 0.2/(1 + cut) > restric_a[j]) {
//mrb   range_inter_a = (double) 0.2/(1 + cut);
//mrb    valpha_opt[j] = expalpha_opt[j] - range_inter_a; 
//mrb }
//mrb    else
//mrb  if (expalpha_opt[j] - (double) 0.1/(1 + cut) > restric_a[j]) {
//mrb   range_inter_a = (double) 0.1/(1 + cut);
//mrb    valpha_opt[j] = expalpha_opt[j] - range_inter_a;
//mrb }
//mrb    else {
//mrb     range_inter_a = 1.0;
//mrb     valpha_opt[j] = expalpha_opt[j];
//mrb     }
   
   
//mrb if (expbeta_opt[j]  - (double) 1.f/(1 + cut) > 1.1) {
//mrb   range_inter_b = (double) 1.f/(1 + cut);
//mrb    vbeta_opt[j]  = expbeta_opt[j]  - range_inter_b;
//mrb }
//mrb     else
//mrb if (expbeta_opt[j]  - (double) 0.9/(1 + cut) > 1.1) {
//mrb   range_inter_b = (double) 0.9/(1 + cut);
//mrb    vbeta_opt[j]  = expbeta_opt[j]  - range_inter_b;
//mrb }
//mrb     else
//mrb if (expbeta_opt[j]  - (double) 0.8/(1 + cut) > 1.1) {
//mrb   range_inter_b = (double) 0.8/(1 + cut);
//mrb    vbeta_opt[j]  = expbeta_opt[j]  - range_inter_b;
//mrb }
//mrb     else
//mrb if (expbeta_opt[j]  - (double) 0.7/(1 + cut) > 1.1) {
//mrb   range_inter_b = (double) 0.7/(1 + cut);
//mrb    vbeta_opt[j]  = expbeta_opt[j]  - range_inter_b;
//mrb }
//mrb     else
//mrb if (expbeta_opt[j]  - (double) 0.6/(1 + cut) > 1.1) {
//mrb   range_inter_b = (double) 0.6/(1 + cut);
//mrb    vbeta_opt[j]  = expbeta_opt[j]  - range_inter_b;
//mrb }
//mrb     else
//mrb  if (expbeta_opt[j]  - (double) 0.5/(1 + cut) > 1.1) {
//mrb   range_inter_b = (double) 0.5/(1 + cut);
//mrb     vbeta_opt[j]  = expbeta_opt[j]  - range_inter_b; 
//mrb }
//mrb      else
//mrb  if (expbeta_opt[j]  - (double) 0.4/(1 + cut) > 1.1) {
//mrb   range_inter_b = (double) 0.4/(1 + cut);
//mrb     vbeta_opt[j]  = expbeta_opt[j]  - range_inter_b;
//mrb }
//mrb      else
//mrb  if (expbeta_opt[j]  - (double) 0.3/(1 + cut) > 1.1) {
//mrb   range_inter_b = (double) 0.3/(1 + cut);
//mrb     vbeta_opt[j]  = expbeta_opt[j]  - range_inter_b;
//mrb }
//mrb      else
//mrb  if (expbeta_opt[j]  - (double) 0.2/(1 + cut) > 1.1) {
//mrb   range_inter_b = (double) 0.2/(1 + cut);
//mrb     vbeta_opt[j]  = expbeta_opt[j]  - range_inter_b;
//mrb }
//mrb      else
//mrb  if (expbeta_opt[j]  - (double) 0.1/(1 + cut) > 1.1) {
//mrb   range_inter_b = (double) 0.1/(1 + cut);
//mrb     vbeta_opt[j]  = expbeta_opt[j]  - range_inter_b;
//mrb }
//mrb      else {
//mrb        range_inter_b = 1.0;
//mrb        vbeta_opt[j]  = expbeta_opt[j];
//mrb       }
   
    //mrb rango_a[j] = expalpha_opt[j] + (double) 1.0f/(1 + cut);
    //mrb rango_b[j] = expbeta_opt[j] + (double) 1.0f/(1 + cut);
//mrb     rango_a[j] = expalpha_opt[j] + range_inter_a;
//mrb     rango_b[j] = expbeta_opt[j] + range_inter_b;
   
//mrb     total_a[j] = (int) pasos*rango_a[j];
//mrb     total_b[j] = (int) pasos*rango_b[j];
   
//mrb     incr_a   = (double) rango_a[0]/total_a[0];
//mrb     incr_b   = (double) rango_b[0]/total_b[0];
 
   //mrb1 printf("\nHOLA2 %f %f %f %f %f %f %d %d %f\n", expalpha_opt[j], expbeta_opt[j], valpha_opt[j], vbeta_opt[j], rango_a[j], rango_b[j], total_a[j], total_b[j], incr);
//mrb  }
//mrb  iteraciones = 200; //iteraciones/4 + iteraciones/2;
 

// } //segundo for w

 for (j = 0; j < 36; j++) {
   if (k1[j] != 0) {
     printf(" %d FINAL alpha %lf  beta %lf\n", j, expalpha_opt[j], expbeta_opt[j]);
   }
 }

 for (j = 0; j < nt; j++)
     expo[j] = expo_opt[j];

 // mike
   if (strcmp(basis,"GTOs") == 0) {
      sprintf(bound,"%s", bound_pol);
  }



 return (1);
}
