#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

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
                  double  gamma_nicp,
                  int     imprime,
                  int     plasma,
                  int     ident) {

 int i, j, k;
 int index, todos, ang_test, deg, test_scf;
 double step_work, energia;
 int bt;

 int expo_no_correcto,
     checking_expo;

     checking_expo    = 0;
     expo_no_correcto = 0;

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
                double  gamma_nicp,
                int     imprime,
                int plasma);

 int onsa, twsa, thsa, fosa, fisa, sisa, sesa, eisa;
 int twpa, thpa, fopa, fipa, sipa, sepa, eipa;
 int thda, foda, fida, sida, seda, eida;
 int fota, fita, sita, seta, eita;
 int figa, siga, sega, eiga;
 int siha, seha, eiha;
 int seia, eiia;
 int eija;

 int onsb, twsb, thsb, fosb, fisb, sisb, sesb, eisb;
 int twpb, thpb, fopb, fipb, sipb, sepb, eipb;
 int thdb, fodb, fidb, sidb, sedb, eidb;
 int fotb, fitb, sitb, setb, eitb;
 int figb, sigb, segb, eigb;
 int sihb, sehb, eihb;
 int seib, eiib;
 int eijb;

 double expo_opt[500];
 int    k1[40];
 double factor1, factor2;
 double expalpha_temp[60], expbeta_temp[60];
 double expalpha_opt[60], expbeta_opt[60];
 double expalpha[60], expbeta[60];

 double energia_0;
 int    step_work2;
 double gamma_nicp; // mike, nicp --> nonideal classical plasma

  for (j = 0; j < nt; j++)
     expo_opt[j] = 0.f;

  for (j = 0; j < 60; j++){
      expalpha_temp[j]    = 0.f;
      expbeta_temp[j]     = 0.f;
      expalpha_opt[j]     = 0.1;
      expbeta_opt[j]      = 1.1;
      expalpha[j]         = 0.f;
      expbeta[j]          = 0.f;
     }
    int steplim;
    int index_t;


     steplim = 100;

 if (ident == 1) {
     step_work        =  0.05;
     step_work2       =  0.025;
    do {
    index_t = 0;
    energia_0 = 1e9;
     for (j = 0; j < 60; j++) {
       expalpha[j]  = expalpha_opt[j];
       expbeta[j]  = expbeta_opt[j];
   //   printf("expoalpha %lf expobeta %lf  || %lf     %lf\n",(double) expalpha_opt[j], expbeta_opt[j], expalpha[j], expbeta[j]);
     }
     k         = 0;
     todos     = 0;
     factor1   = 0.f;
     factor2   = 0.f;
     printf("\nENTRA\n");

       for (onsa = 0 ; onsa < steplim; onsa++) {//1s
            expbeta[0]  = expbeta_opt[0];

        for (onsb = 0 ; onsb < steplim; onsb++) {//1s
                         k1[0] = 0;
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
                        if (k == 0)
                          factor1 =(double) expalpha[index_t];
                          else
                              factor1 =(double) expalpha[index_t]*pow(expbeta[index_t],(double) k);

                         printf("expoalpha %lf expobeta %lf total %lf (k %d|||deg %d) expo_opt %lf %lf\n",(double) expalpha[index_t],
                                                                                   (double) expbeta[index_t],
                                                                                   factor1, k, deg, expalpha_opt[index_t], expbeta_opt[index_t]);
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


                   if (expo_no_correcto == 1)
                     todos = nt;
                     else
                       todos = index + 1;

                  }while(todos < nt);


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
				   gamma_nicp,
                                   0,
                                   plasma);
                     else
                         test_scf = 1;

                   if (test_scf != 0)
                     energia = 1.e10;

                   if (energia < energia_0){
                        energia_0 = energia;
                        printf("\n HOLA1 energia %lf\n", energia);

                           expalpha_temp[0] = expalpha[0];
                           expbeta_temp[0] = expbeta[0];

                        for (i = 0; i < nt; i++)
                            expo_opt[i] = expo[i];
                     //mrb  for (i = 0; i < 2; i++)
                     //mrb     printf("e %lf  %lf", expalpha_temp, expbeta_temp);
                     }

              expbeta[0]  = expbeta[0] + step_work2;
     }//1s
              expalpha[0] = expalpha[0] + step_work;
     }//1s


    if ( expalpha_temp[0] - step_work < 0.f)
       expalpha_opt[0] = expalpha_temp[0];
      else
       expalpha_opt[0] = expalpha_temp[0] - step_work;
       

    if (expbeta_temp[0]  - step_work2 < 1.0)
       expbeta_opt[0]  = expbeta_temp[0];
      else
       expbeta_opt[0]  = expbeta_temp[0]  - step_work2;

    step_work = (double) step_work/10.f;
    step_work2 = (double) step_work2/10.f;


  }while (step_work > 0.001);
  


   for (i = 0; i < nt; i++)
      expo[i] = expo_opt[i];

 }

 if (ident == 2) {
  steplim = 10;
  step_work        =  1.0;
  step_work2       =  0.5;

 do {

      index_t = 0;
      energia_0 = 1e9;
     for (j = 0; j < 60; j++) {
       expalpha[j]  = expalpha_opt[j];
       expbeta[j]  = expbeta_opt[j];
   //   printf("expoalpha %lf expobeta %lf  || %lf     %lf\n",(double) expalpha_opt[j], expbeta_opt[j], expalpha[j], expbeta[j]);
     }
       k         = 0;
       todos     = 0;
       factor1   = 0.f;
       factor2   = 0.f;
       printf("\nENTRA\n");
   


        for (twsa = 0 ; twsa < steplim; twsa++) {//1s
           expbeta[1]  = expbeta_opt[1];

        for (twsb = 0 ; twsb < steplim; twsb++) {//1s
           expbeta[0]  = expbeta_opt[0];
           expalpha[0]  = expalpha_opt[0];

        for (onsa = 0 ; onsa < steplim; onsa++) {//1s
           expbeta[0]  = expbeta_opt[0];
        for (onsb = 0 ; onsb < steplim; onsb++) {//1s

                  for(j = 0; j <= 40; j++)
                         k1[j] = 0;
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
                                                                                   factor1, k, deg, expalpha_opt[index_t], expbeta_opt[index_t]);
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

                  }while(todos < nt);


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
				   gamma_nicp,
                                   0,
                                   plasma);
                     else
                         test_scf = 1;

                   if (test_scf != 0)
                     energia = 1.e10;

                   if (energia < energia_0){
                        energia_0 = energia;
                        printf("\n HOLA1 energia %lf\n", energia);

                        for (j = 0; j < 60; j++) {
                           expalpha_temp[j] = expalpha[j];
                           expbeta_temp[j] = expbeta[j];
                        }
                        for (i = 0; i < nt; i++)
                            expo_opt[i] = expo[i];
                     //mrb  for (i = 0; i < 2; i++)
                     //mrb     printf("e %lf  %lf", expalpha_temp, expbeta_temp);
                     }

              expbeta[0]  = expbeta[0] + step_work;
     }//1s
              expalpha[0] = expalpha[0] + step_work;
     }//1s
              expbeta[1]  = expbeta[1] + step_work;
     }//1s
              expalpha[1] = expalpha[1] + step_work;
     }//1s

    step_work = (double) step_work/2.f;
    step_work2 = (double) step_work2/2.f;
   
   for (i = 0; i < 60; i++) {

    if ( expalpha_temp[i] - step_work < 0.f)
       expalpha_opt[i] = expalpha_temp[i];
      else
       expalpha_opt[i] = expalpha_temp[i] - step_work;


    if (expbeta_temp[i]  - step_work2 < 1.0)
       expbeta_opt[i]  = expbeta_temp[i];
      else
       expbeta_opt[i]  = expbeta_temp[i]  - step_work2;
   }
  
 steplim = 20;

 }while(step_work > 0.001);


   for (i = 0; i < nt; i++)
      expo[i] = expo_opt[i];

 }

}
