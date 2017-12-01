#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
       

 
int ep2_cpu(char   *espin,
                   int     nt,
                   int     elecalfa,
                   int     elecbeta,
                   double *arreglo,
                   double *vectalfa,
                   double *vectbeta,
                   double *valalfa,
                   double *valbeta,
                   int     orbital,
                   double *arreglo_mo,
                   double *arreglo_mo_ab,
                   double *total_energy)
 {//label 1

 extern double ao_to_mo_cpu(int     nt,
                            int     i,
                            int     a,
                            int     j,
                            int     b,
                            double *vectalfa,
                            double *vectbeta,
                            double *integrales_bie);


  int cuad, cuart, cubo;
  int i, j, a, b, cont;

  double  time_1, time_4, suma,  suma1, factor1, factor2, diff;
  double  dfactor1, dsuma1, die, func, dfunc, omega_tmp, omega;
  int     ncis;
  double  ie;
  int     interno;

  double au_to_eV = (double) 27.21138;

  time_1 = time (NULL);

  cuad = nt*nt;
  cubo = nt*cuad;
  cuart = cuad*cuad;

  for (j = 0; j < cuart; j++) {
    arreglo_mo[j] = (double) 0.;
    arreglo_mo_ab[j] = (double) 0.;
     }

  if (strcmp(espin,"alpha") == 0) {//label alpha

    for (a = elecalfa; a < nt; a++)
      for (i = 0; i < elecalfa; i++)
        for (j = i + 1; j < elecalfa; j++) {//label 2
   
            suma             = ao_to_mo_cpu(nt,
                                            orbital,
                                            j,
                                            a,
                                            i,
                                            vectalfa,
                                            vectbeta,
                                            arreglo);
   
            ncis             = orbital + j*nt + a*cuad + i*cubo;
            arreglo_mo[ncis] = suma;
   
            suma             = 0.f;
   
            suma             = ao_to_mo_cpu(nt,
                                            orbital,
                                            i,
                                            a,
                                            j,
                                            vectalfa,
                                            vectbeta,
                                            arreglo);
            
            ncis             = orbital + i*nt + a*cuad + j*cubo;
            arreglo_mo[ncis] = suma;
        }//label 2
           
    suma = 0.f;
   
    for (a = elecbeta; a < nt; a++)
      for (i = 0; i < elecalfa; i++)
        for (j = 0; j < elecbeta; j++) {//label 3
   
            suma                = ao_to_mo_cpu(nt,
                                               orbital,
                                               i,   
                                               a,   
                                               j,   
                                               vectalfa,
                                               vectbeta,
                                               arreglo);
   
            ncis                = orbital + i*nt + a*cuad + j*cubo;
            arreglo_mo_ab[ncis] = suma;
           }//label 3
   
    suma = 0.f;
   
     for (i = 0; i < elecalfa; i++) 
       for (a = elecalfa ; a < nt; a++) 
         for (b = a + 1; b < nt; b++) { //label 4
   
            suma             = ao_to_mo_cpu(nt,
                                            orbital,
                                            a,   
                                            i,   
                                            b,   
                                            vectalfa,
                                            vectbeta,
                                            arreglo);
   
            ncis             = orbital + a*nt + i*cuad + b*cubo;
            arreglo_mo[ncis] = suma;
   
            suma             = 0.f;
   
            suma             = ao_to_mo_cpu(nt,
                                            orbital,
                                            b,
                                            i,
                                            a,
                                            vectalfa,
                                            vectbeta,
                                            arreglo);
   
            ncis             = orbital + b*nt + i*cuad + a*cubo;
            arreglo_mo[ncis] = suma;
   
            }//label 4
   
    suma = 0.f;
   
     for (i = 0; i < elecbeta; i++) 
       for (a = elecalfa ; a < nt; a++) 
         for (b = elecbeta; b < nt; b++) {//label 5
   
            suma                = ao_to_mo_cpu(nt,
                                               orbital,
                                               a,
                                               i,
                                               b,
                                               vectalfa,
                                               vectbeta,
                                               arreglo);
            
            ncis                = orbital + a*nt + i*cuad + b*cubo;
            arreglo_mo_ab[ncis] = suma;
   
          
            }//label 5
   
   
    omega = valalfa[orbital];
    cont = 0;
    do {
       omega_tmp = omega;
       suma1 = 0.f;
       dsuma1 = 0.f;
       for (a = elecalfa; a < nt; a++)
         for (i = 0; i < elecalfa; i++)
           for (j = i + 1; j < elecalfa; j++) {
             ncis = orbital + j*nt + a*nt*nt + i*nt*nt*nt;
             factor1 = arreglo_mo[ncis];
             ncis = orbital + i*nt + a*nt*nt + j*nt*nt*nt;
             factor2 = arreglo_mo[ncis];
             factor1 = (factor2-factor1)/(omega_tmp + valalfa[a] - valalfa[i] - valalfa[j]);
             factor1 = factor1*factor1*(omega_tmp + valalfa[a] - valalfa[i] - valalfa[j]);
             dfactor1 = factor1/(omega_tmp + valalfa[a] - valalfa[i] - valalfa[j]);
             suma1 = suma1 + factor1;
             dsuma1 = dsuma1 - dfactor1;
           }
   
       for (a = elecbeta; a < nt; a++)
         for (i = 0; i < elecalfa; i++)
           for (j = 0; j < elecbeta; j++) {
             ncis = orbital + i*nt + a*nt*nt + j*nt*nt*nt;
             factor1 = arreglo_mo_ab[ncis];
             factor1 = factor1*factor1/(omega_tmp + valbeta[a] - valalfa[i] - valbeta[j]);
             dfactor1 = factor1*(1/(omega_tmp + valbeta[a] - valalfa[i] - valbeta[j]));
             suma1 = suma1 + factor1;
             dsuma1 = dsuma1 - dfactor1;
           }
   
       ie = suma1;
       die = dsuma1;
   
       suma1 = 0.f;
       dsuma1 = 0.f;
       for (i = 0; i < elecalfa; i++) 
         for (a = elecalfa; a < nt; a++) 
           for (b = a + 1; b < nt; b++) {
             ncis = orbital + b*nt + i*nt*nt + a*nt*nt*nt;
             factor1 = arreglo_mo[ncis];
             ncis = orbital + a*nt + i*nt*nt + b*nt*nt*nt;
             factor2 = arreglo_mo[ncis];
             factor1 = (factor2-factor1)/(omega_tmp + valalfa[i] - valalfa[a] - valalfa[b]);
             factor1 = factor1*factor1*(omega_tmp + valalfa[i] - valalfa[a] - valalfa[b]);
             dfactor1 = factor1/(omega_tmp + valalfa[i] - valalfa[a] - valalfa[b]);
             suma1 = suma1 + factor1;
             dsuma1 = dsuma1 - dfactor1;
        }
   
       for (i = 0; i < elecbeta; i++) 
        for (a = elecalfa ; a < nt; a++) 
          for (b = elecbeta; b < nt; b++) {
             ncis = orbital + a*nt + i*nt*nt + b*nt*nt*nt;
             factor1 = arreglo_mo_ab[ncis];
             factor1 = factor1*factor1/(omega_tmp + valbeta[i] - valalfa[a] - valbeta[b]);
             dfactor1 = factor1*(1/(omega_tmp + valbeta[i] - valalfa[a] - valbeta[b]));
             suma1 = suma1 + factor1;
             dsuma1 = dsuma1 - dfactor1;
        }
   
       ie = ie + suma1;
       die = die + dsuma1;
   
       func = omega_tmp - valalfa[orbital] - ie;
       dfunc = 1.f - die;
   
       omega = omega_tmp - func/dfunc;
   
       diff = fabs(omega - omega_tmp);
       cont++;
       printf("IP, iter  = %f, %d\n", -omega, cont);
   
    } while (diff > 1.e-6 && cont < 100);
   
    printf("Koopmans theorem     = %10.4f (%10.4f eV)\n", -valalfa[orbital], -valalfa[orbital]*au_to_eV);
    printf("Ionization potential = %10.4f (%10.4f eV) \n", -omega, -omega*au_to_eV);
   
    *total_energy = omega;
   
    time_4 = time (NULL);
    printf("EP2 time            = %5.2f s\n", time_4 - time_1);
   
   }//label alpha
  else {//label beta


      for (a = elecbeta; a < nt; a++)
        for (i = 0; i < elecbeta; i++)
          for (j = i + 1; j < elecbeta; j++) {//label 6
     
              suma             = ao_to_mo_cpu(nt,
                                              orbital,
                                              j,
                                              a,
                                              i,
                                              vectalfa,
                                              vectbeta,
                                              arreglo);
     
              ncis             = orbital + j*nt + a*cuad + i*cubo;
              arreglo_mo[ncis] = suma;
     
              suma             = 0.f;
              suma             = ao_to_mo_cpu(nt,
                                              orbital,
                                              i,
                                              a,
                                              j,
                                              vectalfa,
                                              vectbeta,
                                              arreglo);
              
            
              ncis             = orbital + i*nt + a*cuad + j*cubo;
              arreglo_mo[ncis] = suma;
          }//label 6
     
      suma = 0.f;
      for (a = elecalfa; a < nt; a++)
        for (i = 0; i < elecbeta; i++)
          for (j = 0; j < elecalfa; j++) {//label 7
              suma                 = ao_to_mo_cpu(nt,
                                                  orbital,
                                                  i,
                                                  a,
                                                  j,
                                                  vectalfa,
                                                  vectbeta,
                                                  arreglo);
     
               ncis                = orbital + i*nt + a*cuad + j*cubo;
               arreglo_mo_ab[ncis] = suma;
          }//label 7
       
       suma = 0.f;
       for (i = 0; i < elecbeta; i++) 
         for (a = elecbeta ; a < nt; a++) 
           for (b = a + 1; b < nt; b++) {//label 8
     
              suma             = ao_to_mo_cpu(nt,
                                              orbital,
                                              a,
                                              i,
                                              b,
                                              vectalfa,
                                              vectbeta,
                                              arreglo);
     
              ncis             = orbital + a*nt + i*cuad + b*cubo;
              arreglo_mo[ncis] = suma;
     
              suma = 0.f;
              suma             = ao_to_mo_cpu(nt,
                                              orbital,
                                              b,
                                              i,
                                              a,
                                              vectalfa,
                                              vectbeta,
                                              arreglo);
            
     
              ncis             = orbital + b*nt + i*cuad + a*cubo;
              arreglo_mo[ncis] = suma;
        } //label 8
     
       suma = 0.f;
       for (i = 0; i < elecalfa; i++) 
         for (a = elecbeta ; a < nt; a++) 
           for (b = elecalfa; b < nt; b++) {//label 9
     
     
              suma                = ao_to_mo_cpu(nt,
                                                 orbital,
                                                 a,
                                                 i,
                                                 b,
                                                 vectalfa,
                                                 vectbeta,
                                                 arreglo);
     
     
              ncis                = orbital + a*nt + i*cuad + b*cubo;
              arreglo_mo_ab[ncis] = suma;
     
        }//label 9        

      omega = valbeta[orbital];
      cont = 0;
      do {
         omega_tmp = omega;
         suma1 = 0.f;
         dsuma1 = 0.f;
         for (a = elecbeta; a < nt; a++)
           for (i = 0; i < elecbeta; i++)
             for (j = i + 1; j < elecbeta; j++) {
               ncis = orbital + j*nt + a*nt*nt + i*nt*nt*nt;
               factor1 = arreglo_mo[ncis];
               ncis = orbital + i*nt + a*nt*nt + j*nt*nt*nt;
               factor2 = arreglo_mo[ncis];
               factor1 = (factor2-factor1)/(omega_tmp + valbeta[a] - valbeta[i] - valbeta[j]);
               factor1 = factor1*factor1*(omega_tmp + valbeta[a] - valbeta[i] - valbeta[j]);
               dfactor1 = factor1/(omega_tmp + valbeta[a] - valbeta[i] - valbeta[j]);
               suma1 = suma1 + factor1;
               dsuma1 = dsuma1 - dfactor1;
             }
     
         for (a = elecalfa; a < nt; a++)
           for (i = 0; i < elecbeta; i++)
             for (j = 0; j < elecalfa; j++) {
               ncis = orbital + i*nt + a*nt*nt + j*nt*nt*nt;
               factor1 = arreglo_mo_ab[ncis];
               factor1 = factor1*factor1/(omega_tmp + valalfa[a] - valbeta[i] - valalfa[j]);
               dfactor1 = factor1*(1/(omega_tmp + valalfa[a] - valbeta[i] - valalfa[j]));
               suma1 = suma1 + factor1;
               dsuma1 = dsuma1 - dfactor1;
             }
     
         ie = suma1;
         die = dsuma1;
     
         suma1 = 0.f;
         dsuma1 = 0.f;
         for (i = 0; i < elecbeta; i++) 
           for (a = elecbeta; a < nt; a++) 
             for (b = a + 1; b < nt; b++) {
               ncis = orbital + b*nt + i*nt*nt + a*nt*nt*nt;
               factor1 = arreglo_mo[ncis];
               ncis = orbital + a*nt + i*nt*nt + b*nt*nt*nt;
               factor2 = arreglo_mo[ncis];
               factor1 = (factor2-factor1)/(omega_tmp + valbeta[i] - valbeta[a] - valbeta[b]);
               factor1 = factor1*factor1*(omega_tmp + valbeta[i] - valbeta[a] - valbeta[b]);
               dfactor1 = factor1/(omega_tmp + valbeta[i] - valbeta[a] - valbeta[b]);
               suma1 = suma1 + factor1;
               dsuma1 = dsuma1 - dfactor1;
          }
     
         for (i = 0; i < elecalfa; i++) 
          for (a = elecbeta ; a < nt; a++) 
            for (b = elecalfa; b < nt; b++) {
               ncis = orbital + a*nt + i*nt*nt + b*nt*nt*nt;
               factor1 = arreglo_mo_ab[ncis];
               factor1 = factor1*factor1/(omega_tmp + valalfa[i] - valbeta[a] - valalfa[b]);
               dfactor1 = factor1*(1/(omega_tmp + valalfa[i] - valbeta[a] - valalfa[b]));
               suma1 = suma1 + factor1;
               dsuma1 = dsuma1 - dfactor1;
          }
     
         ie = ie + suma1;
         die = die + dsuma1;
     
         func = omega_tmp - valbeta[orbital] - ie;
         dfunc = 1.f - die;
     
         omega = omega_tmp - func/dfunc;
     
         diff = fabs(omega - omega_tmp);
         cont++;
         printf("IP, iter  = %f, %d\n", -omega, cont);
     
      } while (diff > 1.e-6 && cont < 100);
     
      printf("Koopmans theorem     = %10.4f (%10.4f eV)\n", -valbeta[orbital], -valbeta[orbital]*au_to_eV);
      printf("Ionization potential = %10.4f (%10.4f eV) \n", -omega, -omega*au_to_eV);
     
      *total_energy = omega;
     
     
      time_4 = time (NULL);
      printf("EP2 time            = %5.2f s\n", time_4 - time_1);

      }//label beta


  return 0;
 

 }//label 1

