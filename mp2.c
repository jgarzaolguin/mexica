#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#ifdef _OPENMP
 #include <omp.h>
#else
 #define omp_get_thread_num() 0
#endif

double ao_to_mo_cpu(int     nt, 
                    int     i, 
                    int     a, 
                    int     j, 
                    int     b, 
                    double *vectalfa, 
                    double *vectbeta,
                    double *integrales_bie)
{
  int cont, indice_i, indice_j, indice_k, indice_l, interno, cuad, cubo, cuart;
  double orbital_a, orbital_s, orbital_b, orbital_r, coef, suma;

  cuad = nt*nt;
  cubo = cuad*nt;
  cuart = cubo*nt;

  #pragma omp parallel shared(nt, cuad, cubo, cuart, i, a, j, b, vectalfa, vectbeta, integrales_bie, suma) private(cont, indice_i, indice_j, indice_k, indice_l, coef, orbital_a, orbital_s, orbital_b, orbital_r, interno)
  {
  suma = 0.f;
  #pragma omp for reduction(+:suma)
  for (cont = 0; cont < cuart; cont++) {
    indice_i = cont/cubo;
    indice_j = cont - indice_i*cubo;
    indice_j = indice_j/cuad;
    indice_k = cont - indice_i*cubo - indice_j*cuad;
    indice_k = indice_k/nt;
    indice_l = cont - indice_i*cubo - indice_j*cuad - indice_k*nt;
    orbital_a = vectalfa[nt*indice_i + i];
    orbital_s = vectalfa[nt*indice_j + a];
    orbital_b = vectbeta[nt*indice_k + j];
    orbital_r = vectbeta[nt*indice_l + b];
    coef = orbital_a*orbital_s*orbital_b*orbital_r;
    interno = indice_i + indice_j*nt + indice_k*cuad + indice_l*cubo;
    suma = suma + coef*integrales_bie[interno];
  }
 } //OpenMP done
  return(suma);
 }

int mp2_cpu(int     nt, 
            int     compara, 
            int     elecalfa, 
            int     elecbeta,
            double *vectalfa, 
            double *vectbeta,
            double *valalfa, 
            double *valbeta,
            double *integrales_bie, 
            double *arreglo_mo, 
            double *total_energy)
{
  int i, j, a, b, cont;
  int indice;
  int indice_i, indice_j, indice_k, indice_l;
  int cuad, cuart, cubo;
  int ncis, ncis_eq;

  double coef, orbital_a, orbital_b, orbital_s, orbital_r;
  double time_1, time_4, suma,  suma1, factor1, factor2, denom, suma_mp2;

 extern double ao_to_mo_cpu(int     nt, 
                            int     i, 
                            int     a, 
                            int     j, 
                            int     b, 
                            double *vectalfa, 
                            double *vectbeta,
                            double *integrales_bie);

  time_1 = time (NULL);

  cuad = nt*nt;
  cubo = nt*cuad;
  cuart = cuad*cuad;

  for (j = 0; j < cuart; j++) 
    arreglo_mo[j] = (double) 0.;

// Viene la suma sobre alfas

  suma_mp2 = (double) 0.;



  int interno;
  int integrales = 0;
  for (i = 0; i < elecalfa; i++)
    for (a = elecalfa; a < nt; a++)
      for (j = 0; j < elecalfa; j++)
        for (b = elecalfa; b < nt; b++) {
          if ((i == j && a <= b) || (i < j)) {
           suma = ao_to_mo_cpu(nt, i, a, j, b, vectalfa, vectalfa, integrales_bie);
            ncis = i + a*nt + j*cuad + b*cubo;
            arreglo_mo[ncis] = suma;
          } else {
           ncis_eq = j + b*nt + i*cuad + a*cubo;
           ncis = i + a*nt + j*cuad + b*cubo;
           arreglo_mo[ncis] = arreglo_mo[ncis_eq];
         }
        }
// Eval'uo (a_alpha j_alpha|b_alpha i_alpha)
  for (i = 0; i < elecalfa; i++)
    for (a = elecalfa; a < nt; a++)
      for (j = 0; j < elecalfa; j++)
        for (b = elecalfa; b < nt; b++) {
            ncis = a + j*nt + b*cuad + i*cubo;
            ncis_eq = j + a*nt + i*cuad + b*cubo;
            arreglo_mo[ncis] = arreglo_mo[ncis_eq];
      }


  if (compara == 0) {
// Si es de capa cerrada solamente har'a esta parte y terminar'a
    suma1 = 0.f;
    for (i = 0; i < elecalfa; i++)
      for (a = elecalfa; a < nt; a++)
        for (j = 0; j < elecalfa; j++)
          for (b = elecalfa; b < nt; b++) {
            ncis = i + a*nt + j*cuad + b*cubo;
            factor1 = arreglo_mo[ncis];
            ncis = a + j*nt + b*cuad + i*cubo;
            factor2 = arreglo_mo[ncis];
            denom = valalfa[i] + valalfa[j] - valalfa[a] - valalfa[b];
            suma1 = suma1 + factor1*(2.f*factor1 - factor2)/denom;
          }
    suma_mp2 = suma1;
  } else {
// Esta parte es para capa abierta
//    Primero hago la suma sobre alfas para MP2
      suma1 = 0.f;
      for (i = 0; i < elecalfa; i++)
        for (a = elecalfa; a < nt; a++)
          for (j = 0; j < elecalfa; j++)
            for (b = elecalfa; b < nt; b++) {
              ncis = i + a*nt + j*cuad + b*cubo;
              factor1 = arreglo_mo[ncis];
              ncis = a + j*nt + b*cuad + i*cubo;
              factor2 = arreglo_mo[ncis];
              denom = valalfa[i] + valalfa[j] - valalfa[a] - valalfa[b];
              suma1 = suma1 + factor1*(factor1 - factor2)/denom;
            }  
      suma1 = (double) 0.5 * suma1;
      printf("Alpha contribution       = %f\n", suma1);
      suma_mp2 = suma_mp2 + suma1;
// Eval'uo (i_beta a_beta|j_beta b_beta)
  for (j = 0; j < cuart; j++) arreglo_mo[j] = (double) 0.;
      for (i = 0; i < elecbeta; i++)
        for (a = elecbeta; a < nt; a++)
          for (j = 0; j < elecbeta; j++)
            for (b = elecbeta; b < nt; b++) {
              if ((i == j && a <= b) || i < j) {
                suma = ao_to_mo_cpu(nt, i, a, j, b, vectbeta, vectbeta, integrales_bie);
                ncis = i + a*nt + j*cuad + b*cubo;
                arreglo_mo[ncis] = suma;
              } else {
                  ncis_eq = j + b*nt + i*cuad + a*cubo;
                  ncis = i + a*nt + j*cuad + b*cubo;
                  arreglo_mo[ncis] = arreglo_mo[ncis_eq];
              }
            }
// Eval'uo (a_beta j_beta|b_beta i_beta)
      for (i = 0; i < elecbeta; i++)
        for (a = elecbeta; a < nt; a++)
          for (j = 0; j < elecbeta; j++)
            for (b = elecbeta; b < nt; b++) {
              ncis = a + j*nt + b*cuad + i*cubo;
              ncis_eq = j + a*nt + i*cuad + b*cubo;
              arreglo_mo[ncis] = arreglo_mo[ncis_eq];
            }
//    Hago la suma sobre betas para MP2
      suma1 = 0.f;
      for (i = 0; i < elecbeta; i++)
        for (a = elecbeta; a < nt; a++)
          for (j = 0; j < elecbeta; j++)
            for (b = elecbeta; b < nt; b++) {
              ncis = i + a*nt + j*cuad + b*cubo;
              factor1 = arreglo_mo[ncis];
              ncis = a + j*nt + b*cuad + i*cubo;
              factor2 = arreglo_mo[ncis];
              denom = valbeta[i] + valbeta[j] - valbeta[a] - valbeta[b];
              suma1 = suma1 + factor1*(factor1 - factor2)/denom;
            }
      suma1 = (double) 0.5 * suma1;
      printf("Beta contribution        = %f\n", suma1);
      suma_mp2 = suma_mp2 + suma1;
// Evaluaci'on de integrales alfa-beta
// Eval'uo (i_alfa a_alfa|j_beta b_beta)
  for (j = 0; j < cuart; j++) arreglo_mo[j] = (double) 0.;
      for (i = 0; i < elecalfa; i++)
        for (a = elecalfa; a < nt; a++)
          for (j = 0; j < elecbeta; j++)
            for (b = elecbeta; b < nt; b++) {
                suma = ao_to_mo_cpu(nt, i, a, j, b, vectalfa, vectbeta, integrales_bie);
                ncis = i + a*nt + j*cuad + b*cubo;
                arreglo_mo[ncis] = suma;
          }
//    Hago la suma sobre alfa-beta para MP2
      suma1 = 0.f;
      for (i = 0; i < elecalfa; i++)
        for (a = elecalfa; a < nt; a++)
          for (j = 0; j < elecbeta; j++)
            for (b = elecbeta; b < nt; b++) {
              ncis = i + a*nt + j*cuad + b*cubo;
              factor1 = arreglo_mo[ncis];
              denom = valalfa[i] + valbeta[j] - valalfa[a] - valbeta[b];
              suma1 = suma1 + factor1*factor1/denom;
            }
      suma1 = (double) 0.5 * suma1;
      printf("Alpha-beta contribution  = %f\n", suma1);
      suma_mp2 = suma_mp2 + suma1;
// Evaluaci'on de integrales beta-alfa
// Eval'uo (i_beta a_beta|j_alfa b_alfa)
  for (j = 0; j < cuart; j++) arreglo_mo[j] = (double) 0.;
      for (i = 0; i < elecbeta; i++)
        for (a = elecbeta; a < nt; a++)
          for (j = 0; j < elecalfa; j++)
            for (b = elecalfa; b < nt; b++) {
                suma = ao_to_mo_cpu(nt, i, a, j, b, vectalfa, vectbeta, integrales_bie);
                ncis = i + a*nt + j*cuad + b*cubo;
                arreglo_mo[ncis] = suma;
          }
//    Hago la suma sobre beta-alfa para MP2
      suma1 = 0.f;
      for (i = 0; i < elecbeta; i++)
        for (a = elecbeta; a < nt; a++)
          for (j = 0; j < elecalfa; j++)
            for (b = elecalfa; b < nt; b++) {
              ncis = i + a*nt + j*cuad + b*cubo;
              factor1 = arreglo_mo[ncis];
              denom = valbeta[i] + valalfa[j] - valbeta[a] - valalfa[b];
              suma1 = suma1 + factor1*factor1/denom;
            }
      suma1 = (double) 0.5 * suma1;
      printf("Beta-alpha contribution  = %f\n", suma1);
      suma_mp2 = suma_mp2 + suma1;

  }  // Termina else para sistema de capa abierta

  time_4 = time (NULL);
  printf("MP2 correlation energy = %15.7f \n", suma_mp2);
  *total_energy = suma_mp2;

  printf("MP2 time in seconds    = %15.2f \n", time_4 - time_1);


  return 0;
  }

