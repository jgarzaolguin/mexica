/* Matrices para el Hartree-Fock atomico
   Jorge Garza, Junio del 2006 */
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


void imprimemat(int nt, double* mat)
{
 int i, j, k;

 k = -1;
 for (i = 0; i < nt; i++) {
    for (j = 0; j < nt; j++) {
      k++;
      printf("%8.6f  ", mat[k]);
    }
    printf("\n");
 }

 }

/* void multmat(int nt, double** mata, double** matb, double** matc)
This is used for the multiplication between two two-dimension matrixs 
{
 int i, j, k;
 double suma;

 for (i = 1; i <= nt; i++)
   for (j = 1; j <= nt; j++) {
     suma = 0.f;
     for (k = 1; k <= nt; k++)
       suma = suma + mata[i][k]*matb[k][j];
     matc[i][j] = suma;
   }
 } */

void multmat(int nt, double* mata, double* matb, double* matc)
{
/* C = A*B, where all these quantities are matrixs */
 int i, j, k, m, elemento1, elemento2;
 double suma;

 k = 0;
 for (i = 0; i < nt; i++)
   for (j = 0; j < nt; j++) {
     suma = (double)0;
     for (m = 0; m < nt; m++) {
       elemento1 = m*nt + i;
       elemento2 = j*nt + m;
       suma += mata[elemento1]*matb[elemento2];
     }
     matc[k] = suma;
     k++;
   }
 }

void summat(int tamanio, double alfa, double* mata, double beta, double* matb, double* matc)
{
/* Multiplicacion entre dos matrices en arreglos unidimensionales */
 int k;
 
 for (k = 0; k < tamanio; k++) matc[k] = alfa*mata[k] + beta*matb[k];

 }

void copiamat(int tamanio, double* mata, double* matb)
{
 int k;

 for (k = 0; k < tamanio; k++) matb[k] = mata[k];
 }

void eigensystem_matdens(int     nt, 
                         char   *tipo, 
                         int     elec, 
                         double *matx, 
                         double *matfock,
                         double *vals, 
                         double *vectsfin, 
                         double *matp,
                         int    *revisa_1)
{
 int dim, sizebi;
 double *vects, *newmatfock;

 int i, k, j, tempo;

 extern void newfock(int, double*, double*, double*);
// extern void valorespropios(int, double*, double*, double*);
 extern void valorespropios(int     nt,
                            double *mata,
                            double *valores,
                            double *vectores,
                            int    *revisa);

 extern void multmat(int, double*, double*, double*);
 extern void matdens(int, int, double*, double*, char*);
//mrb extern void memoria_double_uni(int, double**, char*);
 extern void memoria_double_uni(int      tamanio,
                                double **arreglo,
                                char    *titulo);
 extern void imprimemat(int , double* );

 dim = nt*nt;
 sizebi = dim*sizeof(double);

// Creo dos arreglos
 newmatfock =NULL;
 memoria_double_uni(sizebi, &newmatfock, "Newmatfock");

 vects =NULL;
 memoria_double_uni(sizebi, &vects, "vects");

 newfock(nt, matx, matfock, newmatfock);
 valorespropios(nt, newmatfock, vals, vects, &tempo);
 *revisa_1 = tempo;
// printf("PRO1 (%d)   ", tempo); 
 multmat(nt, matx, vects, vectsfin);
// Obtiene la matriz de Fock
 matdens(nt, elec, vectsfin, matp, tipo);
//Libero los dos arreglos

 free(vects);
 vects = 0;
 free(newmatfock);
 newmatfock = 0;
 }

void transforma(int     nt, 
                double *mats, 
                double *matx)
{
 int i, j, k, tempo;
 double *w, *test2;

// extern void valorespropios(int, double*, double*, double*);
 extern void valorespropios(int     nt,
                            double *mata,
                            double *valores,
                            double *vectores,
                            int    *revisa);


 w = (double *)malloc(nt*sizeof(double));
 if (w == NULL)
 {
  free(w);
  w = 0;
  printf("No hay suficiente memoria para w\n");
  exit(1);
  }

 test2 = (double *)malloc(nt*nt*sizeof(double));
 if (test2 == NULL)
 {
  free(test2);
  test2 = 0;
  printf("No hay suficiente memoria para test2\n");
  exit(1);
  }

 valorespropios(nt, mats, w, test2, &tempo);
// printf("\nPRO2 %d\n", tempo); 
 
 k = 0;
 for (i = 0; i < nt; i++)
   for (j = 0; j < nt; j++) {
     matx[k] = test2[k]/sqrt(w[i]);
     k++;
   }

 free(test2);
 test2 = 0;
 free(w);
 w = 0;
 }


void newfock(int nt, double* matx, double* matfock, double* newmatfock)
{
 int i, j, k, m, elemento1, elemento2, cont;
 double suma1, suma2 ;

 cont = 0;
 for (i = 0; i < nt; i++)
    for (j = 0; j < nt; j++) {
       suma1 = (double)0;
       for (k = 0; k < nt; k++) {
          suma2 = (double)0;
          for (m = 0; m < nt; m++) {
            elemento1 = j*nt + m;
            elemento2 = k*nt + m;
            suma2 += matfock[elemento2]*matx[elemento1];
          }
          elemento1 = i*nt + k;
          suma1 += matx[elemento1]*suma2; 
       }
       cont = i*nt + j;
       newmatfock[cont] = suma1;
    }
 }

void matdens(int nt, int elec, double* vectores, double* matp, char* tipo)
{
 int i, j, m, k, compara, peso, elemento1, elemento2;
 double suma;

 if (strcmp(tipo,"rhf") == 0 || strcmp(tipo,"rks") == 0) compara = 0;
 else compara = 1;

 if (compara == 0) peso = 2;
 else              peso = 1;

 k = 0;
 for (i = 0; i < nt; i++)
   for (j = 0; j < nt; j++)  {
     suma = (double)0;
     for (m = 0; m < elec/peso; m++) {
       elemento1 = i*nt + m;
       elemento2 = j*nt + m;
       suma += vectores[elemento1]*vectores[elemento2];
     }
     matp[k] = (double)peso*suma;
     k++;
   }
 }

void matdensespin(int nt, int elec, double vectores[][10], double matp[][10])
{
 int i, j, a;
 double suma;
                                                                                 
 
 for (i = 1; i <= nt; i++)
    for (j = i; j <= nt; j++)  {
       suma = 0.0;
       for (a = 1; a <= elec; a++)
          suma = suma + vectores[i][a]*vectores[j][a];
       matp[i][j] = suma;
    }
  for (i = 1; i <= nt; i++)
    for (j = 1; j < i; j++)
       matp[i][j] = matp[j][i];
 }

void matrixg(int nt, double* matp, double* matg, int* np, double* expo)
{
 int i, j, k, m, element_k, total_elements, elemento1, elemento2;

 int k2, cont, index_i, index_j, proc_come;
 int bloque, restante, do_ini, do_fin, total_size, sizebi;
 double part_real;
 double *array_temp, *array_temp_2;

 double sum1, sum2, num1, num2;

 extern int indexes(int, int, int*, int*);
 extern double bielectronic_integral_free(int, int, int, int, int, int*, double*);
 extern void memoria_double_uni(int, double**, char*);


 total_elements = nt*nt;

 cont = -1;
 for (k2 = 0; k2 < total_elements; k2++) {
   indexes(nt, k, &index_i, &index_j);
   index_i = index_i - 1;
   index_j = index_j - 1;
//   i = k2/nt;
//   part_real = fmod(k2,nt);
//   j = (int) part_real;
   i = index_i;
   j = index_j;
   sum1 = (double)0;
   for (k = 0; k < nt; k++) {
     sum2 = (double)0;
     for (m = 0; m < nt; m++)  {
       num1 =     bielectronic_integral_free(0, i, j, k, m, np, expo);
       num2 = (double)0.5*bielectronic_integral_free(0, i, m, k, j, np, expo);
       elemento1 = m*nt + k;
       sum2 += matp[elemento1]*(num1 - num2);
     }
     sum1 = sum1 + sum2;
   }
   cont++;
   matg[cont] = sum1;
 } 

 }

void matrixgespin(int nt, double matp[][10], double matpespin[][10], double matg[][10], int np[10], double expo[10])
{
 int i, j, k, m;
 double sum1, sum2, num1, num2;
 extern double bielectronic_integral_free(int, int, int, int, int, int [], double []);
                                                                                 
 for (i = 1; i <= nt; i++)
    for (j = i; j <= nt; j++) {
        sum1 = 0.0;
        for (k = 1; k <= nt; k++) {
           sum2 = 0.0;
           for (m = 1; m <= nt; m++)  {
              num1 = matp[m][k]*bielectronic_integral_free(0, i, j, k, m, np, expo);
              num2 = matpespin[m][k]*bielectronic_integral_free(0, i, m, k, j, np, expo);
              sum2 = sum2 + num1 - num2;
           }
           sum1 = sum1 + sum2;
        }
        matg[i][j] = sum1;
    }
                                                                                 
 for (i = 1; i <= nt; i++)
    for (j = 1; j < i; j++)
       matg[i][j] = matg[j][i];
 }

void valorespropios(int     nt, 
                    double *mata, 
                    double *valores, 
                    double *vectores,
                    int    *revisa)
{
 int i, j, k, lwork, info;
 double workopt;
 double *work;
 double *vectorestemp;

 extern void dsyev(char*, char*, int*, double*, int*, double*, double*, int*, int*);

 vectorestemp = (double *)malloc(nt*nt*sizeof(double));
 if (vectorestemp == NULL) {
   free(vectorestemp);
   vectorestemp = 0;
   printf("No hay suficiente memoria para vectorestemp\n");
   exit(1);
 }

/* Copio matriz A a vectorestemp */

//mrb k = 0;
//mrb for (i = 1; i <= nt; i++)
//mrb   for (j = 1; j <= nt;j++) {
//mrb     vectorestemp[k] = mata[k];
//mrb     k++;
//mrb   }
 
 for (i = 0; i < nt*nt; i++)
    vectorestemp[i] = mata[i];

 lwork =  -1;
 dsyev_ ("V", "U", &nt, vectorestemp, &nt, valores, &workopt, &lwork, &info);
 
 lwork = (int)workopt;
 work = (double *)malloc(lwork*sizeof(double));
 if (work == NULL) {
   free(work);
   work = 0;
   printf("No hay suficiente memoria para work\n");
   exit(1);
 }

 dsyev_ ("V", "U", &nt, vectorestemp, &nt, valores, work, &lwork, &info);

 *revisa = 0;
 if (info != 0) {
   printf("Problema en valores propios\n");
   *revisa = 1;
    //exit(1);
 }

 k = 0;
 for (i = 1; i <= nt; i++)
   for (j = 1; j <= nt;j++) {
     vectores[k] = vectorestemp[k];
     k++;
   }

 free(vectorestemp);
 vectorestemp = 0;
 free(work);
 work = 0;
}

void fase(int nt, double** vector1, double** vector2)
{
 int i, j;
 for (j = 1; j <= nt; j++)
 {
   if (vector1[1][j]*vector2[1][j] < (double)0)
     for (i = 1; i <= nt; i++)
       vector2[i][j] = -vector2[i][j];
  }
 }

void memoria_double_uni(int      tamanio, 
                        double **arreglo, 
                        char    *titulo)
{
 *arreglo = (double *)malloc(tamanio);
 if (arreglo == NULL) {
   free(arreglo);
   arreglo = 0;
   printf("No hay suficiente memoria para %s \n", titulo);
   exit(1);
 }
 }

void memoria_integer_uni(int tamanio, int** arreglo, char* titulo)
{
 *arreglo = (int *)malloc(tamanio);
 if (arreglo == NULL) {
   free(arreglo);
   arreglo = 0;
   printf("No hay suficiente memoria para %s \n", titulo);
   exit(1);
 }
 }

int indexes(int size_N, int elemento, int *indice_i, int *indice_j)
{
 int  i, j;

 i = elemento/size_N;
 j = elemento - i*size_N;

 *indice_i = i + 1;
 *indice_j = j + 1;

 return 0;
 }

void guarda_bielec(char *using_gamma, int nt, double *integrales_bie,
                   int* np, int* mang, int* ncm, double* expo, double Rc,
                   char* bound)
{
 int i, j, k, m, high_l, cont, k2, index_i, index_j, sizebi, element_k, proc_come;
 int total_elements, bloque, restante, do_ini, do_fin, total_size;
 
 double pi, *two_l_plus_1;

 extern double doselec(int, int, int, int, double, char*, double*, int*,
                       int*, int*, double*);

 pi = 4.f*atan(1.f);
 high_l = 50;
 two_l_plus_1 = (double *)malloc(high_l*sizeof(double));

 for (i = 0; i < high_l; i++)
   two_l_plus_1[i] = 4.f*pi/((double) 2*i + 1.f);

 total_elements = nt*nt;

 cont = -1;
 for(k2 = 0; k2 < total_elements; k2++) {
   index_i = k2/nt;
   index_j = k2 - index_i*nt;
   i = index_i;
   j = index_j;
   for(k = 0; k < nt; k++)
     for(m = 0; m < nt; m++) {
       cont++;
//       printf("proc= %d, %d, k2: %d,indices = %d, %d, %d, %d\n", myproc, cont, k2, i, j, k, m);
       integrales_bie[cont] = doselec(i, j, k, m, Rc, bound, expo, np, mang, ncm, two_l_plus_1);
//      printf("Ya tengo temp\n");
     }
  }

  free(two_l_plus_1);

 }
